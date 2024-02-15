"""
Module for managing API queries both for eROSITA and for other databases for cross referencing.
"""
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from itertools import repeat

import numpy as np
import pandas as pd
import requests.exceptions
import sqlalchemy
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.ipac.ned import Ned
from astroquery.simbad import Simbad
from tqdm.auto import tqdm

from pyROS.utils import EHalo, devLogger, mylog, split,cgparams


Simbad.add_votable_fields("otype")

_included_erosita_columns = cgparams["REFDB"]['eRASS1_Catalog']['schema']
_included_erosita_columns_mod = [v[0] for v in _included_erosita_columns.values()]


class Database:
    """
    Abstract representation class for a standard reference database (NED, SIMBAD, etc.)
    """

    # -- backend -- #
    astroquery_obj = None          #: Corresponding astroquery object.
    max_call_frequency = 0         #: The maximum frequency of calls (calls/sec)
    max_threads = 0                #: Maximum number of allowed threads.
    used_columns = {}              #: The columns to include for this database reference class.
    retries = 5                    #: number of retries to give.

    def __new__(cls, name, *args, **kwargs):
        """Construct a database object from a given name."""
        _subcls_dict = {u.__name__: u for u in cls.__subclasses__()}

        if name in _subcls_dict:
            obj = object.__new__(_subcls_dict[name])
            return obj
        else:
            raise ValueError(
                f"Database {name} is not a recognized reference database."
            )

    @classmethod
    def query_radius(
        cls,
        objects,
        radii,
        maxworkers=10,
        connection=None,
        group_size=20,
    ):
        """
        Query this database at a given radius around a standardized eROSITA object.

        Parameters
        ----------
        objects: list
            A list of ``pyROS`` objects to query over.
        radii: list
            A list (shape matching that of ``objects``) containing the search radii for each object.
        maxworkers: int
            The maximum number of allowed threads. This is superceded by the maximum allowed workers for the database.
        connection: str
            The database URL to which data should be written.
        group_size: int
            The minimum target group size for query batching. In general, each thread is given objects in batches of this size.

        Returns
        -------
        None
        """
        # Processing the sizes and setting up variables
        #----------------------------------------------#
        if len(radii) == 1:
            radii = [radii]*len(objects)

        assert len(radii) == len(objects), f"There must be the same number of radii {len(radii)} as supplied objects {len(objects)}."

        _obj_coords = [SkyCoord(ra=o.RA*u.deg,dec=o.DEC*u.deg,frame="icrs") for o in objects]
        _obj_data = [[getattr(o,col) for col in _included_erosita_columns] for o in objects]

        # Preparing threading
        #--------------------#
        if cls.max_threads is not None:
            _mtps = np.amin([cls.max_threads, maxworkers])
        else:
            _mtps = maxworkers

        _threading_timeout_threshold = (1/cls.max_call_frequency)*_mtps

        devLogger.info(
            f"Query to {cls.__name__} running on max threads {_mtps} with max call frequency {cls.max_call_frequency}. Timeout is {_threading_timeout_threshold} seconds."
        )

        _staged_object_coordinates = split(_obj_coords,(len(_obj_coords)//group_size) + 1)
        _staged_radii              = split(radii,(len(_obj_coords)//group_size) + 1)
        _staged_object_data        = split(_obj_data,(len(_obj_data)//group_size) + 1)

        progress_bar = tqdm(total=len(objects),desc=f"XREF from {cls.__name__}",leave=True,disable=(not cgparams['system']['display']['progress_bars']))

        # Passing to threadpool
        #----------------------#
        with ThreadPoolExecutor(max_workers=_mtps) as executor:
            results = executor.map(
                cls._mtqrad,
                _staged_object_coordinates,
                _staged_radii,
                _staged_object_data,
                repeat(_threading_timeout_threshold),
                repeat(connection),
                repeat(progress_bar)
            )

        output_errors = 0
        for res in results:
            output_errors += len(res)

        mylog.info(
            f"Completed query. Found {output_errors} errors from {len(objects)} queries."
        )

    @classmethod
    def _mtqrad(cls, coordinates, radii, object_data, timeout, connection, progress_bar):
        devLogger.debug(f"Querying {cls.__name__} for {len(coordinates)} objects [Thread={threading.current_thread().name}].")

        _thread_return_array = []
        _thread_error_array = []
        # Coordinate Cycling
        #-------------------#
        for o_coord,o_rad,o_data in zip(coordinates,radii,object_data):
            _safe_proceed_flag = True # hold for a valid result.
            _safe_proceed_counter = 0 # the counter to check for infinite loop.
            while _safe_proceed_flag:
                _safe_clock_start = time.perf_counter()
                _safe_proceed_flag = False

                with warnings.catch_warnings():
                    warnings.filterwarnings(action="ignore", module=".*")

                    ## -- Call Procedure -- #
                    try:
                        output_table = cls.astroquery_obj.query_region(o_coord,radius=o_rad)
                    except requests.exceptions.ConnectTimeout as timeout_exception:
                        _safe_proceed_counter += 1

                        if _safe_proceed_counter <= cls.retries:
                            _safe_proceed_flag = True
                            continue
                        else:
                            _thread_error_array.append(o_data)
                            continue

                    ## -- Table Manipulation -- #
                    if output_table is None:
                        continue #--> there was no response, we just proceed as normal.

                    o = output_table.to_pandas()

                    for header,value in zip(_included_erosita_columns_mod,o_data):
                        # add the eRASS data.
                        o[header] = [value]*len(o)

                    _thread_return_array.append(o)

                while time.perf_counter() - _safe_clock_start < timeout:
                    # WAIT FOR TIMEOUT REASONS IF NECESSARY.
                    time.sleep(timeout/5)

        # Post Processing
        #----------------#
        if len(_thread_return_array) != 0:
            # there are output arrays to manage.
            output_dataframe = pd.concat(_thread_return_array,ignore_index=True)

            # DataFrame Manipulations

            output_dataframe = output_dataframe.loc[:, list(cls.used_columns.keys()) + _included_erosita_columns_mod]
            output_dataframe.rename(columns={k:v[0] for k,v in cls.used_columns.items()}, inplace=True)

            # Saving the output to the SQL database.
            with threading.Lock():
                engine = sqlalchemy.create_engine(f"sqlite:///{connection}")
                with engine.connect() as conn:
                    devLogger.info(
                        f"Writing {len(output_dataframe)} to {connection}. [Thread={threading.current_thread().name}]"
                    )
                    output_dataframe.to_sql(
                        f"XREF_{cls.__name__}",
                        con=conn,
                        if_exists="append",
                        index=False,
                        dtype={
                            ** {k:v[1] for k,v in cls.used_columns.items()},
                            ** {v[0]:v[1] for k,v in _included_erosita_columns.items()}
                        },
                    )
            devLogger.info(f"Thread: {threading.current_thread()} -- Completed.")

        with threading.Lock():
            progress_bar.update(n=len(coordinates))

        return _thread_error_array



class NED(Database):
    # -- backend -- #
    astroquery_obj     = Ned                                             #: Corresponding astroquery object.
    max_call_frequency = cgparams['REFDB']['NED']['max_call_frequency']  #: The maximum frequency of calls (calls/sec)
    max_threads        = cgparams['REFDB']['NED']['max_threads']         #: Maximum number of allowed threads.
    used_columns       = cgparams['REFDB']['NED']['schema']              #: The columns to include for this database reference class.
    retries            = cgparams['REFDB']['NED']['retries']             #: number of retries to give.


    def __init__(self, *args, **kwargs):
        pass

    def __new__(cls, *args, **kwargs):
        obj = object.__new__(cls)
        return obj


class SIMBAD(Database):
    # -- backend -- #
    astroquery_obj     = Simbad                                             #: Corresponding astroquery object.
    max_call_frequency = cgparams['REFDB']['SIMBAD']['max_call_frequency']  #: The maximum frequency of calls (calls/sec)
    max_threads        = cgparams['REFDB']['SIMBAD']['max_threads']         #: Maximum number of allowed threads.
    used_columns       = cgparams['REFDB']['SIMBAD']['schema']              #: The columns to include for this database reference class.
    retries            = cgparams['REFDB']['SIMBAD']['retries']             #: number of retries to give.

    def __init__(self, *args, **kwargs):
        pass

    def __new__(cls, *args, **kwargs):
        obj = object.__new__(cls)
        return obj

if __name__ == "__main__":
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    db = Database("SIMBAD")

    data = db.query_radius(
        [SkyCoord(ra=0, dec=0, unit=(u.deg, u.deg))],
        [15 * u.arcmin],
        ["A"],
        ["A"],
        maxworkers=1,
        maxgroup_size=1,
    )
