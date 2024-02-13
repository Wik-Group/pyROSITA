import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from itertools import repeat

import numpy as np
import pandas as pd
import sqlalchemy
from astroquery.ipac.ned import Ned
from astroquery.simbad import Simbad
from sqlalchemy import types as sql_type

from pyROSITA.utils import EHalo, devLogger, mylog, split

_dtypes = {
    "Object": sql_type.TEXT,
    "RA": sql_type.TEXT,
    "DEC": sql_type.TEXT,
    "Type": sql_type.TEXT,
    "UID": sql_type.BIGINT,
    "EXT": sql_type.FLOAT,
    "ERASS_RA": sql_type.TEXT,
    "ERASS_DEC": sql_type.TEXT
}

eras_columns = ["UID","EXT","ERASS_RA","ERASS_DEC"]


class Database:
    """Standard class for implementing different database types."""

    _astroquery_obj = None  # The object reference for building astroquery calls.
    _timeout_req = 0  #
    _max_threads_per_second = None
    name = "abstract"

    columns = None

    def __new__(cls, name, *args, **kwargs):
        """creates a new database given the standardized name."""
        _subcls_dict = {u.__name__: u for u in cls.__subclasses__()}
        if name in _subcls_dict:
            obj = object.__new__(_subcls_dict[name])
            return obj
        else:
            raise ValueError(
                f"The database {name} is not recognized as a valid database"
            )

    @classmethod
    def query_radius(cls, coordinates, radii, uid, ext, maxworkers=10, connection=None,maxgroup_size=20):
        """

        Parameters
        ----------
        coordinates
        radii

        Returns
        -------

        """
        assert len(coordinates) == len(
            radii
        ), f"Coordinates {len(coordinates)} must be same lenth as radii {len(radii)}."

        # sorting out the thread count.
        if cls._max_threads_per_second is not None:
            # a maximum thread rate is implemented to prevent overcalling the database.
            _mtps = np.amin([cls._max_threads_per_second, maxworkers])
        else:
            _mtps = maxworkers

        devLogger.info(
            f"Query to {cls.name} running on max threads {_mtps} with timeout {cls._timeout_req}. Call frequency approx. {_mtps*cls._timeout_req} / sec."
        )

        # split for the threads
        _coordinates, _radii, _uid, _ext = (
            split(coordinates, len(coordinates)//maxgroup_size),
            split(radii, len(radii)//maxgroup_size),
            split(uid, len(uid)//maxgroup_size),
            split(ext, len(ext)//maxgroup_size),
        )

        # manage the connection
        #
        # If a connection is provided, it needs to be passed to the multithread. Otherwise, we need to build an output system.
        #
        if connection is None:
            _output_array = []
            _conn = None
        else:
            _conn = connection
            _output_array = []

        with EHalo(
            text=f"Querying {cls.name} on {_mtps} threads with {len(coordinates)} calls [groups = {len(coordinates)//maxgroup_size}]."
        ):
            with ThreadPoolExecutor(max_workers=_mtps) as executor:
                results = executor.map(
                    cls._mtqrad,
                    _coordinates,
                    _radii,
                    _uid,
                    _ext,
                    repeat(cls._timeout_req),
                    repeat(_conn),
                )

            errors = []
            for result in results:
                _output_array.append(result[0])
                errors += result[1]

        mylog.info(
            f"Completed query. Found {len(errors)} errors from {len(coordinates)} queries."
        )
        if _conn is None:
            return pd.concat(_output_array)

    @classmethod
    def _mtqrad(cls, coordinates, radii,UID,EXT, timeout, connection):
        result = []
        errors = []
        devLogger.info(
            f"Running {len(coordinates)} calls on {threading.current_thread()}."
        )
        for _c, _r,_uid,_ext in zip(coordinates, radii,UID,EXT):
            timer = time.perf_counter()
            _exit_checker = False
            while not _exit_checker:
                _exit_checker = True
                with warnings.catch_warnings():
                    warnings.filterwarnings(action="ignore", module=".*")
                    try:
                        output_table = cls._astroquery_obj.query_region(_c, radius=_r)

                        if output_table is not None:
                            o = output_table.to_pandas()
                            o["UID"] = len(o) * [_uid]
                            o["EXT"] = len(o) * [_ext]
                            o["ERASS_RA"] = len(o) * [_c.ra.to_value('deg')]
                            o["ERASS_DEC"] = len(o) * [_c.dec.to_value('deg')]
                        else:
                            o = pd.DataFrame({k: [] for k in list(cls.columns.keys()) + eras_columns})


                        result.append(o)

                    except Exception as exp:
                        devLogger.error(exp.__str__())
                        errors.append(exp)

                while time.perf_counter() - timer < timeout:
                    time.sleep(timeout / 10)

        # reconstructing the outputs.
        if len(result):
            ret = pd.concat(result, ignore_index=True)
            ret = ret.loc[:, list(cls.columns.keys())+eras_columns]

            ret.rename(columns=cls.columns, inplace=True)
            if connection is not None:
                with threading.Lock():
                    engine = sqlalchemy.create_engine(f"sqlite:///{connection}")
                    with engine.connect() as conn:
                        devLogger.info(
                            f"Thread: {threading.current_thread()} -- Writing {len(ret)} lines to DB."
                        )
                        ret.to_sql(
                            f"XREF_{cls.name}",
                            con=conn,
                            if_exists="append",
                            index=False,
                            dtype={
                                k: v for k, v in _dtypes.items() if k in ret.columns
                            },
                        )
                devLogger.info(f"Thread: {threading.current_thread()} -- Completed.")

                return None, errors
            else:
                devLogger.info(
                    f"Thread: {threading.current_thread()} -- Returning {len(ret)} lines."
                )
                return ret, errors
        else:
            return None, errors


class NED(Database):
    _astroquery_obj = Ned  # The object reference for building astroquery calls.
    _timeout_req = 1  #
    _max_threads_per_second = 15
    name = "NED"

    columns = {"Object Name": "Object", "RA": "RA", "DEC": "DEC", "Type": "Type"}

    def __init__(self, *args, **kwargs):
        pass

    def __new__(cls, *args, **kwargs):
        obj = object.__new__(cls)
        return obj


class SIMBAD(Database):
    _astroquery_obj = Simbad  # The object reference for building astroquery calls.
    _timeout_req = 1  #
    _max_threads_per_second = 5
    name = "SIMBAD"

    columns = {"MAIN_ID": "Object", "RA": "RA", "DEC": "DEC"}

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
        [SkyCoord(ra=0, dec=0, unit=(u.deg, u.deg))], [15 * u.arcmin]
    )
    print(data.columns)
    print(data)
