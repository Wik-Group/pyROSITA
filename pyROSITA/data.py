"""
Scripts for data management procedures
"""
import os
import pathlib as pt
import threading
from itertools import repeat
from pyROSITA.queries import pull_NED_radius,columns,pull_SIMBAD_radius
import astropy.units as u
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from sqlalchemy import create_engine

from pyROSITA.utils import EHalo, mylog, split


def download_data_product(url, output_location):
    """
    Download the data product from the specified url and save it to the output location (directory).

    Parameters
    ----------
    url
    output_location

    Returns
    -------

    """
    # the download must take place.
    mylog.info(f"Downloading from {url}.")
    with EHalo(text="Downloading..."):
        response = requests.get(url)
    # check run
    if not response.ok:
        raise ValueError(f"The download of {url} failed.")
    with open(os.path.join(output_location, pt.Path(url).name), mode="wb") as _f:
        _f.write(response.content)
    # unzip
    import tarfile

    with EHalo(text="Unzipping..."):
        with tarfile.open(os.path.join(output_location, pt.Path(url).name)) as _f:
            _f.extractall(output_location)


class eROSITACatalog:
    def __init__(self, filename, format="fits", _db_file=":memory:"):
        from astropy.table import Table

        self.data = Table.read(filename, format=format).to_pandas()
        self.engine = create_engine(f"sqlite:///{_db_file}")
        error_radii = self._coord_radius(factor=1)

        self.data["ERR_RAD"] = error_radii

    def __len__(self):
        return len(self.data)

    def columns(self):
        return self.data.columns

    def _coord_radius(self, factor=1):
        # fetching coordinates and errors from the catalog.
        ra, dec, erau, eral, edecu, edecl = tuple(
            [
                np.array(self.data[key])
                for key in [
                    "RA",
                    "DEC",
                    "RA_UPERR",
                    "RA_LOWERR",
                    "DEC_UPERR",
                    "DEC_LOWERR",
                ]
            ]
        )

        era = np.array([erau, eral])
        edec = np.array([edecu, edecl])
        err_radius = np.sqrt(
            np.amax(np.abs(edec), axis=0) ** 2
            + (np.sin(dec) * np.amax(np.abs(era), axis=0)) ** 2
        )

        return factor * err_radius

    def cross_reference(self, max_workers=1,databases = None,maxsize=20):
        """

        Parameters
        ----------
        search_locations

        Returns
        -------

        """
        mylog.info(f"Generating a cross reference for {len(self)} eROSITA records.")

        if databases is None:
            databases = ["NED"]
        # --> passing to the thread manager.
        indices = split(np.arange(len(self.data)), len(self.data)//maxsize)
        from concurrent.futures import ThreadPoolExecutor


        for database in databases:
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                executor.map(self._mt_cross_reference, repeat(self.data), indices,repeat(database))

    def _mt_cross_reference(self, dataframe, indices,database):
        mylog.info(
            f"Cross referencing {len(indices)} objects on {threading.current_thread()}."
        )

        out = []

        for _i, _p in enumerate(
            zip(
                list(dataframe.iloc[indices]["DEC"]),
                list(dataframe.iloc[indices]["RA"]),
                list(dataframe.iloc[indices]["ERR_RAD"]),
            )
        ):

            _ra, _dec, _err_radius = _p
            # pulling the data
            pull = None
            if database == "NED":
                pull = pull_NED_radius(SkyCoord(ra=_ra, dec=_dec, unit=(u.deg, u.deg)),5 * _err_radius * u.arcsec)
            elif database == "SIMBAD":
                pull = pull_SIMBAD_radius(SkyCoord(ra=_ra, dec=_dec, unit=(u.deg, u.deg)),5 * _err_radius * u.arcsec)
            else:
                mylog.error(f"The database {database} doesn't exist.")
                continue

            out.append(pull)


        for i in range(len(out)):
            out[i]["UID"] = dataframe.iloc[indices[i]]["UID"]
            out[i]["IAUNAME"] = dataframe.iloc[indices[i]]["IAUNAME"]
            out[i]["offset"] = 5 * dataframe.iloc[indices[i]]["ERR_RAD"]

        # now we can join all of the databases.
        ret = pd.concat(out, ignore_index=True)

        ret = ret.loc[:,columns[database]]

        with threading.Lock():
            mylog.info(
                f"Thread: {threading.current_thread()} -- Writing {len(ret)} lines to DB."
            )
            ret.to_sql(name=f"XREF_{database}", con=self.engine, if_exists="append")
        mylog.info(f"Thread: {threading.current_thread()} -- Completed.")
