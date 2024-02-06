"""
Scripts for data management procedures
"""
import os
import pathlib as pt
import threading
from itertools import repeat
from pyROSITA.queries import pull_NED_radius
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

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            executor.map(self._mt_cross_reference, repeat(self.data), indices,repeat(databases))

    def _mt_cross_reference(self, dataframe, indices,databases):
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
            for database in databases:
                if database == "NED":
                    subpull = pull_NED_radius(SkyCoord(ra=_ra, dec=_dec, unit=(u.deg, u.deg)),5 * _err_radius * u.arcsec)
                else:
                    mylog.error(f"The database {database} doesn't exist.")
                    continue
                if pull is None:
                    pull = subpull.loc[:,:]
                else:
                    pull = pd.concat([pull,subpull])
            out.append(pull)


        for i in range(len(out)):
            out[i]["UID"] = dataframe.iloc[indices[i]]["UID"]
            out[i]["IAUNAME"] = dataframe.iloc[indices[i]]["IAUNAME"]
            out[i]["offset"] = 5 * dataframe.iloc[indices[i]]["ERR_RAD"]

        # now we can join all of the databases.
        ret = pd.concat(out, ignore_index=True)
        ret = ret.loc[
            :,
            [
                "Object Name",
                "RA",
                "DEC",
                "Type",
                "Velocity",
                "Redshift",
                "References",
                "Associations",
                "UID",
                "IAUNAME",
                "offset",
                "source_db"
            ],
        ]

        with threading.Lock():
            mylog.info(
                f"Thread: {threading.current_thread()} -- Writing {len(ret)} lines to DB."
            )
            try:
                ret.to_sql(name="XREF_BASE", con=self.engine, if_exists="append")
            except Exception as exp:
                mylog.error(exp.__str__())
        mylog.info(f"Thread: {threading.current_thread()} -- Completed.")
