"""
Scripts for data management procedures
"""
import os
import pathlib as pt

import astropy.units as u
import numpy as np
import requests
from astropy.coordinates import SkyCoord

from pyROSITA.queries import Database, _dtypes, eras_columns
from pyROSITA.utils import EHalo, mylog, devLogger


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
        self.db_file = _db_file
        error_radii = self._coord_radius(factor=1)

        self.data["ERR_RAD"] = error_radii

    def __len__(self):
        return len(self.data)

    def columns(self):
        return self.data.columns

    @property
    def coordinates(self):
        return [
            SkyCoord(ra=kra * u.deg, dec=kdec * u.deg)
            for kra, kdec in zip(self.data.iloc[:]["RA"], self.data.iloc[:]["DEC"])
        ]

    @property
    def ext(self):
        return np.array(self.data.iloc[:]["EXT"].values)

    @property
    def uid(self):
        return np.array(self.data.iloc[:]["UID"].values)

    def _coord_radius(self, factor=1):
        # fetching coordinates and errors from the catalog.
        ra, dec, erau, eral, edecu, edecl, ext = tuple(
            [
                np.array(self.data[key])
                for key in [
                    "RA",
                    "DEC",
                    "RA_UPERR",
                    "RA_LOWERR",
                    "DEC_UPERR",
                    "DEC_LOWERR",
                    "EXT",
                ]
            ]
        )

        era = np.array([erau, eral])
        edec = np.array([edecu, edecl])
        err_radius = (
            np.sqrt(
                np.amax(np.abs(edec), axis=0) ** 2
                + (np.sin(dec) * np.amax(np.abs(era), axis=0)) ** 2
            )
            + ext
        )

        return factor * err_radius * u.arcsec

    def cross_reference(self, databases=None):
        """

        Parameters
        ----------
        search_locations

        Returns
        -------

        """
        import sqlalchemy as sql

        mylog.info(f"Generating a cross reference for {len(self)} eROSITA records.")

        if databases is None:
            databases = ["NED", "SIMBAD"]
        # --> passing to the thread manager.
        outs = []
        for database in databases:
            mylog.info(f"Generating XREF[{database}]...")

            db = Database(database)

            mylog.info(f"Generating table XREF_{database} in {self.db_file}.")

            meta_data = sql.MetaData()
            cols = [
                    sql.Column(k, v)
                    for k, v in _dtypes.items()
                    if (k in list(db.columns.values())) or (k in eras_columns)
                ]
            _ = sql.Table(
                f"XREF_{database}",
                meta_data,
                *cols,
            )
            devLogger.debug(f"Table XREF_{database} has {cols}.")
            _temp_engine = sql.create_engine(f"sqlite:///{self.db_file}")
            meta_data.create_all(_temp_engine)

            outs.append(
                db.query_radius(
                    self.coordinates,
                    self._coord_radius(factor=1),
                    self.uid,
                    self.ext,
                    connection=self.db_file,
                    maxworkers=10,
                )
            )

        return outs


if __name__ == "__main__":
    pass
