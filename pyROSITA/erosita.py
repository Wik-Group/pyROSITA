"""
Class definitions for defining eROSITA specific data types.
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
    download a data product from eROSITA from a url.

    Parameters
    ----------
    url: str
        The URL to the data product to download.
    output_location: str
        The location at which to save the data product.

    Returns
    -------
    None

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
    """
    Class representation of eROSITA data catalogs.
    """
    def __init__(self, filename, format="fits"):
        """
        Initializes the :py:class:`eROSITACatalog` instance.

        Parameters
        ----------
        filename: str
            The name of the file in which the catalog is stored on disk.
        format: str, optional
            The format of the storage file. Default is ``fits``.
        """
        from astropy.table import Table
        self.filename = filename
        self.data = Table.read(filename, format=format).to_pandas()

        # property backbones
        self._coordinates = None

        error_radii = self.create_search_radii()
        self.data["ERR_RAD"] = error_radii

    def __str__(self):
        return f"<eROSITA Catalog @ {self.filename}>"
    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return len(self.data)

    def __getattr__(self, item):
        if hasattr(super(),item):
            return super().__getattribute__(self,item)
        elif item in self.columns:
            return np.array(self.data.loc[:,item])
        else:
            raise AttributeError(f"No such attribute {item}.")


    @property
    def columns(self):
        """The available columns"""
        return list(self.data.columns)

    @property
    def coordinates(self):
        """Sky coordinates for the objects in the database."""
        if self._coordinates is None:
            self._coordinates = self._construct_coordinates()

        return self._coordinates

    def create_search_radii(self,
                            account_for_extended_sources=True,
                            ):
        """
        creates the relevant search radii for the cross-catalog creation.
        Parameters
        ----------
        account_for_extended_sources: bool, default=``True``.

        Returns
        -------

        """
        # calculating the sky radius to consider in arcseconds^2.
        _ra_dec_error_radius = np.sqrt(
                np.amax(np.abs([self.RA_LOWERR,self.RA_UPERR]), axis=0) ** 2
                + (np.sin(np.deg2rad(self.DEC)) * np.amax(np.abs([self.DEC_LOWERR,self.DEC_UPERR]), axis=0)) ** 2)/(60)

        if account_for_extended_sources:
            error_radius = np.amax([_ra_dec_error_radius,self.EXT/60],axis=0)
        else:
            error_radius = _ra_dec_error_radius

        return error_radius

    def _construct_coordinates(self):
        return [
            SkyCoord(ra=kra * u.deg, dec=kdec * u.deg)
            for kra, kdec in zip(self.data.iloc[:]["RA"], self.data.iloc[:]["DEC"])
        ]


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
                    maxworkers=5,
                    maxgroup_size=5
                )
            )

        return outs

    def _add_table_to_xref(self,database):
        import sqlalchemy as sql
        eng = sql.create_engine(f"sqlite:///{database}")
        with eng.connect() as conn:
            self.data.to_sql("eROSITA",con=conn,if_exists='replace')




if __name__ == "__main__":
    database_directory = "/home/ediggins/pyROSITA_test"
    catalog_path = os.path.join(database_directory, "eRASS1_Hard.v1.0.fits")

    q = eROSITACatalog(catalog_path)
    q._add_table_to_xref(os.path.join(database_directory,"XREF.db"))
