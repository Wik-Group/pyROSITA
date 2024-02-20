"""
Class definitions for defining eROSITA specific data types.
"""
import os
import pathlib as pt

import astropy.units as u
import numpy as np
import requests
from astropy.coordinates import SkyCoord

from pyROS.erosita.sources import eRASS1Source
from pyROS.queries import Database, _included_erosita_columns
from pyROS.utils import EHalo, _enforce_style, mylog

_object_map = {"eRASS1": eRASS1Source}


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
    Class representation of an eROSITA catalog.
    """

    def __init__(self, filename, format="fits", catalog_type="eRASS1"):
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

        # properties initialization #
        self._coordinates = None
        self._objects = None

        # Loading the source type object
        self.source_type = _object_map[catalog_type]

    def __str__(self):
        return f"<eROSITA Catalog @ {self.filename}>"

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return len(self.data)

    @property
    def source_objects(self):
        if self._objects is None:
            self._objects = self.source_type.from_pandas(self.data)

        return self._objects

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

    def _xref_radii(self):
        """
        Calculates matched XREF radii for passing through the XREF generation procedures.
        """
        # calculating the sky radius to consider in arcmin^2.
        _ra_dec_error_radius = np.sqrt(
            np.amax(np.abs([self.RA_LOWERR, self.RA_UPERR]), axis=0) ** 2
            + (
                np.sin(np.deg2rad(self.DEC))
                * np.amax(np.abs([self.DEC_LOWERR, self.DEC_UPERR]), axis=0)
            )
            ** 2
        ) / (60)

        error_radius = 3 * _ra_dec_error_radius + (self.EXT / 60)

        # fixing nans
        error_radius[np.where(np.isnan(error_radius))] = 1.0

        return error_radius * u.arcmin

    def _construct_coordinates(self):
        return [
            SkyCoord(ra=kra * u.deg, dec=kdec * u.deg)
            for kra, kdec in zip(self.data.iloc[:]["RA"], self.data.iloc[:]["DEC"])
        ]

    def __getattr__(self, item):
        if hasattr(super(), item):
            return super().__getattribute__(self, item)
        elif item in self.columns:
            return np.array(self.data.loc[:, item])
        else:
            raise AttributeError(f"No such attribute {item}.")

    def xref(
        self,
        filename,
        groupsize=20,
        maxthreads=10,
        included_databases="all",
        overwrite=False,
        **kwargs,
    ):
        """

        Parameters
        ----------
        filename
        groupsize
        maxthreads
        kwargs

        Returns
        -------

        """
        import sqlalchemy as sql

        mylog.info(f"Generating XREF database for {len(self)} eROSITA records.")

        # Managing IO
        # ------------#
        if included_databases == "all":
            databases = [k.__name__ for k in Database.__subclasses__()]
        else:
            databases = included_databases

        if os.path.exists(filename):
            # the database file already exists; we are going to delete it.
            if overwrite:
                mylog.info(
                    f"Found existing XREF database at {filename}. Overwrite={overwrite}."
                )
                os.remove(filename)
            else:
                raise ValueError(
                    f"Found existing XREF database at {filename}. Overwrite={overwrite}."
                )
        else:
            # generate the necessary pathway.
            pt.Path(filename).parents[0].mkdir(parents=True, exist_ok=True)

        # Loop over DB
        # -------------#
        for database in databases:
            mylog.info(f"Generating XREF for [{database}]...")

            # create reference to the database
            db = Database(database)

            # Creating table in database
            # ---------------------------#
            _db_meta_data = sql.MetaData()
            cols = [
                sql.Column(k, v)
                for k, v in {
                    **{v[0]: v[1] for k, v in db.used_columns.items()},
                    **{v[0]: v[1] for k, v in _included_erosita_columns.items()},
                    "DELTA": sql.FLOAT,
                }.items()
            ]
            _ = sql.Table(f"XREF_{db.__class__.__name__}", _db_meta_data, *cols)
            _temp_engine = sql.create_engine(f"sqlite:///{filename}")
            _db_meta_data.create_all(_temp_engine)

            # Cross Referencing
            # ------------------#

            db.query_radius(
                self.source_objects,
                self._xref_radii(),
                maxworkers=maxthreads,
                group_size=groupsize,
                connection=filename,
            )

    def add_table_to_xref(self, database):
        import sqlalchemy as sql

        eng = sql.create_engine(f"sqlite:///{database}")
        with eng.connect() as conn:
            self.data.to_sql("eROSITA", con=conn, if_exists="replace")

    @_enforce_style
    def plot_field(
        self, field, ax=None, fig=None, y_scale="log", x_scale="linear", *args, **kwargs
    ):
        """
        Plot a field vs radius from this model using Matplotlib.

        Parameters
        ----------
        field : string
            The field to plot.
        fig : Matplotlib Figure
            The figure to plot in. Default; None, in which case
            one will be generated.
        ax : Matplotlib Axes
            The axes to plot in. Default: None, in which case
            one will be generated.
        y_scale: str
            The scaling on the y-axis.
        x_scale: str
            The scaling on the x-axis.
        Returns
        -------
        The Figure, Axes tuple used for the plot.

        """
        import matplotlib.pyplot as plt

        if fig is None:
            fig = plt.figure(figsize=(10, 10))
        if ax is None:
            ax = fig.add_subplot(111)

        ax.hist(getattr(self, field), *args, **kwargs)
        ax.set_yscale(y_scale)
        ax.set_xscale(x_scale)

        ax.set_ylabel("Catalog Source Count / [N]")

        ax.set_xlabel(f"{field}")

        return fig, ax


if __name__ == "__main__":
    cat = eROSITACatalog("/home/ediggins/pyROSITA_test/eRASS1_Hard.v1.0.fits")
    cat.xref(
        "./xref.db",
        groupsize=5,
        maxthreads=1,
        overwrite=True,
        included_databases=["NED", "SIMBAD"],
    )
