# ==============================================================#
# Generates the XREF database from NED references              #
# ==============================================================#
import os
import pathlib as pt

from pyROSITA.erosita import eROSITACatalog
from pyROSITA.utils import mylog

# --------------------------------------------------------------#
# Settings
database_directory = "/scratch/general/vast/u1281896/EROSITA/XREF_FEB"
catalog_path = os.path.join(database_directory, "eRASS1_Hard.v1.0.fits")
first_n = None  # for debugging, only catalog the first N entries.
databases = "all"  # NYI
# --------------------------------------------------------------#
# Setup
if not os.path.exists(database_directory):
    mylog.info(f"{database_directory} doesn't currently exist. Creating it.")
    pt.Path(database_directory).mkdir(parents=True)

if catalog_path is None:
    from pyROSITA.erosita import download_data_product

    mylog.info("The catalog path was not specified. Downloading.")
    download_data_product(
        "https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Hard.tar.gz",
        database_directory,
    )
    catalog_path = os.path.join(database_directory, "eRASS1_Hard.v1.0.fits")

q = eROSITACatalog(catalog_path, _db_file=os.path.join(database_directory, "XREF.db"))

if first_n is not None:
    q.data = q.data.iloc[:first_n]

q.cross_reference(databases=["NED", "SIMBAD"])
