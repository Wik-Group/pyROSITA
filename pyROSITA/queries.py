import pandas as pd
from astroquery.ipac.ned import Ned
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from pyROSITA.utils import mylog
import warnings
columns = {
      "NED":[
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
    "SIMBAD":['MAIN_ID', 'RA', 'DEC', 'RA_PREC', 'DEC_PREC', 'COO_ERR_MAJA',
       'COO_ERR_MINA', 'COO_ERR_ANGLE', 'COO_QUAL', 'COO_WAVELENGTH',
       'COO_BIBCODE', 'SCRIPT_NUMBER_ID', "UID",
          "IAUNAME",
          "offset", 'source_db']
}
def pull_NED_radius(coords,radius):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pull = Ned.query_region(
            coords,
            radius=radius,
        )
        if pull is not None:
            o = pull.to_pandas()
        else:
            o = pd.DataFrame({k: [] for k in columns["NED"]})
        o['source_db'] = ["NED"] * len(o)
        return o
def pull_SIMBAD_radius(coords,radius):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                pull = Simbad.query_region(
                    coords,
                    radius=radius,
                )
            except UserWarning:
                pass
            if pull is not None:
                o = pull.to_pandas()
            else:
                o = pd.DataFrame({k:[] for k in columns["SIMBAD"]})
            o['source_db'] = ["NED"] * len(o)
            return o
    except Exception as exp:
        mylog.critical(exp)

if __name__ == '__main__':
    import astropy.units as u
    data = pull_SIMBAD_radius(SkyCoord(ra=0,dec=0,unit=(u.deg, u.deg)),radius=15*u.arcmin)
    print(data.columns)
    print(data)

