from astroquery.ipac.ned import Ned
from astropy.coordinates import SkyCoord
def pull_NED_radius(coords,radius):
    pull = Ned.query_region(
        coords,
        radius=radius,
    )
    o = pull.to_pandas()
    o['source_db'] = ["NED"] * len(o)
    return o

