# Herpich F.R. -- 2019-12-13
# Plots the S-PLUS footprint on the sky
# using the tiles_nc.csv file
# and the 2MASS all sky catalog
#
# Usage: python plot_footprint.py
#
# Output: splus_footprint.png
#
# Notes: The 2MASS catalog is available at
# http://cdsarc.u-strasbg.fr/viz-bin/Cat?II/281
#
# Last modified: 2023-07-19
# Modified by: Herpich F.R. CASU/IoA Cambridge
# email: fabio.herpich@ast.cam.ac.uk
# --------------------------------------------------

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier

query2mass = False
if query2mass:
    # TODO: test and implement querying to get the allsky 2MASS pointsource catalog
    print('Querying 2MASS catalog...')
    Vizier.ROW_LIMIT = 10000000
    vizier = Vizier(columns=['RAJ2000', 'DEJ2000', 'Hmag'])
    cat = vizier.query_constraints(catalog='II/246', Hmag=(0.0, 16.0))

    cat = Vizier.query_region(
        SkyCoord(ra=0, dec=0, unit=(u.deg, u.deg), frame='icrs'),
        radius=180 * u.deg, catalog='II/246')[0]

    c0 = SkyCoord(ra=cat['RAJ2000'], dec=cat['DEJ2000'],
                  unit=(u.deg, u.deg), frame='icrs')
    ra_rad0 = c0.ra.wrap_at(180 * u.deg).radian
    dec_rad0 = c0.dec.radian
    h, ex, ey = np.histogram2d(ra_rad0, dec_rad0, bins=(500, 500))
    xpos, ypos = np.meshgrid(ex[:-1], ey[:-1])
    h[h < 1] = 1.

print('reading splus table...')
t = ascii.read('tiles_nc.csv')

c = SkyCoord(ra=t['RA'], dec=t['DEC'], unit=(u.hour, u.deg), frame='icrs')
ra_rad = c.ra.wrap_at(180 * u.deg).radian
dec_rad = c.dec.radian

plt.figure(figsize=(16, 8.4))
ax = plt.subplot(111, projection="aitoff")
plt.grid(True)

if query2mass:
    print('ploting 2MASS...')
    ax.scatter(xpos, ypos, c=np.log10(h.T), s=10,
               marker='H', edgecolor='None', cmap='Greys')

mask_obs = (t['STATUS'] == 1) | (t['STATUS'] == 2) | (
    t['STATUS'] == 4) | (t['STATUS'] == 5) | (t['STATUS'] == 6)
lbl = r'$\mathrm{Observed:\ %i}$' % mask_obs.sum()
ax.scatter(ra_rad[mask_obs], dec_rad[mask_obs],
           marker='H', s=8, color='limegreen', label=lbl)
mask_foo = (t['STATUS'] == -1) | (t['STATUS'] == -
                                  2) | (t['STATUS'] == 0) | (t['STATUS'] == 3)
lbl = r'$\mathrm{To\ be\ observed:\ %i}$' % mask_foo.sum()
ax.scatter(ra_rad[mask_foo], dec_rad[mask_foo],
           marker='H', s=8, color='r', label=lbl)

plt.setp(ax.get_xticklabels(), fontsize=18)
plt.setp(ax.get_yticklabels(), fontsize=18)
plt.legend(loc='upper right', scatterpoints=1, markerscale=3, shadow=True,
           bbox_to_anchor=[1.02, 1.07], fancybox=True, fontsize=20)
plt.subplots_adjust(top=0.95, bottom=0.05, right=0.9, left=0.05)
plt.grid(True)

pathtosave = './splus_footprint.png'
print('saving fig %s...' % pathtosave)
plt.savefig(pathtosave, format='png', dpi=300)
plt.show()
