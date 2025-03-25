#!/bin/python

from astropy.utils import iers
import sys
from astropy.coordinates import EarthLocation, AltAz, get_moon, get_sun
from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
iers.conf.auto_download = False


if ('--nightstarts' in sys.argv) | ('-ns' in sys.argv):
    if ('--nightstarts' in sys.argv):
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--nightstarts']
    else:
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '-ns']
    night_starts = np.array(sys.argv)[index + 1].item()
    if len(night_starts.split('-')) < 3:
        raise IOError('wrong format for date')
    else:
        night_starts = Time(night_starts).value.split()[0]
else:
    raise IOError('night start not given')

if ('--nightends' in sys.argv) | ('-ne' in sys.argv):
    if ('--nightstarts' in sys.argv):
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '--nightends']
    else:
        index = np.arange(len(sys.argv))[np.array(sys.argv) == '-ne']
    night_ends = np.array(sys.argv)[index + 1].item()
    if len(night_ends.split('-')) < 3:
        raise IOError('wrong format for date')
    else:
        night_ends = Time(night_ends).value.split()[0]
else:
    raise IOError('final night not given')

f = open('moon_brightness.csv', 'w')

while Time('%s 23:59:59' % night_starts) < Time(night_ends):

    initime = Time('%s %s' % (night_starts, '23:59:59'))
    midnight = initime + 1*u.day
    delta_midnight = np.linspace(-12, 12, 500)*u.hour
    mysite = EarthLocation(lat=-30.2*u.deg, lon=-70.8*u.deg, height=2200*u.m)
    frame_tonight = AltAz(obstime=midnight+delta_midnight,
                          location=mysite)

    times_time_overnight = midnight + delta_midnight
    frame_time_overnight = AltAz(obstime=times_time_overnight, location=mysite)
    sunaltazs_time_overnight = get_sun(
        times_time_overnight).transform_to(frame_time_overnight)

    moon_time_overnight = get_moon(times_time_overnight)
    moonaltazs_time_overnight = moon_time_overnight.transform_to(
        frame_time_overnight)

    # moon phase
    sun_pos = get_sun(initime)
    moon_pos = get_moon(initime)
    elongation = sun_pos.separation(moon_pos)
    moon_phase = np.arctan2(sun_pos.distance*np.sin(elongation),
                            moon_pos.distance - sun_pos.distance*np.cos(elongation))
    moon_brightness = (1. + np.cos(moon_phase))/2.0

    unt = delta_midnight.value[sunaltazs_time_overnight.alt <= -18*u.deg].max()
    unt -= delta_midnight.value[sunaltazs_time_overnight.alt <= -18*u.deg].min()

    print(night_starts, '%.2f' % moon_brightness.value, '%.2f' % unt)
    f.write('%s,,%.2f,,,%.2f\n' % (night_starts, unt, moon_brightness.value))

    night_starts = (Time('%s 23:59:59' % night_starts) +
                    1*u.day).value.split()[0]

f.close()
