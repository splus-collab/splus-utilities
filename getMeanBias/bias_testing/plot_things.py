# F R H 2022-03-22 herpich@usp.br

from astropy.io import ascii
from astropy.time import Time
import matplotlib.pyplot as plt

t = ascii.read('bias_data_20170901-20220320.csv')
nt = ascii.read('bias_data_20220324.csv')

dates = Time([repr(d)[:4] + '-' + repr(d)[4:6] + '-' + repr(d)[6:] + 'T' +
              ('%06i'%j)[:2] + ':' + ('%06i'%j)[2:4] + ':' + ('%06i'%j)[4:]
              for d,j in zip(t['DATE'], t['TIME'])])
ndates = Time([repr(d)[:4] + '-' + repr(d)[4:6] + '-' + repr(d)[6:] + 'T' +
               ('%06i'%j)[:2] + ':' + ('%06i'%j)[2:4] + ':' + ('%06i'%j)[4:]
               for d,j in zip(nt['DATE'], nt['TIME'])])

mask = t['RMS'] > 50
per1 = dates.datetime < Time('2017-12-14T12:00:00').datetime
per2 = (dates.datetime > Time('2017-12-14T12:00:00').datetime)
per2 &= (dates.datetime < Time('2020-11-01T12:00:00').datetime)
per3 = dates.datetime < Time('2020-11-01T12:00:00').datetime

fig = plt.figure(figsize=(6,8))

ax = fig.add_subplot(511)
ax.scatter(dates.datetime, t['LEVEL'], marker='*')
ax.scatter(dates.datetime[mask], t['LEVEL'][mask], marker='*', color='r')
ax.scatter(ndates.datetime, nt['LEVEL'], marker='o', color='m')
ymin = min(t['LEVEL']) - abs(min(t['LEVEL']))*0.1
ymax = max(t['LEVEL']) + abs(max(t['LEVEL']))*0.1
ax.set_ylim(ymin, ymax)
ax.plot([Time('2017-12-14T12:00:00').datetime, Time('2017-12-14T12:00:00').datetime],
        [ymin, ymax], '-', c='r', lw=2)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid()
ax.set_ylabel('mean(C - mean(C))')


ax = fig.add_subplot(512)
ax.scatter(dates.datetime, t['RMS'], marker='*')
ax.scatter(dates.datetime[mask], t['RMS'][mask], marker='*', color='r')
ax.scatter(ndates.datetime, nt['RMS'], marker='o', color='m')
ymin = min(t['RMS']) - abs(min(t['RMS']))*0.1
ymax = max(t['RMS']) + abs(max(t['RMS']))*0.1
ax.set_ylim(ymin, ymax)
ax.plot([Time('2017-12-14T12:00:00').datetime, Time('2017-12-14T12:00:00').datetime],
        [ymin, ymax], '-', c='r', lw=2)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid()
ax.set_ylabel('rms')

ax = fig.add_subplot(513)
ax.scatter(dates.datetime, t['PG0_50'], marker='*')
ax.scatter(dates.datetime[mask], t['PG0_50'][mask], marker='*', color='r')
ax.scatter(ndates.datetime, nt['PG0_50'], marker='o', color='m')
ymin = min(t['PG0_50']) - abs(min(t['PG0_50']))*0.1
ymax = max(t['PG0_50']) + abs(max(t['PG0_50']))*0.1
ax.set_ylim(ymin, ymax)
ax.plot([Time('2017-12-14T12:00:00').datetime, Time('2017-12-14T12:00:00').datetime],
        [ymin, ymax], '-', c='r', lw=2)
plt.grid()
ax.set_ylabel('Gain S (V)')

ax = fig.add_subplot(514)
ax.scatter(dates.datetime, t['PG0_60'], marker='*')
ax.scatter(dates.datetime[mask], t['PG0_60'][mask], marker='*', color='r')
ax.scatter(ndates.datetime, nt['PG0_60'], marker='o', color='m')
ymin = min(t['PG0_60']) - abs(min(t['PG0_60']))*0.1
ymax = max(t['PG0_60']) + abs(max(t['PG0_60']))*0.1
ax.set_ylim(ymin, ymax)
ax.plot([Time('2017-12-14T12:00:00').datetime, Time('2017-12-14T12:00:00').datetime],
        [ymin, ymax], '-', c='r', lw=2)
plt.grid()
ax.set_ylabel('A2T (V)')

ax = fig.add_subplot(515)
ax.scatter(dates.datetime, t['PG5_40'], marker='*')
ax.scatter(dates.datetime[mask], t['PG5_40'][mask], marker='*', color='r')
ax.scatter(ndates.datetime, nt['PG5_40'], marker='o', color='m')
ymin = min(t['PG5_40']) - abs(min(t['PG5_40']))*0.1
ymax = max(t['PG5_40']) + abs(max(t['PG5_40']))*0.1
ax.set_ylim(ymin, ymax)
ax.plot([Time('2017-12-14T12:00:00').datetime, Time('2017-12-14T12:00:00').datetime],
        [ymin, ymax], '-', c='r', lw=2)
plt.grid()
ax.set_xlabel('Time')
ax.set_ylabel('Port 13 ADC Offset (V)')

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.05)
plt.show()
