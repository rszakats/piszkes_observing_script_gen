#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
A simple script to create ccdsh script(s) for the 1m RCC
telescope at Piszkéstető Stn.
Example ccdsh sequence:

rcc autofocus
slew Herba
sequence -n Herba-%N-%F -V -x -j name=Herba 1000*([r,180,delay=5]) -u 22:20
"""

import os
import sys

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.visualization import LinearStretch, ZScaleInterval, imshow_norm
from astropy.wcs import WCS
from astroquery.jplhorizons import Horizons
from astroquery.skyview import SkyView


def getexptime(mag):
    if mag > 15.:
        exp = 180
    if mag < 15. and mag > 14.:
        exp = 120
    if mag < 14. and mag > 13.:
        exp = 90
    if mag < 13.:
        exp = 60
    return exp

def create_finder(date, start, end, astnum, targetname):
    if int(str(start).split(":")[0]) < 12.:
        jds = Time(f"{date} {start}:00").jd + 1
    else:
        jds = Time(f"{date} {start}:00").jd
    if int(str(end).split(":")[0]) < 12.:
        jde = Time(f"{date} {end}:00").jd + 1
    else:
        jde = Time(f"{date} {end}:00").jd
    jdm = (jds+jde)/2.
    obj = Horizons(id=f"{astnum}:", location='561', epochs=jds)
    eph = obj.ephemerides()
    if eph['targetname'][0].split()[1] != targetname:
        print(f"Error! Targetname mismatch!\n{eph['targetname'][0].split()[1]} {targetname}")
    ra = eph['RA'][0]
    dec = eph['DEC'][0]
    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    target = FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg),
                         name=eph['targetname'][0])
    stre = LinearStretch()
    inter = ZScaleInterval()
    hdu = SkyView.get_images(position=coord, coordinates='icrs',
                             survey='DSS2 Red', radius=15*u.arcmin,
                             grid=False)[0][0]
    wcs = WCS(hdu.header)
    coord1 = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    pixels = wcs.world_to_pixel(coord1)
    ax = plt.axes(projection=wcs)
    im, norm = imshow_norm(np.log(hdu.data), ax, origin='lower',
                           interval=inter, cmap='gray_r',
                           stretch=stre)

    ax.set(xlabel='RA', ylabel='DEC')
    plt.title(target.name+str(" ")+str(Time(jds, format='jd').isot))
    ax.grid(False, alpha=0)
    ax.scatter(pixels[0], pixels[1], s=50, edgecolor='blue', marker='o',
                facecolor='red')
    plt.savefig(f"scripts/{d}_{astnum}_{targetname}_finder_plot.png", dpi=150)
    plt.close()
    return ra, dec

start = f"set mount track on\n"\
        f"set dome slit open\n"\
        f"rcc dome auto\n"\
        f"rcc tubecover open\n"\
        f"!sleep 15\n"\
        f"rcc mirrorcover open\n"\
        f"!sleep 15\n"

close = f"set mount track off\n"\
        f"rcc dome manual\n"\
        f"rcc mirrorcover close\n"\
        f"!sleep 15\n"\
        f"rcc tubecover close\n"\
        f"!sleep 15\n"\
        f"rcc slew 0 47.8\n"\
        f"set dome slit close\n"\
        f"!sleep 90\n"\
        f"set dome azimuth=156.3\n"\
        f"!sleep 90\n"\
        f"!echo Done!"

infile = './input.dat'

if os.path.exists(infile):
    inp = Table.read(infile, format='ascii.csv')
else:
    print("No input file!")
    sys.exit(0)

if os.path.exists('scripts') is False:
    os.mkdir('scripts')

dates = np.unique(inp['date'])

for d in dates:
    dat = inp[inp['date'] == d]
    for dd in dat:
        exptime = getexptime(dd['mag'])
        target = dd['astname'].split()[1]
        astnum = dd['astname'].split()[0]
        ra, dec = create_finder(dd['date'], dd['ut_beg'], dd['ut_end'], astnum, target)
        if os.path.exists(f'scripts/{d}.ccdsh') is False:
            with open(f'scripts/{d}.ccdsh', "w") as ofile:
                towrite = f"# date    ut_beg ut_end  astname   mag  type\n"\
                          f"# {dd['date']} {dd['ut_beg']} {dd['ut_end']} {dd['astname']} {dd['mag']} {dd['type']}\n"\
                          f"rcc autofocus apply\n"\
                          f"{start}"\
                          f"!echo Observing {dd['astname']} until {dd['ut_end']} UT.\n"\
                          f"slew {target}\n"\
                          f"# slew ra={ra} dec={dec}  # in case name resolution doesn't work\n"\
                          f"wait\n"\
                          f"sequence -n {target}-%N-%F -V -x -j name={target} 1000*([r,{exptime},delay=3]) -u {dd['ut_end']}\n"
                ofile.write(towrite)
        else:
            with open(f'scripts/{d}.ccdsh', "a") as ofile:
                towrite = f"# date    ut_beg ut_end  astname   mag  type\n"\
                          f"# {dd['date']} {dd['ut_beg']} {dd['ut_end']} {dd['astname']} {dd['mag']} {dd['type']}\n"\
                          f"rcc autofocus apply\n"\
                          f"!echo sleep 15\n"\
                          f"!echo Observing {dd['astname']} until {dd['ut_end']} UT.\n"\
                          f"slew {target}\n"\
                          f"# slew ra={ra} dec={dec}  # in case name resolution doesn't work\n"\
                          f"wait\n"\
                          f"sequence -n {target}-%N-%F -V -x -j name={target} 1000*([r,{exptime},delay=3]) -u {dd['ut_end']}\n"\
                          f"{close}"
                ofile.write(towrite)
                