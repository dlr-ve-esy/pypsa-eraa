#!/usr/bin/env python
# coding: utf-8

import glob
import pandas as pd
import numpy as np
import xarray as xr
import pypsa
import os

TY = '2030'
odir = '../results/timeseries/'
BASEDIR = '../results/'

countries = ['DE', 'DK', 'ES']

idirs_pilot = sorted(glob.glob(BASEDIR+f'pilot_*TY{TY}_????.nc'))

for icy, pilot_f in enumerate(idirs_pilot):

    pilot = pypsa.Network(pilot_f)
    CY = pilot.snapshots.year[0]
    print(f'CY = {CY}')

    idat_pilot = pilot.generators_t.p

    alt_files = sorted(glob.glob(BASEDIR+f'pilot_*TY{TY}_{CY}_*.nc'))

    odat = {country: pd.DataFrame(index=pilot.snapshots, dtype=float) for country in countries}

    for alt_f in alt_files:

        alt = alt_f.split(f'_{CY}_')[-1].split('.')[0]

        idat_alt = pypsa.Network(alt_f).generators_t.p
        
        for country in countries:
    
            gen_name = pilot.generators.index[(pilot.generators.carrier=='Load') & (pilot.generators.bus==country)]
    
            odat[country]['PECD'] = idat_pilot[gen_name]
            odat[country][alt] = idat_alt[gen_name]
            
    for k,v in odat.items():
        ofile = f'{odir}/{k}_shedding_TY{TY}.csv'
        if not os.path.exists(ofile):
            v.to_csv(ofile)
        else:
            v.to_csv(ofile, mode='a', header=False)