#!/usr/bin/env python
# coding: utf-8

import glob
import pandas as pd
import numpy as np
import xarray as xr
import pypsa
import os

TY = '2030'
odir = '../results/'
BASEDIR = '../results/'

idirs_pilot = sorted(glob.glob(BASEDIR+f'pilot_*TY{TY}_????.nc'))

for icy, pilot_f in enumerate(idirs_pilot):

    pilot = pypsa.Network(pilot_f)
    CY = pilot.snapshots.year[0]
    print(f'CY = {CY}')

    Reading_f = glob.glob(BASEDIR+f'pilot_*TY{TY}_{CY}_Reading.nc')[0]
    ninja_f = glob.glob(BASEDIR+f'pilot_*TY{TY}_{CY}_ninja.nc')[0]

    shedding_gens = pilot.generators.index[pilot.generators.carrier=='Load']

    idat_pilot = pilot.generators_t.p[shedding_gens]
    idat_Reading = pypsa.Network(Reading_f).generators_t.p[shedding_gens]
    idat_ninja = pypsa.Network(ninja_f).generators_t.p[shedding_gens]
    
    for gen_name in shedding_gens:
        country = pilot.generators.bus[gen_name]
        ofile = f'{odir}/{country}_shedding_TY{TY}.csv'

        df = pd.DataFrame(index=pilot.snapshots, dtype=float)
        df['PECD'] = idat_pilot[gen_name]
        df['Reading'] = idat_Reading[gen_name]
        df['ninja'] = idat_ninja[gen_name]
        
        if not os.path.exists(ofile):
            df.to_csv(ofile)
        else:
            df.to_csv(ofile, mode='a', header=False)