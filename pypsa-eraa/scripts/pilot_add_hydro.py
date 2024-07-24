#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pypsa
import yaml
import glob

from pypsa.descriptors import get_switchable_as_dense as as_dense

def add_inflows(network, config):

    short_names = config['carriers_short_names']
    basedir = snakemake.params.basedir

    inflow = pd.read_csv(snakemake.input.hydro_inflow_timeseries, header=0, index_col=0)

    hydro_gens = pd.DataFrame(index=inflow.columns)
    hydro_gens['carrier'] = [c.split(' - ')[1] for c in inflow.columns]
    hydro_gens['carrier_short'] = hydro_gens['carrier'].map(short_names)
    hydro_gens['node'] = [c.split(' - ')[0] for c in inflow.columns]
    hydro_gens['name_short'] = hydro_gens[['node', 'carrier_short']].agg(' '.join, axis=1)
    
    inflow.columns = inflow.columns.map(hydro_gens['name_short'])
    inflow.index = network.snapshots
    
    inflow = inflow.loc[:, inflow.columns[inflow.columns.isin(network.stores.index)]]
    buses = inflow.columns
    inflow.columns = [c+' inflow' for c in inflow.columns]
    
    network.add('Carrier', 'surface_runoff')
    network.madd(
        'Generator',
        inflow.columns,
        bus = buses,
        carrier = 'surface_runoff',
        p_nom = inflow.max(),
        p_max_pu = inflow / inflow.max(),
    )

def constraints_from_rule(network, rule, bounds):

    hydro_rules = pd.read_csv('pypsa-eraa/resources/hydro_rules.csv', header=0, index_col=0)

    if rule not in hydro_rules.index:
        print(f'Rule {rule} not defined, skipping')
        
    else:
        c = hydro_rules.component[rule]
        attr = hydro_rules.attribute[rule]
        direction = hydro_rules.direction[rule]
        sus = bounds.columns
    
        df = network.df(c)
        pnl = network.pnl(c)
    
        if direction == 'dispatch':
            idx = df.index[df.bus0.isin(sus)]
            set_attr = 'bus0'
        elif direction == 'store':
            idx = df.index[df.bus1.isin(sus)]
            set_attr = 'bus1'
        else:
            idx = sus
    
        if direction and hydro_rules.compute_pu[rule]:
            bounds /= df.set_index(df[set_attr]).p_nom[sus]
            bounds.columns = [f'{c} {direction}' for c in bounds.columns]
            bounds = bounds.abs().clip(upper=1.)

        if attr == 'e':
            pnl['e_min_pu'].loc[:,bounds.columns] = bounds.fillna(0.)-0.01
            pnl['e_max_pu'].loc[:,bounds.columns] = bounds.fillna(1.)+0.01
        else:
            # if component is already constraint, set bound to minimum of both constraints.
            pnl = pnl[attr]
            existing_bounds = pnl.reindex(columns=idx).fillna(bounds)
            bounds = bounds.where(
                bounds < existing_bounds[bounds.columns], 
                existing_bounds[bounds.columns]
            )
            pnl.loc[:,bounds.columns] = bounds
        
        print((rule, bounds.mean()))

def add_hydro_constraints_from_rules(network, config):

    short_names = config['carriers_short_names']
    basedir = snakemake.params.basedir

    hydro_files = sorted(glob.glob(basedir+f'hydro_uniform_*.csv'))

    for hydro_f in hydro_files[:17]:
        limits = pd.read_csv(hydro_f, index_col=0, header=0)
    
        hydro_sus = pd.DataFrame(index=limits.columns)
        hydro_sus['carrier'] = [c.split(' - ')[1] for c in limits.columns]
        hydro_sus['carrier_short'] = hydro_sus['carrier'].map(short_names)
        hydro_sus['node'] = [c.split(' - ')[0] for c in limits.columns]
        hydro_sus['name_short'] = hydro_sus[['node', 'carrier_short']].agg(' '.join, axis=1)
    
        limits.columns = limits.columns.map(hydro_sus['name_short'])
        limits.index = network.snapshots
        limits.drop(limits.columns[~limits.columns.isin(network.stores.index)], axis=1, inplace=True)
    
        constraint = hydro_f.split('_')[-1][:-4]

        constraints_from_rule(network, constraint, limits)

def add_hydro_constraints(network, config):

    short_names = config['carriers_short_names']
    basedir = snakemake.params.basedir

    hydro_files = sorted(glob.glob(basedir+f'hydro_uniform_*.csv'))

    for hydro_f in hydro_files:
        limits = pd.read_csv(hydro_f, index_col=0, header=0)
    
        hydro_sus = pd.DataFrame(index=limits.columns)
        hydro_sus['carrier'] = [c.split(' - ')[1] for c in limits.columns]
        hydro_sus['carrier_short'] = hydro_sus['carrier'].map(short_names)
        hydro_sus['node'] = [c.split(' - ')[0] for c in limits.columns]
        hydro_sus['name_short'] = hydro_sus[['node', 'carrier_short']].agg(' '.join, axis=1)
    
        limits.columns = limits.columns.map(hydro_sus['name_short'])
        limits.index = network.snapshots
        limits.drop(limits.columns[~limits.columns.isin(network.stores.index)], axis=1, inplace=True)
    
        constraint = hydro_f.split('_')[-1][:-4]
    
        if constraint == 'Maximum Generated energy MWh per week':
            ### Limits have been equally distributed to hours within week
            ### set p_max_pu of dispatching link
    
            disp_links = network.links.index[network.links.bus0.isin(hydro_sus.name_short)]
            max_disp_pu = limits / network.links.set_index(network.links.bus0).p_nom[limits.columns]
            max_disp_pu.columns = [c+' dispatch' for c in max_disp_pu.columns]
            network.links_t.p_max_pu.loc[:,max_disp_pu.columns] = max_disp_pu.clip(upper=1.)
    
        elif constraint == 'Maximum Generation MW':
            ### set p_max_pu of dispatching link
    
            disp_links = network.links.index[network.links.bus0.isin(hydro_sus.name_short)]
            max_disp_pu = limits / network.links.set_index(network.links.bus0).p_nom[limits.columns]
            max_disp_pu.columns = [c+' dispatch' for c in max_disp_pu.columns]

            # if link is already in p_max_pu, set p_max_pu to minimum of both constraints.
            existing_p_max_pu = network.links_t.p_max_pu.reindex(columns=disp_links).fillna(max_disp_pu)
            max_disp_pu = max_disp_pu.where(
                max_disp_pu < existing_p_max_pu[max_disp_pu.columns], 
                existing_p_max_pu[max_disp_pu.columns]
            )
            
            network.links_t.p_max_pu.loc[:,max_disp_pu.columns] = max_disp_pu.clip(upper=1.)

        elif constraint == 'Maximum Pumping MW':
            ### set p_max_pu of storing link
    
            store_links = network.links.index[network.links.bus1.isin(hydro_sus.name_short)]
            max_store_pu = limits / network.links.set_index(network.links.bus1).p_nom[limits.columns]
            max_store_pu.columns = [c+' store' for c in max_store_pu.columns]

            # if link is already in p_max_pu, set p_max_pu to minimum of both constraints.
            existing_p_max_pu = network.links_t.p_max_pu.reindex(columns=store_links).fillna(max_store_pu)
            max_store_pu = max_store_pu.where(
                max_store_pu < existing_p_max_pu[max_store_pu.columns],
                existing_p_max_pu[max_store_pu.columns]
            )

            network.links_t.p_max_pu.loc[:,max_store_pu.columns] = max_store_pu.clip(upper=1.)
    
        elif constraint == 'Maximum Reservoir level, historical (ratio)':
            ### set e_max_pu of store
            network.stores_t.e_max_pu[limits.columns] = limits
    
        elif constraint == 'Minimum Generated energy MWh per week':
            ### Limits have been equally distributed to hours within week
            ### set p_min_pu of dispatching link
    
            disp_links = network.links.index[network.links.bus0.isin(hydro_sus.name_short)]
            min_disp_pu = limits / network.links.set_index(network.links.bus0).p_nom[limits.columns]
            min_disp_pu.columns = [c+' dispatch' for c in min_disp_pu.columns]

            network.links_t.p_min_pu.loc[:,min_disp_pu.columns] = min_disp_pu.clip(upper=1.)
    
        elif constraint == 'Minimum Generation MW':
            ### set p_min_pu of dispatching link
    
            disp_links = network.links.index[network.links.bus0.isin(hydro_sus.name_short)]
            min_disp_pu = limits / network.links.set_index(network.links.bus0).p_nom[limits.columns]
            min_disp_pu.columns = [c+' dispatch' for c in min_disp_pu.columns]

            # if link is already in p_min_pu, set p_min_pu to maximum of both constraints.
            existing_p_min_pu = network.links_t.p_min_pu.reindex(columns=disp_links).fillna(min_disp_pu)
            min_disp_pu = min_disp_pu.where(
                min_disp_pu > existing_p_min_pu[min_disp_pu.columns], 
                existing_p_min_pu[min_disp_pu.columns]
            )
            
            network.links_t.p_min_pu.loc[:,min_disp_pu.columns] = min_disp_pu.clip(upper=1.)

        elif constraint == 'Minimum Reservoir level historical (ratio)':
            ### set e_min_pu of store
            network.stores_t.e_min_pu[limits.columns] = limits

#            e_initial = limits.iloc[0,:] * network.stores.e_nom[limits.columns]
#            network.stores.loc[limits.columns, 'e_initial'] = e_initial
#            network.stores.loc[limits.columns, 'e_initial_per_period'] = True
#            network.stores.loc[limits.columns, 'e_cyclic_per_period'] = False
    
        elif constraint == 'Minimum Reservoir level technical (ratio)':
            ### set e_min_pu of store
            network.stores_t.e_min_pu[limits.columns] = limits

#            e_initial = limits.iloc[0,:] * network.stores.e_nom[limits.columns]
#            network.stores.loc[limits.columns, 'e_initial'] = e_initial
#            network.stores.loc[limits.columns, 'e_initial_per_period'] = True
#            network.stores.loc[limits.columns, 'e_cyclic_per_period'] = False
    
        elif constraint == 'Reservoir level at beginning of week (ratio)':
            ### e_max_pu = e_min_pu at beginning of week
            e_min_pu = limits.fillna(0.)
            e_max_pu = limits.fillna(1.)
            network.stores_t.e_min_pu[e_min_pu.columns] = e_min_pu
            network.stores_t.e_max_pu[e_max_pu.columns] = e_max_pu

#            e_initial = limits.iloc[0,:] * network.stores.e_nom[limits.columns]
#            network.stores.loc[limits.columns, 'e_initial'] = e_initial
#            network.stores.loc[limits.columns, 'e_initial_per_period'] = True
#            network.stores.loc[limits.columns, 'e_cyclic_per_period'] = False

        else:
            print(f'Warning from add_hydro.py: hydro constraint not specified for constraint = {constraint}')

        # The following constraints are not modeled in ERAA
        if False:
            if constraint == 'Maximum Pumped Energy MWh per week':
                ### Limits have been equally distributed to hours within week
                ### set p_max_pu of storing link
        
                store_links = network.links.index[network.links.bus1.isin(hydro_sus.name_short)]
                max_store_pu = limits / network.links.set_index(network.links.bus1).p_nom[limits.columns]
                max_store_pu.columns = [c+' store' for c in max_store_pu.columns]
                network.links_t.p_max_pu[store_links] = max_store_pu.clip(upper=1.)
            
            elif constraint == 'Minimum Pumped energy MWh per week':
                ### Limits have been equally distributed to hours within week
                ### set p_min_pu of storing link
        
                store_links = network.links.index[network.links.bus1.isin(hydro_sus.name_short)]
                min_store_pu = limits / network.links.set_index(network.links.bus1).p_nom[limits.columns]
                min_store_pu.columns = [c+' store' for c in min_store_pu.columns]
                network.links_t.p_min_pu[store_links] = min_store_pu.clip(upper=1.)
    
            elif constraint == 'Minimum Pumping MW':
                ### set p_min_pu of storing link
        
                store_links = network.links.index[network.links.bus1.isin(hydro_sus.name_short)]
                min_store_pu = limits / network.links.set_index(network.links.bus1).p_nom[limits.columns]
                min_store_pu.columns = [c+' store' for c in min_store_pu.columns]
    
                # if link is already in p_min_pu, set p_min_pu to maximum of both constraints.
                existing_p_min_pu = network.links_t.p_min_pu
                min_store_pu = min_store_pu.where(
                    min_store_pu.reindex(columns=existing_p_min_pu.columns, fill_value=0.) > existing_p_min_pu, 
                    existing_p_min_pu
                ).fillna(min_store_pu)
                
                network.links_t.p_min_pu[store_links] = min_store_pu.clip(upper=1.)

def check_consistency(network, strategy='min_over_max'):
    ### check consistency of hydro constraints (partly taken from pypsa.components.consistency_check)
    max_pu = as_dense(network, 'Link', "p_max_pu", network.snapshots)
    min_pu = as_dense(network, 'Link', "p_min_pu", network.snapshots)

    diff = max_pu - min_pu
    diff = diff[diff < 0].dropna(axis=1, how="all")
    for col in diff.columns:
        if strategy == 'min_over_max':
            network.links_t.p_max_pu.loc[diff[col].dropna().index, col] = min_pu.loc[diff[col].dropna().index, col]
        elif strategy == 'max_over_min':
            network.links_t.p_min_pu.loc[diff[col].dropna().index, col] = max_pu.loc[diff[col].dropna().index, col]
            

if __name__ == "__main__":

    ### add hydro inflows and constraints and export network
    network = pypsa.Network(snakemake.input.network)
    add_inflows(network, snakemake.config)
#    add_hydro_constraints(network, snakemake.config)
    add_hydro_constraints_from_rules(network, snakemake.config)
    check_consistency(network)
    
    network.name = 'DestinE Pilot with VRE profiles and hydro'
    network.export_to_netcdf(snakemake.output[0])

