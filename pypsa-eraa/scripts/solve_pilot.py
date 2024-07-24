#!/usr/bin/env python
# coding: utf-8

import pypsa
import pandas as pd
import os
import numpy as np

def add_shedding(network):
    loads_i = network.loads.index
    buses = network.loads.bus.values
    names = [c+' shedding' for c in loads_i]
    p_nom = network.loads_t.p_set.max().values
    p_max_pu = (network.loads_t.p_set / p_nom).values
    
    network.add('Carrier', 'Load')
    network.madd('Generator', names, p_nom=p_nom, p_max_pu=p_max_pu, marginal_cost=1e3, bus=buses, carrier='Load')

def speed_up(network):

    network.loads_t.p_set = network.loads_t.p_set.round(0)
    network.generators_t.p_max_pu = network.generators_t.p_max_pu.round(3)

    HVAC_links = network.links.index[network.links.carrier=='HVAC']
    HVDC_links = network.links.index[network.links.carrier=='HVDC']
    network.links.efficiency.loc[HVAC_links] = 0.998
    network.links.efficiency.loc[HVDC_links] = 0.998
    network.links_t.p_max_pu = network.links_t.p_max_pu.round(3)
    network.links_t.p_min_pu = network.links_t.p_min_pu.round(3)

    ### noisy cost
    network.generators.marginal_cost[network.generators.index[network.generators.carrier=='solar']] = 0.1
    network.generators.marginal_cost[network.generators.index[network.generators.carrier=='onwind']] = 0.2
    network.generators.marginal_cost[network.generators.index[network.generators.carrier=='offwind']] = 0.25
    network.generators.marginal_cost += 1e-2 + 2e-3 * (
        np.random.random(len(network.generators.index)) - 0.5
    )
    network.generators.marginal_cost = network.generators.marginal_cost.clip(lower=0.)

def solve_network(network, snapshots, solver_name='glpk', solver_options={}):

    ### p_min_pu of dispatch links may cause infeasibilities when inflow is insufficient to serve the minimum dispatch requirement and stores energy is bounded by e_min_pu
    disp_links = network.links.index[network.links.carrier.isin(['CLPHS dispatch'])] #'hydro dispatch', 'OLPHS dispatch'])]
    disp_links_lower = disp_links[disp_links.isin(network.links_t.p_min_pu.columns)]
    network.links_t.p_min_pu[disp_links_lower] = 0

    ### allow for flexibility in state of charge at beginning of week to avoid infeasibilities -> consider adding a penalty term to objective
    network.stores_t.e_min_pu *= 0.8
    network.stores_t.e_min_pu = (network.stores_t.e_min_pu * 1.1).clip(upper=1.)
    
    ### fill NaNs at end of time series for NO stores
#    network.links_t.p_max_pu['NOM1 OLPHS dispatch'].fillna(method='ffill', inplace=True)
#    network.links_t.p_max_pu['NON1 OLPHS dispatch'].fillna(method='ffill', inplace=True)
#    network.links_t.p_max_pu['NOS0 OLPHS dispatch'].fillna(method='ffill', inplace=True)
    
    sol = network.optimize(snapshots, solver_name=solver_name, solver_options=solver_options, assign_all_duals=True)

    return((sol, network))

if __name__ == "__main__":

    status = pd.DataFrame()
    status['target_year'] = [snakemake.wildcards['target_year']]
    status['simulation_year'] = [snakemake.wildcards['simulation_year']]
    status['climatic_data_source'] = ['PECD']
    
    nf = snakemake.input[0]
    network = pypsa.Network(nf)    
    add_shedding(network)
    speed_up(network)
    
    status['status'] = ['running']

    if os.path.exists('pypsa_experiment_status.csv'):
        status.to_csv('pypsa_experiment_status.csv', mode='a', header=False, index=False)
    else:
        print('Creating status csv file')
        status.to_csv('pypsa_experiment_status.csv', index=False)
    
    sol, network_solved = solve_network(
        network, 
        network.snapshots[:-48], 
        solver_name=snakemake.config['solve']['solver_name'],
        solver_options=snakemake.config['solve']['solver_options']
    )
    network_solved.name = f'DestinE Pilot {sol[1]}'

    if sol[1] != 'optimal':
        network_solved.export_to_netcdf(snakemake.output[0].replace('.nc', '_infeasible.nc'))
        dummy = pypsa.Network()
        dummy.export_to_netcdf(snakemake.output[0])
    else:
        network_solved.export_to_netcdf(snakemake.output[0])

    status['status'] = sol[1]
    status.to_csv('pypsa_experiment_status.csv', mode='a', header=False, index=False)

    
