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
    disp_links = network.links.index[network.links.carrier.isin(['hydro dispatch', 'OLPHS dispatch', 'CLPHS dispatch'])]
    disp_links = disp_links[disp_links.isin(network.links_t.p_min_pu.columns)]
    network.links_t.p_min_pu[disp_links] = 0 ### solves!

    ### fill NaNs at end of time series for NO stores
    network.links_t.p_max_pu['NOM1 OLPHS dispatch'].fillna(method='ffill', inplace=True)
    network.links_t.p_max_pu['NON1 OLPHS dispatch'].fillna(method='ffill', inplace=True)
    network.links_t.p_max_pu['NOS0 OLPHS dispatch'].fillna(method='ffill', inplace=True)
    
    network.stores_t.e_max_pu['NO OLPHS OLPHS'][-48:] = network.stores_t.e_max_pu['NO OLPHS OLPHS'][-50]

    sol = network.optimize(snapshots, solver_name=solver_name, solver_options=solver_options, assign_all_duals=True)

    return((sol, network))

if __name__ == "__main__":

    nf = snakemake.input['network']
    ALT = snakemake.wildcards['ALT']
    
    network = pypsa.Network(nf)    
    add_shedding(network)
    network.remove('Link', 'LUB1-BE00 HVAC')
    speed_up(network)    
    
    sol, network_solved = solve_network(
        network, 
        network.snapshots[:-49], 
        solver_name=snakemake.config['solve']['solver_name'],
        solver_options=snakemake.config['solve']['solver_options']
    )
    network_solved.name = f'DestinE Pilot with {ALT} {sol[1]}'

    network_solved.export_to_netcdf(snakemake.output[0])
