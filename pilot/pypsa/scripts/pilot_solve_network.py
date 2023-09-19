#!/usr/bin/env python
# coding: utf-8

import pypsa

def add_shedding(network):
    loads_i = network.loads.index
    buses = network.loads.bus.values
    names = [c+' shedding' for c in loads_i]
    p_nom = network.loads_t.p_set.max().values
    p_max_pu = (network.loads_t.p_set / p_nom).values
    
    network.add('Carrier', 'Load')
    network.madd('Generator', names, p_nom=p_nom, p_max_pu=p_max_pu, marginal_cost=1e3, bus=buses, carrier='Load')

def solve_network(network, snapshots, config):

    ### p_min_pu of dispatch links may cause infeasibilities when inflow is insufficient to serve the minimum dispatch requirement and stores energy is bounded by e_min_pu
    disp_links = network.links.index[network.links.carrier.isin(['hydro dispatch', 'OLPHS dispatch', 'CLPHS dispatch'])]
    disp_links = disp_links[disp_links.isin(network.links_t.p_min_pu.columns)]
    network.links_t.p_min_pu[disp_links] = 0 ### solves!

    sol = network.optimize(snapshots, **config['solve']['lopf_kwargs'])#, solver_options=config['solve']['solver_options'])

    return((sol, network))

if __name__ == "__main__":

    network = pypsa.Network(snakemake.input[0])
    add_shedding(network)

    network.remove('Link', 'LUB1-BE00 HVAC')

    sol, network_solved = solve_network(network, network.snapshots, snakemake.config)
    network_solved.name = f'DestinE Pilot {sol[1]}'

    network_solved.export_to_netcdf(snakemake.output[0])