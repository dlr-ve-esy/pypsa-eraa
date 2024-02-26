#!/usr/bin/env python
# coding: utf-8

import pypsa
import pandas as pd
import numpy as np

def solve_network(network, snapshots, solver_name='glpk', solver_options={}):

    ### p_min_pu of dispatch links may cause infeasibilities when inflow is insufficient to serve the minimum dispatch requirement and stores energy is bounded by e_min_pu
    disp_links = network.links.index[network.links.carrier.isin(['hydro dispatch', 'OLPHS dispatch', 'CLPHS dispatch'])]
    disp_links = disp_links[disp_links.isin(network.links_t.p_min_pu.columns)]
    network.links_t.p_min_pu[disp_links] = 0 ### solves!

    sol = network.optimize(snapshots, solver_name=solver_name, solver_options=solver_options, assign_all_duals=True)

    return((sol, network))

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

if __name__ == "__main__":

    network_f = '/home/brschy/repos/destine/pilot/pypsa/results/pilot_elec-vre-hydro_simpl_TY2030_2002.nc'

    network = pypsa.Network(network_f)

    speed_up(network)

    sol, network_solved = solve_network(
        network, 
        network.snapshots[:-49], 
        solver_name='gurobi',
        solver_options={
            'method': 2,
            'crossover': 0,
            'threads': 4,
            'BarConvTol': 1.e-6,
            'Seed': 123,
            'AggFill': 0,
            'PreDual': 0,
            'GURO_PAR_BARDENSETHRESH': 200
            }
    )
    network_solved.name = f'DestinE Pilot {sol[1]}'

    network_solved.export_to_netcdf(network_f)
