#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import pypsa
import yaml

### define loads
def define_loads(network, dem_f):
    
    loads_p_set = pd.read_csv(dem_f, header=0, index_col=0)
    loads_p_set.index = network.snapshots

    network.madd('Load', loads_p_set.columns, bus=[c.split(' ')[0] for c in loads_p_set.columns], carrier=[c.split(' ')[-1] for c in loads_p_set.columns], p_set=loads_p_set)

### define generators
def define_generators(network, config, capacities=None):

    carriers = config['generation']['generator_types']
    short_names = config['carriers_short_names']
    
    gen_specs = capacities.loc[capacities.index[capacities.technology.isin(carriers)],:]

    gen_specs.technology = gen_specs.technology.map(short_names).astype(str)
    gen_specs['name'] = gen_specs[['node', 'technology']].agg(' '.join, axis=1)
    gen_specs.set_index('name', inplace=True)

    network.madd('Carrier', gen_specs.technology.unique())
    network.madd('Generator', gen_specs.index, bus=gen_specs.node, carrier=gen_specs.technology, p_nom=gen_specs.nominal_capacity, marginal_cost=gen_specs.marginal_cost)
    
### define storage_units
def define_stores(network, config, capacities=None, sus_energy=None):

    generator_types = config['generation']['generator_types']
    carriers = config['storage']['storage_types']
    storage_efficiencies = config['storage']['storage_efficiencies']
    storage_soc_initial = config['storage']['soc_initial']
    
    sus_id = pd.Index([c if (capacities.technology[c].split(' (')[0] in carriers) & (capacities.technology[c] not in generator_types) else np.nan
                   for c in capacities.index]).dropna()
    sus_specs = capacities.loc[sus_id,:]
    sus_specs['technology_short'] = sus_specs.technology.map(config['carriers_short_names']).astype(str)
    sus_specs['name'] = sus_specs[['node', 'technology']].agg(' '.join, axis=1)
    sus_specs.set_index('name', inplace=True)
    
    sus_energy['name'] = sus_energy[['node', 'technology']].agg(' '.join, axis=1)
    sus_energy.set_index('name', inplace=True)
    
    sus_p_nom_off = sus_specs.nominal_capacity[[sus_specs.technology[c].split('(')[-1][:-1] in ['Offtake', 'Pumping'] for c in sus_specs.index]]
    sus_p_nom_off.index = [id.split(' (')[0] for id in sus_p_nom_off.index]
    sus_p_nom_disp = sus_specs.nominal_capacity[[sus_specs.technology[c].split('(')[-1][:-1] in ['Injection', 'Turbine'] for c in sus_specs.index]]
    sus_p_nom_disp.index = [id.split(' (')[0] for id in sus_p_nom_disp.index]
    
    sus_efficiencies = pd.DataFrame(storage_efficiencies).T
    
    network.madd('Carrier', sus_specs.technology_short.unique())
    
    for su in sus_specs.index:
        carrier = sus_specs.technology[su].split(' (')[0]
        carrier_short = sus_specs.technology_short[su]
        node = sus_specs.node[su]
        name = f'{node} {carrier}'
        name_short = f'{node} {carrier_short}'
        e_nom = sus_energy.nominal_capacity[(sus_energy.technology==carrier)&(sus_energy.node==node)]
    
        if (not e_nom.empty) & (name_short not in network.stores.index) & (node in network.buses.index):
    
            network.add('Bus', name_short, x=network.buses.x[node], y=network.buses.y[node])
    
            if carrier in storage_soc_initial.keys():
                e_initial_per_period = True
                e_initial = storage_soc_initial[carrier] / 100 * e_nom
                e_cyclic = False
            else:
                e_initial_per_period = False
                e_initial = 0.
                e_cyclic = True
            
            network.add(
                'Store',
                name_short,
                bus = name_short,
                carrier = carrier_short,
                e_nom = e_nom,
                e_initial = e_initial,
                e_initial_per_period = e_initial_per_period,
                e_cyclic = e_cyclic
            )
    
            ### dispatch_link
            if name not in sus_p_nom_disp.index:
                print(f'Warning store {name} has turbining power of 0.')
                p_nom_disp = 0.
            else:
                p_nom_disp = sus_p_nom_disp[name]
                
            network.add(
                'Link',
                f'{name_short} dispatch',
                bus0 = name_short,
                bus1 = node,
                efficiency = sus_efficiencies.dispatch[carrier],
                p_nom = p_nom_disp,
                carrier = f'{carrier_short} dispatch'
            )
    
            ### store link
            if name in sus_p_nom_off.index:
                network.add(
                    'Link',
                    f'{name_short} store',
                    bus0 = node,
                    bus1 = name_short,
                    efficiency = sus_efficiencies.store[carrier],
                    p_nom = sus_p_nom_off[name],
                    carrier = f'{carrier_short} store'
                )

def add_electricity(network, config):

    ### config from configfile
#    TY = snakemake.wildcards[0] #config['scenario']['target_year']    
#    simulation_year = int(snakemake.wildcards[1]) #config['scenario']['simulation_year']
#    scenario = f'National estimates {TY}'
#    basedir_static = snakemake.params.basedir_static
#    basedir_climatic = snakemake.params.basedir_climatic
    
    ### read input files
    capacities = pd.read_csv(snakemake.input.generation_technical, header=0, index_col=0)#pd.read_csv(basedir_static+f'generation_technical.csv', header=0, index_col=0)
    capacities["technology"] = capacities["technology"].str.strip()
    sus_energy = pd.read_csv(snakemake.input.storage_technical, header=0, index_col=0)#pd.read_csv(basedir_static+f'storage_technical.csv', header=0, index_col=0)

    ### add network components
    define_loads(network, dem_f=snakemake.input.demand_timeseries)
    define_generators(network, config, capacities)
    define_stores(network, config, **{
        'capacities': capacities,
        'sus_energy': sus_energy
    })

    network.name = 'DestinE Pilot with electricity'
    return(network)

if __name__ == "__main__":

    ### add electricity and export network
    network = pypsa.Network(snakemake.input.network)
    network_with_elec = add_electricity(network, snakemake.config)
    network_with_elec.export_to_netcdf(snakemake.output[0])



