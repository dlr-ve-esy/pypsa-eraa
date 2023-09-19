#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pypsa
import yaml


### read net transfer capacities, specify links
def initialise_with_links(config):

    network = pypsa.Network(name='DestinE Pilot base')
    TY = snakemake.wildcards[0]
    simulation_year = snakemake.wildcards[1]
    
    basedir = config['io']['databasedir']+f'scenario_{TY}_{simulation_year}/'
    
    bus_locations = pd.read_csv(basedir+'bus_locations.csv', header=0, index_col=0)
    network.madd('Bus', bus_locations.index, x=bus_locations.x, y=bus_locations.y)
    
    for link_type in config['transmission']['link_types']:
    
        network.add('Carrier', link_type)
    
        link_specs = pd.read_csv(basedir+f'links_technical_{link_type}_{TY}.csv', header=0, index_col=0)
    
        link_p_max_pu = pd.read_csv(basedir+f'links_cf_timeseries_{link_type}_{TY}.csv', header=0, index_col=0).reindex(columns=link_specs.index, fill_value=1.)
        link_p_min_pu = -link_p_max_pu[link_specs.index[link_specs.bidirectional]]
        
        snapshots = pd.DatetimeIndex(pd.date_range(start=f'{simulation_year}-01-01 00:00:00', freq='H', periods=link_p_max_pu.shape[0], name='snapshots'))
        link_p_max_pu.index = snapshots
        link_p_min_pu.index = snapshots
    
        link_specs.index = [f'{ll} {link_type}' for ll in link_specs.index]
        link_p_max_pu.columns = [f'{ll} {link_type}' for ll in link_p_max_pu.columns]
        link_p_min_pu.columns = [f'{ll} {link_type}' for ll in link_p_min_pu.columns]
    
        if link_type == config['transmission']['link_types'][0]:
            network.set_snapshots(snapshots)
    
        network.madd('Link', 
             link_specs.index, 
             bus0=link_specs['start'], 
             bus1=link_specs['end'], 
             p_min_pu=link_p_min_pu, 
             p_max_pu=link_p_max_pu, 
             p_nom=link_specs['nominal_capacity'], 
             carrier=link_type
        )

    return(network)

def drop_constants(network):
    ### identify links with constant p_min_pu/p_max_pu, drop from time series dataframe
    constant_p_min_pu = network.links.index[network.links_t.p_min_pu.min() == network.links_t.p_min_pu.max()]
    constant_p_max_pu = network.links.index[network.links_t.p_max_pu.min() == network.links_t.p_max_pu.max()]
    
    network.links.loc[constant_p_min_pu, 'p_min_pu'] = network.links_t.p_min_pu.min()[constant_p_min_pu]
    network.links_t.p_min_pu.drop(constant_p_min_pu, axis=1, inplace=True)
    network.links.loc[constant_p_max_pu, 'p_max_pu'] = network.links_t.p_max_pu.max()[constant_p_max_pu]
    network.links_t.p_max_pu.drop(constant_p_max_pu, axis=1, inplace=True)
    
    ### identify duplicated links (those from 'hidden columns' in excel sheet) and drop
    duplicated_links = network.links.index[[len(l.split('.'))==2 for l in network.links.index]]
    network.mremove('Link', duplicated_links)


if __name__ == "__main__":

    network = initialise_with_links(snakemake.config)
#    drop_constants(network)

    ### export network
    network.export_to_netcdf(snakemake.output[0])

