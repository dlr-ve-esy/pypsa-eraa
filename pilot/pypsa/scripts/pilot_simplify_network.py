#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import matplotlib as mpl
import pypsa
import glob
import yaml

from pypsa.clustering.spatial import get_clustering_from_busmap
from pypsa.descriptors import get_switchable_as_dense as get_as_dense

def aggregate_links(network, remove_zero_capacity_links=False):
    links_list = network.links.index
    p_max_pu = get_as_dense(network, 'Link', 'p_max_pu', network.snapshots)
    p_min_pu = get_as_dense(network, 'Link', 'p_min_pu', network.snapshots)
    new_p_min_pu = pd.DataFrame()
    new_p_max_pu = pd.DataFrame()
    new_p_nom = pd.Series()
    droplinks = []
    for link in links_list:

        if link in droplinks: continue
        
        bus0, bus1 = network.links.loc[link, ['bus0', 'bus1']]
        forward_p_max_pu = p_max_pu[link]
        bidirectional = p_min_pu[link].min() < 0.
        if not bidirectional:
            try:
                reverse_direction_link = network.links.index[(network.links.bus0==bus1) & (network.links.bus1==bus0)][0]
            except:
                print(f'Warning: Link {link} is not bidirectional but seems not to have a reverse direction link neither...')
                continue
                
            reverse_p_max_pu = p_max_pu[reverse_direction_link]
            p_nom = np.max([network.links.p_nom[link], network.links.p_nom[reverse_direction_link]])

            new_p_max_pu[link] = forward_p_max_pu.mul(network.links.p_nom[link]) / p_nom
            new_p_min_pu[link] = -reverse_p_max_pu.mul(network.links.p_nom[reverse_direction_link]) / p_nom
            new_p_nom[link] = p_nom
            droplinks.append(reverse_direction_link)

            if (remove_zero_capacity_links) & (p_nom == 0.):
                print('yes')
                droplinks.append(link)
                new_p_min_pu.drop(link, inplace=True)
                new_p_max_pu.drop(link, inplace=True)

    return((droplinks, new_p_nom, new_p_max_pu, new_p_min_pu))

def simplify_network(network):
    ### aggregate links
#    droplinks, new_p_nom, new_p_max_pu, new_p_min_pu = aggregate_links(network, remove_zero_capacity_links=True)

#    network.mremove('Link', droplinks)
#    network.links_t.p_max_pu.loc[:,new_p_max_pu.columns] = new_p_max_pu
#    network.links_t.p_min_pu.loc[:,new_p_min_pu.columns] = new_p_min_pu

    #### cluster network
    stores_strategies = {
        'e_max_pu': "capacity_weighted_average",
        'e_min_pu': "capacity_weighted_average",
        'e_initial': "sum",
        'e_nom': 'sum'    
    }

    network.loads.carrier = 'total_demand'
    
    busmap = {bus: bus[:2] for bus in network.buses.index}
    busmap['DEKF'] = 'KF' # Kriegers Flak
    busmap['DKKF'] = 'KF'
    busmap['LUV1'] = 'DE' # pumped storage Vianden (located in Luxenbourg but connected to the German grid)

    stores = network.stores
    for store in stores.index:
        busmap[stores.bus[store]] = f'{stores.bus[store][:2]} {stores.carrier[store]}'

    Clustering = get_clustering_from_busmap(network, busmap, with_time=True, aggregate_one_ports={'Generator', 'Load', 'Store'}, one_port_strategies={'Store': stores_strategies})

    n_clustered = Clustering.network
    
    n_clustered.links_t.p_max_pu.drop(n_clustered.links_t.p_max_pu.columns[~n_clustered.links_t.p_max_pu.columns.isin(n_clustered.links.index)], axis=1, inplace=True)
    n_clustered.links_t.p_min_pu.drop(n_clustered.links_t.p_min_pu.columns[~n_clustered.links_t.p_min_pu.columns.isin(n_clustered.links.index)], axis=1, inplace=True)

    return(n_clustered)

if __name__ == "__main__":

    network = pypsa.Network(snakemake.input[0])

    n_clustered = simplify_network(network)
    n_clustered.name = 'Simple DestinE Pilot with VRE profiles and hydro'
    n_clustered.export_to_netcdf(snakemake.output[0])
