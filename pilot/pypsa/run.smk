configfile: "pypsa/pypsa_config.yaml"

rule prepare_networks:
    input:
        link_specs_HVAC = PEMMDB_PATH + 'TY{target_year}/links_technical_HVAC.csv',
        link_specs_HVDC = PEMMDB_PATH + 'TY{target_year}/links_technical_HVDC.csv',
        link_p_max_pu_HVAC = PEMMDB_PATH + 'TY{target_year}/links_timeseries_HVAC.csv',
        link_p_max_pu_HVDC = PEMMDB_PATH + 'TY{target_year}/links_timeseries_HVDC.csv',
        bus_locations = RESOURCE_PATH + 'bus_locations.csv'
    output:
        'pypsa/networks/pilot_base_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_prepare_network.py'

rule add_electricity:
    input:
        network = 'pypsa/networks/pilot_base_TY{target_year}_{simulation_year}.nc',
        generation_technical = PEMMDB_PATH + 'TY{target_year}/generation_technical.csv',
        storage_technical = PEMMDB_PATH + 'TY{target_year}/storage_technical.csv',
        demand_timeseries = PECD_PATH + 'TY{target_year}/CY{simulation_year}/demand_timeseries.csv'
    output:
        'pypsa/networks/pilot_elec_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_add_electricity.py'

rule add_vre:
    input:
        network = 'pypsa/networks/pilot_elec_TY{target_year}_{simulation_year}.nc',
        generation_vre_timeseries = PECD_PATH + 'TY{target_year}/CY{simulation_year}/generation_vre_timeseries.csv',
        hydro_inflow_timeseries = PECD_PATH + 'TY{target_year}/CY{simulation_year}/hydro_inflow_timeseries.csv',
    output:
        'pypsa/networks/pilot_elec-vre_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_add_vre_profiles.py'

rule add_hydro:
    input:
        network = 'pypsa/networks/pilot_elec-vre_TY{target_year}_{simulation_year}.nc',
        hydro_inflow_timeseries = PECD_PATH + 'TY{target_year}/CY{simulation_year}/hydro_inflow_timeseries.csv',
    params:
        basedir = PECD_PATH + 'TY{target_year}/'
    output:
        'pypsa/networks/pilot_elec-vre-hydro_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_add_hydro.py'

rule simplify_network:
    input:
        'pypsa/networks/pilot_elec-vre-hydro_TY{target_year}_{simulation_year}.nc'
    output:
        'pypsa/networks/pilot_elec-vre-hydro_simpl_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_simplify_network.py'

rule solve_network:
    input:
        'pypsa/networks/pilot_elec-vre-hydro_simpl_TY{target_year}_{simulation_year}.nc'
    output:
        'pypsa/results/pilot_elec-vre-hydro_simpl_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_solve_network.py'