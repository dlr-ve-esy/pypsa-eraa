configfile: "pypsa_config.yaml"

TARGETYEARS = config['scenario']['target_years']
SIMULATIONYEARS = config['scenario']['simulation_years']

rule evaluate_adequacy:
    input: expand('results/pilot_elec-vre-hydro_simpl_TY{TY}_{SY}.nc', TY=TARGETYEARS, SY=SIMULATIONYEARS)

rule prepare_networks:
    output:
        'networks/pilot_base_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_prepare_network.py'

rule add_electricity:
    input:
        'networks/pilot_base_TY{target_year}_{simulation_year}.nc'
    output:
        'networks/pilot_elec_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_add_electricity.py'

rule add_vre:
    input:
        'networks/pilot_elec_TY{target_year}_{simulation_year}.nc'
    output:
        'networks/pilot_elec-vre_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_add_vre_profiles.py'

rule add_hydro:
    input:
        'networks/pilot_elec-vre_TY{target_year}_{simulation_year}.nc'
    output:
        'networks/pilot_elec-vre-hydro_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_add_hydro.py'

rule simplify_network:
    input:
        'networks/pilot_elec-vre-hydro_TY{target_year}_{simulation_year}.nc'
    output:
        'networks/pilot_elec-vre-hydro_simpl_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_simplify_network.py'

rule solve_network:
    input:
        'networks/pilot_elec-vre-hydro_simpl_TY{target_year}_{simulation_year}.nc'
    output:
        'results/pilot_elec-vre-hydro_simpl_TY{target_year}_{simulation_year}.nc'
    script:
        'scripts/pilot_solve_network.py'