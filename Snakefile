import glob

configfile: "main_config.yaml"
inputs = []

TARGETYEARS = config['scenario']['target_years']
CLIMATICDATASOURCES = config['climatic_data_sources']
PECDALTERNATIVES = config['PECD_alternatives']['include']
ALL_CLIMATICYEARS = range(config['scenario']['climatic_year_from'], config['scenario']['climatic_year_to'])
SIMULATIONYEARS = config['scenario']['climatic_years']

if SIMULATIONYEARS == 'all':
    try:
        SIMULATIONYEARS = range(config['scenario']['climatic_year_from'], config['scenario']['climatic_year_to'])
    except:
        print('When climatic_years = "all", climatic_year_from and climatic_year_to must be specified')

RESOURCE_PATH = "resources/"
MODEL_PATH = config['model'] + "/"
DOWNLOAD_PATH = RESOURCE_PATH + "raw/"
PEMMDB_PATH = RESOURCE_PATH + "PEMMDB/"
PECD_PATH = RESOURCE_PATH + 'climatic/PECD/'

wildcard_constraints:
    cds = CLIMATICDATASOURCES,
    target_year = "(\d{4})",
    simulation_year = "(\d{4})"

include: PEMMDB_PATH + "prepare.smk"
include: PECD_PATH + 'prepare.smk'
include: MODEL_PATH + "run.smk"
include: RESOURCE_PATH + 'climatic/retrieve_climate_data.smk'

### build pilot(s) with PEEMDB + PECD and collect additional data sources (if applicable)
if CLIMATICDATASOURCES:
    inputs.extend(expand(RESOURCE_PATH + "climatic/{CDS}/CY{SY}/generation_vre_timeseries.csv", CDS=CLIMATICDATASOURCES, SY=SIMULATIONYEARS))

if config['model'] == 'pypsa':
    if config['solve']['solve']:
        odir = 'pypsa/results/'
    else:
        odir = 'pypsa/networks/'
    if CLIMATICDATASOURCES:
        inputs.extend(expand(odir+'pilot_elec-vre-hydro_simpl_TY{TY}_{SY}_{cds}.nc', TY=TARGETYEARS, SY=SIMULATIONYEARS, cds=CLIMATICDATASOURCES))
    if PECDALTERNATIVES:
#        idirs = sorted(glob.glob(config['PECD_alternatives']['dir']+'PECD2021*'))
#        ALTERNATIVES = [idir.split('PECD2021_')[-1] for idir in idirs]
        inputs.extend(expand(odir+'pilot_elec-vre-hydro_simpl_TY{TY}_{SY}_{alt}.nc', TY=TARGETYEARS, SY=SIMULATIONYEARS, alt=PECDALTERNATIVES))
    if not CLIMATICDATASOURCES and not PECDALTERNATIVES:
        inputs.extend(expand(odir+'pilot_elec-vre-hydro_simpl_TY{TY}_{SY}.nc', TY=TARGETYEARS, SY=SIMULATIONYEARS))

if config['model'] == 'remix':
    inputs.extend(expand(MODEL_PATH + "results/TY{TY}/CY{CY}/commodity_balance_annual.csv", TY=TARGETYEARS, CY=SIMULATIONYEARS)),
    inputs.extend(expand(MODEL_PATH + "results/TY{TY}/CY{CY}/slack_commodity_balance.csv", TY=TARGETYEARS, CY=SIMULATIONYEARS)),

rule run_all:
    input: inputs

rule preprocess_input_data:
    input:
        inputs,
        expand(PECD_PATH + "TY{TY}/CY{CY}/hydro_inflow_timeseries.csv", TY=TARGETYEARS, CY=SIMULATIONYEARS),
        expand(PECD_PATH + "TY{TY}/CY{CY}/demand_timeseries.csv", TY=TARGETYEARS, CY=SIMULATIONYEARS),
        expand(PECD_PATH + "TY{TY}/CY{CY}/generation_vre_timeseries.csv", TY=TARGETYEARS, CY=SIMULATIONYEARS),

rule prepare_remix_model:
    input:
        inputs,
        expand(MODEL_PATH + "modeldata/TY{TY}/CY{CY}/final_model/converter_coefficient.csv", TY=TARGETYEARS, CY=SIMULATIONYEARS),

rule remix_hydro_constraints:
    input:
        inputs,
        expand(MODEL_PATH + "modeldata/TY{TY}/CY{CY}/hyrdo_constraints/storage_levelprofile.csv", TY=TARGETYEARS, CY=SIMULATIONYEARS),