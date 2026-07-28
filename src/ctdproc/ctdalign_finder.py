def ctdalign_finder(yaml_file):


    import yaml
    import warnings

    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        module="seabirdscientific"
    )

    import datetime as datetime
    import numpy as np
    import pandas as pd
#    import gsw

# # Import modified Sea-Bird 
    import seabirdscientific.conversion as conv
    import seabirdscientific.instrument_data as id_sbs
    import seabirdscientific.processing as proc

    from .load import load_raw_data, convert_data
    from .parsers import CTDCoefficients, load_xmlcon
    from .processing import (
        filt_interp,
        crosshigh,
    )
    from .utilities import extract_head_hex, add_timestamp


    with open(yaml_file, "r") as file:
        config = yaml.safe_load(file)
    if type(config['constants']['sample_interval'])==str:
        numi,denomi=config['constants']['sample_interval'].split('/')
        sample_interval=float(numi)/float(denomi)
        
    offset_alignment_values = []
    step = config.get("alignment_settings",{})
    #Here a loop is needed
    #lower example to track ['pressure','temperature','conductivity'] --> conductivity will lag temperature
    for j in list(range(config['constants']['ind'][0],config['constants']['ind'][1]+1)):

        raw_data = load_raw_data("../data/"+config['constants']['name']+'_'+str(j).zfill(3)+'.hex')
        data = convert_data(raw_data,"../data/"+config['constants']['name']+'_'+str(j).zfill(3)+'.xmlcon',sample_interval)

        for var_lead, var_lag, max_lag, pre_up,pre_lo in zip(

                step.get("variable_lead", []),
                step.get("variable_lag", []),
                step.get("max_lag",[]),
                step.get("pressure_upper",[]),
                step.get("pressure_lower",[])
                ):
            print(var_lead,var_lag,max_lag,pre_up,pre_lo)

            func_name = step.get('name')
            _,lg_dn,_,lg_up=crosshigh(data,['pressure',var_lead,var_lag],
                                      max_lag,pre_up,pre_lo,sample_interval,False)

                #allocate lg_up,lg_dn

            offset_alignment_values.append({
             'cast_no': j,
             'variable_lead': var_lead,
             'variable_lag': var_lag,
             'lag_down_cast': lg_dn,
             'lag_up_cast': lg_up
             })


    df_lags = pd.DataFrame(offset_alignment_values)



    df_means = df_lags.groupby(['variable_lead', 'variable_lag'])[['lag_down_cast', 'lag_up_cast']].mean().reset_index()
    df_means['cast_no'] = 'AVERAGE'  # Etiqueta clara para identificar el promedio

    df_means = df_means[['cast_no', 'variable_lead', 'variable_lag', 'lag_down_cast', 'lag_up_cast']]


    with open('alignment_report.csv', 'w', encoding='utf-8') as f:
        f.write("# === DETAILED LOG PER CAST ===\n")
        df_lags.to_csv(f, index=False)
    
        f.write("\n# === CAMPAIGN SUMMARY AVERAGES ===\n")
        df_means.to_csv(f, index=False)

import sys
def main():
    if len(sys.argv) != 2:

        print("Usage: python -m ctdproc.ctdalign_finder config_align.yaml")
        sys.exit(1)
        
    yaml_file = sys.argv[1]
    ctdalign_finder(yaml_file)
if __name__ == "__main__":
    main()
