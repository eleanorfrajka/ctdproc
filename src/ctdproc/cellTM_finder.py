def cellTM_finder(yaml_file):


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
    import matplotlib.pyplot as plt

# # Import modified Sea-Bird 
    import seabirdscientific.conversion as conv
    import seabirdscientific.instrument_data as id_sbs
    import seabirdscientific.processing as proc

    from .load import load_raw_data, convert_data
    from .parsers import CTDCoefficients, load_xmlcon
    from .processing import (
        filt_interp,
        crosshigh,
        find_opt_alp_tat_fast,

    )
    from .utilities import extract_head_hex, add_timestamp
    from .pipeline import (
    apply_align,
    apply_low_filter,
    apply_wild_edit,
    apply_celltm,
   # apply_loop_edit,
   # apply_bin,
    )

    with open(yaml_file, "r") as file:
        config = yaml.safe_load(file)
    if type(config['constants']['sample_interval'])==str:
        numi,denomi=config['constants']['sample_interval'].split('/')
        sample_interval=float(numi)/float(denomi)
    BAD_FLAG_VALUE = np.float64(-9.99e-29)    
    tau_alpha_values = []
    calib_step = config.get("celltm_settings",{})
    #Here a loop is needed
    #lower example to track ['pressure','temperature','conductivity'] --> conductivity will lag temperature
    for j in list(range(config['constants']['ind'][0],config['constants']['ind'][1]+1)):

        raw_data = load_raw_data("../data/"+config['constants']['name']+'_'+str(j).zfill(3)+'.hex')
        data = convert_data(raw_data,"../data/"+config['constants']['name']+'_'+str(j).zfill(3)+'.xmlcon',sample_interval)

        for feat in config.get('features') or []:
            func_name = feat.get('name')

            # print(vars_to_process)
            if func_name == 'align':

                data, log = apply_align(data,feat,sample_interval,BAD_FLAG_VALUE,)
            #    text_process.extend(log)

            elif func_name == 'low_filter':
                data, log = apply_low_filter(data,feat,sample_interval,)
            #   text_process.extend(log)
            
            elif func_name == 'wild_edit':

                data, log = apply_wild_edit(data,feat,BAD_FLAG_VALUE,)
            #    text_process.extend(log)


        for pre_up, pre_lo, bin_T,T_to_use,C_to_use in zip(

                calib_step.get("pressure_upper", []),
                calib_step.get("pressure_lower", []),
                calib_step.get("bin_temperature",[]),
                calib_step.get("temperature_touse",[]),
                calib_step.get("conductivity_touse",[])
                ):
            print(j,pre_up,pre_lo,bin_T)

            alpha_range = np.arange(0.01, 0.09, 0.001) #FIXED
            tau_range = np.arange(1, 16, 0.5) #FIXED
            MAT = np.zeros((len(alpha_range), len(tau_range)))

            #CHECK HERE TO MODIFIY THE NAME OF THIS FUNCTION
            MAT=find_opt_alp_tat_fast(alpha_range,tau_range,data,
                                      ['pressure',T_to_use,
                                       C_to_use,'corrected_sal'],pre_up,pre_lo,bin_T)

            row_idx, col_idx = np.unravel_index(np.argmin(MAT), MAT.shape)
            best_alpha = alpha_range[row_idx]
            best_tau = tau_range[col_idx]
            min_val = MAT[row_idx, col_idx]
            #allocate tau and alpha

            tau_alpha_values.append({

                'cast_no': j,
                'temperature': T_to_use,
                'conductivity': C_to_use,
                'tau': best_tau,
                'alpha': best_alpha,
                'abs(S-S_corrected)':min_val
                })

            plt.figure(figsize=(7, 5))
            plt.contourf(tau_range, alpha_range, MAT, cmap='gray_r', levels=50)
            plt.colorbar(label='Average Salinity Error')
            plt.plot(best_tau, best_alpha, 'ro', label=f'Optimum (tau={best_tau:.1f}, alpha={best_alpha:.4f})')
            plt.ylabel('Alpha')
            plt.xlabel('Tau')
            plt.title(f'Cell Thermal Mass Optimization - Cast {str(j).zfill(3)} ({C_to_use})')
            plt.legend()
            
            plt.savefig(f"celltm_plot_cast_{str(j).zfill(3)}_{C_to_use}.png", dpi=150)
            plt.close()
    df_tau_alpha = pd.DataFrame(tau_alpha_values)

    

    df_means = df_tau_alpha.groupby(['temperature', 'conductivity'])[['tau', 'alpha', 'abs(S-S_corrected)']].mean().reset_index()
    df_means['cast_no'] = 'AVERAGE'  
    df_means = df_means[['cast_no', 'temperature', 'conductivity', 'tau', 'alpha', 'abs(S-S_corrected)']]
    #df_means = df_means[['cast_no', 'temperature', 'conductivity', 'tau', 'alpha', 'abs(S-S_corrected)']]

    with open('alpha_tau_report.csv', 'w', encoding='utf-8') as f:
        f.write("# === DETAILED LOG PER CAST ===\n")
        df_tau_alpha.to_csv(f, index=False)

        f.write("\n# === CAMPAIGN SUMMARY AVERAGES ===\n")
        df_means.to_csv(f, index=False)



import sys
def main():
    if len(sys.argv) != 2:

        print("Usage: python -m ctdproc.cellTM_finder config_align.yaml")
        sys.exit(1)
        
    yaml_file = sys.argv[1]
    cellTM_finder(yaml_file)
if __name__ == "__main__":
    main()
