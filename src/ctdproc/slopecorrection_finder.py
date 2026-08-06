def slopecorrection_finder(yaml_file):


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
    import gsw
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
        slope_for_correction,
        bottle_avg,
        
        

    )
    from .utilities import extract_head_hex, add_timestamp, bl_reader
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
    path_dir_data = '../data/'
    a_b_stations_m_values = []
    detailed_points = [] 
    slope_step = config.get("slope_settings",{})
    #Here a loop is needed
    #lower example to track ['pressure','temperature','conductivity'] --> conductivity will lag temperature
    data_salt_name = slope_step.get("bottle_sal_file",[])
    data_salt = pd.read_csv(data_salt_name)
    data_salt['No_bottle'] = data_salt['No_bottle'].astype(int)
    for offset, duration, temp, cond in zip(
            

            slope_step.get("scan_offset",[]),
            slope_step.get("scan_duration",[]),
            slope_step.get("temperature_sensors",[]),
            slope_step.get("conductivity_sensors",[])


            ):
        
        print(data_salt,temp,cond,offset,duration)
        a=np.zeros(data_salt.shape[0]) #Create vector for a and b
            #b=np.zeros(data_salt.shape[0])

        processed_stations = {}
    
        unique_stations = data_salt['Station'].unique()


        for station in unique_stations:


            hex_file = path_dir_data+station+".hex"
            xmlconfile = path_dir_data+station+".XMLCON"
            bl_file = path_dir_data+station+'.bl'

            raw_data = load_raw_data(hex_file)
            data=convert_data(raw_data,xmlconfile,sample_interval)
            data['scan']=data.index
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
                
            bl_data = bl_reader(bl_file)
            from_bl_data_example=bottle_avg(data,bl_data,offset,duration,sample_interval)
            processed_stations[station] = from_bl_data_example.set_index('bottle')
                
        temp_vector = np.zeros(data_salt.shape[0])
        press_vector = np.zeros(data_salt.shape[0])
   


        for k in range(data_salt.shape[0]):

            station = data_salt['Station'].iloc[k]
            bottle_no = data_salt['No_bottle'].iloc[k]
            station_data = processed_stations[station]
            a[k] = station_data.loc[bottle_no, cond]
            temp_vector[k] = station_data.loc[bottle_no, temp]
            press_vector[k] = station_data.loc[bottle_no, 'pressure']
            

        b = gsw.C_from_SP(data_salt['Sal_avg'].values, temp_vector, press_vector) / 10 # S/m
        stations = data_salt['Station'].values
        bottles = data_salt['No_bottle'].values
        m = slope_for_correction(a,b)


        a_b_stations_m_values.append({
            'temperature_sensor':temp,
            'conductivity_sensor':cond,
        # 'station_no':stations,
        # 'CTD_conductivity': a,
        # 'Bottle_conductivity': b,
            'slope_coefficient': m
            })
            
        for st, bt, c_ctd, c_btl in zip(stations, bottles, a, b):
            detailed_points.append({
                'temperature_sensor': temp,
                'conductivity_sensor': cond,
                'Station': st,
                'No_bottle': bt,
                'CTD_conductivity_a': c_ctd,
                'Bottle_conductivity_b': c_btl
                })



        fig, ax = plt.subplots(1, 1, figsize=(15, 5))
        
        clean_labels = pd.Series(stations).str.replace(config['constants']['name']+'_', '')
        
    
        ax.plot(clean_labels, a / b, '.', markersize=8, color='blue', label='Bottle Ratios (a/b)')
        
    
        ax.axhline(y=1.0, color='black', linestyle='-', linewidth=1)
        ax.axhline(y=m, color='green', linestyle='-.', linewidth=1.5, label=f'Calculated Slope m ({m:.6f})')

        
    
        ax.tick_params(axis='x', rotation=90, labelsize=8)
        ax.set_xlabel('Stations')
        ax.set_ylabel('Ratio CTD / Salinometer')
        ax.set_title(f'Conductivity Ratio Validation & Slope Report\nSensor Configuration: Temp={temp} | Cond={cond}')
      #  ax.set_ylim(0.994, 1.006)
        ax.legend(loc='upper right')
        ax.grid(True, linestyle=':', alpha=0.5)
        
        plt.tight_layout()
        plt.savefig(f"Conductivity_ratio_report_{cond}.png", dpi=150)
        plt.close()


    

    df_summary = pd.DataFrame(a_b_stations_m_values)
    df_detailed = pd.DataFrame(detailed_points)


    with open('slope_report.csv', 'w', encoding='utf-8') as f:
        f.write("# === CAMPAIGN SUMMARY AVERAGES (COMPUTED PARAMETERS) ===\n")
        df_summary.to_csv(f, index=False)

        f.write("\n# === DETAILED LOG PER STATION / BOTTLE / SENSOR ===\n")
        df_detailed.to_csv(f, index=False)





import sys
def main():
    if len(sys.argv) != 2:

        print("Usage: python -m ctdproc.slopecorrection_finder config_saltcorr.yaml")
        sys.exit(1)
        
    yaml_file = sys.argv[1]
    slopecorrection_finder(yaml_file)
if __name__ == "__main__":
    main()
