def CTD_process(yaml_file):
    import yaml
#    import os
#    import sys
    import warnings

    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        module="seabirdscientific"
    )
#    sys.path.append(os.path.abspath("../"))

    import datetime as datetime
    import numpy as np
    import pandas as pd
    import gsw

# # Import modified Sea-Bird 
    import seabirdscientific.conversion as conv
    import seabirdscientific.instrument_data as id_sbs
    import seabirdscientific.processing as proc


# Import local src
 #   from src.parsers import CTDCoefficients, load_xmlcon
 #   from src.processing import filt_interp, crosshigh, alp_tau, find_opt_alp_tat
 #   from src.utilities import extract_head_hex
    # Import local ctdproc package
    from .load import load_raw_data, convert_data
    from .parsers import CTDCoefficients, load_xmlcon
    from .processing import (
        filt_interp,
        crosshigh,
        alp_tau,
        find_opt_alp_tat,
    )
    from .utilities import extract_head_hex, add_timestamp
    from .pipeline import (
    apply_align,
    apply_low_filter,
    apply_wild_edit,
    apply_celltm,
    apply_loop_edit,
    apply_bin,
    apply_slope_correction,
    )
    from pathlib import Path

    with open(yaml_file, "r") as file:
        config = yaml.safe_load(file)
    if type(config['constants']['sample_interval'])==str:
        numi,denomi=config['constants']['sample_interval'].split('/')
        sample_interval=float(numi)/float(denomi)

    data_dir = Path(config['constants'].get('data_dir', '../data'))
    output_dir = Path(config['constants'].get('output_dir', '.'))
    output_dir.mkdir(parents=True, exist_ok=True)


    BAD_FLAG_VALUE = np.float64(-9.99e-29)

    decimales = {
        
        'pressure': 3,
        'temperature': 4,
        'conductivity': 6,
        'temperature2': 4,
        'conductivity2': 6,
        'latitude': 5,
        'longitude': 5,
        'salinity': 4,
        'salinity2': 4,
        'elapsed_time': 3,
        'oxygen': 3,
        'chlorophyll': 4,
        'turbidity': 4}
    
    for j in list(range(config['constants']['ind'][0],config['constants']['ind'][1]+1)):
        
    

        
        base = config['constants']['name'] + '_' + str(j).zfill(3)
        raw_data = load_raw_data(data_dir / (base + '.hex'))
        data = convert_data(raw_data, data_dir / (base + '.xmlcon'), sample_interval)
        head=[]
        head_coeff=[]
        head=extract_head_hex(data_dir / (base + '.hex'))
        head_coeff=extract_head_hex(data_dir / (base + '.XMLCON'))
        text_content =[]
        text_content.extend(head)
        text_content.append("* \n")
        text_content.extend(["* " + line for line in head_coeff[19:]]) 
        text_content.append("* \n")
        # text_content
        text_process=[]
    
        # for step in config.get('features', []):
        for step in config.get('features') or []:
            func_name = step.get('name')
            print(func_name)
        
            # print(vars_to_process)
            if func_name == 'align':
                
                data, log = apply_align(data,step,sample_interval,BAD_FLAG_VALUE,)
                text_process.extend(log)

            elif func_name == 'low_filter':
                data, log = apply_low_filter(data,step,sample_interval,)
                text_process.extend(log)

            elif func_name == 'wild_edit':
                
                data, log = apply_wild_edit(data,step,BAD_FLAG_VALUE,)
                text_process.extend(log)
      
            elif func_name == 'celltm':
                data, log = apply_celltm(data,step,sample_interval,)
                text_process.extend(log)

            elif func_name == 'slope_corr':
                data, log = apply_slope_correction(data,step,)
                text_process.extend(log)
    
            elif func_name == 'loop_edit':
                data, log = apply_loop_edit(data,step,sample_interval,BAD_FLAG_VALUE,)
                text_process.extend(log)

            elif func_name == 'bin':
                bin_results, bin_logs, filenames = apply_bin(data,step, BAD_FLAG_VALUE, base,)
                filenames = [str(output_dir / f) for f in filenames]
            
                # text_process.extend(log)
            

            
                # text_content.extend(text_process)
                for bin_data, bin_log, filename in zip(
                    
                    bin_results,
                    bin_logs,
                    filenames
                    ):
                    
                    file_text = []
                
                        # all processing before binning
                    file_text.extend(text_content)
                    file_text.extend(text_process)
                
                        # only this bin information
                    file_text.extend(bin_log)
                
                    with open(filename, "w", encoding="utf-8") as f:
                        f.writelines(file_text)
                
                    bin_data.round(decimales).to_csv(
                        filename,
                        mode="a",
                        index=False
                    )


        if config['constants']['convert'] == 'yes':
            nominal_name = str(output_dir / (base + "_24hz.cnv"))
        
            if config and config.get('constants', {}).get('time_ref') is not None:
                
                ref_time = config['constants']['time_ref']
        
                scans = np.arange(data.shape[0])
        
                data['timeJ'] = add_timestamp(ref_time,head,scans)
        
            data['elapsed_time'] = (data['elapsed_time'] - data['elapsed_time'].iloc[0])
        
            data = data.drop(columns=['flag'])
            file_text = []
            file_text.extend(text_content)
            file_text.extend(text_process)
            file_text.append("* Nominal resolution\n")
        
            with open(nominal_name, "w", encoding="utf-8") as f:
                f.writelines(file_text)
                
            data.round(decimales).to_csv(
                nominal_name,
                mode="a",
                index=False
            )           

        print('processed_'+str(j).zfill(3))

