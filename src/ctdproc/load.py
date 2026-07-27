def load_raw_data(hex_file):
    import seabirdscientific.instrument_data as id_sbs

    raw_data = id_sbs.read_hex_file(
         filepath=hex_file,
         instrument_type=id_sbs.InstrumentType.SBE911Plus,
         enabled_sensors=[
             id_sbs.Sensors.Temperature,
             id_sbs.Sensors.Conductivity,
             id_sbs.Sensors.Pressure,
             id_sbs.Sensors.SecondaryTemperature,
             id_sbs.Sensors.SecondaryConductivity,
             id_sbs.Sensors.ExtVolt0,
             id_sbs.Sensors.ExtVolt1,
             id_sbs.Sensors.ExtVolt2,
             id_sbs.Sensors.ExtVolt3,
             id_sbs.Sensors.ExtVolt4,
             id_sbs.Sensors.ExtVolt5,
             id_sbs.Sensors.ExtVolt6,
             id_sbs.Sensors.ExtVolt7,
             # id_sbs.Sensors.SystemTime,
             id_sbs.Sensors.nmeaLocation,
             # id_sbs.Sensors.nmeaTime,
             #It is really important to have the correct sensors that is manually edited but it can be automatized
          ],
          frequency_channels_suppressed=0,
          voltage_words_suppressed=0,
      )

    return raw_data




def convert_data(raw_data,xmlconfile,sample_interval):
    from .parsers import CTDCoefficients, load_xmlcon
    import seabirdscientific.conversion as conv
    import numpy as np
    import pandas as pd
    import gsw
    coeffs = load_xmlcon(xmlconfile)
    temperature = conv.convert_temperature_frequency(
        frequency=raw_data["temperature"],
        coefs=coeffs.temperature_primary,
        standard="ITS90",
        units="C",
         )
    pressure = conv.convert_pressure_digiquartz(
        pressure_count=raw_data["digiquartz pressure"],
        compensation_voltage=raw_data["temperature compensation"],
        coefs=coeffs.pressure,
        units="dbar",
        sample_interval=sample_interval,
        )


    temperature2 = conv.convert_temperature_frequency(
        frequency=raw_data["secondary temperature"],
        coefs=coeffs.temperature_secondary,
        standard="ITS90",
        units="C",
     )

    conductivity = conv.convert_conductivity(
        conductivity_count=raw_data["conductivity"],
        temperature=temperature,
        pressure=pressure,
        coefs=coeffs.conductivity_primary,
        scalar=1/10,
       )
    conductivity2 = conv.convert_conductivity(
        conductivity_count=raw_data["secondary conductivity"],
        temperature=temperature2,
        pressure=pressure,
        coefs=coeffs.conductivity_secondary,
        scalar=1/10,
       )

    chlorophyll = conv.convert_eco(
        raw=raw_data['volt 7'],
        coefs=coeffs.chlorophyll
    )

    turbidity = conv.convert_eco(
        raw=raw_data['volt 6'],
        coefs=coeffs.turbidity
    )

    salinity=gsw.SP_from_C(
        C=conductivity*10,
        t=temperature,
        p=pressure
    )

    oxygen = conv.convert_sbe43_oxygen(
        voltage=raw_data['volt 2'],
        temperature=temperature,
        pressure=pressure,
        salinity=salinity,
        coefs=coeffs.oxygen_primary,
        apply_tau_correction=True,
        apply_hysteresis_correction=True,
        window_size=1,
        sample_interval=sample_interval

    )
#data_forbin[['pre','tem','con','tem2','con2','lat','lon','elap_t','oxyg','chlo','turb','sal','sal2','flag']]
    flag = np.zeros(len(temperature))
    scans = np.arange(len(temperature))
    elapsed_time = (np.cumsum((np.ones(len(temperature))*sample_interval))-sample_interval).astype(float)
    # data['elapse_time']=elapsed_time
    data = pd.DataFrame(
         {
             "pressure": pressure,
             "temperature": temperature,
             "conductivity": conductivity,
             "temperature2": temperature2,
             "conductivity2": conductivity2,
             "latitude": raw_data["NMEA Latitude"],
             "longitude": raw_data["NMEA Longitude"],
             "elapsed_time": elapsed_time,
             # "scans": scans,
             "oxygen": oxygen,
             "chlorophyll": chlorophyll,
             "turbidity": turbidity,
             # "system_time": raw_data["system time"],
             "flag": flag,
               }
     )

    return data
