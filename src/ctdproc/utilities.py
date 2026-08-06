# src/utilities.py

import logging
logger = logging.getLogger(__name__)
import datetime as datetime
import pandas as pd


def extract_head_hex(
    hex_file: str
) -> list:
    """ Extract information from .hex and .XMLCON
    """
    for_heading = []

    with open(hex_file, "r", encoding="latin-1") as file:
        for line in file:

            if line.startswith("*END*"):
                break


            if not line.strip():
                continue


            for_heading.append(line)
    return for_heading


def add_timestamp(ref_time,header,scans):
    for linea in header:
        if ref_time in linea:
            string_date = linea.split('=')[1].strip()
            init_date = datetime.datetime.strptime(string_date, '%b %d %Y %H:%M:%S')
            break
    base_tojul=(init_date-datetime.datetime(init_date.year,1,1)).total_seconds()
    julian_stamp_time=((base_tojul + (scans / 24.0)) / 86400.0) + 1.0
    return julian_stamp_time


def bl_reader(bl_file):
    df_bl= pd.read_csv(
        bl_file,
        sep = ',',
        comment = '*',
        header = 1,
        skipinitialspace = True,
        names = ['sequence','position','datetime_utc','start_scan','end_scan'],
    )
    return df_bl