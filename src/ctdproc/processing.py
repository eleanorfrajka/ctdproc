# src/processing.py

import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt
import gsw

# Import modified seabirdscientific
from seabirdscientific import processing as proc

logger = logging.getLogger(__name__)

# BAD_FLAG
BAD_FLAG_VALUE = np.float64(-9.99e-29)


def filt_interp(
    data: pd.DataFrame, param: str, tm: float, si: float, flag: bool
) -> np.ndarray:
    """Filter with low pass and linear interpolate if bad flags
    """
    data_i = data[param].copy()
    indx = data_i == BAD_FLAG_VALUE
    data_i[indx] = np.nan

    data_i = data_i.interpolate(limit_direction="both")

    data_interp = proc.low_pass_filter(
        x=data_i.values, time_constant=tm, sample_interval=si
    )

    if flag:
        data_interp[indx] = BAD_FLAG_VALUE

    return data_interp


def crosshigh(
    data: pd.DataFrame,
    param: list,
    maxL: int,
    pi: float,
    pf: float,
    si: float,
    high: bool,
) -> list:
    """Compute correlation of high frequency values for up- down profiles
    """
    data_i = data.copy()
    indx1 = data_i[param[1]] == BAD_FLAG_VALUE
    indx2 = data_i[param[2]] == BAD_FLAG_VALUE

    data_i.loc[indx1, param[1]] = np.nan
    data_i.loc[indx2, param[2]] = np.nan
    data_i[param[1]] = data_i[param[1]].interpolate(limit_direction="both")
    data_i[param[2]] = data_i[param[2]].interpolate(limit_direction="both")

    peak_pressure_idx = np.argmax(data_i[param[0]])

    fpar1_dn = data_i[param[1]][0:peak_pressure_idx]
    fpar1_up = data_i[param[1]][peak_pressure_idx + 1 :]
    fpar2_dn = data_i[param[2]][0:peak_pressure_idx]
    fpar2_up = data_i[param[2]][peak_pressure_idx + 1 :]
    fpar0_dn = data_i[param[0]][0:peak_pressure_idx]
    fpar0_up = data_i[param[0]][peak_pressure_idx + 1 :]

    if high:
        b, a = butter(3, 0.005)
        fpar1_dn_smoothed = filtfilt(b, a, fpar1_dn)
        fpar1_up_smoothed = filtfilt(b, a, fpar1_up)
        fpar2_dn_smoothed = filtfilt(b, a, fpar2_dn)
        fpar2_up_smoothed = filtfilt(b, a, fpar2_up)

        fpar1_dn = fpar1_dn - fpar1_dn_smoothed
        fpar1_up = fpar1_up - fpar1_up_smoothed
        fpar2_dn = fpar2_dn - fpar2_dn_smoothed
        fpar2_up = fpar2_up - fpar2_up_smoothed

    mask_dn = (fpar0_dn >= pi) & (fpar0_dn <= pf)
    mask_up = (fpar0_up >= pi) & (fpar0_up <= pf)

    lag_dn, corre_dn, _, _ = plt.xcorr(
        fpar1_dn[mask_dn], fpar2_dn[mask_dn], maxlags=maxL
    )
    lag_up, corre_up, _, _ = plt.xcorr(
        fpar1_up[mask_up], fpar2_up[mask_up], maxlags=maxL
    )

    best_corr = [
        corre_up[np.argmax(np.abs(corre_up))],
        si * lag_up[np.argmax(np.abs(corre_up))],
        corre_dn[np.argmax(np.abs(corre_dn))],
        si * lag_dn[np.argmax(np.abs(corre_dn))],
    ]
    return best_corr


def alp_tau(
        data: pd.DataFrame, param: list, pi: float, pf: float, tbin: float, figure: bool = False
) -> float:
    """ For a given alpha and Ta compute the average absolute error of salinity between up- down- profile
    """
    data_i = data.copy()
    indx1 = data_i[param[1]] == BAD_FLAG_VALUE
    indx2 = data_i[param[2]] == BAD_FLAG_VALUE

    data_i.loc[indx1, param[1]] = np.nan
    data_i.loc[indx2, param[2]] = np.nan
    data_i[param[1]] = data_i[param[1]].interpolate(limit_direction="both")
    data_i[param[2]] = data_i[param[2]].interpolate(limit_direction="both")
    data_i[param[2]] = gsw.SP_from_C(10*data_i[param[2]].values, t=data_i[param[1]].values, p=data_i[param[0]].values)
    
    peak_pressure_idx = np.argmax(data_i[param[0]])

    fpar1_dn = data_i[param[1]][0:peak_pressure_idx]
    fpar1_up = data_i[param[1]][peak_pressure_idx + 1 :]
    fpar2_dn = data_i[param[2]][0:peak_pressure_idx]
    fpar2_up = data_i[param[2]][peak_pressure_idx + 1 :]
    fpar0_dn = data_i[param[0]][0:peak_pressure_idx]
    fpar0_up = data_i[param[0]][peak_pressure_idx + 1 :]

    mask_dn = (fpar0_dn >= pi) & (fpar0_dn <= pf)
    mask_up = (fpar0_up >= pi) & (fpar0_up <= pf)

    min_max = np.min([fpar1_dn[mask_dn].max(), fpar1_up[mask_up].max()])
    max_min = np.max([fpar1_dn[mask_dn].min(), fpar1_up[mask_up].min()])
    temp_group = np.arange(max_min, min_max, tbin)

    pro_dn = pd.DataFrame({param[1]: fpar1_dn[mask_dn], param[2]: fpar2_dn[mask_dn]})
    pro_up = pd.DataFrame({param[1]: fpar1_up[mask_up], param[2]: fpar2_up[mask_up]})

    pro_dn[param[1] + "group"] = pd.cut(pro_dn[param[1]], bins=temp_group)
    pro_up[param[1] + "group"] = pd.cut(pro_up[param[1]], bins=temp_group)

    avg_fpar2_dn = pro_dn.groupby(
        param[1] + "group", observed=False
    )[param[2]].mean().reset_index()
    avg_fpar2_up = pro_up.groupby(
        param[1] + "group", observed=False
    )[param[2]].mean().reset_index()

    avg_fpar2_dn["temp_mid"] = (
        avg_fpar2_dn[param[1] + "group"]
        .apply(lambda x: x.mid if pd.notnull(x) else np.nan)
        .astype(float)
    )
    avg_fpar2_up["temp_mid"] = (
        avg_fpar2_up[param[1] + "group"]
        .apply(lambda x: x.mid if pd.notnull(x) else np.nan)
        .astype(float)
    )

    if figure:
        ax = avg_fpar2_dn.plot(
            x=param[2], y="temp_mid", color="red", label="dn", linestyle="-."
        )
        avg_fpar2_up.plot(
            ax=ax, y="temp_mid", x=param[2], color="blue", label="up", linestyle="-."
        )

    return float(np.mean(np.abs(avg_fpar2_dn[param[2]] - avg_fpar2_up[param[2]])))


def find_opt_alp_tat(
    alpha_r: np.ndarray,
    tau_r: np.ndarray,
    data: pd.DataFrame,
    param: list,
    pi: float,
    pf: float,
    tbin: float,
    figure: bool,
) -> np.ndarray:
    """Try to find the best alpha and Tau values.
    """
    data_i = data.copy()
    max_value = np.zeros((len(alpha_r), len(tau_r)))

    for r, i in enumerate(alpha_r):
        for c, j in enumerate(tau_r):
            data_i[param[3]] = proc.cell_thermal_mass(
                temperature_C=data_i[param[1]],
                conductivity_Sm=data_i[param[2]],
                amplitude=i,  # alpha
                time_constant=j,  # 1/beta (tau)
                sample_interval=1 / 24,  # Modified 
            )
            max_value[r, c] = alp_tau(
                data_i, [param[0], param[1], param[3]], pi, pf, tbin, figure=False
            )
    return max_value

