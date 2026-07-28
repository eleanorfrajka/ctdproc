import logging

import gsw
import numpy as np
import seabirdscientific.processing as proc
from .processing import filt_interp

logger = logging.getLogger(__name__)



def apply_align(data, step, dt, bad_flag_value):

    text = []

    for var, offset in zip(
        step.get("variables", []),
        step.get("offset", [])
    ):

        data[var] = proc.align_ctd(
            x=data[var],
            offset=offset,
            sample_interval=dt,
            flag_value=bad_flag_value,
        )

        text.append(f"* alignctd {var} {offset} s\n")

    text.append("\n")

    return data, text


def apply_low_filter(data, step, dt):

    text = []

    for var, time_constant in zip(
        step.get("variables", []),
        step.get("time_constant", [])
    ):

        data[var] = filt_interp(
            data,
            var,
            time_constant,
            dt,
            False
        )

        text.append(f"* filter_low_pass {var} {time_constant} s\n")

    text.append("\n")

    return data, text


def apply_wild_edit(data, step, bad_flag_value):

    text = []

    for var, std1, std2, scans, distance in zip(
        step.get("variables", []),
        step.get("std_pass_1", []),
        step.get("std_pass_2", []),
        step.get("scans_per_block", []),
        step.get("distance_to_mean", []),
    ):

        data[var] = proc.wild_edit(
            data=data[var].values,
            flags=data["flag"].values,
            std_pass_1=std1,
            std_pass_2=std2,
            scans_per_block=scans,
            distance_to_mean=distance,
            exclude_bad_flags=True,
            flag_value=bad_flag_value,
        )

        text.append(f"* wildedit {var}\n")
        text.append(f"* wildedit pass1 nstd {std1}\n")
        text.append(f"* wildedit pass2 nstd {std2}\n")
        text.append(f"* wildedit midelta {distance}\n")
        text.append(f"* wildedit npoint {scans}\n")

    text.append("\n")

    return data, text

def apply_celltm(data, step, sample_interval):


    text = []

    for var, amplitude, time_constant in zip(
        step.get("variables", []),
        step.get("amplitude", []),
        step.get("time_constant", []),
    ):

        if var == "conductivity":
            temp_var = "temperature"
            sal_var = "salinity"
        else:
            temp_var = "temperature2"
            sal_var = "salinity2"

        data[var] = proc.cell_thermal_mass(
            temperature_C=data[temp_var],
            conductivity_Sm=data[var],
            amplitude=amplitude,
            time_constant=time_constant,
            sample_interval=sample_interval,
        )

        data[sal_var] = gsw.SP_from_C(
            10 * data[var].values,
            t=data[temp_var].values,
            p=data["pressure"].values
        )

        text.append(f"* celltm {var}\n")
        text.append(f"* celltm A {amplitude}\n")
        text.append(f"* celltm B {time_constant}\n")

    text.append("\n")


    return data, text


def apply_loop_edit(data, step, sample_interval, bad_flag_value):

    text = []

    def _scalar(val):
        return val[0] if isinstance(val, list) else val

    min_velocity = _scalar(step.get("min_velocity"))
    window_size = _scalar(step.get("window_size"))
    mean_speed_percent = _scalar(step.get("mean_speed_percent"))
    remove_surface_soak = _scalar(step.get("remove_surface_soak"))
    min_soak_depth = _scalar(step.get("min_soak_depth"))
    max_soak_depth = _scalar(step.get("max_soak_depth"))
    use_deck_pressure_offset = _scalar(step.get("use_deck_pressure_offset"))

    data["flag"] = proc.loop_edit_pressure(
        pressure=data["pressure"].values,
        latitude=float(np.nanmedian(data["latitude"].values)),
        flag=data["flag"].values,
        sample_interval=sample_interval,
        min_velocity_type=proc.MinVelocityType.FIXED,
        min_velocity=min_velocity,
        window_size=window_size,
        mean_speed_percent=mean_speed_percent,
        remove_surface_soak=remove_surface_soak,
        min_soak_depth=min_soak_depth,
        max_soak_depth=max_soak_depth,
        use_deck_pressure_offset=use_deck_pressure_offset,
        exclude_flags=False,
        flag_value=bad_flag_value,
    )

    text.append("* loopedit\n")
    text.append(f"* loopedit_minvel {min_velocity}\n")
    text.append(f"* loopedit: min_depth {min_soak_depth}\n")
    text.append(f"* loopedit: max_depth {max_soak_depth}\n")
    text.append("\n")

    return data, text


def apply_bin(data, step, bad_flag_value, base_name):

    bin_results = []
    bin_logs = []
    output_names = []

    bin_types = step.get("bin_type", [])
    bin_sizes = step.get("bin_size", [])
    profiles = step.get("type_profile", [])

    data["elapsed_time"] = (
        data["elapsed_time"] - data["elapsed_time"].iloc[0]
    )

    for bin_type, bin_size, profile in zip(
        bin_types,
        bin_sizes,
        profiles
    ):

        bin_text = []

        if profile == "BOTH":
            cast_type = proc.CastType.BOTH
            cast_name = ""

        elif profile == "UP":
            cast_type = proc.CastType.UPCAST
            cast_name = "up"

        elif profile == "DOWN":
            cast_type = proc.CastType.DOWNCAST
            cast_name = "down"

        else:
            raise ValueError(
                f"Unknown cast type: {profile}"
            )


        if bin_type == "elapsed_time":

            bin_data = proc.bin_average(
                dataset=data,
                bin_variable="elapsed_time",
                bin_size=bin_size,
                cast_type=cast_type,
                flag_value=bad_flag_value,
                interpolate=True,
            )

            bin_label = f"{bin_size:g}sec"

            bin_text.append(
                f"* binsize elapsed_time {bin_size}s\n"
            )


        elif bin_type == "scans":

            bin_data = proc.bin_average(
                dataset=data,
                bin_variable="nScan",
                bin_size=bin_size,
                cast_type=cast_type,
                flag_value=bad_flag_value,
                interpolate=True,
            )

            bin_label = f"{bin_size:g}scans"

            bin_text.append(
                f"* binsize nScan {bin_size}\n"
            )


        elif bin_type == "pressure":

            bin_data = proc.bin_average(
                dataset=data,
                bin_variable="pressure",
                bin_size=bin_size,
                cast_type=cast_type,
                flag_value=bad_flag_value,
                interpolate=True,
            )

            bin_label = f"{bin_size:g}dbar"

            bin_text.append(
                f"* binsize pressure {bin_size}dbar\n"
            )


        else:
            raise ValueError(
                f"Unknown bin type: {bin_type}"
            )


        # remove flag for output
        if "flag" in bin_data.columns:
            bin_data = bin_data.drop(columns=["flag"])


        # filename convention:
        # upMIXSED2_038_1dbar.cnv
        # downMIXSED2_038_1dbar.cnv
        # MIXSED2_038_24scans.cnv

        if cast_name:
            filename = f"{cast_name}{base_name}_{bin_label}.cnv"
        else:
            filename = f"{base_name}_{bin_label}.cnv"


        output_names.append(filename)
        bin_results.append(bin_data)
        bin_logs.append(bin_text)


    return bin_results, bin_logs, output_names
