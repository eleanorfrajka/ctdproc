# src/parsers.py

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Optional

# 'modified seabirdscientific'
from seabirdscientific.cal_coefficients import (
    AltimeterCoefficients,
    ConductivityCoefficients,
    ECOCoefficients,
    Oxygen43Coefficients,
    PressureDigiquartzCoefficients,
    TemperatureFrequencyCoefficients,
)


@dataclass
class CTDCoefficients:
    temperature_primary: Optional[TemperatureFrequencyCoefficients] = None
    temperature_secondary: Optional[TemperatureFrequencyCoefficients] = None
    conductivity_primary: Optional[ConductivityCoefficients] = None
    conductivity_secondary: Optional[ConductivityCoefficients] = None
    pressure: Optional[PressureDigiquartzCoefficients] = None
    oxygen_primary: Optional[Oxygen43Coefficients] = None
    oxygen_secondary: Optional[Oxygen43Coefficients] = None
    chlorophyll: Optional[ECOCoefficients] = None
    turbidity: Optional[ECOCoefficients] = None  # It is missing
    altimeter: Optional[AltimeterCoefficients] = None


def load_xmlcon(xml_file: str) -> CTDCoefficients:
    """Read a .xmlcon to extract calibration coefficients.
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()

    coeffs = CTDCoefficients()

    temp_count = 0
    cond_count = 0
    oxy_count = 0

    for sensor in root.findall(".//Sensor"):

        # --------------------------------------------------
        # TEMPERATURE
        # --------------------------------------------------
        temp = sensor.find("TemperatureSensor")
        if temp is not None:
            obj = TemperatureFrequencyCoefficients(
                g=float(temp.findtext("G")),
                h=float(temp.findtext("H")),
                i=float(temp.findtext("I")),
                j=float(temp.findtext("J")),
                f0=float(temp.findtext("F0")),
            )
            if temp_count == 0:
                coeffs.temperature_primary = obj
            else:
                coeffs.temperature_secondary = obj
            temp_count += 1

        # --------------------------------------------------
        # CONDUCTIVITY
        # --------------------------------------------------
        cond = sensor.find("ConductivitySensor")
        if cond is not None:
            eq1 = cond.find('./Coefficients[@equation="1"]')
            if eq1 is not None:  # CHECK
                obj = ConductivityCoefficients(
                    g=float(eq1.findtext("G")),
                    h=float(eq1.findtext("H")),
                    i=float(eq1.findtext("I")),
                    j=float(eq1.findtext("J")),
                    cpcor=float(eq1.findtext("CPcor")),
                    ctcor=float(eq1.findtext("CTcor")),
                    wbotc=float(eq1.findtext("WBOTC")),
                )
                if cond_count == 0:
                    coeffs.conductivity_primary = obj
                else:
                    coeffs.conductivity_secondary = obj
                cond_count += 1

        # --------------------------------------------------
        # PRESSURE
        # --------------------------------------------------
        press = sensor.find("PressureSensor")
        if press is not None:
            coeffs.pressure = PressureDigiquartzCoefficients(
                c1=float(press.findtext("C1")),
                c2=float(press.findtext("C2")),
                c3=float(press.findtext("C3")),
                d1=float(press.findtext("D1")),
                d2=float(press.findtext("D2")),
                t1=float(press.findtext("T1")),
                t2=float(press.findtext("T2")),
                t3=float(press.findtext("T3")),
                t4=float(press.findtext("T4")),
                t5=float(press.findtext("T5")),
                AD590M=float(press.findtext("AD590M")),
                AD590B=float(press.findtext("AD590B")),
                offset=float(press.findtext("Offset")),
                slope=float(press.findtext("Slope")),
            )

        # --------------------------------------------------
        # OXYGEN
        # --------------------------------------------------
        oxy = sensor.find("OxygenSensor")
        if oxy is not None:
            eq1 = oxy.find('./CalibrationCoefficients[@equation="1"]')
            if eq1 is not None:
                obj = Oxygen43Coefficients(
                    soc=float(eq1.findtext("Soc")),
                    v_offset=float(eq1.findtext("offset")),
                    tau_20=float(eq1.findtext("Tau20")),
                    a=float(eq1.findtext("A")),
                    b=float(eq1.findtext("B")),
                    c=float(eq1.findtext("C")),
                    e=float(eq1.findtext("E")),
                    d0=float(eq1.findtext("D0")),
                    d1=float(eq1.findtext("D1")),
                    d2=float(eq1.findtext("D2")),
                    h1=float(eq1.findtext("H1")),
                    h2=float(eq1.findtext("H2")),
                    h3=float(eq1.findtext("H3")),
                )
                if oxy_count == 0:
                    coeffs.oxygen_primary = obj
                else:
                    coeffs.oxygen_secondary = obj
                oxy_count += 1

        # --------------------------------------------------
        # CHLOROPHYLL
        # --------------------------------------------------
        chl = sensor.find("FluoroWetlabECO_AFL_FL_Sensor")
        if chl is not None:
            coeffs.chlorophyll = ECOCoefficients(
                slope=float(chl.findtext("ScaleFactor")),
                offset=float(chl.findtext("Vblank")),
            )

        # --------------------------------------------------
        # TURBIDITY
        # --------------------------------------------------
        tur = sensor.find("TurbidityMeter")
        if tur is not None:
            coeffs.turbidity = ECOCoefficients(
                slope=float(tur.findtext("ScaleFactor")),
                offset=float(tur.findtext("DarkVoltage")),
            )

        # --------------------------------------------------
        # ALTIMETER
        # --------------------------------------------------
        alt = sensor.find("AltimeterSensor")
        if alt is not None:
            coeffs.altimeter = AltimeterCoefficients(
                slope=float(alt.findtext("ScaleFactor")),
                offset=float(alt.findtext("Offset")),
            )

    return coeffs

