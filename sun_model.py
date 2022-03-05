import math

import numpy as np


def liu_jordan(day_of_year: int,
               hour_of_day: int,
               tilt_degrees: float,
               azimuth_degrees: float,
               latitude: float,
               indirect_radiation: float,
               direct_radiation: float,
               reflectivity: float) -> float:
    """Liu Jordan model for estimating solar radiation on titled surface at the given latitude.

    Args:
        day_of_year (int): Day of year.
        hour_of_day (int): Hour of day in 0-23 convention
        tilt_degrees (float): Angle of tilt in degrees in range 0-90
        azimuth_degrees (float): Angle of azimuth in degrees.
        latitude (float): Latitude in degrees.
        indirect_radiation (float): Indirect solar radiation in W/m2.
        direct_radiation (float): Direct solar radiation in W/m2.
        reflectivity (float): Reflectivity of the surface.

    Returns:
        float: Estimated solar radiation in W/m2
    """
    degrees_to_radians = 3.14159 / 180

    omega = 15 * (12 - hour_of_day)
    som = math.sin(omega * degrees_to_radians)
    com = math.cos(omega * degrees_to_radians)

    deklin = 23.45 * math.sin(360 * (284 + day_of_year) / 365 * degrees_to_radians)
    sin_deklin = math.sin(degrees_to_radians * deklin)
    cos_deklin = math.cos(degrees_to_radians * deklin)

    sin_tilt = math.sin(degrees_to_radians * tilt_degrees)
    cos_tilt = math.cos(degrees_to_radians * tilt_degrees)

    sin_azimuth = math.sin(degrees_to_radians * azimuth_degrees)
    cos_azimuth = math.cos(degrees_to_radians * azimuth_degrees)

    sin_latitude = math.sin(degrees_to_radians * latitude)
    cos_latitude = math.cos(degrees_to_radians * latitude)

    R0 = (1 - cos_tilt) / 2
    Rd = (1 + cos_tilt) / 2

    nominator_1 = sin_deklin * (sin_latitude * cos_tilt -
                                cos_latitude * sin_tilt * cos_azimuth)
    nominator_2 = cos_deklin * (cos_latitude * cos_tilt * com +
                                sin_latitude * sin_tilt * cos_azimuth * com +
                                sin_tilt * sin_azimuth * som)
    denominator = sin_deklin * sin_latitude + cos_deklin * cos_latitude * com
    Rb = (nominator_1 + nominator_2) / denominator
    Rb = np.clip(Rb, 0, 3)

    total_radiation = Rb * indirect_radiation + \
                      Rd * direct_radiation + \
                      (indirect_radiation + direct_radiation) * reflectivity * R0

    return total_radiation
