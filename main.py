import matplotlib.pyplot as plt
import pandas as pd

import sun_model as sm


def get_solar_radiation_for_year(dataset_path: str,
                                 tilt_degrees: float,
                                 azimuth_degrees: float,
                                 latitude: float,
                                 reflectivity: float) -> pd.DataFrame:
    dataset = pd.read_csv(dataset_path)
    solar_radiation = []

    for _, row in dataset.iterrows():
        solar_radiation.append(
            sm.liu_jordan(day_of_year=row['Day'],
                          hour_of_day=row['Hour'],
                          tilt_degrees=tilt_degrees,
                          azimuth_degrees=azimuth_degrees,
                          latitude=latitude,
                          indirect_radiation=row['IDH'], direct_radiation=row['ISH'],
                          reflectivity=reflectivity))

    dataset['LJ_model'] = solar_radiation
    return dataset


if __name__ == '__main__':
    simulation_df = get_solar_radiation_for_year(
        dataset_path='LiuJordanData.csv',
        tilt_degrees=90.,
        azimuth_degrees=-45.,
        latitude=52.1,
        reflectivity=0.1)
    print(simulation_df)

    simulation_df[['IDH', 'ISH', 'LJ_model']].plot()
    plt.show()
