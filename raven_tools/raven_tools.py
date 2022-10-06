import glob
from pathlib import Path
import os


def set_filenames(year):
    """

    :param year:
    """
    forcing_paths[0] = f"data_obs/RhiresD_ch01h.swiss.lv95_{year}01010000_{year}12310000_clipped.nc"
    forcing_paths[1] = f"data_obs/TabsD_ch01r.swiss.lv95_{year}01010000_{year}12310000_pet.nc"
    forcing_paths[2] = f"data_obs/TmaxD_ch01r.swiss.lv95_{year}01010000_{year}12310000_clipped.nc"
    forcing_paths[3] = f"data_obs/TminD_ch01r.swiss.lv95_{year}01010000_{year}12310000_clipped.nc"


def write_rvt(start, end):
    """

    :param start:
    :param end:
    :type grd: GeoDataFrame
    :type filename: Path
    :param filename: Path to the grid weight text file
    :param grd: Grid derived from the netCDF file
    """
    # Write to GridWeights.txt
    date_list = [*range(start, end + 1)]
    file_type = ":FileType          rvt ASCII Raven 3.5"
    author = ":WrittenBy         Peter Zweifel"
    creation_date = ":CreationDate      April 2022"
    description = "#\n# Emulation of GR4J simulation of Broye\n#------------------------------------------------------------------------\n"
    gauge = [
        ":Gauge PYR2034\n",
        "  :Latitude    46.835913\n",
        "  :Longitude 6.9360708\n",
        "  :Elevation  441.0\n",
        ":EndGauge\n"
    ]

    flow_observation = [
        "# observed streamflow\n",
        ":RedirectToFile data_obs/BroPay_Q_2034_daily.rvt"
    ]
    with open((os.path.join(model_dir, model_name)), 'w') as ff:
        ff.write(f"#########################################################################\n"
                 f"{file_type}\n"
                 f"{author}\n"
                 f"{creation_date}\n")
        ff.write(f"{description}\n")
        ff.write(f"# meteorological forcings\n")

        for d in date_list:
            ff.write(f"# Year {d}\n")
            # set_filenames(d)
            forcing_type['Rainfall'][
                2] = f"    :FileNameNC           data_obs/RhiresD_v2.0_swiss.lv95/out/RhiresD_ch01h.swiss.lv95_{d}01010000_{d}12310000_clipped.nc"
            forcing_type['Average Temperature'][
                2] = f"    :FileNameNC           data_obs/TabsD_v2.0_swiss.lv95/out/TabsD_ch01r.swiss.lv95_{d}01010000_{d}12310000_clipped.nc"
            forcing_type['Maximum Temperature'][
                2] = f"    :FileNameNC           data_obs/TmaxD_v2.0_swiss.lv95/out/TmaxD_ch01r.swiss.lv95_{d}01010000_{d}12310000_clipped.nc"
            forcing_type['Minimum Temperature'][
                2] = f"    :FileNameNC           data_obs/TminD_v2.0_swiss.lv95/out/TminD_ch01r.swiss.lv95_{d}01010000_{d}12310000_clipped.nc"
            for f in forcing_type.values():
                for t in f:
                    ff.write(f"{t}\n")
            ff.write("\n")

        ff.writelines(gauge)
        ff.writelines(flow_observation)


if __name__ == '__main__':
    model_dir = Path("/media/mainman/Data/RAVEN/testmodels/GR4J/model")
    model_name = Path("raven_broye_gr4j.rvt")
    forcing_paths = [
        Path(model_dir, "data_obs/RhiresD_v2.0_swiss.lv95/out"),
        Path(model_dir, "data_obs/TabsD_v2.0_swiss.lv95/out"),
        Path(model_dir, "data_obs/TmaxD_v2.0_swiss.lv95/out"),
        Path(model_dir, "data_obs/TminD_v2.0_swiss.lv95/out"),
    ]

    tabsd_files = [f for f in glob.glob(str(forcing_paths[0]) + "/*.nc")]
    start_year = 2000
    end_year = 2018
    forcing_type = {
        'Rainfall': [
            ":GriddedForcing           Rainfall",
            "    :ForcingType          RAINFALL",
            "    :FileNameNC           data_obs/RhiresD_ch01h.swiss.lv95_{d}01010000_{d}12310000_clipped.nc",
            "    :VarNameNC            RhiresD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            "    :RedirectToFile       data_obs/GridWeights.txt ",
            ":EndGriddedForcing"],
        'Average Temperature': [
            ":GriddedForcing           Average Temperature",
            "    :ForcingType          TEMP_AVE",
            "    :FileNameNC           data_obs/TabsD_ch01r.swiss.lv95_{d}01010000_{d}12310000_clipped.nc",
            "    :VarNameNC            TabsD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            "    :RedirectToFile       data_obs/GridWeights.txt ",
            ":EndGriddedForcing"],
        'Maximum Temperature': [
            ":GriddedForcing           Maximum Temperature",
            "    :ForcingType          TEMP_MAX",
            "    :FileNameNC           data_obs/TmaxD_ch01r.swiss.lv95_{d}01010000_{d}12310000_clipped.nc",
            "    :VarNameNC            TmaxD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            "    :RedirectToFile       data_obs/GridWeights.txt ",
            ":EndGriddedForcing"],
        'Minimum Temperature': [
            ":GriddedForcing           Minimum Temperature",
            "    :ForcingType          TEMP_MIN",
            "    :FileNameNC           data_obs/TminD_ch01r.swiss.lv95_{d}01010000_{d}12310000_clipped.nc",
            "    :VarNameNC            TminD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            "    :RedirectToFile       data_obs/GridWeights.txt ",
            ":EndGriddedForcing"
        ]
    }

    for f in forcing_type.values():
        print(f)

    write_rvt(1981, 1982)