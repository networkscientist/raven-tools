# import pandas as pd
# import geopandas as gpd
# data_path = "/home/sirian/Applications/Hydrology/BasinMaker/broye_gr4j/data/"
# out_path = "/home/sirian/Applications/Hydrology/BasinMaker/broye_gr4j/output/"
# shapes = gpd.read_file(data_path + "/Catchment/Broye_Payerne.shp")
# df = gpd.read_file("/home/sirian/Applications/Hydrology/RAVEN/testmodels/GR4J/catchment/CH-0057_attributes.geojson")
# required_columns = [
#     "id",
#     "area_ch1903plus",
#     "lab_x",
#     "lab_y",
#     "water_name",
#     "place",
#     "a0401_eu_dem_v11_e40n20crp_chv1_0",
#     "a0404_eu_dem_v11_e40n20_slp8v1_0",
#     "a0407_eu_dem_v11_asp8sm_maskv1_0"
# ]
# df_new=df[required_columns]
# cmt = gpd.read_file("/home/sirian/Applications/Hydrology/RAVEN/testmodels/GR4J/catchment/Broye_Payerne.shp")
# cmt.plot()
# cmt.to_crs(epsg=2056).plot()
# cmt_new = cmt.merge(df_new, on='id', )
# cmt_new.to_file("result.shp")