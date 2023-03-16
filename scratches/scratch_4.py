import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from geopandas import GeoDataFrame
from netCDF4 import Dataset
from shapely.geometry import Polygon

# Read in the bounding box shape file from the clipping folder as a GeoPandas DataFrame
bbox = gpd.read_file(bounding_box_filename)
bbox.set_crs(2056, allow_override=True)
# bbox = bbox.to_crs("EPSG:2056")

# Read in the clipped netCDF file into a xarray Dataset
ds: Dataset = xr.open_dataset(netcdf_filepath, engine="netcdf4")

# Since we're only interested in the grid (not the actual cell values), only use 1 day
# Try to use 0-1 instead of strings to avoid having to use concrete date values.
# ds = ds.sel(time=slice('1962-01-01', '1962-01-02'))
start_date = f"01-01-{start_year}"
start_date = datetime.strptime(str(start_date), "%d-%m-%Y")
ds = ds.sel(time=slice(start_date, (start_date + timedelta(days=1))))
# Select the actual data value column
column_name = re.findall("[A-Za-z]+", forcing_name)[0]
xarr: xr.DataArray = ds[column_name]
# Convert Dataset to DataFrame and reset the index
df: pd.DataFrame = xarr.to_dataframe().reset_index()
# Convert the DataFrame to a GeoDataFrame, using the northing and easting provided by the netCDF. Furthermore,
# change the projection to LV95
data_column = pd.Series(df[column_name])
combi_catchment_grid: GeoDataFrame = gpd.GeoDataFrame(data_column, geometry=gpd.points_from_xy(df.E, df.N),
                                                      crs="EPSG:2056")

# From the GeoDataFrame that contains the netCDF grid, get the extent as points
xmin, ymin, xmax, ymax = combi_catchment_grid.total_bounds
# Since the netCDF has a cell size of 1000m, set this here
length: int = 1000
wide: int = 1000

# Create the northing and easting coordinates of the bounding vertex for each cell
cols: list = list(np.arange(xmin, xmax + wide, wide))
rows: list = list(np.arange(ymin, ymax + length, length))

# initialize the Polygon list
polygons: list[Polygon] = []
cell_id: list[str] = []

# Create the GeoDataFrame with the grid, providing column names plus polygons for the geometry
grid_cols = ['row', 'col', 'cell_id', 'polygons', 'area', 'area_rel']
grid: GeoDataFrame = gpd.GeoDataFrame(columns=grid_cols, geometry='polygons', crs="EPSG:2056")
# Loop over each cell and create the corresponding Polygon, then append it to the Polygon list
# Also write the cell id from cell_id=i_row*n_columns + i_column
for ix, x in enumerate(cols):
    for iy, y in enumerate(rows):
        polygons.append(Polygon([(x, y), (x + wide, y), (x + wide, y + length), (x, y + length)]))
        cid = int(iy * len(cols) + ix)
        cell_id.append(f"{str(cid)}")

# Use the polygon list in the GeoDataFrame and set the projection accordingly
grid["polygons"] = polygons
grid["cell_id"] = cell_id

# Every cell that is not within the catchment area will have area_rel set to zero
grid["area_rel"] = 0
grid = grid.set_crs(2056, allow_override=True)

if export_shp:
    # Export the grid to a shape file
    grid.to_file(str(str(out_path) + ".shp"))
return grid
