import rasterio
import numpy as np
import folium
from rasterio.features import shapes
import geopandas as gpd
from shapely.geometry import shape
import streamlit as st
from streamlit_folium import st_folium, folium_static

# Define the path to your DTM file and the elevation threshold
dtm_file_path = r'./data/raw_dtm_merge_10m.tiff'
slp_file_path = r'./data/slope_dtm_merge_10m_resample.tif'
asp_file_path = r'./data/aspect_dtm_merge_10m_resample.tif'

st.title('Grape grower')

########################################################################################
# add user selectors to the app
########################################################################################

st.markdown('## Selection Criteria')
st.markdown('#### Elevation')
# add elevation slider
elevation_threshold = st.slider('Select max elevation (m)', min_value=0, max_value=1000, 
                          value=300, step=50) 

st.markdown('#### Slope')
# add slope sliders for min and max
min_slope = st.slider('Select min slope (degrees)', min_value=0, max_value=88, 
                          value=25, step=1) 
max_slope = st.slider('Select max slope (degrees)', min_value=0, max_value=88, 
                          value=45, step=1)
st.markdown('#### Ideal Aspect')
st.write('(Could allow to put another _ok aspect_ option later)')
st.write('Aspect is in compass degrees where 0 is North, 180 is South')
min_ideal_asp = st.slider('Select min ideal aspect (degrees)', min_value=0, max_value=360, 
                          value=180, step=1) 
max_ideal_asp = st.slider('Select max ideal aspect (degrees)', min_value=0, max_value=360, 
                          value=220, step=1)

min_area = st.slider('Select minimum area (in hectares)', min_value=0.1, max_value=50.1, value=1.0, step=0.5)

########################################################################################
# process dtm
########################################################################################

@st.cache_data
def process_dtm(elevation_threshold):
    # Read the DTM data from the .tiff file
    with rasterio.open(dtm_file_path, 'r+') as src:
        src.crs= rasterio.crs.CRS({"init": "epsg:2180"})
        dtm = src.read(1)
        transform = src.transform

    # Create a mask for areas below the given elevation
    below_dtm_threshold = dtm < elevation_threshold

    # Get the shapes of areas below the threshold
    mask_shapes = shapes(below_dtm_threshold.astype(np.uint8), transform=transform)

    # Convert shapes to a GeoDataFrame
    geometries = [shape(geom) for geom, value in mask_shapes if value == 1]
    dtm_select = gpd.GeoDataFrame(geometry=geometries, crs="EPSG:2180")
    return dtm_select

########################################################################################
# process slope
########################################################################################

@st.cache_data
def process_slp(min_slope, max_slope):
    # Read the slope data from the .tiff file
    with rasterio.open(slp_file_path, 'r+') as src:
        src.crs= rasterio.crs.CRS({"init": "epsg:2180"})
        slp = src.read(1)
        transform = src.transform

    # Create a mask for areas within the selected slope range
    within_slp_threshold = (slp >= min_slope) & (slp <= max_slope)

    # Get the shapes of areas within the threshold
    mask_shapes = shapes(within_slp_threshold.astype(np.uint8), transform=transform)

    # Convert shapes to a GeoDataFrame
    geometries = [shape(geom) for geom, value in mask_shapes if value == 1]
    slp_select = gpd.GeoDataFrame(geometry=geometries, crs="EPSG:2180")
    return slp_select

# ########################################################################################
# # process aspect
# ########################################################################################

@st.cache_data
def process_aspect(min_ideal_asp, max_ideal_asp):
    # Read the slope data from the .tiff file
    with rasterio.open(asp_file_path, 'r+') as src:
        src.crs= rasterio.crs.CRS({"init": "epsg:2180"})
        asp = src.read(1)
        transform = src.transform

    # Create a mask for areas within the selected slope range
    within_asp_threshold = (asp >= min_ideal_asp) & (asp <= max_ideal_asp)

    # Get the shapes of areas within the threshold
    mask_shapes = shapes(within_asp_threshold.astype(np.uint8), transform=transform)

    # Convert shapes to a GeoDataFrame
    geometries = [shape(geom) for geom, value in mask_shapes if value == 1]
    asp_select = gpd.GeoDataFrame(geometry=geometries, crs="EPSG:2180")
    return asp_select

@st.cache_data
def all_criteria_met(_elev, _slp, _asp, smallest):
            # intersect all areas that meet selected criteria
            intersect_all_criteria = _elev.overlay(_slp).overlay(_asp)
            # calculate area of intersected areas in hectares
            intersect_all_criteria['area_ha'] = intersect_all_criteria['geometry'].area/10**4
            # select the areas that meet the defined area criteria into a new geodataframe
            all_criteria = intersect_all_criteria[intersect_all_criteria['area_ha']>= smallest]
            # select the areas that do not meet the defined area criteria into a new geodataframe
            all_criteria_small = intersect_all_criteria[intersect_all_criteria['area_ha']< smallest]
            return all_criteria, all_criteria_small
     

if 'run_analysis' not in st.session_state:
        st.session_state.run_analysis = False
def click_button():
            st.session_state.run_analysis = True

st.button('Click to run first analysis', on_click=click_button)

if st.session_state.run_analysis:
            # run the processing functions
            dtm_select = process_dtm(elevation_threshold)
            st.write(f'Processed DTM selection giving {dtm_select.shape[0]} features')
            slp_select = process_slp(min_slope, max_slope)
            st.write(f'Processed slope selection giving {slp_select.shape[0]} features')
            asp_select = process_aspect(min_ideal_asp, max_ideal_asp)
            st.write(f'Processed aspect selection giving {asp_select.shape[0]} features')
            all_criteria, all_criteria_small = all_criteria_met(dtm_select, slp_select, asp_select, min_area)

            st.write(f'Ran intersect of all criteria of elevation, slope, aspect and minimum area. {all_criteria.shape[0]} features meet all criteria. Two layers have been added to the map - **All criteria met**, and **All criteria met (small)**.')

#############################################################################
# Create a Folium map
#############################################################################

latitude, longitude = 50.9601, 15.9751  # Replace with appropriate values
m = folium.Map(location=[latitude, longitude], tiles='OpenStreetMap', zoom_start=10)
folium.TileLayer('cartodbpositron').add_to(m)
#folium.TileLayer('cartodbdark_matter').add_to(m)
folium.TileLayer(
        tiles = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        attr = 'Esri',
        name = 'Esri Satellite',
        overlay = False,
        control = True
       ).add_to(m)

# criteria met styles
blue = {'fillColor': '#021076', 'color': '#021076'}
yellow = {'fillColor': '#eed959', 'color': '#eed959'}

popup = folium.GeoJsonPopup(fields=["area_ha"])

folium.GeoJson(all_criteria, name='All criteria met', 
               #popup=popup, 
               #zoom_on_click=True, 
               style_function=lambda x:blue).add_to(m)
folium.GeoJson(all_criteria_small, name='All criteria met (small)', 
               #popup=popup, 
               #zoom_on_click=True, 
               style_function=lambda x:yellow).add_to(m)

folium.LayerControl().add_to(m)
st_folium(m, width=700, height=500)