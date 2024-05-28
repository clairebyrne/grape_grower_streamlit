import rasterio
import numpy as np
import folium
from rasterio.features import shapes
import geopandas as gpd
from shapely.geometry import shape
import streamlit as st
from streamlit_folium import st_folium

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


########################################################################################
# process dtm
########################################################################################

@st.cache_data
def process_dtm():
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
def process_slp():
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
def process_aspect():
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


# run the processing functions
dtm_select = process_dtm()
st.write(f'Processed DTM selection giving {dtm_select.shape[0]} features')
slp_select = process_slp()
st.write(f'Processed slope selection giving {slp_select.shape[0]} features')
asp_select = process_aspect()
st.write(f'Processed aspect selection giving {asp_select.shape[0]} features')


#############################################################################
# Create a Folium map
#############################################################################

# Set the initial location to the center of the area of interest
latitude, longitude = 50.9601, 15.9751  # Replace with appropriate values
m = folium.Map(location=[latitude, longitude], zoom_start=10)

orange = {'fillColor': '#ff9302', 'color': '#ff9302'}
pink = {'fillColor': '#fc417c', 'color': '#fc417c'}
green = {'fillColor': '#5dd4a2', 'color': '#5dd4a2'}

col1, col2, col3 = st.columns(3)

with col1:
    if 'dtm_to_map' not in st.session_state:
        st.session_state.dtm_to_map = False
    def click_button():
        st.session_state.dtm_to_map = True

    st.button('Add DTM to map?', on_click=click_button)

    if st.session_state.dtm_to_map:
        st.write('Added DTM layer...')
        folium.GeoJson(dtm_select, name='DTM selection', style_function=lambda x:orange).add_to(m)

with col2:
    if 'slp_to_map' not in st.session_state:
        st.session_state.slp_to_map = False
    def click_button():
        st.session_state.slp_to_map = True

    st.button('Add Slope to map?', on_click=click_button)

    if st.session_state.slp_to_map:
        st.write('Added slope layer...')
        folium.GeoJson(slp_select, name='Slope selection', style_function=lambda x:pink).add_to(m)

with col3:
    if 'asp_to_map' not in st.session_state:
        st.session_state.asp_to_map = False
    def click_button():
        st.session_state.asp_to_map = True

    st.button('Add Aspect to map?', on_click=click_button)

    if st.session_state.asp_to_map:
        st.write('Added aspect layer...')
        folium.GeoJson(asp_select, name='Ideal Aspect', style_function=lambda x:green).add_to(m)

folium.LayerControl().add_to(m)

# Display the map in Streamlit

st_folium(m, width=700, height=500)
