import obspy
import numpy as np
import pathlib
from pyproj import Proj,transform,Geod
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import rasterio
from rasterio.plot import show
from rasterio.warp import calculate_default_transform, reproject, Resampling
from datetime import datetime
from datetime import timedelta
from matplotlib.dates import date2num, DateFormatter
from location.compute_backazimuths import get_station_coordinates
from location.compute_backazimuths import get_station_grid_locations
import geopandas as gpd
import cartopy
import cartopy.crs as ccrs
from shapely import geometry
from collections import namedtuple



def transform_imagery(file,dst_crs):
    with rasterio.open(file) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rasterio.open('data/imagery/' + dst_crs + '.TIF', 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)
            
            
            
def get_station_coordinates(xml_path):
    stat_coords = []
    inv = obspy.read_inventory(xml_path + "*XML")
    for s in inv.get_contents()['stations']:
        channel = inv.get_contents()['networks'][0] + "." + s.split(' ')[0].split('.')[1] + "..HHZ"
        lat = inv.get_coordinates(channel)["latitude"]
        lon = inv.get_coordinates(channel)["longitude"]
        stat_coords.append([lon,lat])
    _, idx = np.unique(stat_coords,axis=0,return_index=True)
    stat_coords = np.array(stat_coords)[np.sort(idx)]
    return stat_coords



def get_station_grid_locations(station_lon_lat_coords,crs):
    # convert station coordinates to x and y and take average station location
    p2 = Proj(crs,preserve_units=False)
    p1 = Proj(proj='latlong',preserve_units=False)
    [stat_x,stat_y] = transform(p1,p2,station_lon_lat_coords[:,0],station_lon_lat_coords[:,1])
    return np.stack((stat_x,stat_y),axis=1)            
            
    
    
def plot_catalog_and_big_event_backazimuths(backazimuths,big_event_backazimuth,array_centroid,station_grid_coords,colors):

    # open LANDSAT imagery file and plot as it is, in EPSG:3031
    original_file = "data/imagery/epsg:3245.TIF"

    sat_imagery = rasterio.open(original_file)
    sat_data = sat_imagery.read(1)

    # Construct figure and axis to plot on
    fig,ax = plt.subplots(figsize=(12,12))
    axes_coords = np.array([0, 0, 1, 1])
    ax_image = fig.add_axes(axes_coords)

    # get corners of imagery extent
    p2 = Proj("EPSG:3245",preserve_units=False)
    p1 = Proj(proj='latlong',preserve_units=False)

    # plot imagery
    bounds = sat_imagery.bounds
    horz_len = bounds[2]-bounds[0]
    vert_len = bounds[3]-bounds[1]
    plot_bounds = [bounds[0]+0.35*horz_len,bounds[2]-0.25*horz_len,bounds[1]+0.25*vert_len,bounds[3]-0.4*vert_len]
    c_map = [np.min(sat_data),np.max(sat_data)]
    c_map_range = c_map[1]-c_map[0]
    ax_image.imshow(sat_data,cmap='gray',extent=[bounds[0],bounds[2],bounds[1],bounds[3]],interpolation='none', vmin=np.min(sat_data), vmax=np.max(sat_data)-0.2*c_map_range)
    #ax_image.imshow(sat_data,cmap='gray',extent=[bounds[0],bounds[2],bounds[1],bounds[3]])

    # define, transform, and plot lat/lon grid
    lat = [-74,-74.5,-75,-75.5]
    lon = [-98,-100,-102,-104]
    x_lab_pos=[]
    y_lab_pos=[]
    line = np.linspace(-110,-90,100)
    for i in lat:
        line_x,line_y = transform(p1,p2,line,np.linspace(i,i,100))
        ax_image.plot(line_x,line_y,linestyle='--',linewidth=2,c='w',alpha=0.75)
        y_lab_pos.append(line_y[np.argmin(np.abs(line_x-plot_bounds[0]))])
    line = np.linspace(-80,-70,100)
    for i in lon:
        line_x,line_y = transform(p1,p2,np.linspace(i,i,100),line)
        ax_image.plot(line_x,line_y,linestyle='--',linewidth=2,c='w',alpha=0.75)
        x_lab_pos.append(line_x[np.argmin(np.abs(line_y-plot_bounds[2]))])

    # set ticks and labels for lat/lon grid
    ax_image.set_xticks(x_lab_pos)
    ax_image.set_xticklabels(labels=[str(lon[0]) + '$^\circ$',str(lon[1]) + '$^\circ$',str(lon[2]) + '$^\circ$',str(lon[3]) + '$^\circ$'],fontsize=25)
    ax_image.set_xlabel("Longitude",fontsize=25)
    ax_image.set_yticks(y_lab_pos)
    ax_image.set_yticklabels(labels=[str(lat[0]) + '$^\circ$',str(lat[1]) + '$^\circ$',str(lat[2]) + '$^\circ$',str(lat[3]) + '$^\circ$'],fontsize=25)
    ax_image.set_ylabel("Latitude",fontsize=25)
    ax_image.set_xlim([plot_bounds[0],plot_bounds[1]])
    ax_image.set_ylim([plot_bounds[2],plot_bounds[3]])
    #ax_image.set_title("Event Backazimuths",fontsize=25)
    
    # properly center the polar plot on the array centroid
    x_pos = (array_centroid[0]-plot_bounds[0])/(plot_bounds[1]-plot_bounds[0])
    y_pos = (array_centroid[1]-plot_bounds[2])/(plot_bounds[3]-plot_bounds[2])
    width = 0.3

    # make polar plot centered at array centroid
    ax_polar = fig.add_axes([x_pos-width/2,y_pos-width/2,width,width], projection = 'polar')
    ax_polar.set_theta_zero_location('N')
    ax_polar.set_theta_direction(-1)

    radius,bins = np.histogram(backazimuths[~np.isnan(backazimuths)]*np.pi/180,bins=np.linspace(0,2*np.pi,37))
    patches = ax_polar.bar(bins[:-1], radius, zorder=1, align='edge', width=np.diff(bins),
                     edgecolor='black',facecolor=colors[0],fill=True, linewidth=1,alpha = .5)
    radius,bins = np.histogram(big_event_backazimuth*np.pi/180,bins=np.linspace(0,2*np.pi,37))
    patches = ax_polar.bar(bins[:-1], radius*3, zorder=1, align='edge', width=np.diff(bins),
                     edgecolor='black',facecolor=colors[1],fill=True, linewidth=1,alpha = .5)

    # Remove ylabels for area plots (they are mostly obstructive)
    ax_polar.set_yticks([])
    ax_polar.axis('off')

    # plot station locations
    ax_stats = fig.add_axes(axes_coords)
    ax_stats.scatter(station_grid_coords[:,0],station_grid_coords[:,1],marker="^",c='black',s=100)
    ax_stats.set_xlim([plot_bounds[0],plot_bounds[1]])
    ax_stats.set_ylim([plot_bounds[2],plot_bounds[3]])
    ax_stats.axis('off')
    
    # add legend for spatial groups
    x,y = transform(p1,p2,-102.75,-75.175)
    box = matplotlib.patches.Rectangle((x-1000, y-5000), 20000, 7000, linewidth=1, edgecolor='k', facecolor='w')
    ax_stats.add_patch(box)
    ax_stats.text(x, y-500, "Events",c=colors[0],fontsize=20)
    ax_stats.text(x, y-3000, "Big event",c=colors[1],fontsize=20)
   
    # add North arrow
    line_x,line_y = transform(p1,p2,np.linspace(-102.2,-102.2,100),np.linspace(-74.65,-74.6,100))
    ax_stats.plot(line_x,line_y,color='w',linewidth = 5)
    ax_stats.scatter(line_x[-1],line_y[-1],marker=(3,0,4),c='w',s=400)
    ax_stats.text(line_x[-1]-3500,line_y[-1]-2500,"N",color='w',fontsize=25)
    
    # add scale bar
    ax_stats.plot([plot_bounds[0]+(plot_bounds[1]-plot_bounds[0])/2-5000,plot_bounds[0]+(plot_bounds[1]-plot_bounds[0])/2+15000],[plot_bounds[2]+12000,plot_bounds[2]+12000],color='w',linewidth = 5)
    ax_stats.text(plot_bounds[0]+(plot_bounds[1]-plot_bounds[0])/2,plot_bounds[2]+9000,"20 km",color='w',fontsize=25)

    # add inset figure of antarctica
    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    ax_inset = fig.add_axes([0.765,0.7,0.275,0.275],projection = ccrs.SouthPolarStereo())
    ax_inset.set_extent([-180, 180, -90, -65], crs=ccrs.PlateCarree())
    geom = geometry.box(minx=-103,maxx=-99,miny=-75.5,maxy=-74.5)
    ax_inset.add_geometries([geom], crs=ccrs.PlateCarree(), edgecolor='r',facecolor='none', linewidth=1)
    ax_inset.add_feature(cartopy.feature.OCEAN, facecolor='#A8C5DD', edgecolor='none')
    
    plt.savefig("outputs/catalog_and_big_event_backazimuths.png",bbox_inches="tight",dpi=100)

    
    
def plot_big_event_backazimuths(backazimuths,freq,array_centroid,station_grid_coords,colors):

    # open LANDSAT imagery file and plot as it is, in EPSG:3031
    original_file = "data/imagery/epsg:3245.TIF"

    sat_imagery = rasterio.open(original_file)
    sat_data = sat_imagery.read(1)

    # Construct figure and axis to plot on
    fig,ax = plt.subplots(figsize=(12,12))
    axes_coords = np.array([0, 0, 1, 1])
    ax_image = fig.add_axes(axes_coords)

    # get corners of imagery extent
    p2 = Proj("EPSG:3245",preserve_units=False)
    p1 = Proj(proj='latlong',preserve_units=False)

    # plot imagery
    bounds = sat_imagery.bounds
    horz_len = bounds[2]-bounds[0]
    vert_len = bounds[3]-bounds[1]
    plot_bounds = [bounds[0]+0.35*horz_len,bounds[2]-0.25*horz_len,bounds[1]+0.25*vert_len,bounds[3]-0.4*vert_len]
    c_map = [np.min(sat_data),np.max(sat_data)]
    c_map_range = c_map[1]-c_map[0]
    ax_image.imshow(sat_data,cmap='gray',extent=[bounds[0],bounds[2],bounds[1],bounds[3]],interpolation='none', vmin=np.min(sat_data), vmax=np.max(sat_data)-0.2*c_map_range)
    #ax_image.imshow(sat_data,cmap='gray',extent=[bounds[0],bounds[2],bounds[1],bounds[3]])

    # define, transform, and plot lat/lon grid
    lat = [-74,-74.5,-75,-75.5]
    lon = [-98,-100,-102,-104]
    x_lab_pos=[]
    y_lab_pos=[]
    line = np.linspace(-110,-90,100)
    for i in lat:
        line_x,line_y = transform(p1,p2,line,np.linspace(i,i,100))
        ax_image.plot(line_x,line_y,linestyle='--',linewidth=2,c='w',alpha=0.75)
        y_lab_pos.append(line_y[np.argmin(np.abs(line_x-plot_bounds[0]))])
    line = np.linspace(-80,-70,100)
    for i in lon:
        line_x,line_y = transform(p1,p2,np.linspace(i,i,100),line)
        ax_image.plot(line_x,line_y,linestyle='--',linewidth=2,c='w',alpha=0.75)
        x_lab_pos.append(line_x[np.argmin(np.abs(line_y-plot_bounds[2]))])

    # set ticks and labels for lat/lon grid
    ax_image.set_xticks(x_lab_pos)
    ax_image.set_xticklabels(labels=[str(lon[0]) + '$^\circ$',str(lon[1]) + '$^\circ$',str(lon[2]) + '$^\circ$',str(lon[3]) + '$^\circ$'],fontsize=25)
    ax_image.set_xlabel("Longitude",fontsize=25)
    ax_image.set_yticks(y_lab_pos)
    ax_image.set_yticklabels(labels=[str(lat[0]) + '$^\circ$',str(lat[1]) + '$^\circ$',str(lat[2]) + '$^\circ$',str(lat[3]) + '$^\circ$'],fontsize=25)
    ax_image.set_ylabel("Latitude",fontsize=25)
    ax_image.set_xlim([plot_bounds[0],plot_bounds[1]])
    ax_image.set_ylim([plot_bounds[2],plot_bounds[3]])
    #ax_image.set_title("Event Backazimuths",fontsize=25)
    
    # properly center the polar plot on the array centroid
    x_pos = (array_centroid[0]-plot_bounds[0])/(plot_bounds[1]-plot_bounds[0])
    y_pos = (array_centroid[1]-plot_bounds[2])/(plot_bounds[3]-plot_bounds[2])
    width = 0.3

    # make polar plot centered at array centroid
    ax_polar = fig.add_axes([x_pos-width/2,y_pos-width/2,width,width], projection = 'polar')
    ax_polar.set_theta_zero_location('N')
    ax_polar.set_theta_direction(-1)

    for i in range(len(backazimuths)):
        radius,bins = np.histogram(backazimuths[i]*np.pi/180,bins=np.linspace(0,2*np.pi,int(360)))
        patches = ax_polar.bar(bins[:-1], radius*3, zorder=1, align='edge', width=np.diff(bins),
                         edgecolor='black',facecolor=colors[i],fill=True, linewidth=1,alpha = .5)

    # Remove ylabels for area plots (they are mostly obstructive)
    ax_polar.set_yticks([])
    ax_polar.axis('off')

    # plot station locations
    ax_stats = fig.add_axes(axes_coords)
    ax_stats.scatter(station_grid_coords[:,0],station_grid_coords[:,1],marker="^",c='black',s=100)
    ax_stats.set_xlim([plot_bounds[0],plot_bounds[1]])
    ax_stats.set_ylim([plot_bounds[2],plot_bounds[3]])
    ax_stats.axis('off')
    
    # add legend for spatial groups
    x,y = transform(p1,p2,-102.75,-75.175)
    box = matplotlib.patches.Rectangle((x-1000, y-5000), 20000, 7000, linewidth=1, edgecolor='k', facecolor='w')
    ax_stats.add_patch(box)
    for i in range(len(backazimuths)):
        ax_stats.text(x, y-500-3000*i, str(freq[i]) + " Hz",c=colors[i],fontsize=20)
   
    # add North arrow
    line_x,line_y = transform(p1,p2,np.linspace(-102.2,-102.2,100),np.linspace(-74.65,-74.6,100))
    ax_stats.plot(line_x,line_y,color='w',linewidth = 5)
    ax_stats.scatter(line_x[-1],line_y[-1],marker=(3,0,4),c='w',s=400)
    ax_stats.text(line_x[-1]-3500,line_y[-1]-2500,"N",color='w',fontsize=25)
    
    # add scale bar
    ax_stats.plot([plot_bounds[0]+(plot_bounds[1]-plot_bounds[0])/2-5000,plot_bounds[0]+(plot_bounds[1]-plot_bounds[0])/2+15000],[plot_bounds[2]+12000,plot_bounds[2]+12000],color='w',linewidth = 5)
    ax_stats.text(plot_bounds[0]+(plot_bounds[1]-plot_bounds[0])/2,plot_bounds[2]+9000,"20 km",color='w',fontsize=25)

    # add inset figure of antarctica
    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    ax_inset = fig.add_axes([0.765,0.7,0.275,0.275],projection = ccrs.SouthPolarStereo())
    ax_inset.set_extent([-180, 180, -90, -65], crs=ccrs.PlateCarree())
    geom = geometry.box(minx=-103,maxx=-99,miny=-75.5,maxy=-74.5)
    ax_inset.add_geometries([geom], crs=ccrs.PlateCarree(), edgecolor='r',facecolor='none', linewidth=1)
    ax_inset.add_feature(cartopy.feature.OCEAN, facecolor='#A8C5DD', edgecolor='none')
    
    plt.savefig("outputs/big_event_backazimuths.png",bbox_inches="tight",dpi=100)