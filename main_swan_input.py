import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr
import datetime

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.offsetbox import AnchoredText

#import oceanwaves


def swan_grid_reg(lon_w, lon_e, lat_s, lat_n, elemente_size):

    """

    Generate Swan regular grid

    Parameters:
    -----------

    lon_e (float): longitude east point
    lon_w (float): longitude west point

    lat_s (float): latitude south point
    lon_n (float): latitude north point

    elemente_size (float): ERA5 grid size resolution in degree


    Returns:
    --------

    grid plot shape


    """

    # GRID GENERATION FOR SWAN

    ele_ = 1 / elemente_size

    lon_ = np.arange(lon_w, lon_e, 1 / ele_)
    lat_ = np.arange(lat_s, lat_n, 1 / ele_)

    x, y = np.meshgrid(lon_, lat_)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([lon_e, lon_w, lat_s, lat_n], crs=ccrs.PlateCarree())

    plt.scatter(x, y)
    # Put a background image on for nice sea rendering.
    ax.stock_img()

    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

    SOURCE = 'Natural Earth'
    LICENSE = 'public domain'

    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(states_provinces, edgecolor='gray')

    # Add a text annotation for the license information to the
    # the bottom right corner.
    # text = AnchoredText('\u00A9 {}; license: {}'
    #                    ''.format(SOURCE, LICENSE),
    #                    loc=4, prop={'size': 12}, frameon=True)
    # ax.add_artist(text)

    plt.show()

    # GENERATE COORDINATE

    # Extract boundary points from the grid

    vert_N = [x[-1, :], y[-1, :]]
    vert_E = [x[:, -1], y[:, -1]]
    vert_S = [x[0, :], y[0, :]]
    vert_W = [x[:, 0], y[:, 0]]

    # GENERATE COORDINATE grid for ERA5

    ele_era5 = int((ele_) / (1 / 0.5))
    coord_N = np.column_stack((vert_N[0], vert_N[1]))[:][::ele_era5]
    coord_S = np.column_stack((vert_S[0], vert_S[1]))[:][::ele_era5]
    coord_E = np.column_stack((vert_E[0], vert_E[1]))[:][::ele_era5]
    coord_W = np.column_stack((vert_W[0], vert_W[1]))[:][::ele_era5]

    return vert_N, vert_E, vert_S, vert_W, coord_N, coord_S, coord_E, coord_W


def swan_grid_curv(vert_S, vert_E, filename):

    """
    Write SWAN GRID for COAWST system.

    Parameters:
    -----------

    - Vert_S (array): Vert_S generated from reg_grid_swan

    - Vert_E (array): Vert_E generated from reg_grid_swan

    - filename (str): Name the file grid output


    Returns:
    --------

    - swan grid (grid.txt)
    - write line for grid in swan.swn

    """

    x = vert_S[0]
    y = vert_E[1]

    # Grid resolution

    reso_x = np.round(np.abs(x[0] - x[1]), 3)
    reso_y = np.round(np.abs(x[0] - x[1]), 3)

    # Adapting to grid to SWAN

    vert_east = np.append(y, y[-1] + reso_x)
    vert_south = np.append(x, x[-1] + reso_y)

    # Generating grid

    xi, yi = np.meshgrid(vert_south, vert_east)

    grid = np.hstack([xi.flatten(), yi.flatten()])
    grid_ = np.round(grid, 3)
    df = pd.DataFrame({'grid': grid_})

    # Save grid

    filename_grid = filename
    df.to_csv(filename_grid, header=False, index=False)

    # Write swan line to grid input

    shape = xi.shape
    shape_x = shape[1]
    shape_y = shape[0]
    file = open('swan_grid_input.txt', 'w')
    file.write('&&----SWAN GRID----&& \n\n')
    file.write(f'CGRID CURV {shape_x - 1} {shape_y - 1} CIR 36 0.05 1 24 \n')
    file.write(f"READGRID COORDINATES 1 '{filename_grid}' 4 0 0 FREE \n")
    file.close()

    return xi, yi


def boundary(era5_data, coords):
    bound_wave = []

    for b in range(len(coords)):
        bound = era5_data.sel(

            dict(longitude=(coords[b][0]), latitude=coords[b][1]),
            method='nearest',

        )

        bound_wave.append(bound)

    boundary = xr.concat(bound_wave, dim="bound_wave")

    return boundary


def wave_era52swan_TPAR(boundary_swan, location_name, output_paths=''):

    """
    Generate boundary condition to SWAN from ERA5 WAVE

    PARAMETERS:
    -----------

    boundary (netcdf): generate by boundary function

    location (str): Which location is the input wave parameters
                    e.g. N = North, E = East, S = South and W = West

    RETURNS
    --------

    boundary file to run in SWAN

    """

    for i in range(len(boundary_swan.bound_wave)):
        df = pd.DataFrame()
        df['time'] = time_model
        df['swh'] = boundary_swan.swh.values[i]
        df['tp'] = boundary_swan.pp1d.values[i]
        df['dir'] = boundary_swan.mwd.values[i]
        df['spread'] = np.ones(len(time_model)) * 20

        file = open(output_paths + 'boundary_wave_' + location_name + '_' + str(i) + '.txt', 'w')
        file.write('TPAR \n')
        file.write(df.to_string(header=False, index=False))
        file.close()


def write_namelist_wave_boundary(location, vert, boundary_swan, grid_size_x, grid_size_y):

    """
    Location (str): use for coordinates: east, west, south, north

    vert (str): boundary position (N, S, E or W)

    boundary_swan: Boundary location respectively to position (generated from boundary function!!)

    grid_size_x (int): grid size in longitude

    grid_size_y (int): grid size in latitude

    """

    resolution = int(len(vert[0]) / len(boundary_swan.bound_wave))
    size = len(boundary_swan.bound_wave)

    x = grid_size_x
    y = grid_size_y

    if location == 'east':

        file = open('boundary_namelist_wave_E.txt', 'w')
        file.write('&&----wave_boundary----&& \n\n')
        file.write('BOUND SHAPESPEC JONSWAP MEAN DSPR DEGREES \n')

        name = 'E'
        resol = resolution
        bound_count = 0
        lat_0 = 1

        while resol <= resolution * size:
            file.write(
                f"BOUNDSPEC SEGMENT IJ {grid_size_x} {lat_0} {grid_size_x} {resol} VARIABLE FILE 0 'boundary_wave_" + name + "_" + str(
                    bound_count) + ".txt \n\n")
            lat_0 += resolution
            resol += resolution
            bound_count += 1

        file.close()


    elif location == 'west':

        file = open('boundary_namelist_wave_W.txt', 'w')
        file.write('&&----wave_boundary----&& \n\n')
        file.write('BOUND SHAPESPEC JONSWAP MEAN DSPR DEGREES \n')

        name = 'W'
        resol = resolution
        bound_count = 0
        lat_0 = 1

        while resol <= resolution * size:
            file.write(
                f"BOUNDSPEC SEGMENT IJ 1 {lat_0} 1 {resol} VARIABLE FILE 0 'boundary_wave_" + name + "_" + str(
                    bound_count) + ".txt \n\n")
            lat_0 += resolution
            resol += resolution
            bound_count += 1

        file.close()

    elif location == 'south':

        file = open('boundary_namelist_wave_S.txt', 'w')
        file.write('&&----wave_boundary----&& \n\n')
        file.write('BOUND SHAPESPEC JONSWAP MEAN DSPR DEGREES \n')

        name = 'S'
        resol = resolution
        bound_count = 0
        lon_0 = 1

        while resol <= resolution * size:
            file.write(
                f"BOUNDSPEC SEGMENT IJ {lon_0} 1 {resol} 1 VARIABLE FILE 0 'boundary_wave_" + name + "_" + str(
                    bound_count) + ".txt \n\n")
            lon_0 += resolution
            resol += resolution
            bound_count += 1

    elif location == 'north':

        file = open('boundary_namelist_wave_N.txt', 'w')
        file.write('&&----wave_boundary----&& \n\n')
        file.write('BOUND SHAPESPEC JONSWAP MEAN DSPR DEGREES \n')

        name = 'N'
        resol = resolution
        bound_count = 0
        lon_0 = 1

        while resol <= resolution * size:
            file.write(
                f"BOUNDSPEC SEGMENT IJ {lon_0} {grid_size_y} {resol} {grid_size_y} VARIABLE FILE 0 '"
                f"boundary_wave_" + name + "_" + str(
                    bound_count) + ".txt \n\n")
            lon_0 += resolution
            resol += resolution
            bound_count += 1

    file.close()


def bathy2swan_reg(filename, gebco, lon_west, lon_east, lat_north, lat_south):
    """
    Transform GEBCO DATA to SWAN bottom input

    PARAMETERS:
    -----------

    filename (str): The name that the file will be save.

    gebco (netcdf): The GEBCO Database

    lon_west (float): Longitude point in WEST

    lon_east (float): Longitude point in EAST

    lat_north (float): Latitude point in NORTH

    lat_south (float): Latitude point in SOUTH

    RETURN:
    -------

    Bathymetric file and the LAT and LON origin.
    The GEBCO file begging from bottom left.

    """

    import xarray as xr
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    data = xr.open_dataset(gebco)

    # GRID size

    bathy = data.sel(lat=slice(lat_north, lat_south),
                     lon=slice(lon_west, lon_east))

    df_bathy = pd.DataFrame(bathy.elevation)
    df_bathy.to_csv(filename + '.bot', sep=' ', header=False, index=False, float_format='%7.3f')

    lat_origin = bathy.lat[0].values
    lon_origin = bathy.lon[0].values

    lon_len = np.shape(df_bathy)[1]
    lat_len = np.shape(df_bathy)[0]

    print(f"The GEBCO file begging from TOP left: Latitude_0 = {lat_origin} e Longitude_0 = {lon_origin} ")

    # WRITE NAME LIST BATHY

    write_boundary = open('bathy_namelist_input', 'w')
    write_boundary.write('&& BATHY INPUT &&\n\n')
    write_boundary.write(
        f"INPGRID  BOTTOM REGular {lon_origin} {lat_origin} 0 {lon_len - 1} {lat_len - 1} {0.00416667} {0.00416667} \n")
    write_boundary.write(f"READINP BOTTOM -1 '{filename + '.bot'}' 2 0 FREE \n\n")
    write_boundary.close()

    plt.contourf(bathy.lon, bathy.lat, bathy.elevation, cmap='RdBu', vmin=-1500, vmax=0)
    plt.plot(bathy.lon[0], bathy.lat[0], marker='o', color='r')
    plt.colorbar()


def bathy2swan_curv(gebco, filename, vert_S, vert_E):

    x = vert_S[0]
    y = vert_E[1]

    # Grid resolution
    reso_x = np.round(np.abs(x[0] - x[1]), decimals=3)
    reso_y = np.round(np.abs(x[0] - x[1]), decimals=3)

    # Adapting to grid to SWAN
    vert_east = np.append(y, y[-1] + reso_x)
    vert_south = np.append(x, x[-1] + reso_y)

    # Generating grid
    xi, yi = np.meshgrid(vert_south, vert_east)

    gebco_ = xr.open_dataset(gebco)
    gebco = gebco_.sel(lon=slice(xi.min(), xi.max()), lat=slice(yi.max(), yi.min()))

    gebco_swan = gebco.sel(lon=xi[0], lat=yi[:, 0], method='nearest')
    bot = gebco_swan.elevation.values
    df_bot = pd.DataFrame(bot)
    df_bot.to_csv(filename, index=False, header=False, sep=' ')

    lon_len = np.shape(xi)[1]
    lat_len = np.shape(yi)[0]

    # WRITE NAME LIST BATHY

    write_boundary = open('bathy_namelist_input', 'w')
    write_boundary.write('&& BATHY INPUT &&\n\n')
    write_boundary.write(f"INPGRID  BOTTOM CURV 0 0 {lon_len - 1} {lat_len - 1} EXC 9.999000e+003 \n")
    write_boundary.write(f"READINP BOTTOM -1 '{filename + '.bot'}' 4 0 FREE \n\n")
    write_boundary.close()

    plt.contourf(gebco_swan.lon, gebco_swan.lat, gebco_swan.elevation, cmap='RdBu', vmin=-1500, vmax=0)
    plt.plot(gebco_swan.lon[0], gebco_swan.lat[0], marker='o', color='r')
    plt.colorbar()

    return df_bot

# SPEC 2D FROM WAVE PARTITION
# USED FOR BOUNDARY CONDITION


def spec2d(Hs, Tp, pdir, spread=30, units='deg', normalize=True):

    import numpy as np



    """ This function returns the 2d-spectrum from Hs, Tp, pdir and spread

    Edite by Nícolas and Marília to be used in the master project 

    -----------------------------------------------------

    REFERENCE:
    - https://www.orcina.com/webhelp/OrcaFlex/Content/html/Waves,Wavespectra.htm
    - http://research.dnv.com/hci/ocean/bk/c/a28/s3.htm
    - https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/applications/waves/

    -----------------------------------------------------

    Parameters: 

    Hs = []; % Sig. higth 
    Tp = []; % Peak Period  
    pdir = [];  % Peak Direction  
    spread = []; % Spread in degrees. If no value is passed, default is used (30°)
    units = [deg] or [rad]; # Direction unit. Default is 'deg'


     """

    # Setup parameters:

    freq = np.array([0.0500, 0.0566, 0.0642, 0.0727, 0.0824, 0.0933, 0.1057, 0.1198,
                     0.1357, 0.1538, 0.1742, 0.1974, 0.2236, 0.2533, 0.2870, 0.3252,
                     0.3684, 0.4174, 0.4729, 0.5357, 0.6070, 0.6877, 0.7791, 0.8827,
                     1.0000])  # frequency bin in swan

    directions = np.arange(0., 360., 10)  # Directional bin swan

    if Tp == 0:
        spec_2d = np.zeros((len(freq), len(directions)))
    else:
        fp = 1 / Tp  # peak frequency

        gam = 3.3  # default value according to SWAN manual

        g = 9.81  # gravity

        if spread > 31.5:
            ms = 1
        elif 31.5 >= spread > 27.6:
            ms = 2
        elif 27.5 >= spread > 24.9:
            ms = 3
        elif 24.9 >= spread > 22.9:
            ms = 4
        elif 22.9 >= spread > 21.2:
            ms = 5
        elif 21.2 >= spread > 19.9:
            ms = 6
        elif 19.9 >= spread > 18.8:
            ms = 7
        elif 18.8 >= spread > 17.9:
            ms = 8
        elif 17.9 >= spread > 17.1:
            ms = 9
        elif 17.1 >= spread > 14.2:
            ms = 10
        elif 14.2 >= spread > 12.4:
            ms = 15
        elif 12.4 >= spread > 10.2:
            ms = 20
        elif 10.2 >= spread > 8.9:
            ms = 30
        elif 8.9 >= spread > 8:
            ms = 40
        elif 8 >= spread > 7.3:
            ms = 50
        elif 7.3 >= spread > 6.8:
            ms = 60
        elif 6.8 >= spread > 6.4:
            ms = 70
        elif 6.4 >= spread > 6.0:
            ms = 80
        elif 6.0 >= spread > 5.7:
            ms = 90
        elif 5.7 >= spread > 4:
            ms = 100
        elif 4 >= spread > 2.9:
            ms = 400
        elif 7.3 >= spread >= 2.0:
            ms = 800
        else:
            ms = 1000

        sigma = freq * 0
        sigma[freq < fp] = 0.07
        sigma[freq >= fp] = 0.09
        sigma = np.array(sigma)

        # Pierson-Moskowitz Spectrum

        alpha = 1 / (0.06533 * gam ** 0.8015 + 0.13467) / 16;  # Yamaguchi (1984), used in SWAN

        pm = alpha * Hs ** 2 * Tp ** -4 * freq ** -5 * np.exp(-1.25 * (Tp * freq) ** -4)

        # apply JONSWAP shape

        jon = pm * gam ** np.exp(-0.5 * (Tp * freq - 1) ** 2. / (sigma ** 2.))

        jon[np.isnan(jon)] = 0

        # Optionally correct total energy of user-discretized spectrum to match Hm0,
        # as above methods are only an approximation

        eps = np.finfo(float).eps

        if normalize is True:
            corr = Hs ** 2 / (16 * trapz_and_repeat(jon, freq))
            jon = jon * corr

        # Directional Spreading

        '''Generate wave spreading'''

        from math import gamma

        directions = np.asarray(directions, dtype=float)

        # convert units to radians
        if units.lower().startswith('deg'):
            directions = np.radians(directions)
            pdir = np.radians(pdir)
        elif units.lower().startswith('rad'):
            pass
        else:
            raise ValueError('Unknown units: %s')

        # compute directional spreading
        A1 = (2. ** ms) * (gamma(ms / 2 + 1)) ** 2. / (np.pi * gamma(ms + 1))
        cdir = A1 * np.maximum(0., np.cos(directions - pdir))
        # cdir = np.maximum(0., np.cos(directions - pdir))**s

        # convert to original units
        if units.lower().startswith('deg'):
            directions = np.degrees(directions)
            pdir = np.degrees(pdir)
            cdir = np.degrees(cdir)

        # normalize directional spreading
        if normalize:
            cdir /= trapz_and_repeat(cdir, directions - pdir, axis=-1)

        cdir = np.array([list(cdir)] * len(jon))

        jon_list = list(jon)

        jon = np.array([ele for ele in jon_list for i in range(len(directions))]
                       ).reshape(len(freq), len(directions))

        jon2 = jon * cdir
        jon2[np.isnan(jon2)] = 0
        spec_2d = jon2

    return spec_2d


def write_specs(time, coords, specs, output_filename, ext='.bnd'):

    import pandas as pd

    ''' 
    WHAT: 

    Writes Swan spectral file from spec numpy array

    INPUT:
    time: time vector with same length as second dimension os specs
    coords: list of lists of [[lon1, lat1], [lon2, lat2], ...]. same 
    length as first dimension of specs
    specs: 4d numpy array with dims [location, time, freq, dir]
    output_filename: str to be inside the filename 
    ext: str with the extension of the file. default option is '.bnd'
    '''

    datetime = pd.to_datetime(time)

    for loc, coord in enumerate(coords.bound_wave):
        fin = open('specs_' + output_filename + str(loc) + str(ext), 'w')
        fin.write(
            '''SWAN   1                                Swan standard spectral file, version
$   Data produced by SWAN version 40.51AB             
$   Project:                 ;  run number:     
TIME                                    time-dependent data
     1                                  time coding option
LONLAT                                  locations in spherical coordinates
1
     {coord[0]:.6f}      {coord[1]:.6f}
AFREQ                                   absolute frequencies in Hz
    25
         0.0500
         0.0566
         0.0642
         0.0727
         0.0824
         0.0933
         0.1057
         0.1198
         0.1357
         0.1538
         0.1742
         0.1974
         0.2236
         0.2533
         0.2870
         0.3252
         0.3684
         0.4174
         0.4729
         0.5357
         0.6070
         0.6877
         0.7791
         0.8827
         1.0000
NDIR                                   spectral nautical directions in degr
36
       0
       10
       20
       30
       40
       50
       60
       70
       80
       90
       100
       110
       120
       130
       140
       150
       160
       170
       180
       190
       200
       210
       220
       230
       240
       250
       260
       270
       280
       290
       300
       310
       320
       330
       340
       350
QUANT
     1                                  number of quantities in table
VaDens                                  variance densities in m2/Hz/degr
m2/Hz/degr                              unit
   -0.9999                              exception value
            ''')
        fin.close()
        fin = open('specs_' + output_filename + str(loc) + str(ext), "a+")
        for line in range(len(datetime)):
            fin.write(
   '''{:%Y%m%d.%H%M}
FACTOR
1
{:>10}'''.format(
   	datetime[line], 
   	pd.DataFrame(specs[loc][line][:][:]).to_csv(index=False, header=False,
                                                                sep=',',
                                                                float_format='%7.5f',
                                                                na_rep=-0.9999,
                                                                lineterminator='\n')))
        
        fin = open('specs_' + output_filename + str(loc) + str(ext), "rt")
        data = fin.read()
        data = data.replace(',', ' ')
        fin.close()

        fin = open('specs_' + output_filename + str(loc) + str(ext), "wt")
        fin.write(data)
        fin.close()


def spec2d_era5(bound, hs_wave, tp_wave, pdir_wave, spread_wave):
    size = bound.bound_wave.shape[0]

    dfs = []

    for i in range(size):
        Hs = bound[hs_wave][i].values
        spread = np.rad2deg(bound[spread_wave][i].values)
        pdir = bound[pdir_wave][i].values
        Tp = bound[tp_wave][i].values
        time = bound.time

        df = pd.DataFrame({
            'time': time,
            'hs': Hs,
            'tp': Tp,
            'pdir': pdir,
            'spread': spread,
        })

        dfs.append(df.fillna(0))

    ##
    specs = []

    for df in dfs:
        loc = []
        aa = []

        for row in df.itertuples():
            try:
                spec = spec2d(row.hs, row.tp, row.pdir, row.spread)
                loc.append(spec)

            except:
                loc.append(spec)

        specs.append(loc)

    specs = np.array(specs)
    print(specs.shape)

    return specs


def spec_bound(bound):
    sea = spec2d_era5(bound, hs_wave='shww', tp_wave='mpww', pdir_wave='mdww', spread_wave='dwww')
    swell_1 = spec2d_era5(bound, hs_wave='p140121', tp_wave='p140123', pdir_wave='p140122', spread_wave='dwps')
    swell_2 = spec2d_era5(bound, hs_wave='p140124', tp_wave='p140126', pdir_wave='p140125', spread_wave='dwps')
    swell_3 = spec2d_era5(bound, hs_wave='p140127', tp_wave='p140129', pdir_wave='p140128', spread_wave='dwps')

    wave_cond = sea + swell_1 + swell_2 + swell_3

    return wave_cond

# Atmospheric and ocean boundary condition for SWAN

# WIND

def wind_era52swan(era5_database, name, time_init, time_end, lon_e, lon_w, lat_s, lat_n):
    """
    Generate SWAN wind input

    Parameters:
    -----------

    era5_data (dataset): ERA5 database

    name (srt): filename output

    time_init (datetime): Initial time (e.g. format '2017-09-01T00:00:00.000000000')

    time_end (datetime): Initial time (e.g. format '2017-09-01T00:00:00.000000000' )

    Returns:
    --------

    SWAN wind input


    """

    era5_database = xr.open_dataset(era5_database)
    era5_ = era5_database.sel(time=slice(time_init, time_end))
    era5_data_ = era5_.sel(longitude=slice(lon_w, lon_e), latitude=slice(lat_n, lat_s))
    era5_data = era5_data_ #.sel(expver=1)

    # GRID WIND

    u10 = era5_data.u10.values
    v10 = era5_data.v10.values

    # TIME ARRAY

    time_0 = era5_data.time.values[0]
    time_end = era5_data.time.values[-1]
    time = pd.date_range(time_0, time_end, freq='1H')
    time_model = time.format(formatter=lambda x: x.strftime('%Y%m%d.%H%M%S'))

    wind_filename_output = 'wind_era5_' + name + '.wnd'
    file = open(wind_filename_output, 'w')

    for t in range(len(time_model)):
        file.write(time_model[t])
        file.write('''
    ''')
        file.close()
        file = open(wind_filename_output, 'a+')
        file.write(pd.DataFrame(u10[t]).to_csv(index=False, header=False, na_rep=0, float_format='%7.3f'))
        file.write(pd.DataFrame(v10[t]).to_csv(index=False, header=False, na_rep=0, float_format='%7.3f'))

    file = open(wind_filename_output, "rt")
    data = file.read()
    data = data.replace(',', ' ')
    file.close()
    file = open(wind_filename_output, "wt")
    file.write(data)
    file.close()

    # WRITE INPUT WIND FORMAT TO SWAN

    lat_0 = era5_data.latitude[-1].values
    lon_0 = era5_data.longitude[0].values
    lon_len = len(era5_data.longitude)
    lat_len = len(era5_data.latitude)
    lon_ele_size = (era5_data.longitude[1] - era5_data.longitude[0]).values
    lat_ele_size = (era5_data.latitude[1] - era5_data.latitude[0]).values

    write_boundary = open('wind_namelist_input', 'w')
    write_boundary.write('&& WIND INPUT &&\n\n')
    write_boundary.write(
        f"INPGRID WIND REGular {lon_0} {lat_0} 0 {lon_len - 1} {lat_len - 1} {abs(lon_ele_size)} {abs(lat_ele_size)} NONSTAT {time_model[0]} 1 HR {time_model[-1]}\n")
    write_boundary.write(f"READINP WIND 1 '{wind_filename_output}' 4 0 FREE \n\n")
    write_boundary.close()

    # PLOT WIND SPEED MAP

    wspd = np.sqrt(u10 ** 2 + v10 ** 2)

    plt.contourf(era5_data.longitude, era5_data.latitude, wspd[0])
    plt.plot(era5_data.longitude[0], era5_data.latitude[0], marker='o')
    plt.colorbar()

    return era5_data

# CURRENT

def current_cmds2swan(cmds_database, name, time_init, time_end, lon_e, lon_w, lat_s, lat_n):


    """
    Generate SWAN current input

    Parameters:
    -----------

    cmds_data (dataset): Copernicus Marine Data Store

    name (srt): filename output

    time_init (datetime): Initial time (e.g. format '2017-09-01T00:00:00.000000000')

    time_end (datetime): Initial time (e.g. format '2017-09-01T00:00:00.000000000' )

    Returns:
    --------

    SWAN current input


    """

    # LOAD DATA
    cmds_database = xr.open_dataset(cmds_database)
    cmds_ = cmds_database.sel(time=slice(time_init, time_end))
    cmds_sl = cmds_.sel(longitude=slice(lon_w, lon_e), latitude=slice(lat_s, lat_n))
    cmds_data = cmds_sl.isel(depth=0)

    # CURRENT GRID

    uo = cmds_data.uo.values
    vo = cmds_data.vo.values

    # TIME ARRAY

    time_0 = cmds_data.time.values[0]
    time_end = cmds_data.time.values[-1]
    time = pd.date_range(time_0, time_end, freq='1D')
    time_model = time.format(formatter=lambda x: x.strftime('%Y%m%d.%H%M%S'))

    current_filename_output = 'current_cmds_' + name + '.dat'
    file = open(current_filename_output, 'w')

    for t in range(len(time_model)):
        file.write(time_model[t])
        file.write('''
    ''')
        file.close()
        file = open(current_filename_output, 'a+')
        file.write(pd.DataFrame(uo[t]).to_csv(index=False, header=False, na_rep=0, float_format='%7.3f'))
        file.write(pd.DataFrame(vo[t]).to_csv(index=False, header=False, na_rep=0, float_format='%7.3f'))

    file = open(current_filename_output, "rt")
    data = file.read()
    data = data.replace(',', ' ')
    file.close()
    file = open(current_filename_output, "wt")
    file.write(data)
    file.close()

    # WRITE CURRENT INPUT FORMAT TO SWAN

    lat_0 = cmds_data.latitude[0].values
    lon_0 = cmds_data.longitude[0].values
    lon_len = len(cmds_data.longitude)
    lat_len = len(cmds_data.latitude)
    lon_ele_size = (cmds_data.longitude[1] - cmds_data.longitude[0]).values
    lat_ele_size = (cmds_data.latitude[1] - cmds_data.latitude[0]).values

    write_boundary = open('current_namelist_input', 'w')
    write_boundary.write('&& CURRENT INPUT &&\n\n')
    write_boundary.write(
        f"INPGRID CURRENT REGular {lon_0} {lat_0} 0 {lon_len - 1} {lat_len - 1} {abs(lon_ele_size)} {abs(lat_ele_size)} NONSTAT {time_model[0]} 1 DAY {time_model[-1]}\n")
    write_boundary.write(f"READINP CURRENT 1 '{current_filename_output}' 4 0 FREE \n\n")
    write_boundary.close()

    # PLOT CURRENT SPEED MAP

    cspd = np.sqrt(uo ** 2 + vo ** 2)

    plt.contourf(cmds_data.longitude, cmds_data.latitude, cspd[0])
    plt.colorbar()

    return cmds_data

# SEA LEVEL SURFACE

def ssh_cmds2swan(cmds_database, name, time_init, time_end, lon_e, lon_w, lat_s, lat_n):

    """
    Generate SWAN current input

    Parameters:
    -----------

    cmds_data (dataset): Copernicus Marine Data Store

    name (srt): filename output

    time_init (datetime): Initial time (e.g. format '2017-09-01T00:00:00.000000000')

    time_end (datetime): Initial time (e.g. format '2017-09-01T00:00:00.000000000' )

    Returns:
    --------

    SWAN current input


    """

    # LOAD DATA

    cmds_database = xr.open_dataset(cmds_database)
    cmds_ = cmds_database.sel(time=slice(time_init, time_end))
    cmds_sl = cmds_.sel(longitude=slice(lon_w, lon_e), latitude=slice(lat_s, lat_n))
    cmds_data = cmds_sl.isel(depth=0)

    # Sea surface high GRID

    ssh = cmds_data.zos.values

    # TIME ARRAY

    time_0 = cmds_data.time.values[0]
    time_end = cmds_data.time.values[-1]
    time = pd.date_range(time_0, time_end, freq='1D')
    time_model = time.format(formatter=lambda x: x.strftime('%Y%m%d.%H%M%S'))

    ssh_filename_output = 'ssh_cmds_' + name + '.dat'
    file = open(ssh_filename_output, 'w')

    for t in range(len(time_model)):
        file.write(time_model[t])
        file.write('''
    ''')
        file.close()
        file = open(ssh_filename_output, 'a+')
        file.write(pd.DataFrame(ssh[t]).to_csv(index=False, header=False, na_rep=0, float_format='%7.3f'))

    file = open(ssh_filename_output, "rt")
    data = file.read()
    data = data.replace(',',' ')
    file.close()
    file = open(ssh_filename_output, "wt")
    file.write(data)
    file.close()

    # WRITE CURRENT INPUT FORMAT TO SWAN

    lat_0 = cmds_data.latitude[0].values
    lon_0 = cmds_data.longitude[0].values
    lon_len = len(cmds_data.longitude)
    lat_len = len(cmds_data.latitude)
    lon_ele_size = (cmds_data.longitude[1] - cmds_data.longitude[0]).values
    lat_ele_size = (cmds_data.latitude[1] - cmds_data.latitude[0]).values

    write_boundary = open('ssh_namelist_input', 'w')
    write_boundary.write('&& SSH INPUT &&\n\n')
    write_boundary.write(
        f"INPGRID WLEV REGular {lon_0} {lat_0} 0 {lon_len - 1} {lat_len - 1} {abs(lon_ele_size)} {abs(lat_ele_size)} NONSTAT {time_model[0]} 1 DAY {time_model[-1]}\n")
    write_boundary.write(f"READINP WLEV 1 '{ssh_filename_output}' 4 0 FREE \n\n")
    write_boundary.close()

    ## PLOT SSH SPEED MAP

    plt.contourf(cmds_data.longitude, cmds_data.latitude, cmds_data.zos[0])
    plt.colorbar()

    return cmds_data
    
  
  
def iterable(arr):
    '''Returns an iterable'''
    
    try:
        iter(arr)
        return arr
    except:
        return (arr,)
        
        
        
def expand_and_repeat(mtx, shape=None, repeat=None,
                      exist_dims=None, expand_dims=None):
    '''Expands matrix and repeats matrix contents along new dimensions

    Provide ``shape`` and ``exist_dims`` or ``expand_dims``, or
    ``repeat`` and ``expand_dims``.

    Parameters
    ----------
    mtx : numpy.ndarray
      Input matrix
    shape : tuple, optional
      Target shape of output matrix
    repeat : tuple or int, optional
      Repititions along new dimensions
    exist_dims : tuple or int, optional
      Indices of dimensions in target shape that are present in input matrix
    expand_dims : tuple or int, optional
      Indices of dimensions in target shape that are not present in input matrix

    Returns
    -------
    numpy.ndarray
      Matrix with target shape

    Examples
    --------
    >>> expand_and_repeat([[1,2,3],[4,5,6]], shape=(2,3,4), exist_dims=(0,1))
    >>> expand_and_repeat([[1,2,3],[4,5,6]], shape=(2,3,4), expand_dims=(2,))
    >>> expand_and_repeat([[1,2,3],[4,5,6]], shape=(2,3,4), expand_dims=2)
    >>> expand_and_repeat([[1,2,3],[4,5,6]], repeat=(4,), expand_dims=(2,))
    >>> expand_and_repeat([[1,2,3],[4,5,6]], repeat=4, expand_dims=2)

    '''
    
    mtx = np.asarray(mtx)
    
    if shape is not None:
        shape = iterable(shape)
        
        if mtx.ndim > len(shape):
            raise ValueError('Nothing to expand. Number of matrix '
                             'dimensions (%d) is larger than the '
                             'dimensionality of the target shape '
                             '(%d).' % (mtx.ndim, len(shape)))
        
        if exist_dims is not None:
            exist_dims = iterable(exist_dims)

            if len(exist_dims) != len(set(exist_dims)):
                raise ValueError('Existing dimensions should be unique.')
            
            if mtx.ndim != len(exist_dims):
                raise ValueError('Number of matrix dimensions (%d) '
                                 'should match the number of existing '
                                 'dimensions (%d).' % (mtx.ndim, len(exist_dims)))

            expand_dims = [i
                           for i in range(len(shape))
                           if i not in exist_dims]
                             
        elif expand_dims is not None:
            expand_dims = iterable(expand_dims)
            
            if len(expand_dims) != len(set(expand_dims)):
                raise ValueError('Expanding dimensions should be unique.')
            
            if len(shape) - mtx.ndim != len(expand_dims):
                raise ValueError('Dimensionality of the target shape '
                                 'minus the number of matrix dimensions '
                                 '(%d) should match the number of expanding '
                                 'dimensions (%d).' % (len(shape) - mtx.ndim, len(expand_dims)))
            
            exist_dims = [i
                          for i in range(len(shape))
                          if i not in expand_dims]
            
        else:
            raise ValueError('Target shape undetermined. Provide '
                             '``exist_dims`` or ``expand_dims``.')

        repeat = [n
                  for i, n in enumerate(shape)
                  if i in expand_dims]

        for i1, i2 in enumerate(exist_dims):
            if shape[i2] != mtx.shape[i1]:
                raise ValueError('Current matrix dimension (%d = %d) '
                                 'should match target shape (%d = %d).' % (i1, mtx.shape[i1], i2, shape[i2]))

    elif repeat is not None and expand_dims is not None:
        repeat = iterable(repeat)
        expand_dims = iterable(expand_dims)

        if len(expand_dims) != len(set(expand_dims)):
            raise ValueError('Expanding dimensions should be unique.')

        if len(repeat) != len(expand_dims):
            raise ValueError('Number of repititions (%d) should '
                             'match the number of expanding '
                             'dimensions (%d).' % (len(repeat), len(expand_dims)))
            
    else:
        raise ValueError('Target shape undetermined. Provide '
                         '``shape`` and ``exist_dims`` or '
                         '``expand_dims``, or ``repeat`` and ``expand_dims``.')

    for i, n in zip(expand_dims, repeat):
        mtx = np.expand_dims(mtx, i).repeat(n, axis=i)

    return mtx



def trapz_and_repeat(mtx, x, axis=-1):

    if axis < 0:
        axis += len(mtx.shape)

    return expand_and_repeat(np.trapz(mtx, x, axis=axis),
                             shape=mtx.shape, expand_dims=axis)
