"""
Routines to extract requested parameters in CCHDO's netcdf bottle parameters
"""

# Set some module level variables that can either be manually set or will be set during other routines
from datetime import datetime, timedelta
verbose = False
really_verbose = False
datapath_glob = ''
exponames_glob = []
# This is WOCE date offset (as specified in the standard) in python's datetime type
woce_ref_time = datetime(1980,1,1)

import pandas

def set_exponames( exponames ):
    global exponames_glob
    exponames_glob = exponames

def set_verbosity( verbosity = 0 ):
    global verbose
    global really_verbose
    if verbosity == 0:
        verbose = False
        really_verbose = False
    elif verbosity == 1:
        verbose = True
        really_verbose = False
    elif verbosity >  1:
        verbose = True
        really_verbose = True

def update_expo_list( datapath='', outfile = "" ):
    """
    Returns all unique expocodes and optionally write them to a file for later reference
    Inputs:
        datapath: Root path where CCHDO files are stored
        outfile:  The name of the file to save the list of expocodes
    """
    from glob import glob
    from netCDF4 import Dataset
    from os.path import join

    global datapath_glob
    global exponames_glob

    if len(datapath)>0:
        datapath_glob = datapath
    else:
        datapath = datapath_glob

    filelist = glob(join(datapath,'*.nc'))
    # Retrieve a list of all the netcdf files in the given directory
    # Loop through every file to find all possible expocodes
    exponames = set()
    for file in filelist:
        data = Dataset(file)
        exponames.add(getattr(data,'EXPOCODE'))
        data.close()

    # Write the exponames as a text file if requested
    if len(outfile) > 0:
        wfile = open(outfile,'w')
        [ wfile.write("%s\n" % line) for line in exponames ]
        wfile.close()

    if verbose:
        print("Read %d expocodes from %s" % (len(exponames),datapath))

    exponames_glob = list(exponames)

    return exponames
def read_expo_list( readfile ):
    """
    Read the list of expocodes in readfile, return it, and set exponames_glob
    Inputs:
        readfile: Name of the file containing a list of expocodes
        datapath: Root path where CCHDO files are stored
    """
    from os.path import join
    global exponames_glob

    expofile = open(readfile,'r')
    exponames = expofile.readlines()
    exponames = [ string.rstrip('\n') for string in exponames ]
    exponames = list(set(exponames))
    exponames_glob = exponames

    return exponames

def extract_expo_fields( exponame, qc_flags = [2], fields_in = set(), fields_0_in = set(), datapath = datapath_glob, aux_fields = True ):
    """
    Extract fields from the bottle data associated with the given exponame
    Inputs:
        exponame: The name of the expedition to read data from
        qc_flags: A list of 'acceptable' water sample QC codes based on WOCE standard
        fields:   A list of the fields to extract from the bottle data, by default
                  pressure, temperature, and salinity are extracted
        fields_0: Any additional scalar values in each data file, by default
                  time, latitude, and longitude
        datapath: Root path where CCHDO files are stored
        aux_fields: If true, add common auxiliary fields, in-situ density, sigma0, sigma2,
                    conservative temperature, surface potential temperature, and absolute
                    salinity
    Output:
        expodata: A dictionary with all the extracted fields
    """
    from os.path import join
    from glob import glob
    from netCDF4 import Dataset, chartostring
    import numpy as np
    import gsw

    global datapath_glob
    global exponames_glob

    if len(datapath)>0:
        datapath_glob = datapath
    else:
        datapath = datapath_glob

    # Set some default fields to extract and then append the input ones
    fields = set(['pressure','temperature','bottle_salinity'])
    [ fields.add(field) for field in fields_in ]
    fields_0 = set(['time','latitude','longitude'])
    [ fields_0.add(field) for field in fields_0_in ]
    # Convert fields and field_0 back into list so we can iterate over them
    fields = list(fields)
    fields_0 = list(fields_0)

    filelist = glob(join(datapath,'*%s*.nc' % exponame))
    if verbose:
        print("Expocode %s has %d files" % (exponame, len(filelist) ) )

    # Initialize the dictionary used to store the extracted fields
    botdata = {}
    stadata = {}
    numskip = {}
    for field in fields:
        botdata[field] = np.array([])
        numskip[field] = 0
    for field in fields_0:
        botdata[field] = np.array([])
        stadata[field] = np.array([])

    if len(filelist) == 0:
        print("WARNING: No input files associated with %s" % exponame)
        return botdata

    for file in filelist:
        # Read in the scalar fields and replicate it so that every bottle measurement
        # has a corresponding lat, lon, time
        expo_ncvars = Dataset(join(datapath,file)).variables
        nk = expo_ncvars['pressure'][:].size
        if really_verbose:
            print("File %s has %d measurements" % (file, nk))
        for field in fields_0:
            stadata[field] = np.append(stadata[field],expo_ncvars[field][:])
            botdata[field] = np.append(botdata[field],expo_ncvars[field][:]*np.ones(nk))
        # Read in the bottle fields and also the QC field if it exists
        for field in fields:
            # Try to read the field in the open netCDF
            try:
                # Might need to try to find the CTD salinity instead of bottle
                if field == 'bottle_salinity':
                    try:
                        tempvar = expo_ncvars[field][:]
                    except:
                        tempvar = expo_ncvars['salinity'][:]
                else:
                    tempvar = expo_ncvars[field][:]
                # Attempt to read the QC code associated with the field. If it exists, apply the
                # requested quality control flags
                try:
                    varflags = expo_ncvars[field+'_QC'][:]
                    # Loop over all requested qc flags
                    flagvar = np.zeros( varflags.shape, dtype = bool )
                    for flag in qc_flags:
                        flagvar[ varflags == flag ] = True
                    # Apply the flags
                    tempvar[ np.logical_not(flagvar) ] = np.nan
                except:
                    if really_verbose and (field != 'temperature' and field != 'pressure'):
                        print("QC for field %s does not exist" % field)
                botdata[field] = np.append(botdata[field],tempvar)
            # If the variable can't be found in the field, append NaNs to main consistency in array shapes
            # Record that the a variable was skipped
            except:
                numskip[field] += 1
                botdata[field] = np.append(botdata[field],np.zeros(nk)*np.nan)
                if verbose:
                    print("Could not find %s in file %s" % (field,file) )
        # Closeout the file
        Dataset(join(datapath,file)).close()

    # Use the reference time to convert to python datetime format
    botdata['time'] = [ woce_ref_time + timedelta(minutes=time) for time in botdata['time'] ]
    stadata['time'] = [ woce_ref_time + timedelta(minutes=time) for time in stadata['time'] ]
    # Convert longitudes to 0 to 360
    botdata['longitude'] = np.mod(botdata['longitude'],360)
    stadata['longitude'] = np.mod(stadata['longitude'],360)

    # Print out how many files did not have the requested field in them
    if verbose:
        [print("Field %s could not be found in %d files" % (field, numskip[field])) for field in fields if numskip[field]>0 ]

    if aux_fields:
        botdata['abssal'] = gsw.conversions.SA_from_SP( botdata['bottle_salinity'], botdata['pressure'], botdata['longitude'], botdata['latitude'])
        botdata['ctemp' ] = gsw.conversions.CT_from_t(  botdata['abssal'],botdata['temperature'],botdata['pressure'])
        botdata['ptemp0'] = gsw.conversions.pt0_from_t( botdata['abssal'],botdata['temperature'],botdata['pressure'])
        botdata['sigma0'] = gsw.density.sigma0( botdata['abssal'],botdata['temperature'])
        botdata['sigma2'] = gsw.density.sigma2( botdata['abssal'],botdata['temperature'])
        botdata['rhoinsitu'] = gsw.density.rho( botdata['abssal'],botdata['temperature'],botdata['pressure'])


    botdata = pandas.DataFrame.from_dict(botdata)
    stadata = pandas.DataFrame.from_dict(stadata)
    # Calculate distance along transect
    botdata.sort_values('time',inplace=True)
    stadata.sort_values('time',inplace=True)
    botdata['distance'] = np.append(0,gsw.distance(np.array(botdata.longitude), np.array(botdata.latitude),p=np.zeros(botdata.latitude.shape)).cumsum())
    stadata['distance'] = np.append(0,gsw.distance(np.array(stadata.longitude), np.array(stadata.latitude),p=np.zeros(stadata.latitude.shape)).cumsum())

    expodata = {}
    expodata['station'] = stadata
    expodata['bottle'] =  botdata
    expodata['expocode'] = exponame
    return expodata

def extract_all_expos( qc_flags = [2], fields_in = set(), fields_0_in = set(), datapath = datapath_glob, aux_fields = True ):
    """
    Extracts the requested fields from every expocode listed in the datapath
        qc_flags: A list of 'acceptable' water sample QC codes based on WOCE standard
        fields:   A list of the fields to extract from the bottle data, by default
                  pressure, temperature, and salinity are extracted
        fields_0: Any additional scalar values in each data file, by default
                  time, latitude, and longitude
        datapath: Root path where CCHDO files are stored
        aux_fields: If true, add common auxiliary fields, in-situ density, sigma0, sigma2,
                    conservative temperature, surface potential temperature, and absolute
                    salinity
    Output:
        expodata: A two-level dictionary with all the cruises and all the extracted fields
    """

    global datapath_glob
    global exponames_glob

    if len(datapath)>0:
        datapath_glob = datapath
    else:
        datapath = datapath_glob

    # In the case where exponames has yet to be defined, get the expocode list
    if len(exponames_glob) == 0:
        exponames = update_expo_list(datapath)
        exponames_glob = exponames
    else:
        exponames = exponames_glob

    expodata = {}

    for exponame in exponames:
        expodata[exponame] = extract_expo_fields( exponame, qc_flags, fields_in, fields_0_in, datapath, aux_fields )

    return expodata

def plot_expo_transect( expodata, proj=None, figsize = (12,6) ):
    """
    Plot the transect of a given expedition
    Inputs:
        expodata: Dictionary containing the latitude and longitude fields
        proj:     Specify a particular map projection to use
        figuse:   Specify the size of a figure
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt

    if proj is None:
        proj = ccrs.PlateCarree(central_longitude=0.5*(20+380))
    fig = plt.figure(figsize = figsize)
    ax = plt.axes(projection=proj)
    ax.plot(expodata[station]['longitude'],expodata[station]['latitude'], 'bo', transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAND,facecolor='0.25')
    ax.set_global()
    plt.title(expodata['expocode'].upper())
    plt.show()

def calc_weights(xp, yp, x, y, xcorr, ycorr):
    """
    Calculate the weights for the Barnes objective mapping. Essentially a gaussian dropoff as a function of distance
    Note that this returns the unnormalized weights under the assumption that some points may not be valid for
    different fields
    Inputs:
        xp: one set of x-coordinates
        yp: one set of y-coordinates
        x : the other set of x-coordinates
        y:  the other set of y-coordinates
        xcorr: x correlation length scale in units of the x coordinates
        ycorr: y correlation length scale in units of the x coordinates
    """
    import itertools
    import numpy as np

    # If the correlation length scale(s) are provided, we can save some computational cost by scaling immediately
    inv_corr = 1./xcorr
    xp = xp * inv_corr
    x  = x  * inv_corr
    inv_corr = 1./ycorr
    yp = yp * inv_corr
    y  = y  * inv_corr

    # First calculate the x-component of the weights
    dp = np.fromiter([ x1 - x2 for (x1,x2) in itertools.product( xp, x ) ], np.float )
    dp = dp*dp
    w = -dp
    # Calculate the y component of the weights
    dp = np.array([ y1 - y2 for (y1,y2) in itertools.product( yp, y ) ])
    dp = dp*dp
    if xcorr is None:
        inv_corr = 1./(4*np.mean(dp)) # Note the factor of four is because dp has already been squared (2*dp)**2
        w = w - (dp*inv_corr)
    else:
        w = w - dp
    w = np.exp(w)
    return w.reshape(xp.size,-1)

def barnes(x, y, z, xi, yi, xcorr, ycorr, npass = 2, w = None, wp = None):
    """
    Maps scattered data onto the specified grid using a Barnes style interpolation
    """
    import numpy as np
    from multiprocessing import Pool

    z = np.ma.masked_invalid(z)
    if w is None or wp is None:
        x = np.array(x)
        y = np.array(y)

        xp = np.array(xi).reshape(-1)
        yp = np.array(yi).reshape(-1)

        # Weight matrix from output points to input points
        if wp is None:
            wp = calc_weights(xp,yp,x,y,xcorr,ycorr)
        if verbose:
            print("Size of weight matrix, output to input: %d %d" % wp.shape)

        # Weight matrix from input points to input poitns
        if w is None:
            w = calc_weights(x,y,x,y,xcorr,ycorr)
        if verbose:
            print("Size of weight matrix, input to input: %d %d" % w.shape)

    # Make w and wp a masked array
    w  = np.ma.MaskedArray(w, copy=True, mask = np.zeros(w.shape))
    wp = np.ma.MaskedArray(wp,copy=True, mask = np.zeros(wp.shape))
    # Mask nan points of input
    inmask = np.isnan(z)
    w.mask[:,inmask] = True
    w.mask[inmask,:] = True
    wp.mask[:,inmask] = True

    w  = w/w.sum(axis=-1)[:,np.newaxis]
    wp = wp/wp.sum(axis=-1)[:,np.newaxis]

    zpp = zp = np.sum(wp*z,axis=-1)
    zpb = np.zeros(z.shape)
    if verbose:
        print("Standard deviation of initial guess: %f" % np.std(zp))
    for iter in range(2,npass+1):
      if np.mod(iter,2)==0:
        zpa = zpb + np.sum((z-zpb)*w, axis=-1)
        zp  = zp  + np.sum((z-zpa)*wp,axis=-1)
        zpn = zp

      else:
        zpb = zpa + np.sum((z-zpa)*w, axis=-1)
        zp  = zp  + np.sum((z-zpb)*wp,axis=-1)
        zpp = zp
      if verbose:
          print('RMS adjustment after pass %d: %f' % (iter, np.std(zpn-zpp)))

    return zp

def grid_transect_variables( botdata, fields,  xgrid = None, ygrid = None, xvar = 'distance', yvar = 'pressure', nx = None, ny = None, xcorr=None, ycorr = 50, npass = 2 ):
    """
    Grid the scattered bottle data to the requested points using a Barnes objective mapping
    Inputs:
        botdata: Dictionary containing the data from the cruise
        fields:   A list of fields to be gridded
    Optional inputs:
        xgrid:    The x-value of points of the output grid
        ygrid:    The y-value of points of the output grid
        xvar:     The name of the x variable (default: longitude)
        yvar:     The name of the y variable (default: pressure)
        nx:       Number of points to grid along the x-axis (overrides xgrid)
        ny:       Number of points to grid along the y-axis (overrides ygrid)
        xcorr:    Correlation scale in units of x-variable along x
        ycorr:    Correlation scale in units of y-variable along y
        npass:    Number of passes to use in the objective mapping

    """
    import numpy as np
    import itertools
    from pandas import DataFrame
    from matplotlib.path import Path
    from scipy.interpolate import griddata
    if not (nx is None) and (ny is None):
        print("nx and ny must either both be none or valid integers")

    # Filter the input fields for valid points
    if (type(botdata[xvar])) is not type(np.array([])):
        x = np.array(botdata[xvar])
        y = np.array(botdata[yvar])
    notnan = np.logical_not(np.logical_or(np.isnan(x), np.isnan(y) ))
    x = x[notnan]
    y = y[notnan]

    # Initialize an empty dictionary containing all the remapped fields
    remap = {}
    # Make the grid if only a number of points is specified
    if (nx is not None) and (ny is not None):
        xgrid = np.linspace(x.min(), x.max(), nx)
        ygrid = np.linspace(y.min(), y.max(), ny)
    xgrid, ygrid = np.meshgrid(xgrid,ygrid)
    remap[xvar] = xgrid
    remap[yvar] = ygrid

    # Mask the output grid based on the outer edges of the input data
    # First loop to get the 'top' boundary (in pressure space this would be the surfacemost points
    bound_poly = []
    for time in sorted(botdata.time.unique()):
        bound_poly.append( (botdata.loc[botdata['time']==time][xvar].mean(),
                            botdata.loc[botdata['time']==time][yvar].min()) )
    # Next loop in reverse time to get the bottom boundary (e.g. biggest pressure)
    for time in sorted(botdata.time.unique(),reverse=True):
        bound_poly.append( (botdata.loc[botdata['time']==time][xvar].mean(),
                            botdata.loc[botdata['time']==time][yvar].max()) )
    # Make an array with the points
    remap['bound_poly'] = Path(bound_poly,closed=True)

    # Do only the points within the original data extent
    points = np.vstack( (xgrid.reshape(-1),ygrid.reshape(-1) ) ).T
    inpoints = remap['bound_poly'].contains_points(points)
    inpoints_grid = inpoints.reshape(xgrid.shape)

    xp = xgrid.reshape(-1)
    yp = ygrid.reshape(-1)

    xp = xp[inpoints]
    yp = yp[inpoints]
    # Need to calculate w and wp for the first field of interpolation
    # If the correlation length scales are not specified, calculate it based on twice the average distance between points
    # Note that the work array dp is set to zero to free RAM after use
    if xcorr is None:
        dp = np.fromiter([ x1 - x2 for (x1,x2) in itertools.product( xp, xp ) ], np.float )
        xcorr = 2.*np.mean(np.abs(dp))
        dp = 0.
    if ycorr is None:
        dp = np.fromiter([ x1 - x2 for (x1,x2) in itertools.product( yp, yp ) ], np.float )
        ycorr = 2.*np.mean(np.abs(dp))
        dp = 0.

    if verbose:
        print("Begin objective mapping of fields with %d passes, %f x length scale, %f y length scale" % (npass, xcorr, ycorr) )
    if verbose:
        print("Calculating weights from %d inputs and %d outputs" % (x.size, xp.size))
    wp = calc_weights(xp, yp, x, y, xcorr, ycorr)
    if verbose:
        print("Calculating internal weights from %d inputs" % x.size)
    w =  calc_weights(x, y, x, y, xcorr, ycorr)
    # Make wp and w masked arrays so that we can modify just the masks when there's missing bottle data
    wp = np.ma.MaskedArray(wp)
    w  = np.ma.MaskedArray(w )

    for zvar in fields:

        z = botdata[zvar][notnan]
        if verbose:
            print("Number of measurements for %s: %d" % (zvar,z.size))
        if zvar == 'pressure':
            # Pressure is treated as a special case because it should be relatively accurate and smooth already
            zopt = griddata((x,y),z,(xp,yp),rescale = True).flatten()
        else:
            zopt = barnes(x,y,z,xgrid,ygrid,xcorr,ycorr,npass, w, wp)
        # Initially say that the grid is all masked
        remap[zvar] = np.ma.masked_invalid(np.zeros(xgrid.shape)*np.nan)
        # Unmask the points in the valid dataset
        remap[zvar].mask[inpoints_grid] = False
        # Store the mapped variables
        remap[zvar][inpoints_grid] = zopt

    return remap

def remap_fields( remap, fields, xvar, yvar1, yvar2, ny = None, ynew = None ):
    """
    Remap fields within the remap dictionary or pandas dataframe from one coordinate to another. For example, regridding from density to pressure
    space.
    Input:
        remap:  Dictionary or dataframe from containing relevant data
        fields: Which fields to remap
        xvar:   Name of the x-coordinate to use
        yvar1:  The name of the original coordinate (e.g. sigma2)
        yvar2:  The name of the y-coordinate (e.g. pressure'), note that yvar2 must have also been mapped previously
        ny:     If this is specified, will attempt to regrid onto a linearly spaced grid between min/max of yvar2
        ynew:   The new grid discretized along yvar2 (probably the right way to use this routine)
    """
    import numpy as np
    from scipy.interpolate import griddata
    from matplotlib.path import Path


    if ynew is None and ny is None:
        print("ynew or ny must be specified")
        return None

    if ynew is None:
        ynew = np.linspace(remap[yvar2].min(), remap[yvar2].max(), ny)

    xin = remap[xvar].flatten()
    yin = remap[yvar2].flatten()

    xnew = np.unique(xin)
    xint, yint = np.meshgrid(xnew,ynew)

    # Setup the output dictionary
    regrid = {}
    regrid[xvar] = xint
    regrid[yvar2] = yint

    # Make a polygon defining the input extent of the data, note that we're assuming a regular grid
    # First traverse along the x-direction finding the smallest valid value
    nx, ny = remap[yvar2].shape
    bound_poly = []
    for i in range(0,nx):
        bound_poly.append( (remap[xvar][i,0], remap[yvar2][i,:].min() ) )
    # Do the same, but in reverse to get the largest valid value
    for i in reversed(range(0,nx)):
        bound_poly.append( (remap[xvar][i,0], remap[yvar2][i,:].max() ) )
    # Make an array with the points
    regrid['bound_poly'] = Path(bound_poly,closed=True)

    # Make mask based on the valid input points
    points = np.vstack( (xint.flatten(),yint.flatten() ) ).T
    inpoints = regrid['bound_poly'].contains_points(points)
    inpoints = inpoints.reshape(xint.shape)

    for field in fields:
        zin = remap[field].flatten()
        notnan = np.logical_not( np.logical_or( np.logical_or( np.isnan(xin), np.isnan(yin) ), np.isnan(zin) ) )
        # Do that actual regridding
        zout = griddata( (xin[notnan], yin[notnan]), zin[notnan], (xint.flatten(), yint.flatten()), rescale = True)
        zout = zout.reshape( xint.shape )
        # Now store the results of the regridding and mask appropraitely
        regrid[field] = np.ma.masked_where( inpoints, zout )

    return regrid

def map_gridded_fields_to_expo(stadata, griddeddata, gridfields, tidx = None, grid_xvar = 'nav_lon', grid_yvar = 'nav_lat', grid_zvar = 'deptht'):
    """
    Interpolate gridded data to the transect latitude and longitude. Data is first interpolated horizontally and then interpolated vertically
    Inputs:
        stadata:        dictionary or pandas dataframe containing station information
        griddeddata:    dictionary describing gridded data
        gridfields:     name of the fields in the gridded data to be mapped onto the transect
    Inputs: (optional)
        tidx:           use a specific time index or (default) average over the time axis
        grid_xvar:      name of the x (longitude) variable in the griddeddata
        grid_yvar:      name of the y (latitude) variable in the griddeddata
        grid_zvar:      name of the z (depth) variable in the griddeddata
    """

    from scipy.interpolate import griddata
    import numpy as np

    # How much of the model data to include outside of the transect (in degrees lat/lon)
    lonwindow = 2
    latwindow = 2
    # Make some shorthand variables
    gridx = np.mod(griddeddata[grid_xvar],360)
    gridy = griddeddata[grid_yvar]
    gridz = griddeddata[grid_zvar]

    # Mask the data in space to make the gridding faster
    mask = np.ones(gridx.shape, dtype = bool)
    mask = np.logical_and( mask, gridx > (stadata['longitude'].min() - lonwindow) )
    mask = np.logical_and( mask, gridx < (stadata['longitude'].max() + lonwindow) )
    mask = np.logical_and( mask, gridy > (stadata['latitude'].min() - latwindow) )
    mask = np.logical_and( mask, gridy < (stadata['latitude'].max() + latwindow) )

    # Initialize storage arrays
    nk = gridz.size
    nsta = len(stadata.index)
    remapped = {}
    remapped['latitude']  = np.zeros((nk,nsta))
    remapped['longitude'] = np.zeros((nk,nsta))
    remapped['distance']  = np.zeros((nk,nsta))
    remapped['pressure']  = np.zeros((nk,nsta))

    for field in gridfields:
        # Use either a single timestep if tidx is provided, otherwise assume that we want to average over the time dimension
        if tidx is None:
            data = griddeddata[field][:,:,:,:].mean(axis = 0)
        else:
            data = griddeddata[field][tidx,:,:,:]
        remapped[field] = np.zeros((nk,nsta))
        if verbose:
            print("Gridding %s with shape" % field, data.shape)
        for k in range(0,nk):
            subdata = data[k,:,:]
            remapped[field][k,:] = griddata((gridx[mask],gridy[mask]),subdata[mask],(stadata['longitude'],stadata['latitude']))
            remapped['latitude'][k,:] = stadata['latitude']
            remapped['longitude'][k,:] = stadata['longitude']
            remapped['distance'][k,:] = stadata['distance']
            remapped['pressure'][k,:] = gridz[k]

    return remapped


