"""
Routines to extract requested parameters in CCHDO's netcdf bottle parameters
"""

# Set some module level variables that can either be manually set or will be set during other routines
from datetime import datetime, timedelta
verbose = True
really_verbose = False
datapath_glob = ''
exponames_glob = []
# This is WOCE date offset (as specified in the standard) in python's datetime type
woce_ref_time = datetime(1980,1,1)

import pandas

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
    expodata = {}
    numskip = {}
    for field in fields:
        expodata[field] = np.array([])
        numskip[field] = 0
    for field in fields_0:
        expodata[field] = np.array([])

    if len(filelist) == 0:
        print("WARNING: No input files associated with %s" % exponame)
        return expodata

    for file in filelist:
        # Read in the scalar fields and replicate it so that every bottle measurement
        # has a corresponding lat, lon, time
        expo_ncvars = Dataset(join(datapath,file)).variables
        nk = expo_ncvars['pressure'][:].size
        if really_verbose:
            print("File %s has %d measurements" % (file, nk))
        for field in fields_0:
            if field == 'station':
                stid = int(np.asscalar(chartostring(expo_ncvars[field][:])))
                expodata[field] = np.append(expodata[field],stid*np.ones(nk))
            else:
                expodata[field] = np.append(expodata[field],expo_ncvars[field][:]*np.ones(nk))
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
                expodata[field] = np.append(expodata[field],tempvar)
            # If the variable can't be found in the field, append NaNs to main consistency in array shapes
            # Record that the a variable was skipped
            except:
                numskip[field] += 1
                expodata[field] = np.append(expodata[field],np.zeros(nk)*np.nan)
                if verbose:
                    print("Could not find %s in file %s" % (field,file) )
        # Closeout the file
        Dataset(join(datapath,file)).close()

    # Use the reference time to convert to python datetime format
    expodata['time'] = [ woce_ref_time + timedelta(minutes=time) for time in expodata['time'] ]
    # Print out how many files did not have the requested field in them
    if verbose:
        [print("Field %s could not be found in %d files" % (field, numskip[field])) for field in fields if numskip[field]>0 ]

    if aux_fields:
        expodata['abssal'] = gsw.conversions.SA_from_SP( expodata['bottle_salinity'], expodata['pressure'], expodata['longitude'], expodata['latitude'])
        expodata['ctemp' ] = gsw.conversions.CT_from_t(  expodata['abssal'],expodata['temperature'],expodata['pressure'])
        expodata['ptemp0'] = gsw.conversions.pt0_from_t( expodata['abssal'],expodata['temperature'],expodata['pressure'])
        expodata['sigma0'] = gsw.density.sigma0( expodata['abssal'],expodata['temperature'])
        expodata['sigma2'] = gsw.density.sigma2( expodata['abssal'],expodata['temperature'])
        expodata['rhoinsitu'] = gsw.density.rho( expodata['abssal'],expodata['temperature'],expodata['pressure'])

    expodata['expocode'] = exponame

    expodata = pandas.DataFrame.from_dict(expodata)
    # Calculate distance along transect
    expodata.sort_values('time',inplace=True)
    expodata['distance'] = np.append(0,gsw.distance(np.array(expodata.longitude), np.array(expodata.latitude),p=np.zeros(expodata.latitude.shape)).cumsum())

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
    ax.plot(expodata['longitude'],expodata['latitude'], 'bo', transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAND,facecolor='0.25')
    ax.set_global()
    plt.title(expodata['expocode'].upper())
    plt.show()

def grid_expo_variables( expodata, xgrid, ygrid, xvar = 'longitude', yvar = 'pressure' ):
    """
    Grid the scattered bottle data to the requested points using a 'universal kriging' approach
    Inputs:
        expodata: Dictionary containing the data from the cruise
        xgrid:    The x-value of points of the output grid
        ygrid:    The y-value of points of the output grid
    Optional inputs:
        xvar:     The name of the x variable (default: longitude)
        yvar:     The name of the y variable (default: pressure)
    """
    return



