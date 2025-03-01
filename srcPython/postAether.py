#!/usr/bin/env python
""" Standard model visualization routines
"""

from glob import glob
import re
from aetherpy.io import read_routines
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from netCDF4 import Dataset
from aetherpy.utils.time_conversion import datetime_to_epoch
import argparse
import os

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(description = 'Post process Aether files')
    parser.add_argument('-rm', \
                        help='removes processed files', \
                        action="store_true")

    args = parser.parse_args()

    return args



#----------------------------------------------------------------------------
# This returns the core of the filename without the _g????.nc
#----------------------------------------------------------------------------

def get_core_file(filename):
    coreFile = ''
    isEnsemble = False
    ensembleFile = ''
    ensembleNumber = -1
    m = re.match('.*([0123]D.*)(_g\d*)(\..*)',filename)
    if m:
        coreFile = m.group(1)
        # check if file is a member of an ensemble:
        check = re.match('.*([0123]D.*)(_m)(\d*)',coreFile)
        if (check):
            ensembleFile = check.group(1)
            isEnsemble = True
            ensembleNumber = int(check.group(3)) + 1

    fileInfo = {'coreFile': coreFile,
                'isEnsemble': isEnsemble,
                'ensembleFile': ensembleFile,
                'ensembleNumber': ensembleNumber,
                'ensembleMembers': -1}
    return fileInfo

#----------------------------------------------------------------------------
# Add to the list of strings if there isn't already an identical string
#----------------------------------------------------------------------------

def if_unique(list, newItem):
    IsFound = False
    index = -1
    for i, item in enumerate(list):
        if (item == newItem):
            IsFound = True
            index = i
    return IsFound, index

#----------------------------------------------------------------------------
# This looks at all of the netcdf files, figures out the core name, and
# adds them to the list of files to process
#----------------------------------------------------------------------------

def get_base_files():
    filelist = sorted(glob('?????*.nc'))
    files = []
    filesInfo = []
    for file in filelist:
        fileInfo = get_core_file(file)
        coreFile = fileInfo['coreFile']
        if (len(coreFile) > 0):
            IsFound, i = if_unique(files, coreFile)
            if (not IsFound):
                files.append(coreFile)
                filesInfo.append(fileInfo)

    # Figure out ensembles:
    # (1) get list of unique ensemble filess
    # (2) count how many files have this unique ensemble name
    
    ensembleFiles = []
    ensembleCounter = []
    for fileInfo in filesInfo:
        if (len(fileInfo['ensembleFile']) > 0):
            IsFound, i = if_unique(ensembleFiles, fileInfo['ensembleFile'])
            if (IsFound):
                ensembleCounter[-1] += 1
            else:
                ensembleFiles.append(fileInfo['ensembleFile'])
                ensembleCounter.append(1)
    
    # (3) store the number of ensemble members:            
    for i, fileInfo in enumerate(filesInfo):
        if (len(fileInfo['ensembleFile']) > 0):
            IsFound, item = if_unique(ensembleFiles, fileInfo['ensembleFile'])
            if (IsFound):
                filesInfo[i]['ensembleMembers'] = ensembleCounter[item]
    
    return filesInfo


#----------------------------------------------------------------------------
# Simply return the index of the matching string from a list
#----------------------------------------------------------------------------

def find_var_index(allVars, varToFind):
    index = -1
    for i,var in enumerate(allVars):
        if (var == varToFind):
            index = i
            break
    return index

#----------------------------------------------------------------------------
# Find the minimum and maximum of a variable in a bunch of blocks
#----------------------------------------------------------------------------

def determine_min_max(allBlockData, varToPlot, altToPlot):
    mini = 1e32
    maxi = -1e32
    for data in allBlockData:
        iVar = find_var_index(data['vars'], varToPlot)
        if (iVar > -1):
            v = data[iVar][2:-2, 2:-2, altToPlot]
            if (np.min(v) < mini):
                mini = np.min(v)
            if (np.min(v) > maxi):
                maxi = np.max(v)
    if (mini < 0):
        if (np.abs(mini) > maxi):
            maxi = np.abs(mini)
        mini = -maxi
    return mini, maxi
    
#----------------------------------------------------------------------------
# Plot a single block of a variable at a given altitude
#----------------------------------------------------------------------------

def plot_block(data, varToPlot, altToPlot, ax, mini, maxi, i):
    iLon = find_var_index(data['vars'], 'lon')
    iLat = find_var_index(data['vars'], 'lat')
    iAlt = find_var_index(data['vars'], 'z')
    if (iAlt < 0):
        iAlt = find_var_index(data['vars'], 'alt')
    iVar = find_var_index(data['vars'], varToPlot)

    if (mini < 0):
        cmap = cm.bwr
    else:
        cmap = cm.plasma
    
    alts = data[iAlt][0][0] / 1000.0  # Convert from m to km

    # Change to 2d representation:

    #lons = data[iLon][2:-2, 2:-2, 0]
    #lats = data[iLat][2:-2, 2:-2, 0]
    #v = data[iVar][2:-2, 2:-2, altToPlot]
    lons = data[iLon][:, :, 0]
    lats = data[iLat][:, :, 0]
    v = data[iVar][:, :, altToPlot]

    #print("lons : ", lons[:, 5])
    #print("lats : ", lats[:, 5])
    #print("value : ", i, v[:, 5]*180.0/np.pi)
    
    nLons, nLats = np.shape(lons)
    xp = np.zeros((nLons+1, nLats+1))
    yp = np.zeros((nLons+1, nLats+1))

    for iX in np.arange(1, nLons):
        for iY in np.arange(1, nLats):
            xp[iX, iY] = (lons[iX-1, iY] + lons[iX, iY])/2.0
            yp[iX, iY] = (lats[iX, iY-1] + lats[iX, iY])/2.0

        xp[iX, 0] = (lons[iX-1, 0] + lons[iX, 0])/2.0
        xp[iX, nLats] = xp[iX, nLats-1]

        yp[iX, 0] = 2 * lats[iX, 0] - yp[iX, 0]
        yp[iX, nLats] = 2 * yp[iX, nLats-1] - lats[iX, nLats-1]
    
    # more xp
    for iY in np.arange(1, nLats):
        xp[0, iY] = 2 * lons[0, iY] - xp[1, iY]
        xp[nLons, iY] = 2 * xp[nLons-1, iY] - lons[nLons-1, iY]
    xp[0, 0] = xp[0, 1]
    xp[nLons, nLats] = xp[nLons, nLats-1]

    xp[0, nLats] = xp[0, nLats-1]
    xp[nLons, 0] = 2 * xp[nLons-1, 0] - lons[nLons-1, 0]
            
    # more yp
    for iY in np.arange(1, nLats):
        yp[0, iY] = (lats[0, iY-1] + lats[0, iY])/2.0
        yp[nLons, iY] = (lats[nLons-1, iY-1] + lats[nLons-1, iY])/2.0
    yp[0, 0] = yp[1, 0]
    yp[nLons, nLats] = yp[nLons-1, nLats]

    yp[0, nLats] = 2 * yp[0, nLats-1] - lats[0, nLats-1]
    yp[nLons, 0] = yp[nLons-1, 0]
    
    plot = ax.scatter(lons, lats, c = v, cmap=cmap, vmin = mini, vmax = maxi)
    
    return plot
    
#----------------------------------------------------------------------------
# make a figure with all of the block plotted for the specified variable
# and altitude
#----------------------------------------------------------------------------

def plot_all_blocks(allBlockData, varToPlot, altToPlot, plotFile):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    mini, maxi = determine_min_max(allBlockData, varToPlot, altToPlot)
    i = 0
    for data in allBlockData:
        if (i < 100):
            plot = plot_block(data, varToPlot, altToPlot, ax, mini, maxi, i)
        i = i+1
    cbar = fig.colorbar(plot, ax = ax)
    cbar.set_label(varToPlot)
    ax.set_title(varToPlot)
    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('Latitude (deg)')
    print('  Outputting file : ', plotFile)
    fig.savefig(plotFile)
    plt.close(fig)
    
#----------------------------------------------------------------------------
# Read all of the block files in, given the core filename
#----------------------------------------------------------------------------

def read_block_files(coreFile):
    print('Figuring out coreFile : ', coreFile)
    fileList = sorted(glob(coreFile + '_g*.nc'))
    header = read_routines.read_aether_headers(fileList)
    allBlockData = []
    for file in fileList:
        print('  Reading file : ', file)
        data = read_routines.read_aether_file(file, file_vars=header['vars'])
        allBlockData.append(data)
    return allBlockData, fileList

#----------------------------------------------------------------------------
# return the nLons (nX), nLats (nY), and nAlts (nZ) of a block's data
#----------------------------------------------------------------------------

def get_sizes(allBlockData):
    lon3d = np.asarray(allBlockData[0][0])
    s = lon3d.shape
    nX = s[0]
    nY = s[1]
    nZ = s[2]
    return nX, nY, nZ
    
#----------------------------------------------------------------------------
# Write a NetCDF file from the data
#----------------------------------------------------------------------------

def write_netcdf(allBlockData, fileName):

    print('  Outputting file : ', fileName)
    ncfile = Dataset(fileName, 'w')

    nBlocks = len(allBlockData)
    nLons, nLats, nZ = get_sizes(allBlockData)

    lon_dim = ncfile.createDimension('lon', nLons)
    lat_dim = ncfile.createDimension('lat', nLats)
    z_dim = ncfile.createDimension('z', nZ)
    block_dim = ncfile.createDimension('block', None)
    time_dim = ncfile.createDimension('time', None)

    oneBlock = allBlockData[0]

    time_out = ncfile.createVariable('time', np.float64, ('time',))
    time_out[0] = datetime_to_epoch(oneBlock["time"])
        
    allNetCDFVars = []
    # create all of the variables
    varList = []
    for iV, v in enumerate(oneBlock['vars']):
        varList.append(v)
        longName = oneBlock['long_name'][iV]
        unitName = oneBlock['units'][iV]
        allNetCDFVars.append(ncfile.createVariable(v, np.float32, \
                                                   ('block', 'lon', 'lat', 'z')))
        allNetCDFVars[-1].units = unitName
        allNetCDFVars[-1].long_name = longName

        for iB, oneBlock in enumerate(allBlockData):
            tmp = np.asarray(oneBlock[iV])
            allNetCDFVars[-1][iB,:,:,:] = tmp
        
    ncfile.close()

#----------------------------------------------------------------------------
# copy block data in one file
#----------------------------------------------------------------------------

def copy_or_add_block_data(allBlockData,
                           oldBlockData = [],
                           factor = 1.0):

    if (len(oldBlockData) > 0):
        doAdd = True
    else:
        doAdd = False

    newBlockData = []

    for ib, oneBlock in enumerate(allBlockData):
        obCopy = {}
        for key in oneBlock.keys():
            if (type(key) == int):
                if (doAdd):
                    obCopy[key] = oldBlockData[ib][key] + \
                        oneBlock[key] * factor
                else:
                    obCopy[key] = oneBlock[key] * factor
            else:
                obCopy[key] = oneBlock[key]
        newBlockData.append(obCopy)
        
    return newBlockData

#----------------------------------------------------------------------------
# copy block data in one file
#----------------------------------------------------------------------------

iCopy_ = 0
iAdd_ = 1
iSub_ = 2
iMult_ = 3
iPower_ = 4

def do_math_on_block_data(blockData1,
                          blockData2 = [],
                          math = iCopy_,
                          factor = 1.0,
                          iLowestVar = 3):

    newBlockData = []

    for ib, oneBlock in enumerate(blockData1):
        obCopy = {}
        for key in oneBlock.keys():
            if (type(key) == int):
                if (key >= iLowestVar):
                    if (math == iCopy_):
                        obCopy[key] = oneBlock[key] * 1.0
                    if (math == iAdd_):
                        obCopy[key] = oneBlock[key] + blockData2[ib][key]
                    if (math == iSub_):
                        obCopy[key] = oneBlock[key] - blockData2[ib][key]
                    if (math == iMult_):
                        obCopy[key] = oneBlock[key] * factor
                    if (math == iPower_):
                        obCopy[key] = oneBlock[key] ** factor
                else:
                    obCopy[key] = oneBlock[key]
            else:
                obCopy[key] = oneBlock[key]                    
        newBlockData.append(obCopy)
    return newBlockData

#----------------------------------------------------------------------------
# Calc Standard Deviation of Ensemble Run
#----------------------------------------------------------------------------

def calc_std_of_ensembles(filesInfo,
                          ensembleIndexList,
                          ensembleMean):
    
    for i, iF in enumerate(ensembleIndexList):
        cF = filesInfo[iF]['coreFile']
        print('---> Going back through corefiles: ', cF)
        allBlockData, filelist = read_block_files(cF)

        # subtract
        diff = do_math_on_block_data(allBlockData,
                                     blockData2 = ensembleMean,
                                     math = iSub_)
        # square
        diffs = do_math_on_block_data(diff,
                                      factor = 2,
                                      math = iPower_)
                
        if (i == 0):
            sums = do_math_on_block_data(diffs,
                                         math = iCopy_)
        else:
            sums = do_math_on_block_data(sums,
                                         blockData2 = diffs,
                                         math = iAdd_)
    factor = 1.0 / fileInfo['ensembleMembers']
    sumsD = do_math_on_block_data(sums,
                                  factor = factor,
                                  math = iMult_)
    stdData = do_math_on_block_data(sumsD,
                                    factor = 0.5,
                                    math = iPower_)

    return stdData


#----------------------------------------------------------------------------
# write and plot data
#----------------------------------------------------------------------------

def write_and_plot_data(dataToWrite,
                        fileStart,
                        fileAddon,
                        iVar,
                        iAlt):

    netcdfFile = fileStart + fileAddon + '.nc'
    write_netcdf(dataToWrite, netcdfFile)

    plotFile = fileStart + fileAddon + '.png'
    var = dataToWrite[0]['vars'][iVar]
    plot_all_blocks(dataToWrite, var, iAlt, plotFile)

    return


#----------------------------------------------------------------------------
# main code
#----------------------------------------------------------------------------

args = parse_args()

filesInfo = get_base_files()

iVar = 3
iAlt = 10

for iFile, fileInfo in enumerate(filesInfo):
    coreFile = fileInfo['coreFile']
    print(coreFile)
    allBlockData, filelist = read_block_files(coreFile)
    write_and_plot_data(allBlockData, coreFile, '', iVar, iAlt)

    if (fileInfo['isEnsemble']):
        factor = 1.0 / float(fileInfo['ensembleMembers'])
        if (fileInfo['ensembleNumber'] == 1):
            ensembleData = copy_or_add_block_data(allBlockData,
                                                  factor = factor)
            ensembleIndexList = [iFile]
        else:
            ensembleData = copy_or_add_block_data(allBlockData,
                                                  oldBlockData = ensembleData,
                                                  factor = factor)
            ensembleIndexList.append(iFile)
        if (fileInfo['ensembleNumber'] == fileInfo['ensembleMembers']):

            write_and_plot_data(ensembleData, fileInfo['ensembleFile'],
                                '_mean', iVar, iAlt)
            
            stdData = calc_std_of_ensembles(filesInfo,
                                            ensembleIndexList,
                                            ensembleData)
            write_and_plot_data(stdData, fileInfo['ensembleFile'],
                                '_std', iVar, iAlt)
                        
    if (args.rm):
        print('  ---> Removing files...')
        for file in filelist:
            command = 'rm -f '+file
            os.system(command)
    
