#! /usr/bin/env python

import math, re, datetime, sys
import numpy as np
from io_tools import *


# Read background binary files from Polair3D output
def read_bkgd_bin(bkgd_species, indir, Nt, Nx, Ny, Nz):
    # Put data in a dict
    bkgd_data = {}
    for species in bkgd_species:
        infile = indir + species + '.bin'
        bkgd_data[species] = np.memmap(infile, dtype='float32', mode='r',
                                       shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]

    return bkgd_data


# Set background for streets
def set_bkgd_bin(street_list, bkgd_data, current_date, date_min, delta_t,
                  Nt, x_min, y_min, delta_x, delta_y, Nx, Ny):
    # Get index of current date
    c_id = int((current_date - date_min).total_seconds() / delta_t)
    if c_id < 0 or c_id >= Nt:
        sys.exit('ERROR: background data not available for this date.')

    for i in range(len(street_list)):
        # Get cell indices (X, Y) in the Polair3D grid
        # Xid = (street_list[i].lon_cen - x_min) / delta_x
        # Yid = (street_list[i].lat_cen - y_min) / delta_y
        # Xid = int(Xid)-1 if Xid%1 < 0.5 else int(Xid)
        # Yid = int(Yid)-1 if Yid%1 < 0.5 else int(Yid)
        Xid, Yid = get_polair_id(street_list[i].lon_cen, street_list[i].lat_cen,
                                 x_min, y_min, delta_x, delta_y, Nx, Ny)
        for key, value in bkgd_data.items():
            street_list[i].bkgd[key] = value[c_id, Yid, Xid]

# Append background to binary files
def append_bkgd_data(street_list, bkgd_files):
    for i, var in enumerate(street_list[0].bkgd.keys()):
        data = np.zeros((len(street_list)), 'float')
        for j in range(len(street_list)):
            data[j] = street_list[j].bkgd[var]
        append_binary(data, bkgd_files[var])


# Read the background concentration from a text file.
def get_background_concentration(input_file, current_date, street_list):


    hasBackground = False
    input_background = open(input_file)
    header = input_background.readline()
    nstreet = 0
    for line in input_background.readlines():
        line = line.replace('\n','')
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        str_times = line_info[0]
        year = int(str_times[0:4])
        month = int(str_times[5:7])
        day = int(str_times[8:10])
        hour = int(str_times[11:13])

        background_date = datetime.datetime(year,month,day,hour)
        if current_date == background_date:
            print("Background concentration data are available for the date ", current_date)
            o3 = float(line_info[1])
            no2 = float(line_info[2])
            no = float(line_info[3])
            hasBackground = True
            break

    if (hasBackground == False):
        print("Error: background concentration data are not available.")
    else:   
        for s in range(len(street_list)):
            street = street_list[s]
            street.bkgd['O3'] = o3
            street.bkgd['NO2'] = no2
            street.bkgd['NO'] = no

    return 0
        
def get_polair_id(lon, lat, x_min, y_min, dx, dy, Nx, Ny):
    # Get cell indices (X, Y) in the Polair3D grid
    # Xid = (street_list[i].lon_cen - x_min) / delta_x
    # Yid = (street_list[i].lat_cen - y_min) / delta_y
    # Xid = int(Xid)-1 if Xid%1 < 0.5 else int(Xid)
    # Yid = int(Yid)-1 if Yid%1 < 0.5 else int(Yid)

    Xid = max(int((lon - x_min + dx / 2.) / dx), 0)
    Xid = min(Xid, Nx - 1)
    Yid = max(int((lat - y_min + dy / 2.) / dy), 0)
    Yid = min(Yid, Ny - 1)

    return Xid, Yid

def get_polair_ind(polair_lon, polair_lat, street):
    min_distance = 99.
    i = 0
    j = 0
    for lon in polair_lon:
        for lat in polair_lat:
            distance = np.sqrt(pow((lon - street.lon_cen), 2.0) + 
                               pow((lat - street.lat_cen), 2.0))
            if distance < min_distance:
                indx = i
                indy = j
                min_distance = distance
            j = j + 1
        j = 0    
        i = i + 1
    return indx, indy

def get_chimere_background_concentration(current_date, street_list, melchior_spec_list,molar_mass_melchior2,chimout_dir,chimout_lab) :

    import netCDF4
    import os,sys

    str_date = current_date.strftime("%Y%m%d")
    input_file=chimout_dir+'/out.'+str_date+'00_'+chimout_lab+'.nc'
    if not os.path.isfile(input_file) :
       print(('CHIMERE background conditions are requested but the file is not found: '+str(input_file)))
       sys.exit()

    nc = netCDF4.Dataset(input_file, 'r')
    chim_times = nc.variables["Times"][:]
    lons = nc.variables["lon"][:]
    lats = nc.variables["lat"][:]
    new_spec_list=[] 
    # New melchior species list from what is actually present in the CHIMERE out file
    for spec in melchior_spec_list:
        if spec in list(nc.variables.keys()):
           new_spec_list.append(spec)
        else:
           print(('Warning!! '+spec+' not found in CHIMERE output file'))

    # Transform CHIMERE date-time to python datetime

    N=chim_times.shape[0] #number of hours to parse
    times=[]
    for i in range(N): #loop over hours
            s1,s2,s3,s4=chim_times[i,0:4]
            YEAR=s1+s2+s3+s4
            s1,s2=chim_times[i,5:7]
            MONTH=s1+s2
            s1,s2=chim_times[i,8:10]
            DAY=s1+s2
            s1,s2=chim_times[i,11:13]
            HOUR=s1+s2
            times.append(datetime.datetime(int(YEAR),int(MONTH),int(DAY),int(HOUR)))
    times=np.array(times)

    for spec in new_spec_list:
        for s in range(len(street_list)):
            street = street_list[s]
            street.background[spec]=0.0


    hasBackground = False
    for t in range(len(times)):
        background_date = times[t]
        if current_date == background_date:
            print("Background data are available for the date ", current_date)
            ind_t = t
            hasBackground = True

    if (hasBackground == False):
        print("Error: background data are not available")
        sys.exit()

    print((lons.shape))
    nx = lons.shape[1]
    ny = lons.shape[0]

    for s in range(len(street_list)):
       street = street_list[s]
       lat1 = street.lat_cen
       lon1 = street.lon_cen
       init_length = 9999.0
       for i in range(nx):
           for j in range(ny):
               lat2 = lats[j, i]
               lon2 = lons[j, i]
               length = distance_on_unit_sphere(lat1, lon1, lat2, lon2)
               if length < init_length:
                   init_length = length
                   ind_i = i
                   ind_j = j
       print((ind_i,ind_j))
       tem2=nc.variables['tem2'][ind_t,ind_j,ind_i] #Kelvin
       psfc=nc.variables['pres'][ind_t,0,ind_j,ind_i] #Pascal

       if psfc > 0 :
          molecular_volume=22.41 * (tem2/273.)  * (1013*10**2)/psfc
       else :
          molecular_volume=22.41

       for spec in new_spec_list:
           conc_ppb=nc.variables[spec][ind_t,0,ind_j,ind_i]
           molecular_mass=molar_mass_melchior2[spec]
           ppb2ug=molecular_mass / molecular_volume   #
           if spec=='NO2' : print((s,ind_t,0,ind_j,ind_i,conc_ppb,conc_ppb * ppb2ug))
           street.background[spec]=conc_ppb * ppb2ug

    return 0     #street.background in ug/m3
