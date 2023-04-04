import re, sys, os
from compute_distance import *
import datetime
import numpy as np

# Class for a node (intersection)
# -------------------------------
class Node:
    def __init__(self, node_id, lon, lat):
        self.lon = lon
        self.lat = lat
        self.id = node_id
        self.eff_id = 0
        self.connected_street = []
        self.removed = False
        self.wind_dir = 0.0
        self.wind_speed = 0.0
        self.pblh = 0.0
        self.ust = 0.0
        self.lmo = 0.0

# Class for a street
# ------------------
class Street:
    def __init__(self, st_id, begin_node, end_node, length, width, height, lon_cen, lat_cen, emission):
        self.id = st_id
        self.eff_id = st_id
        self.begin = begin_node
        self.end = end_node
        self.eff_begin = begin_node
        self.eff_end = end_node
        self.length = length
        self.width = width
        self.height = height
        self.removed = False
        self.lon_cen = lon_cen
        self.lat_cen = lat_cen
        self.wind_dir = 0.0
        self.wind_speed = 0.0
        self.pblh = 0.0
        self.ust = 0.0
        self.lmo = 0.0
        self.sh = 0.0
        self.psfc = 0.0
        self.t2 = 0.0
        self.attenuation = 0.0
        self.background={}     # in ug/m3
        self.emission = emission
        self.eff_emission = emission
        self.typo = 0
        self.lwc = 0.0
        self.rain = 0.0


# Check if two nodes are same
# ---------------------------
def are_nodes_same(node1, node2):
    if ((node1.lon == node2.lon) and (node1.lat == node2.lat)):
        return True
    else:
        lon1, lat1 = node1.lon, node1.lat
        lon2, lat2 = node2.lon, node2.lat
        length = distance_on_unit_sphere(lat1, lon1, lat2, lon2) / 1000. # in km
        if ((length < 0.01)):
            return True
        else:
            return False

# def are_streets_same(node_list, street1, street2):
    # for li in node_list:
    #     if street1.begin == li.id:
    #         node_begin1 = li
    #     if street1.end == li.id:
    #         node_end1 = li
    #     if street2.begin == li.id:
    #         node_begin2 = li
    #     if street2.end == li.id:
    #         node_end2 = li

def are_streets_same(node_list, i, j):
    
    node_begin1 = node_list[2 * i]
    node_end1 = node_list[2 * i + 1]
    node_begin2 = node_list[2 * j]
    node_end2 = node_list[2 * j +1]

    
    # Node11 : the begin node of the street 1
    # Node12 : the end node of the street 1
    # Node21 : the begin node of the street 2
    # Node22 : the end node of the street 2
    # Check if Node11 == Node21 and Node12 == Node22
    if are_nodes_same(node_begin1, node_begin2) and are_nodes_same(node_end1, node_end2):
        return True
    # Check if Node11 == Node22 and Node12 == Node21
    elif are_nodes_same(node_begin1, node_end2) and are_nodes_same(node_end1, node_begin2):
        return True
    else:
        return False

# Merge streets in the case where there are two emission data sets for both ways. 
# --------------------

def merging_street(output_file, node_list, street_list):
    ntemp = 0
    n_street = len(street_list)
    for i in range(n_street - 1):
      if (street_list[i].removed == False):
        is_street_found = False
        j = i + 1
        while (is_street_found == False and j < n_street):
            # if are_streets_same(node_list, street_list[i], street_list[j]):
            if are_streets_same(node_list, i, j):                
                street_list[j].eff_begin = street_list[i].begin
                street_list[j].eff_end = street_list[i].end
                street_list[j].eff_id = street_list[i].id
                street_list[i].eff_emission = street_list[i].emission + street_list[j].emission
                street_list[j].eff_emission = 0.0
                street_list[j].removed = True
                ntemp = ntemp + 1
                is_street_found = True
            j = j + 1

    f = open(output_file, 'w')
    f.write("# Street id \t Effective street id \n")
    for i in range(n_street):
        street = street_list[i]
        f.write(str(street.id) + "\t" + str(street.eff_id) + "\n")
    f.close()

    return ntemp                

# Manual merging of intersections.
# -------------------------------
def manual_merging_street(street_list):
    input_file = "street-merging.txt"
    input_merging = open(input_file)
    print("Manual merging using the street list in " + input_file)
    header = input_merging.readline()
    ntemp = 0
    for line in input_merging.readlines():
        line = line.replace('\n','')
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        if len(line_info) != 2:
            break
        else:
            removed_street_id, remained_street_id = int(line_info[0]), int(line_info[1])
            for i in range(len(street_list)):
                if (street_list[i].id == remained_street_id):
                    for j in range(len(street_list)):
                        if (street_list[j].id == removed_street_id):
                            street_list[j].eff_id = street_list[i].id
                            street_list[j].eff_begin = street_list[i].begin
                            street_list[j].eff_end = street_list[i].end
                            street_list[i].eff_emission = street_list[i].eff_emission + street_list[j].eff_emission
                            street_list[j].eff_emission = 0.0
                            street_list[j].removed = True
                            ntemp = ntemp + 1
    input_merging.close()
    return ntemp

# Merging using a look-up table.
# -------------------------------
def lut_merging_street(lut_file, street_list):
    n_street = len(street_list)
    ntemp = 0
    print("Read the lookup-table: ", lut_file)    
    input_merging = open(lut_file)
    header = input_merging.readline()
    for line in input_merging.readlines():
        line = line.replace('\n','')
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        if len(line_info) != 2:
            break
        else:
            street_id, effective_street_id = int(line_info[0]), int(line_info[1])
            if street_id != effective_street_id:
                for nst in range(n_street):
                    if (street_list[nst].id == effective_street_id):
                        i = nst 
                    elif (street_list[nst].id == street_id):
                        # the emissions in this street are merged into
                        # the street having the index i.
                        j = nst 
                street_list[j].eff_begin = street_list[i].begin
                street_list[j].eff_end = street_list[i].end
                street_list[j].eff_id = street_list[i].id
                street_list[i].eff_emission = street_list[i].emission + street_list[j].emission
                street_list[j].eff_emission = 0.0
                street_list[j].removed = True
                ntemp = ntemp + 1
    input_merging.close()
    return ntemp

# Remove the same nodes
# --------------------
def merging_node(node_list):
    n_node = 0
    for i in range(len(node_list) - 1):
      if (node_list[i].removed == False):  
        for j in range(i + 1, len(node_list)):
          if (node_list[j].removed == False):  
            if ((node_list[i].lon == node_list[j].lon) and 
                (node_list[i].lat == node_list[j].lat)) :
                node_list[j].eff_id = node_list[i].eff_id
                node_list[j].removed = True
                n_node = n_node + 1
                for st_ in node_list[j].connected_street:
                    node_list[i].connected_street.append(st_)
    return n_node


# Remove the near nodes: 
# the distance between the nodes 
# is smaller than 10 m
# --------------------
def merging_near_node(node_list):
    n_node = 0
    for i in range(len(node_list) - 1):
      if (node_list[i].removed == False):  
        lon1, lat1 = node_list[i].lon, node_list[i].lat
        for j in range(i + 1, len(node_list)):
          if (node_list[j].removed == False):  
            lon2, lat2 = node_list[j].lon, node_list[j].lat
            if ((lon1 != lon2) and (lat1 != lat2)):
                length = distance_on_unit_sphere(lat1, lon1, lat2, lon2) / 1000. # in km
                if ((length < 0.01)):
                    id_old = node_list[j].eff_id
                    node_list[j].eff_id = node_list[i].eff_id
                    id_new = node_list[j].eff_id
                    node_list[j].removed = True
                    n_node = n_node + 1
                    for st_ in node_list[j].connected_street:
                        # check if st_ is an element of node_list[i]
                        if (st_ in node_list[i].connected_street):
                            sys.exit("Error: street " + str(st_) + " is " \
                                     "already in node_list.")
                        else:
                            node_list[i].connected_street.append(st_)

                    # Node A was merged to Node B.
                    # Node B was merged to Node C.
                    # Node A should be merged to Node C.
                    for node_ in node_list:
                        if (node_.eff_id == id_old):
                            node_.eff_id = id_new
    return n_node

# Manual merging of intersections.
# -------------------------------
def manual_merging_node(node_list):
    n_node = 0
    input_file = "intersection-merging.txt"
    print("Manual merging of the nodes using the node list in " + input_file)
    input_merging = open(input_file)
    header = input_merging.readline()
    for line in input_merging.readlines():
        line = line.replace('\n','')
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        if len(line_info) != 2:
            break
        else:
            removed_node_id, remained_node_id = int(line_info[0]), int(line_info[1])
            for i in range(len(node_list)):
                if (node_list[i].id == removed_node_id): 
                    for j in range(len(node_list)):
                        if (node_list[j].id == remained_node_id): 
                            id_old = node_list[i].eff_id
                            node_list[i].eff_id = node_list[j].eff_id
                            id_new = node_list[i].eff_id
                            node_list[i].removed = True
                            n_node = n_node + 1
                            for st_ in node_list[i].connected_street:
                                node_list[j].connected_street.append(st_)
                            #     
                            for node_ in node_list:
                                if (node_.eff_id == id_old):
                                    node_.eff_id = id_new
    input_merging.close()
    return n_node

# Get the street width and the builiding height from the input file.
# ------------------------------------------------------------------
def get_street_geog(input_file, street_list):
    # Read street width, length and builiding height.
    print("Read the geographical informations: ", input_file)
    input_street_geog = open(input_file)
    header = input_street_geog.readline()
    nstreet = 0
    for line in input_street_geog.readlines():
        line = line.replace('\n','')        
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        street_id = int(line_info[0])    
        for i in range(len(street_list)):
            street = street_list[i]
            if street.id == street_id:
                street.width = float(line_info[2])
                street.height = float(line_info[3])
                nstreet = nstreet + 1
    input_street_geog.close()
    if nstreet != len(street_list):
        print("Warning: Missing input data: data given for %d streets, but data are needed for %d streets." % (nstreet, len(street_list)))
    return 0


def get_meteo_data(input_dir, current_date, street_list, node_list, wrfout_prefix):
    ka = 0.4
    g = 9.81
    r = 287.

    from scipy.io import netcdf
    str_date = current_date.strftime("%Y-%m-%d")
    input_file = input_dir + wrfout_prefix + "_" + str_date + "_00:00:00"
    f = netcdf.netcdf_file(input_file, 'r')

    start_date = f.__getattribute__("SIMULATION_START_DATE")
    print(str_date, start_date)
    
    previous_date = current_date - datetime.timedelta(seconds = 3600)
    str_prev_date = previous_date.strftime("%Y-%m-%d")
    previous_file = input_dir + wrfout_prefix + "_" + str_prev_date + "_00:00:00"
     
    times = f.variables["Times"][:]
    lons = f.variables["XLONG"][:]
    lats = f.variables["XLAT"][:]
    u10 = f.variables["U10"][:]
    v10 = f.variables["V10"][:]
    pblh = f.variables["PBLH"][:]
    ust = f.variables["UST"][:]
    theta = f.variables["T"][:] # Potential temperature perturbation
    tsk = f.variables["TSK"][:] # Skin temperature
    psfc = f.variables["PSFC"][:] # Surface pressure
    hfx = f.variables["HFX"][:] # Sensible heat
    lh = f.variables["LH"][:] # Latent heat

    nx = lons.shape[2]
    ny = lons.shape[1]
    
    # For chemistry
    sh = f.variables["QVAPOR"][:] # Specific humidity
    t2 = f.variables["T2"][:] # Surface temperature
    lwc = f.variables["QCLOUD"][:] # Cloud water mixing ratio
    pressure_pert = f.variables["P"][:] # Pressure Perturbation
    base_pres = f.variables["PB"][:] # Base state pressure

    rainc = f.variables["RAINC"][:] # Convective Rain
    rainnc = f.variables["RAINNC"][:] # Non-convective Rain
    rain = rainc + rainnc

    prev_accumulated_rain = True
    print(rain.shape)
    if (os.path.isfile(previous_file) and prev_accumulated_rain):
        f_prev = netcdf.netcdf_file(previous_file, 'r')
        rainc_prev = f_prev.variables["RAINC"][:] # Convective Rain
        rainnc_prev = f_prev.variables["RAINNC"][:] # Non-convective Rain
        rain_prev = (rainc_prev + rainnc_prev)
        for i in range(nx):
            for j in range(ny):
                rain[0,j,i] = max((rain[0,j,i] - rain_prev[-1,j,i]), 0.0)
    else:
        print("File for the previous date is not available.")
        
    solar_radiation = f.variables["SWDOWN"][:] # Solar Radiation
    

    
    # Heights
    ptop = f.variables["P_TOP"][0]
    Tiso = f.variables["TISO"][0]
    p00 = f.variables["P00"][0]
    Ts0 = f.variables["T00"][0]
    A = f.variables["TLP"][0]
    Terrain = f.variables["HGT"][:]
    sigma_h = f.variables["ZNU"][:]

    Piso =  p00 * np.exp((Tiso - Ts0) / A)
    aux = np.log(Piso / p00)
    Ziso = - r * A * aux*aux / (2.0*g) - r * Ts0 * aux / g
    nz = theta.shape[1]
   
    print("WRF output data (nt, nz, ny, nx): ", theta.shape)
    hasMeteo = False
    for t in range(len(times)):
        str_times = b''.join(times[t])
        year = int(str_times[0:4])
        month = int(str_times[5:7])
        day = int(str_times[8:10])
        hour = int(str_times[11:13])
        meteo_date = datetime.datetime(year,month,day,hour)
        if current_date == meteo_date:
            print("Meteo data are available for the date ", current_date)
            ind_t = t
            hasMeteo = True
    if (hasMeteo == False):
        print("Error: meteo data are not available")
        sys.exit()

    # Decumulate rain
    if (ind_t == 0):
        rain = rain[ind_t]
    else:
        rain = rain[ind_t] - rain[ind_t - 1]
   
    
    # Get meteo data for the streets
    for s in range(len(street_list)):
        street = street_list[s]
        lat1 = street.lat_cen
        lon1 = street.lon_cen
        init_length = 9999.0
        for i in range(nx):
            for j in range(ny):
                lat2 = lats[0, j, i]
                lon2 = lons[0, j, i]
                length = distance_on_unit_sphere(lat1, lon1, lat2, lon2)
                if length < init_length:
                    init_length = length
                    ind_i = i
                    ind_j = j
        u10_cell = u10[ind_t, ind_j, ind_i]
        v10_cell = v10[ind_t, ind_j, ind_i]
        pblh_cell = pblh[ind_t, ind_j, ind_i]
        ust_cell = ust[ind_t, ind_j, ind_i]
        sh_cell = sh[ind_t, 0, ind_j, ind_i]
        psfc_cell = psfc[ind_t, ind_j, ind_i]
        t2_cell = t2[ind_t, ind_j, ind_i]

        # For chemistry
        lwc_cell = lwc[ind_t, 0, ind_j, ind_i]
        lwc_colon = lwc[ind_t, :, ind_j, ind_i]
        rain_cell = rain[ind_j, ind_i]
        solar_radiation_cell = solar_radiation[ind_t, ind_j, ind_i]
    
        
        # Compute wind direction
        temp = math.atan2(v10_cell, u10_cell) # -pi < temp < pi
        if temp <= (math.pi / 2.0):
            temp2 = math.pi / 2.0 - temp
        else:
            temp2 = (math.pi / 2.0 - temp) + 2.0 * math.pi
        street.wind_dir = temp2 # in radian, 0 for the wind to north, pi to south.

        # Compute heights
        terrain_cell = Terrain[ind_t, ind_j, ind_i]
        sigma_h_colon = sigma_h[ind_t, :]
        gridz_interf = np.zeros([nz + 1], 'float')
        gridz_in = np.zeros([nz], 'float')
        for k in range(nz):
            ps0 = p00 * np.exp(-Ts0 / A + np.sqrt(Ts0*Ts0/ (A*A) - 2. * g * terrain_cell  / (A*r))) - ptop
            refpress = sigma_h_colon[k] * ps0 + ptop;
            if (refpress >= Piso):
                aux = np.log(refpress / p00);
                gridz_in[k] = - r * A * aux * aux / (2.0 * g) - r * Ts0 * aux / g - terrain_cell
            else:
                aux = np.log(refpress / Piso)
                gridz_in[k] = Ziso - r * Tiso * aux / g - terrain_cell
        gridz_interf[0] = 0.0        
        for k in range(1, nz + 1):
            gridz_interf[k] = 2 * gridz_in[k - 1] - gridz_interf[k - 1]

        # Compute attenuation.
        pressure_colon = pressure_pert[ind_t, :, ind_j, ind_i] + base_pres[ind_t, :, ind_j, ind_i]
        theta_colon = theta[ind_t, :, ind_j, ind_i] + 300.
        sh_colon = sh[ind_t, :, ind_j, ind_i]
        rh_colon = np.zeros([nz], 'float')
        crh_colon = np.zeros([nz], 'float')
        cloud_fraction = np.zeros([nz], 'float')
        for k in range(nz):
            rh_colon[k] = compute_relative_humidity(sh_colon[k], 
                                                    theta_colon[k], 
                                                    pressure_colon[k])
            crh_colon[k] = compute_crh(psfc_cell, pressure_colon[k])
            cloud_fraction[k] = compute_cloud_fraction(rh_colon[k],
                                                       crh_colon[k])

        medium_cloudiness, high_cloudiness = compute_cloudiness(nz, 
                                                                pressure_colon,
                                                                cloud_fraction,
                                                                gridz_interf)

        attenuation = np.zeros([nz], 'float')
        compute_attenuation(nz, gridz_in, rh_colon,
                            medium_cloudiness, high_cloudiness, attenuation)
        
        #
        street.wind_speed = math.sqrt(u10_cell ** 2 + v10_cell ** 2)
        street.pblh = pblh_cell
        street.ust = ust_cell
        street.sh = sh_colon[0]
        street.psfc = psfc_cell
        street.t2 = t2_cell
        street.attenuation = attenuation[0]

        # Compute LMO
        theta_cell = theta[ind_t, 0, ind_j, ind_i] + 300.
        theta0_cell = tsk[ind_t, ind_j, ind_i] * pow((psfc[ind_t, ind_j, ind_i] / 101325.0), (-287.0 / 1005.0))
        theta_mean = 0.5 * (theta_cell + theta0_cell)
        evaporation = lh[ind_t, ind_j, ind_i] / 2.5e9
        hfx_cell = hfx[ind_t, ind_j, ind_i]
        lmo = (-1.0 * pow(ust_cell, 3.0) * theta_mean) / (ka * g * (hfx_cell + 0.608 * theta_mean * evaporation))
        street.lmo = lmo
        
        street.lwc = lwc_cell
        street.rain = rain_cell
        
    # Get meteo data for the intersections.    
    for n in range(len(node_list)):
        node = node_list[n]
        lat1 = node.lat
        lon1 = node.lon
        init_length = 9999.0
        for i in range(nx):
            for j in range(ny):
                lat2 = lats[0, j, i]
                lon2 = lons[0, j, i]
                length = distance_on_unit_sphere(lat1, lon1, lat2, lon2)
                if length < init_length:
                    init_length = length
                    ind_i = i
                    ind_j = j

        u10_cell = u10[ind_t, ind_j, ind_i]
        v10_cell = v10[ind_t, ind_j, ind_i]
        pblh_cell = pblh[ind_t, ind_j, ind_i]
        ust_cell = ust[ind_t, ind_j, ind_i]

        # Compute LMO
        theta_cell = theta[ind_t, 0, ind_j, ind_i] + 300.
        theta0_cell = tsk[ind_t, ind_j, ind_i] * pow((psfc[ind_t, ind_j, ind_i] / 101325.0), (-287.0 / 1005.0))
        theta_mean = 0.5 * (theta_cell + theta0_cell)
        evaporation = lh[ind_t, ind_j, ind_i] / 2.5e9
        hfx_cell = hfx[ind_t, ind_j, ind_i]
        lmo = (-1.0 * pow(ust_cell, 3.0) * theta_mean) / (ka * g * (hfx_cell + 0.608 * theta_mean * evaporation))

        temp = math.atan2(v10_cell, u10_cell) # -pi < temp < pi
        if temp <= (math.pi / 2.0):
            temp2 = math.pi / 2.0 - temp
        else:
            temp2 = (math.pi / 2.0 - temp) + 2.0 * math.pi
        node.wind_dir = temp2 # in radian, 0 for the wind to north, pi to south.
        node.wind_speed = math.sqrt(u10_cell ** 2 + v10_cell ** 2)
        node.pblh = pblh_cell
        node.ust = ust_cell
        node.lmo = lmo

    return 0

def compute_attenuation(Nz, gridz_in, RelativeHumidity, 
                        MediumCloudiness, HighCloudiness, Attenuation):
    a = 0.1
    b = 0.3
    c = 1.5
    B = 0.0
    norm = 0.0
    #// Calculation of B.
    for k in range(Nz):
        if (gridz_in[k] < 1500.):
            dz = gridz_in[k + 1] - gridz_in[k]
            if (RelativeHumidity[k] > 0.7):
                  B += (RelativeHumidity[k] - 0.7) * dz
            norm = norm + (1. - 0.7) * dz;

    # // Normalization.
    B = B / norm

    for k in range(Nz):
        Attenuation[k] = (1. - a * HighCloudiness) * (1. - b * MediumCloudiness)* np.exp(-c * B);

    return Attenuation


def compute_cloudiness(Nz, Pressure, CloudFraction, GridZ_interf):
    #            /*** Low clouds ***/
    cloud_max = 0;
    P_0 = 80000.
    P_1 = 45000.
    LowCloudiness = 0.0
    MediumCloudiness = 0.0
    HighCloudiness = 0.0
    # // The first level is excluded.
    for k in range(1, Nz):
        if (Pressure[k] > P_0):
            cloud_max = max(cloud_max, CloudFraction[k])
    below = True; above = False
    k_base = 0; k_top = 0
    for k in range(1, Nz):
        if (Pressure[k] > P_0 and (not above)):
            below = below and ( CloudFraction[k] < 0.5 * cloud_max
                               or CloudFraction[k] == 0 )
            above = (not below) and CloudFraction[k] < 0.5 * cloud_max
            if ((not below) and k_base == 0):
                  k_base = k
            if (above):
                  k_top = k

    if (above):
        while (k < Nz and Pressure[k] > P_0):
                k = k + 1
    k_max = k - 1
    # // Goes up to P_0.
    if (k_base > k_top):
        k_top = k_max + 1
    LowIndices_0 = k_base
    LowIndices_1 = k_top
    k = k_base
    # // k_top == 0 means no cloud.
    while (k < k_top and k_top != 0):
        LowCloudiness = LowCloudiness + CloudFraction[k] * (GridZ_interf[k + 1] - GridZ_interf[k])
        k = k + 1
   
    if (k_top != 0):
        LowCloudiness = LowCloudiness / (GridZ_interf[k_top] - GridZ_interf[k_base])

    # /*** Medium clouds ***/

    cloud_max = 0;
    # // Starts above low clouds.
    for k in range(k_max + 1, Nz):
        if (Pressure[k] > P_1):
            cloud_max = max(cloud_max, CloudFraction[k])
    below = True; above = False
    k_base = 0; k_top = 0;
    for k in range(k_max + 1, Nz):
        if (Pressure[k] > P_1 and (not above)):
            below = below and ( CloudFraction[k] < 0.5 * cloud_max
                        	   or CloudFraction[k] == 0 )
            above = (not below) and CloudFraction[k] < 0.5 * cloud_max
            if ((not below) and k_base == 0):
                k_base = k
            if (above):
                k_top = k

    if (above):
        while (k < Nz and Pressure[k] > P_1):
            k = k + 1
    k_max = k - 1
    # // Goes up to P_1.
    if (k_base > k_top):
        k_top = k_max + 1
    MediumIndices_0 = k_base
    MediumIndices_1 = k_top
    k = k_base
    # // k_top == 0 means no cloud.
    while (k < k_top and k_top != 0):
        MediumCloudiness = MediumCloudiness + CloudFraction[k] * (GridZ_interf[k + 1] - GridZ_interf[k])
        k = k + 1
    if (k_top != 0):
        MediumCloudiness = MediumCloudiness / (GridZ_interf[k_top] - GridZ_interf[k_base])

    # /*** High clouds ***/

    cloud_max = 0;
    # // Starts above low clouds.
    for k in range(k_max + 1, Nz):
        cloud_max = max(cloud_max, CloudFraction[k])
    below = True; above = False
    k_base = 0; k_top = 0
    for k in range(k_max + 1, Nz):
        if (not above):
            below = below and ( CloudFraction[k] < 0.5 * cloud_max
                        	   or CloudFraction[k] == 0 )
            above = (not below) and CloudFraction[k] < 0.5 * cloud_max
        if ((not below) and k_base == 0):
            k_base = k
        if (above):
            k_top = k
    k_max = k - 1;
    # // Goes up to the top.
    if (k_base > k_top):
        k_top = k_max + 1
    HighIndices_0 = k_base
    HighIndices_1 = k_top
    k = k_base;
    # // k_top == 0 means no cloud.
    while (k < k_top and k_top != 0):
        HighCloudiness = HighCloudiness + CloudFraction(h, k, j, i) * (GridZ_interf[k + 1] - GridZ_interf[k])
        k = k + 1
    # if (k_top != 0):
    if (k_top != 0 and k_top != k_base): # YK
        HighCloudiness = HighCloudiness / (GridZ_interf[k_top] - GridZ_interf[k_base])

    return MediumCloudiness, HighCloudiness

def compute_cloud_fraction(rh, crh):
    tmp = rh - crh
    if (tmp < 0. and crh == 1.):
        CloudFraction = 0.;
    else:
        tmp = tmp / (1. - crh)
        CloudFraction = tmp * tmp
    return CloudFraction
 

def compute_crh(SurfacePressure, Pressure):
    coeff0 = 2.
    coeff1 = np.sqrt(3.)
    sig = Pressure / SurfacePressure
    crh = 1.0 - coeff0 * sig * (1.0 - sig) * (1.0 + (sig - 0.5) * coeff1)
    return crh

def compute_relative_humidity(SpecificHumidity, Temperature, Pressure):
    P_sat = 611.2 * np.exp(17.67 * (Temperature - 273.15)
                        / (Temperature - 29.65))
    RelativeHumidity = SpecificHumidity * Pressure / ( (0.62197 * (1.0 - SpecificHumidity) + SpecificHumidity ) * P_sat)
    return RelativeHumidity






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
            street.background['O3'] = o3
            street.background['NO2'] = no2
            street.background['NO'] = no

    return 0

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



def read_traffic_data(input_file, emis_species_list, epsg_code):

    #
    # Define the projections input/output
    # -----------------------------------
    import pyproj
    
    epsg_code = 'epsg:' + epsg_code
    proj_in = pyproj.Proj(init=epsg_code)
    # WGS84
    proj_out = pyproj.Proj(init='epsg:4326') 

    print("=================", input_file)
    node_id = 0
    street_id = 0
    input_street = open(input_file)
    header = input_street.readline()
    header = header.strip("\n")
    header_info = [x for x in re.split('\t| ', header) if x.strip() != '']

    species_ind = []
    for emis_species in emis_species_list:
        species_check = False
        for inds, species_name in enumerate(header_info):
            if emis_species == species_name:
                print(emis_species + " found with the index:", inds) 
                species_ind.append(inds)
                species_check = True
        if (species_check == False):
            sys.exit(emis_species + " is not found in " + input_file)
                
    node_list = []
    street_list = []
    for line in input_street.readlines():
        line = line.replace('\n','')
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        street_id = int(line_info[0])
        node_id = node_id + 1
        id_begin = node_id
        x = float(line_info[3])
        y = float(line_info[4])
        lon1, lat1 = pyproj.transform(proj_in, proj_out, x, y)
        node = Node(node_id, lon1, lat1)
        node.connected_street.append(street_id)
        node_list.append(node)
        node_id = node_id + 1
        id_end = node_id
        x = float(line_info[5])
        y = float(line_info[6])
        lon2, lat2 = pyproj.transform(proj_in, proj_out, x, y)
        node = Node(node_id, lon2, lat2)
        node.connected_street.append(street_id)
        node_list.append(node)

        # Check input data.
        if ((lon1 == lon2) and (lat1 == lat2)):
            sys.exit("Error: a street has two same intersection coordinate " + \
                     "for the node " + str(node_id) + ", lon: " + str(lon1) + \
                     ", lat: " + str(lat1))
        
        # Street length
        length = distance_on_unit_sphere(lat1, lon1, lat2, lon2) # in meter
        width = 0.0
        height = 0.0
        lon_cen = (lon1 + lon2) * 0.5
        lat_cen = (lat1 + lat2) * 0.5

        # Conversion of the unit of the emission input data
        # from ug/km/h to ug/s
        emission = np.zeros([len(species_ind)], 'float')
        for i, ind in enumerate(species_ind):
            emission[i] = float(line_info[species_ind[i]]) * (length / 1000.) / 3600.0

        street = Street(street_id, id_begin, id_end, length, width, height, lon_cen, lat_cen, emission)        

        street_list.append(street)
    print(" --------------------------")
    print(" --- Number of the nodes: ", len(node_list))
    print(" --- Number of the streets: ", len(street_list))
    print(" --------------------------")
    input_street.close()
    return street_list, node_list    

def write_output(node_list, street_list, node_list_eff, street_list_eff, current_date, output_dir, emis_species_list):

    str_date = current_date.strftime("%Y%m%d%H")
    date = str_date[0:8]
    hour = str_date[8:10]

    # Write 
    # street_ID, begin_node, end_node, -----
    file_emission = output_dir + "emission." + date + "-" + hour + ".txt"
    f = open(file_emission, 'w')

    for i in range(len(street_list_eff)):
        street = street_list_eff[i]
        street_id = street.id
        begin_node = street.eff_begin
        end_node = street.eff_end

        nemis = len(street.eff_emission)

        f.write(str(street_id) + "\t" + str(begin_node) + "\t" + 
                str(end_node) + "\t")        
        for s in range(nemis):
            f.write(str(float(street.eff_emission[s])) + "\t")
        f.write("\n")
    f.close()

    # Write
    # node_ID, longitude, latitude, number of connected segments, segment indices
    file_node = output_dir + "intersection.dat"
    f = open(file_node, 'w')
    header = "#id;lon;lat;number_of_streets;1st_street_id;2nd_street_id;...\n"
    f.write(header)
    for i in range(len(node_list_eff)):
        node = node_list_eff[i]
        if node.removed == False:
            lon = node.lon
            lat = node.lat
            node_id = node.id
            ns = len(node.connected_street)
            st_ = ""
            for st in node.connected_street:
                st_ += str(st) + ";"
            f.write(str(node_id) + ";" + str(lon) + ";" + str(lat) + 
                    ";" + str(ns) + ";" + st_ +  "\n")
    f.close()

    # Write
    file_node = output_dir + "street.dat"
    f = open(file_node, 'w')
    header = "#id;begin_inter;end_inter;length;width;height;typo\n"
    f.write(header)
    for i in range(len(street_list_eff)):
        street = street_list_eff[i]
        street_id = street.id
        begin_node = street.eff_begin
        end_node = street.eff_end
        length = street.length
        width = street.width
        height = street.height
        typo = street.typo
        f.write(str(street_id) + ";" + str(begin_node) + ";" + str(end_node) + 
                ";" + str(length) + ";" + str(width) + ";" + str(height) +
                ";" + str(typo) + "\n")
    f.close()
    
    # Write
    import os.path
    
    # Emission data
    emission_array = np.zeros([len(street_list_eff), len(emis_species_list)], 'float')

    wind_dir = np.zeros((len(street_list_eff)), 'float')
    wind_speed = np.zeros((len(street_list_eff)), 'float')
    pblh = np.zeros((len(street_list_eff)), 'float')
    ust = np.zeros((len(street_list_eff)), 'float')
    lmo = np.zeros((len(street_list_eff)), 'float')
    psfc = np.zeros((len(street_list_eff)), 'float')
    t2 = np.zeros((len(street_list_eff)), 'float')
    attenuation = np.zeros((len(street_list_eff)), 'float')
    sh = np.zeros((len(street_list_eff)), 'float')
    lwc = np.zeros((len(street_list_eff)), 'float')
    rain = np.zeros((len(street_list_eff)), 'float')
    background={}
    street0=street_list_eff[0]
    for spec in list(street0.background.keys()):
        background[spec]=np.zeros((len(street_list_eff)), 'float')

    for i in range(len(street_list_eff)):
        street = street_list_eff[i]

        emission_array[i] = street.eff_emission

        wind_dir[i] = street.wind_dir
        wind_speed[i] = street.wind_speed
        pblh[i] = street.pblh
        ust[i] = street.ust
        lmo[i] = street.lmo
        psfc[i] = street.psfc
        t2[i] = street.t2
        sh[i] = street.sh
        attenuation[i] = street.attenuation
        lwc[i] = street.lwc
        rain[i] = street.rain

        for spec in list(background.keys()):
            background[spec][i]=street.background[spec]

    wind_dir_inter = np.zeros((len(node_list_eff)), 'float')
    wind_speed_inter = np.zeros((len(node_list_eff)), 'float')
    pblh_inter = np.zeros((len(node_list_eff)), 'float')
    ust_inter = np.zeros((len(node_list_eff)), 'float')
    lmo_inter = np.zeros((len(node_list_eff)), 'float')
    for i in range(len(node_list_eff)):
        node = node_list_eff[i]
        wind_dir_inter[i] = node.wind_dir
        wind_speed_inter[i] = node.wind_speed
        pblh_inter[i] = node.pblh
        ust_inter[i] = node.ust
        lmo_inter[i] = node.lmo

    return wind_dir, wind_speed, pblh, ust, lmo, psfc, t2, sh, attenuation, background, wind_dir_inter, wind_speed_inter, pblh_inter, ust_inter, lmo_inter, emission_array, lwc, rain


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

def get_polair_ind_v2(lon, lat, x_min, Delta_x, y_min, Delta_y, Nx, Ny):
    # Get cell indices (X, Y) in the Polair3D grid
    # Xid = (street_list[i].lon_cen - x_min) / delta_x
    # Yid = (street_list[i].lat_cen - y_min) / delta_y
    # Xid = int(Xid)-1 if Xid%1 < 0.5 else int(Xid)
    # Yid = int(Yid)-1 if Yid%1 < 0.5 else int(Yid)

    index_x = max(int((lon - x_min + Delta_x / 2.) / Delta_x), 0)
    index_x = min(index_x, Nx - 1)
    index_y = max(int((lat - y_min + Delta_y / 2.) / Delta_y), 0)
    index_y = min(index_y, Ny - 1)

    return index_x, index_y

def is_holiday(date, country_code):
    # Requirement: holidays python library (https://pypi.org/project/holidays/)
    # Install: pip install holidays
    try:
        import holidays
        isCountryFound = False
        for country in holidays.list_supported_countries():
            if country_code == country:
                print(('Found country code "{}"'.format(country_code)))
                country_holidays = holidays.CountryHoliday(country_code)
                isCountryFound = True
                break
        
        if (isCountryFound == False):
            print(('Error: given country code "{}" is not found in supported country codes'.format(country_code)))
            print(holidays.list_supported_countries())
            return False
        
        return date.date() in country_holidays
        
    except ImportError as err:
        print('*** Could not import "holidays" Python module ***')
        print('Please install "holidays", or Holidays are not taken into account.')
        print('To install, type pip install holidays')
        return False
  
    


def utc_to_local(utc, zone = 'Europe/Paris'):
    # Need to install dateutil
    # pip install python-dateutil
    from dateutil import tz
    
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz(zone)

    # Tell the datetime object that it's in UTC time zone since 
    # datetime objects are 'naive' by default
    utc = utc.replace(tzinfo=from_zone)

    # Convert time zone
    return utc.astimezone(to_zone)


# def projection_conversion(epsg_code):
# #
# # Define the projections input/output
# # -----------------------------------
#     import pyproj
#     wgs84 = pyproj.Proj('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
#     lambert93 = pyproj.Proj('+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')


# Read background binary files from Polair3D output
def read_bkgd_bin(bkgd_species, indir, Nt, Nx, Ny, Nz):
    # Put data in a dict
    bkgd_data = {}
    for species in bkgd_species:
        infile = indir + species + '.bin'
        bkgd_data[species] = np.memmap(infile, dtype='float32', mode='r',
                                       shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]

    return bkgd_data


# def get_polair_id(lon, lat, x_min, y_min, dx, dy, Nx, Ny):
#     Xid = max(int((lon - x_min + dx / 2.) / dx), 0)
#     Xid = min(Xid, Nx - 1)
#     Yid = max(int((lat - y_min + dy / 2.) / dy), 0)
#     Yid = min(Yid, Ny - 1)

# Set background for streets
def set_bkgd_bin(street_list, bkgd_data, current_date, date_min, delta_t,
                  Nt, x_min, y_min, delta_x, delta_y, Nx, Ny):
    # Get index of current date
    c_id = int((current_date - date_min).total_seconds() / delta_t)
    if c_id < 0 or c_id >= Nt:
        sys.exit('ERROR: background data not available for this date.')

    for i in range(len(street_list)):

        Xid, Yid = get_polair_ind_v2(street_list[i].lon_cen, street_list[i].lat_cen,
                                     x_min, y_min, delta_x, delta_y, Nx, Ny)
        for key, value in bkgd_data.items():
            street_list[i].background[key] = value[c_id, Yid, Xid]



        
