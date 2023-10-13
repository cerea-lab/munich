#! /usr/bin/env python

import math, datetime, os
import numpy as np
from io_tools import *
from misc import *
from atmopy import talos
from scipy.io import netcdf

    
SB_constant = 5.67 * pow(10,-8) # (Stefan-Boltzmann constant)
ka = 0.4
g = 9.81
r = 287.

""" 
From Soulhac(2011)

# input 
  theta: Potential temperature perturbation
  psfc: Surface pressure
  T0: Skin temperature
  rain : Total precipitation 
  sr : Solar radiation
  cloud_fraction: Cloud fraction (range from 0 to 1)
  ust : Frictionn velocity
# output
  lmo : Monin-Obukhov length
"""
def compute_lmo(theta, psfc, T0,
                rain, sr, cloud_fraction, ust):

    theta = theta + 300.
    theta0 = T0 * pow((psfc / 101325.0), (-287.0 / 1005.0))
   
    theta_mean = 0.5 * (theta + theta0)

    # Fraction ragne from 0 to 1
    cloud_fraction_max = np.max(cloud_fraction[:])
    # Ragne conversion from 0 ~ 1 to 0 ~ 8                
    cloud_fraction_max = cloud_fraction_max * 8.0
    
    if rain >= 0.2:
        alpha = 1.0
    else:
        alpha = 0.0

    # Based on Van Ulden and Holtslag (1985)  
    s_r = math.exp(0.055 * (T0 - 279.0))


    ### Rn : Net radiation (W m-2)
    epsil = 0.88  #(Stull, 1988)
    albedo = 0.18
    Rn_1 = (1 - albedo) * sr 
    Rn_2 = (0.94 * pow(10,-5) * pow(T0, 2.0) - epsil) * \
        SB_constant * pow(T0, 4.0)
    Rn_3 =  60.0 * (cloud_fraction_max / 8.0)
    Rn   = (Rn_1 + Rn_2 + Rn_3) / 1.12
    
    ### H0 : Sensible heat ﬂux between ground and atmosphere (W m-2)
    H0 = (1 + (1 - alpha) * s_r) / (1 + s_r) * 0.7 * Rn - (20. * alpha)
    
    ### Lmo
    density = 1.2 # air density [kg/m3]
    cp = 1005.0   # Specific heat 
    lmo = (-1.0 * density * cp * pow(ust, 3.0) * theta_mean) / (ka * g * H0)

    return lmo

""" 
From the implemented version in Polyphemus

# input 
  theta: Potential temperature perturbation
  psfc: Surface pressure
  tsk: Skin temperature
  lh: Latent heat
  ust: Frictionn velocity
  hfx: Sensible heat
# output
  lmo : Monin-Obukhov length
"""
def compute_lmo_poly(theta, psfc, tsk, lh, ust, hfx):
    
    theta = theta + 300.
    theta0 = tsk * pow((psfc / 101325.0), (-287.0 / 1005.0))
    theta_mean = 0.5 * (theta + theta0)
    evaporation = lh / 2.5e9
    lmo = (-1.0 * pow(ust, 3.0) * theta_mean) / (ka * g * (hfx + 0.608 * theta_mean * evaporation))

    return lmo

""" 
LMO from RMOL of WRF

# input 
  rmol: Inverse of LMO
# output
  lmo : Monin-Obukhov length
"""
def compute_lmo_wrf(rmol):

    if rmol != 0 :
        lmo = 1. / rmol  ### rmol = 1./lmo
    else :
        lmo = 1.
      
    return lmo
    
def wrf_type_1(current_date,
               input_dir,
               wrfout_prefix,
               option_lmo):

    str_date = current_date.strftime("%Y-%m-%d")
    input_file = input_dir + wrfout_prefix + "_" + str_date + "_00:00:00"

    try:
        f = netcdf.netcdf_file(input_file, 'r')
        print(input_file + " is found.")
    except IOError:
        print(input_file + " is not found.")
        try:
            input_file = input_dir + wrfout_prefix + "_" + str_date
            f = netcdf.netcdf_file(input_file, 'r')
            print(input_file + " is found.")
        except IOError:
            print(input_file + " is not found.")
                
    start_date = f.__getattribute__("SIMULATION_START_DATE")
    print(str_date, start_date)

    previous_date = current_date - datetime.timedelta(seconds = 3600)
    str_prev_date = previous_date.strftime("%Y-%m-%d")
    previous_file = input_dir + wrfout_prefix + "_" + str_prev_date + "_00:00:00"
    f_prev = netcdf.netcdf_file(previous_file, 'r')
   

    times = f.variables["Times"][:]
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

    # ===
    # List of the required variables
    #
    # theta = f.variables["T"][:] # Potential temperature perturbation
    # tsk = f.variables["TSK"][:] # Skin temperature
    # psfc = f.variables["PSFC"][:] # Surface pressure
    # hfx = f.variables["HFX"][:] # Sensible heat
    # lh = f.variables["LH"][:] # Latent heat
    # sh = f.variables["QVAPOR"][:] # Specific humidity
    # t2 = f.variables["T2"][:] # Surface temperature
    # lwc = f.variables["QCLOUD"][:] # Cloud water mixing ratio
    # pressure_pert = f.variables["P"][:] # Pressure Perturbation
    # base_pres = f.variables["PB"][:] # Base state pressure
    # rainc = f.variables["RAINC"][:] # Convective Rain
    # rainnc = f.variables["RAINNC"][:] # Non-convective Rain
    # solar_radiation = f.variables["SWDOWN"][:] # Solar Radiation
    # ====
    variable_to_use = ["U10", "V10", "U", "V",
                       "PBLH", "UST", "T",
                       "TSK", "PSFC", "HFX", "LH",
                       "QVAPOR", "T2", "QCLOUD", "P", "PB",
                       "RAINC", "RAINNC", "SWDOWN", "CLDFRA",
                       "HGT", "ZNU"]        

    # Read the required variables and make a dictionary
    dict_data_netcdf_meteo = {}
    for var in variable_to_use:
        dict_data_netcdf_meteo[var] = f.variables[str(var)][ind_t]
   

    # lons = f.variables["XLONG"][:]
    # nx = lons.shape[2]
    # ny = lons.shape[1]
    
    rainc = f.variables["RAINC"][ind_t] # Convective Rain
    rainnc = f.variables["RAINNC"][ind_t] # Non-convective Rain
    rain = rainc + rainnc

    if (ind_t == 0 and os.path.isfile(previous_file)):
        rainc_prev = f_prev.variables["RAINC"][-1] # Convective Rain
        rainnc_prev = f_prev.variables["RAINNC"][-1] # Non-convective Rain
        rain_prev = (rainc_prev + rainnc_prev)
    elif (ind_t == 0 and not (os.path.isfile(previous_file))):
        rain_prev = 0.0
        print("Warning: File for the previous date is not available.")
    else:
        rainc_prev = f.variables["RAINC"][ind_t - 1]
        rainnc_prev = f.variables["RAINNC"][ind_t - 1]
        rain_prev = (rainc_prev + rainnc_prev)

    # Decumulate rain
    rain = np.clip(rain - rain_prev, 0, None)

    if option_lmo == "wrf":
        rmol = f.variables["RMOL"][ind_t]
    else:
        rmol = 0


    
    return input_file, rain, dict_data_netcdf_meteo, rmol

def wrf_type_2(current_date,
               input_file, t_ind,
               delta_t, time_step_wrf,
               option_lmo):

    from scipy.io import netcdf

    try:
        f = netcdf.netcdf_file(input_file, 'r')
        print(input_file + " is found.")
    except IOError:
        print(input_file + " is not found.")

    times = f.variables["time"][:]

    start_date = f.__getattribute__("SIMULATION_START_DATE")

    # Conversion byte to string
    start_date = start_date.decode()
    s_date = datetime.datetime.strptime(str(start_date), "%Y-%m-%d_%H:%M:%S")
    print(start_date, s_date)

    hasMeteo = False
    for t in range(len(times)):
        meteo_date = s_date + datetime.timedelta(seconds = time_step_wrf * t)
        if current_date == meteo_date:
            print("Meteo data are available for the date ", current_date)
            ind_t = t
            hasMeteo = True
    if (hasMeteo == False):
        print("Error: meteo data are not available")
        sys.exit()

    
    variable_to_use = ["U10", "V10",
                       "PBLH", "UST", "TSK", "PSFC", "HFX",
                       "LH", "QVAPOR", "T2", "T",
                       "QCLOUD", "P", "PB", "SWDOWN",
                       "SH2O", "COSZEN", "HGT", "U", "V", "CLDFRA",
                       "ZNU"
    ]
   
    dict_data_netcdf_meteo={}
    for var in variable_to_use:
        print("Read..." + var)
        dict_data_netcdf_meteo[var] = f.variables[str(var)][ind_t]
       

    rainc = f.variables["RAINC"][ind_t] # Convective Rain
    rainnc = f.variables["RAINNC"][ind_t] # Non-convective Rain
    rain = rainc + rainnc
    
    # Decumulate rain
    if (ind_t == 0):
        pass
    else:
        rainc_prev = f.variables["RAINC"][ind_t - 1] # Convective Rain
        rainnc_prev = f.variables["RAINNC"][ind_t - 1] # Non-convective Rain
        rain_prev = rainc_prev + rainnc_prev
        rain = np.clip(rain - rain_prev, 0, None)

    if option_lmo == "wrf":
        rmol = f.variables["RMOL"][ind_t]
    else:
        rmol = 0

    f.close()
    
    return input_file, rain, dict_data_netcdf_meteo, rmol


def get_meteo_data(input_dir, ind_t_munich, delta_t,
                   current_date, street_list, node_list,
                   wrf_config):

    content =  [("file_type", "[type]", "Int"), \
                ("wrfout_prefix", "[type]", "String"), \
                ("filename","[type]","String"), \
                ('time_step_wrf', '[type]', 'Float'), \
                ('option_lmo', '[option]', 'String'), \
                ('option_surface_wind', '[option]', 'String')
    ]

    print("WRF configuration file: ", wrf_config)
    config = talos.Config(wrf_config, content)
    option_lmo = config.option_lmo

    # === WRF Type 1 in meteo.cfg ===
    if config.file_type == 1:
        input_file, rain, dict_data_netcdf_meteo, rmol = \
            wrf_type_1(current_date,
                       input_dir,
                       config.wrfout_prefix,
                       config.option_lmo)
        
        f = netcdf.netcdf_file(input_file, 'r')
        
        lons = f.variables["XLONG"][0]
        lats = f.variables["XLAT"][0]
        
    # === WRF Type 2 in meteo.cfg ===
    elif config.file_type == 2:
        input_file, rain, dict_data_netcdf_meteo, rmol = \
            wrf_type_2(current_date,
                       config.filename, ind_t_munich,
                       delta_t,
                       config.time_step_wrf,
                       config.option_lmo)

        f = netcdf.netcdf_file(input_file, 'r')

        lons = f.variables["XLONG"]
        lats = f.variables["XLAT"]

    # Get the dimensions
    nx = f.dimensions["west_east"]
    ny = f.dimensions["south_north"]
    nz = f.dimensions["bottom_top"]
    print("WRF output data (nz, ny, nx): ", nz, ny, nx)
       
    # Get meteo data for the streets
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

        if config.option_surface_wind == "10m":
            u_cell = dict_data_netcdf_meteo["U10"][ind_j, ind_i]
            v_cell = dict_data_netcdf_meteo["V10"][ind_j, ind_i]
        elif config.file_type == "1st_layer":
            u_cell = dict_data_netcdf_meteo["U"][0, ind_j, ind_i]
            v_cell = dict_data_netcdf_meteo["V"][0, ind_j, ind_i]
            
        pblh_cell = dict_data_netcdf_meteo["PBLH"][ind_j, ind_i]
        ust_cell = dict_data_netcdf_meteo["UST"][ind_j, ind_i]
        sh_cell = dict_data_netcdf_meteo["QVAPOR"][0, ind_j, ind_i]
        psfc_cell = dict_data_netcdf_meteo["PSFC"][ind_j, ind_i]
        t2_cell = dict_data_netcdf_meteo["T2"][ind_j, ind_i]

        # For chemistry
        lwc_cell = dict_data_netcdf_meteo["QCLOUD"][0, ind_j, ind_i]
        lwc_colon = dict_data_netcdf_meteo["QCLOUD"][:, ind_j, ind_i]
        solar_radiation_cell = dict_data_netcdf_meteo["SWDOWN"][ind_j, ind_i]
        cloud_fraction_cell = dict_data_netcdf_meteo["CLDFRA"][:, ind_j, ind_i]
        theta_cell = dict_data_netcdf_meteo["T"][0, ind_j, ind_i]
        tsk_cell = dict_data_netcdf_meteo["TSK"][ind_j, ind_i]
        lh_cell = dict_data_netcdf_meteo["LH"][ind_j, ind_i]
        hfx_cell = dict_data_netcdf_meteo["HFX"][ind_j, ind_i]
      
        rain_cell = rain[ind_j, ind_i]
        
        # Compute wind direction
        temp = math.atan2(v_cell, u_cell) # -pi < temp < pi
        if temp <= (math.pi / 2.0):
            temp2 = math.pi / 2.0 - temp
        else:
            temp2 = (math.pi / 2.0 - temp) + 2.0 * math.pi
        street.wdir = temp2 # in radian, 0 for the wind to north, pi/2 to east, pi to south.

        # Compute heights
        ptop = f.variables["P_TOP"][0]
        Tiso = f.variables["TISO"][0]
        p00 = f.variables["P00"][0]
        Ts0 = f.variables["T00"][0]
        A = f.variables["TLP"][0]
        # Terrain = f.variables["HGT"][:]
        # sigma_h = f.variables["ZNU"][:]

        Piso =  p00 * np.exp((Tiso - Ts0) / A)
        aux = np.log(Piso / p00)
        Ziso = - r * A * aux*aux / (2.0*g) - r * Ts0 * aux / g
        
        terrain_cell = dict_data_netcdf_meteo["HGT"][ind_j, ind_i]
        sigma_h_colon = dict_data_netcdf_meteo["ZNU"][:]
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
        pressure_pert = dict_data_netcdf_meteo["P"][:, ind_j, ind_i]
        base_pres = dict_data_netcdf_meteo["PB"][:, ind_j, ind_i]
        pressure_colon = pressure_pert + base_pres
        theta_colon = dict_data_netcdf_meteo["T"][:, ind_j, ind_i] + 300.
        sh_colon = dict_data_netcdf_meteo["QVAPOR"][:, ind_j, ind_i]
        rh_colon = np.zeros([nz], 'float')
        crh_colon = np.zeros([nz], 'float')
        cloud_fraction_computed = np.zeros([nz], 'float')
        for k in range(nz):
            rh_colon[k] = compute_relative_humidity(sh_colon[k], 
                                                    theta_colon[k], 
                                                    pressure_colon[k])
            crh_colon[k] = compute_crh(psfc_cell, pressure_colon[k])
            cloud_fraction_computed[k] = compute_cloud_fraction(rh_colon[k],
                                                       crh_colon[k])

        medium_cloudiness, high_cloudiness = compute_cloudiness(nz, 
                                                                pressure_colon,
                                                                cloud_fraction_computed,
                                                                gridz_interf)

        attenuation = np.zeros([nz], 'float')
        compute_attenuation(nz, gridz_in, rh_colon,
                            medium_cloudiness, high_cloudiness, attenuation)
        
        #
        street.wspd = math.sqrt(u_cell ** 2 + v_cell ** 2)
        street.pblh = pblh_cell
        street.ust = ust_cell
        street.spechumid = sh_colon[0]
        street.surfpressure = psfc_cell
        street.surftemp = t2_cell
        street.attenuation = attenuation[0]

        # Compute LMO

#        theta_cell = theta[ind_t, 0, ind_j, ind_i] + 300.

        # === LMO from Polyphemus
        if (option_lmo == "polyphemus"):
            lmo = compute_lmo_poly(theta_cell,
                                   psfc_cell,
                                   tsk_cell,
                                   lh_cell,
                                   ust_cell,
                                   hfx_cell)
        
        elif (option_lmo == "sirane"):
            lmo = compute_lmo(theta_cell, psfc_cell, tsk_cell,
                              rain_cell, solar_radiation_cell,
                              cloud_fraction_cell, ust_cell)

        elif (option_lmo == "wrf"):
            rmol_cell = rmol[ind_j, ind_i]
            lmo = compute_lmo_wrf(rmol)
        else:
            sys.exit("Wrong option name for LMO")
    
        street.lmo = lmo
        
        street.liquidwaterc = lwc_cell
        street.rain = rain_cell
        
    # Get meteo data for the intersections.    
    for n in range(len(node_list)):
        node = node_list[n]
        lat1 = node.lat
        lon1 = node.lon
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

        if config.option_surface_wind == "10m":
            u_cell = dict_data_netcdf_meteo["U10"][ind_j, ind_i]
            v_cell = dict_data_netcdf_meteo["V10"][ind_j, ind_i]
        elif config.option_surface_wind == "1st_layer":
            u_cell = dict_data_netcdf_meteo["U"][0, ind_j, ind_i]
            v_cell = dict_data_netcdf_meteo["V"][0, ind_j, ind_i]

        pblh_cell = dict_data_netcdf_meteo["PBLH"][ind_j, ind_i]
        ust_cell = dict_data_netcdf_meteo["UST"][ind_j, ind_i]

        # Compute LMO
        # theta_cell = theta[ind_t, 0, ind_j, ind_i] + 300.
        # theta0_cell = tsk[ind_t, ind_j, ind_i] * pow((psfc[ind_t, ind_j, ind_i] / 101325.0), (-287.0 / 1005.0))
        # theta_mean = 0.5 * (theta_cell + theta0_cell)
        # evaporation = lh[ind_t, ind_j, ind_i] / 2.5e9
        # hfx_cell = hfx[ind_t, ind_j, ind_i]
        # lmo = (-1.0 * pow(ust_cell, 3.0) * theta_mean) / (ka * g * (hfx_cell + 0.608 * theta_mean * evaporation))

        theta_cell = dict_data_netcdf_meteo["T"][0, ind_j, ind_i]
        psfc_cell = dict_data_netcdf_meteo["PSFC"][ind_j, ind_i]
        tsk_cell = dict_data_netcdf_meteo["TSK"][ind_j, ind_i]
        rain_cell = rain[ind_j, ind_i]
        solar_radiation_cell = dict_data_netcdf_meteo["SWDOWN"][ind_j, ind_i]
        cloud_fraction_cell = dict_data_netcdf_meteo["CLDFRA"][:, ind_j, ind_i]
        ust_cell = dict_data_netcdf_meteo["UST"][ind_j, ind_i]
        
        # theta_cell = theta[ind_t, 0, ind_j, ind_i]
        # psfc_cell = psfc[ind_t, ind_j, ind_i]
        # tsk_cell = tsk[ind_t, ind_j, ind_i]
        # rain_cell = rain[ind_j, ind_i]
        # solar_radiation_cell = solar_radiation[ind_t, ind_j, ind_i]
        # cloud_fraction_cell = cloud_fraction[ind_t, :, ind_j, ind_i]
        # ust_cell = ust[ind_t, ind_j, ind_i]

        if (option_lmo == "polyphemus"):
            lmo = compute_lmo_poly(theta_cell,
                                   psfc_cell,
                                   tsk_cell,
                                   lh_cell,
                                   ust_cell,
                                   hfx_cell)
        elif (option_lmo == "sirane"):
            lmo = compute_lmo(theta_cell, psfc_cell, tsk_cell,
                              rain_cell, solar_radiation_cell,
                              cloud_fraction_cell, ust_cell)
        elif (option_lmo == "wrf"):
            rmol_cell = rmol[ind_j, ind_i]
            lmo = compute_lmo_wrf(rmol)
        else:
            sys.exit("Wrong option name for LMO")
        
        temp = math.atan2(v_cell, u_cell) # -pi < temp < pi
        if temp <= (math.pi / 2.0):
            temp2 = math.pi / 2.0 - temp
        else:
            temp2 = (math.pi / 2.0 - temp) + 2.0 * math.pi
        node.wdir = temp2 # in radian, 0 for the wind to north, pi to south.
        node.wspd = math.sqrt(u_cell ** 2 + v_cell ** 2)
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



# Append meteo to binary files
def append_meteo_data(node_list, street_list, street_files, node_files):
    # Street
    street_data = [np.zeros((len(street_list)), 'float') for i in range(len(street_files))]
    for i in range(len(street_list)):
        street_data[0][i] = street_list[i].attenuation
        street_data[1][i] = street_list[i].liquidwaterc
        street_data[2][i] = street_list[i].lmo
        street_data[3][i] = street_list[i].pardb
        street_data[4][i] = street_list[i].pardiff
        street_data[5][i] = street_list[i].pblh
        street_data[6][i] = street_list[i].rain
        street_data[7][i] = street_list[i].soilwater
        street_data[8][i] = street_list[i].solarrad
        street_data[9][i] = street_list[i].spechumid
        street_data[10][i] = street_list[i].surfpressure
        street_data[11][i] = street_list[i].surfrichard
        street_data[12][i] = street_list[i].surftemp
        street_data[13][i] = street_list[i].ust
        street_data[14][i] = street_list[i].wdir
        street_data[15][i] = street_list[i].wspd
    for i in range(len(street_files)):
        append_binary(street_data[i], street_files[i])

    # Node
    node_data = [np.zeros((len(node_list)), 'float') for i in range(len(node_files))]
    for i in range(len(node_list)):
        node_data[0][i] = node_list[i].lmo
        node_data[1][i] = node_list[i].pblh
        node_data[2][i] = node_list[i].ust
        node_data[3][i] = node_list[i].wdir
        node_data[4][i] = node_list[i].wspd
    for i in range(len(node_files)):
        append_binary(node_data[i], node_files[i])


# Set meteo for streets and nodes
def set_meteo_bin(node_list, street_list, meteo_data, current_date, date_min, delta_t,
                  Nt, x_min, y_min, delta_x, delta_y, Nx, Ny):
    # Get index of current date
    c_id = int((current_date - date_min).total_seconds() / delta_t)
    if c_id < 0 or c_id >= Nt:
        sys.exit('ERROR: meteo data not available for this date.')

    # Set meteo for the streets
    for i in range(len(street_list)):
        # Get cell indices (X, Y) in the Polair3D grid
        # Xid = (street_list[i].lon_cen - x_min) / delta_x
        # Yid = (street_list[i].lat_cen - y_min) / delta_y
        # Xid = int(Xid)-1 if Xid%1 < 0.5 else int(Xid)
        # Yid = int(Yid)-1 if Yid%1 < 0.5 else int(Yid)
        Xid, Yid = get_polair_id(street_list[i].lon_cen, street_list[i].lat_cen,
                                 x_min, y_min, delta_x, delta_y, Nx, Ny)
        street_list[i].attenuation = meteo_data['attenuation'][c_id, Yid, Xid]
        street_list[i].liquidwaterc = meteo_data['liquidwaterc'][c_id, Yid, Xid]
        street_list[i].lmo = meteo_data['lmo'][c_id, Yid, Xid]
        street_list[i].pardb = meteo_data['pardb'][c_id, Yid, Xid]
        street_list[i].pardiff = meteo_data['pardiff'][c_id, Yid, Xid]
        street_list[i].pblh = meteo_data['pblh'][c_id, Yid, Xid]
        street_list[i].rain = meteo_data['rain'][c_id, Yid, Xid]
        street_list[i].soilwater = meteo_data['soilwater'][c_id, Yid, Xid]
        street_list[i].solarrad = meteo_data['solarrad'][c_id, Yid, Xid]
        street_list[i].spechumid = meteo_data['spechumid'][c_id, Yid, Xid]
        street_list[i].surfpressure = meteo_data['surfpressure'][c_id, Yid, Xid]
        street_list[i].surfrichard = meteo_data['surfrichard'][c_id, Yid, Xid]
        street_list[i].surftemp = meteo_data['surftemp'][c_id, Yid, Xid]
        street_list[i].ust = meteo_data['ust'][c_id, Yid, Xid]
        u = meteo_data['zonalwind'][c_id, Yid, Xid]
        v = meteo_data['meridiowind'][c_id, Yid, Xid]
        street_list[i].wdir = compute_wdir(u, v)
        street_list[i].wspd = meteo_data['wspd'][c_id, Yid, Xid]

    # Set meteo for the nodes
    for i in range(len(node_list)):
        # Get cell indices (X, Y) in the Polair3D grid
        # Xid = (node_list[i].lon - x_min) / delta_x
        # Yid = (node_list[i].lat - y_min) / delta_y
        # Xid = int(Xid)-1 if Xid%1 < 0.5 else int(Xid)
        # Yid = int(Yid)-1 if Yid%1 < 0.5 else int(Yid)
        Xid, Yid = get_polair_id(node_list[i].lon, node_list[i].lat,
                                 x_min, y_min, delta_x, delta_y, Nx, Ny)
        node_list[i].lmo = meteo_data['lmo'][c_id, Yid, Xid]
        node_list[i].pblh = meteo_data['pblh'][c_id, Yid, Xid]
        node_list[i].ust = meteo_data['ust'][c_id, Yid, Xid]
        u = meteo_data['zonalwind'][c_id, Yid, Xid]
        v = meteo_data['meridiowind'][c_id, Yid, Xid]
        node_list[i].wdir = compute_wdir(u, v)
        node_list[i].wspd = meteo_data['wspd'][c_id, Yid, Xid]


# Compute wind direction from zonal (u) and meridional (v) data
def compute_wdir(u, v):
    wdir = np.arctan2(v, u)
    hpi = np.pi / 2. # Half of Pi
    wdir = np.where(wdir<=hpi, hpi-wdir, (hpi-wdir)+2.*np.pi)
    return wdir

# Read meteo binary files from Polair3D input
def read_meteo_bin(indir, Nt, Nx, Ny, Nz):
    # Put data in a dict
    meteo_data = {}
    infile = indir + 'Attenuation.bin'
    meteo_data['attenuation'] = np.memmap(infile, dtype='float32', mode='r',
                                          shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]
    infile = indir + 'LiquidWaterContent.bin'
    meteo_data['liquidwaterc'] = np.memmap(infile, dtype='float32', mode='r',
                                           shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]
    infile = indir + 'LMO.bin'
    meteo_data['lmo'] = np.memmap(infile, dtype='float32', mode='r',
                                  shape=(Nt, Ny, Nx))
    infile = indir + 'PARdb.bin'
    meteo_data['pardb'] = np.memmap(infile, dtype='float32', mode='r',
                                    shape=(Nt, Ny, Nx))
    infile = indir + 'PARdiff.bin'
    meteo_data['pardiff'] = np.memmap(infile, dtype='float32', mode='r',
                                      shape=(Nt, Ny, Nx))
    infile = indir + 'BoundaryHeight.bin'
    meteo_data['pblh'] = np.memmap(infile, dtype='float32', mode='r',
                                   shape=(Nt, Ny, Nx))
    infile = indir + 'Rain.bin'
    meteo_data['rain'] = np.memmap(infile, dtype='float32', mode='r',
                                   shape=(Nt, Ny, Nx))
    infile = indir + 'SoilWater.bin'
    meteo_data['soilwater'] = np.memmap(infile, dtype='float32', mode='r',
                                        shape=(Nt, Ny, Nx))
    infile = indir + 'SolarRadiation.bin'
    meteo_data['solarrad'] = np.memmap(infile, dtype='float32', mode='r',
                                       shape=(Nt, Ny, Nx))
    infile = indir + 'SpecificHumidity.bin'
    meteo_data['spechumid'] = np.memmap(infile, dtype='float32', mode='r',
                                        shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]
    infile = indir + 'SurfacePressure.bin'
    meteo_data['surfpressure'] = np.memmap(infile, dtype='float32', mode='r',
                                           shape=(Nt, Ny, Nx))
    infile = indir + 'SurfaceRichardson.bin'
    meteo_data['surfrichard'] = np.memmap(infile, dtype='float32', mode='r',
                                          shape=(Nt, Ny, Nx))
    infile = indir + 'SurfaceTemperature.bin'
    meteo_data['surftemp'] = np.memmap(infile, dtype='float32', mode='r',
                                       shape=(Nt, Ny, Nx))
    infile = indir + 'FrictionModule.bin'
    meteo_data['ust'] = np.memmap(infile, dtype='float32', mode='r',
                                  shape=(Nt, Ny, Nx))
    infile = indir + 'WindModule.bin'
    meteo_data['wspd'] = np.memmap(infile, dtype='float32', mode='r',
                                   shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]
    infile = indir + 'ZonalWind.bin'
    meteo_data['zonalwind'] = np.memmap(infile, dtype='float32', mode='r',
                                        shape=(Nt, Nz, Ny, Nx+1))[:, 0, :, :]
    infile = indir + 'MeridionalWind.bin'
    meteo_data['meridiowind'] = np.memmap(infile, dtype='float32', mode='r',
                                          shape=(Nt, Nz, Ny+1, Nx))[:, 0, :, :]

    return meteo_data
