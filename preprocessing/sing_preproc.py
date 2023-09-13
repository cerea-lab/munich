#! /usr/bin/env python

import numpy as np
import sys, datetime, math, os
from atmopy import *
from street_network import *
from optparse import OptionParser
from io_tools import *

parser = OptionParser(usage = "%prog configuration_file")
(options, args) = parser.parse_args()

if not args:
	parser.error("A configuration file is required.")

content = [("emission_dir_weekday", "[input]", "String"), \
           ("emission_dir_weekend", "[input]", "String"), \
           ("epsg_code", "[input]", "String"), \
           ("country_code", "[input]", "String"), \
           ("weekday_file_prefix", "[input]", "String"), \
           ("weekend_file_prefix", "[input]", "String"), \
           ("is_local_hour", "[input]", "Bool"), \
           ("time_zone", "[input]", "String"), \
           ("emission_species", "[input]", "StringList"), \
           ("geog_info", "[input]", "String"), \
           ("background_concentration", "[background]", "String"), \
           ("meteo_dir", "[input]", "String"), \
           ("wrfout_prefix", "[input]", "String"),\
           ("Output_dir", "[output]", "String"), \
           ("Date_min_polair", "[domain]", "DateTime"), \
           ("Delta_t_polair", "[domain]", "Float"), \
           ("Nt_polair", "[domain]", "Int"), \
           ("is_street_merged", "[option]", "Bool"), \
           ("is_street_manually_merged", "[option]", "Bool"), \
           ("is_near_node_merged", "[option]", "Bool"), \
           ("is_node_manually_merged", "[option]", "Bool"), \
           ("is_voc_speciated", "[option]", "Bool"), \
           ("is_nox_speciated", "[option]", "Bool"), \
           ("is_isvoc_speciated", "[option]", "Bool"), \
           ("is_pm10_speciated", "[option]", "Bool"), \
           ("Nsize_sections", "[option]", "Int"), \
           ("Size_dist_ec_om_emis", "[option]", "FloatList"), \
           ("Size_dist_dust_emis", "[option]", "FloatList"), \
           ("om_redist", "[option]", "String"), \
           ("option_background","[background]","Int"), \
           # background from Chimere           
           ("chimere_dir", "[background]", "String"), \
           ("chimout_lab", "[background]", "String"), \
           ("chimere_species","[background]","StringList"), \
           ("melch2molmass_file","[background]","String"), \
           # background from Polair3d
           ('polair3d_dir', '[background]', 'String'),
           ("polair3d_species","[background]","StringList"), \
           ('date_min_bkgd', '[background]', 'DateTime'),
           ('delta_t_bkgd', '[background]', 'Float'),
           ('Nt_bkgd', '[background]', 'Int'),
           ('x_min_bkgd', '[background]', 'Float'),
           ('y_min_bkgd', '[background]', 'Float'),
           ('delta_x_bkgd', '[background]', 'Float'),
           ('delta_y_bkgd', '[background]', 'Float'),
           ('Nx_bkgd', '[background]', 'Int'),
           ('Ny_bkgd', '[background]', 'Int'),
           ('Nz_bkgd', '[background]', 'Int'))
]
config = talos.Config(sys.argv[1], content)

################
# Main program #
################


import shutil
output_dir = config.Output_dir
shutil.rmtree(output_dir, ignore_errors=True)
os.makedirs(output_dir)
emis_dir = output_dir + "/emission/"
os.makedirs(emis_dir)
meteo_dir = output_dir + "/meteo/"
os.makedirs(meteo_dir)
background_dir = output_dir + "/background/"
os.makedirs(background_dir)
textfile_dir = output_dir + "/textfile/"
os.makedirs(textfile_dir)
emis_dir = output_dir + "/grid_emission/"
os.makedirs(emis_dir)
        
emis_species_list = config.emission_species
print(emis_species_list)


# Remove old look-up table.
lut_file = "street-merging-lookup-table.txt"
if (os.path.isfile(lut_file)):
    os.remove(lut_file)

# Number of emitted species
ns_emis = len(emis_species_list)

#LL-------
#Check options speciation
if config.is_voc_speciated:
        isFound = False
        for sp_emis in emis_species_list:
                if sp_emis == "NMHC":
                        isFound = True
                        break
        if isFound == False:
                print("The option is_voc_speciated is activated. To perform the speciation of VOC species, please add NMHC in the emission_species.")
                sys.exit()

if config.is_nox_speciated:
        isFound = False
        for sp_emis in emis_species_list:
                if sp_emis == "NOx":
                        isFound = True
                        break
        if isFound == False:
                print("The option is_nox_speciated is activated. To perform the speciation of NOx species, please add NOx in the emission_species.")
                sys.exit()
                
if config.is_pm10_speciated:
        #PM10
        isFound = False
        for sp_emis in emis_species_list:
                if sp_emis == "PM10":
                        isFound = True
                        break
        if isFound == False:
                print("The option is_pm10_speciated is activated. To perform the speciation of PM10 species, please add PM10 in the emission_species.")
                sys.exit()
        #EC
        isFound = False
        for sp_emis in emis_species_list:
                if sp_emis == "EC":
                        isFound = True
                        break
        if isFound == False:
                print("The option is_pm10_speciated is activated. To perform the speciation of PM10 species, please add EC in the emission_species.")
                sys.exit()
        #OM
        isFound = False
        for sp_emis in emis_species_list:
                if sp_emis == "OM":
                        isFound = True
                        break
        if isFound == False:
                print("The option is_pm10_speciated is activated. To perform the speciation of PM10 species, please add OM in the emission_species.")
                sys.exit()

if config.is_isvoc_speciated:
        isFound = False
        for sp_emis in emis_species_list:
                if sp_emis == "OM":
                        isFound = True
                        break
        if isFound == False:
                print("The option is_isvoc_speciated is activated. To perform the speciation of ISVOC species, please add OM in the emission_species.")
                sys.exit()                
#---------

background_all={}

# Make input files for SinG simulations.
x_min = config.x_min
y_min = config.y_min
Delta_x = config.Delta_x
Delta_y = config.Delta_y
Nx = config.Nx
Ny = config.Ny
x_max = x_min + Delta_x * (Nx - 1) + Delta_x * 0.5
y_max = y_min + Delta_y * (Ny - 1) + Delta_y * 0.5
polair_lon = np.arange(x_min, x_max, Delta_x)
polair_lat = np.arange(y_min, y_max, Delta_y)

emission = np.zeros([ns_emis, 1, Ny, Nx], 'float')

earth_radius_2 = 6371229. * 6371229.
pi = 3.14159265358979323846264


##################################
#                                #
#    IF BACKGROUND FROM Polair3d #
#                                #
##################################
if (config.option_background == 2):

    # Output binary files to append to for streets
    bkgd_files = {x: background_dir + x + '.bin' for x in config.polair3d_species}    

    # Load files using memmap and return a dict
    bkgd_data = read_bkgd_bin(config.polair3d_species, config.polair3d_dir, config.Nt_bkgd,
                              config.Nx_bkgd, config.Ny_bkgd, config.Nz_bkgd)


#################################
#                               #
#    IF BACKGROUND FROM CHIMERE #
#                               #
#################################
elif (config.option_background == 3):

    melchior_spec_list=config.chimere_species
    chimout_dir=config.chimere_dir
    chimout_lab=config.chimout_lab
    
    molar_mass_melchior2={}

    with open(config.melch2molmass_file,'r') as f:
        for ln in f:
            spec=ln.split()[0]
            val=np.float(ln.split()[1])
            molar_mass_melchior2[spec]=val
    f.close()

# Output file names
file_wind_dir = output_dir + "/meteo/WindDirection.bin"
file_wind_speed = output_dir + "/meteo/WindSpeed.bin"
file_pblh = output_dir + "/meteo/PBLH.bin"
file_ust = output_dir + "/meteo/UST.bin"
file_lmo = output_dir + "/meteo/LMO.bin"
file_psfc = output_dir + "/meteo/SurfacePressure.bin"
file_t2 = output_dir + "/meteo/SurfaceTemperature.bin"
file_sh = output_dir + "/meteo/SpecificHumidity.bin"
file_attenuation = output_dir + "/meteo/Attenuation.bin"
file_lwc = output_dir + "/meteo/LiquidWaterContent.bin"
file_rain = output_dir + "/meteo/Rain.bin"


file_wind_dir_inter = output_dir + "/meteo/WindDirectionInter.bin"
file_wind_speed_inter = output_dir + "/meteo/WindSpeedInter.bin"
file_pblh_inter = output_dir + "/meteo/PBLHInter.bin"
file_ust_inter = output_dir + "/meteo/USTInter.bin"
file_lmo_inter = output_dir + "/meteo/LMOInter.bin"


# Get date info
begin_date = config.t_min
delta_t = config.Delta_t # in hr
nt = config.Nt
read_lut = False
date_list = []

for t in range(nt):
    current_date = begin_date + datetime.timedelta(hours = (delta_t * t))
    print("\n=====================================================")
    print("Current date (UTC): ", current_date)
    date_list.append(current_date)

    # Set time index for polair
    print(config.Date_min_polair)
    time_diff = current_date - config.Date_min_polair
    time_diff_seconds = time_diff.days * 24 * 60 * 60 + time_diff.seconds
    indt = int(time_diff_seconds / 3600)
    print("Time index: ", indt) 

    if config.is_local_hour:
            print(('Conversion from {} to UTC'.format(config.time_zone)))
            current_date_local = utc_to_local(current_date, config.time_zone)
    else:
            current_date_local = current_date

    print("Current date (local hour): ", current_date_local)
    str_date = current_date_local.strftime("%Y%m%d%H")
    date = str_date[0:8]
    hour = str_date[8:10]
    # Should take into account the French holidays during the period.
    if (current_date.weekday() >= 0) and (current_date.weekday() <= 4):
            if is_holiday(current_date, config.country_code):
                    print("The current day is a holiday.")
                    input_file = config.emission_dir_weekend + "/" + \
                    config.weekend_file_prefix + hour
            else:
                    print("The current day is a weekday.")
                    input_file = config.emission_dir_weekday + "/" + \
                    config.weekday_file_prefix + hour
    else:
            print("The current day is a weekend.")
            input_file = config.emission_dir_weekend + "/" + \
            config.weekend_file_prefix + hour

    print("Read the input data (segment coordinates and emission rates) from the file --- ",input_file)
    
    street_list, node_list = read_traffic_data(input_file,
                                               emis_species_list, config.epsg_code)

    # Get street geographical informations: length, width, builiding height.
    geog_info = config.geog_info
    get_street_geog(geog_info, street_list)


    # Merging streets if they have same (or very near) nodes.
    if (config.is_street_merged):
            n_street = len(street_list)
            print(" - Initial number of streets: ", n_street)

            # Automatic merging for the separated roads
            lut_file = "street-merging-lookup-table.txt"
            read_lut = os.path.isfile(lut_file)
            if (read_lut):
                    ntemp = lut_merging_street(lut_file, street_list)
            else:
                    ntemp = merging_street(lut_file, node_list, street_list)
                    print(" - Number of streets after merging the same streets: %d" % (n_street - ntemp))

            # Manual merging for the separated roads.
            if (config.is_street_manually_merged):
                    ntemp2 = manual_merging_street(street_list)
                    print(" - Number of streets after manual merging of the streets: %d" % (n_street - ntemp - ntemp2))

    # Make a new node list
    street_list_eff = []
    for ist in range(len(street_list)):
        street = street_list[ist]
        begin_node = street.eff_begin
        end_node = street.eff_end
        if (street.removed == False):
            street_list_eff.append(street)
            for inode in range(len(node_list)):
                node = node_list[inode]
                if node.id == begin_node or node.id == end_node:
                    node.eff_id = node.id

    node_list_temp = []
    for inode in range(len(node_list)):
        node = node_list[inode]
        if node.id == node.eff_id:
            node_list_temp.append(node)

    # Merging nodes if they have same coordinates or they are very near.
    print(" === Merging nodes === ")
    n_node = len(node_list_temp)
    print("Initial number of nodes: ", n_node)
    n_node1 = merging_node(node_list_temp)
    print("Number of nodes after removing the same nodes: ", (n_node - n_node1))
    if config.is_near_node_merged:
            n_node2 = merging_near_node(node_list_temp)
            print("Number of nodes after removing the nearest nodes: ", (n_node - n_node1 - n_node2))
    if config.is_node_manually_merged:
            n_node3 = manual_merging_node(node_list_temp)
            print("Number of nodes after manually removing the nearest nodes: ", (n_node - n_node1 - n_node2 - n_node3))

    node_list_eff = []
    for inode in range(len(node_list_temp)):
        node = node_list_temp[inode]
        if node.id == node.eff_id:
            node_list_eff.append(node)

    for inode in range(len(node_list)):
        node = node_list[inode]
        for j in range(len(street_list_eff)):
            street = street_list_eff[j]
            if street.begin == node.id:
                street.eff_begin = node.eff_id
            if street.end == node.id:
                street.eff_end = node.eff_id

    outfile = config.Output_dir + 'textfile/street_all.csv'
    with open(outfile, 'w') as f:
        f.write('id' + ',' + 'eff_id' + ',' + 'node_begin'\
            + ',' + 'node_eff_begin' + ',' + 'node_end' + ','\
            + 'node_eff_end' + '\n')

        for street in street_list:
            f.write(str(street.id) + ',' + str(street.eff_id) + ',' + str(street.begin)\
                    + ',' + str(street.eff_begin) + ',' + str(street.end) + ','\
                    + str(street.eff_end) + '\n')
    get_meteo_data(config.meteo_dir, current_date, \
                   street_list_eff, node_list_eff, config.wrfout_prefix)

    ### Backgrond concentration ###

    # Text
    if (config.option_background == 1):
        background_concentration_file = config.background_concentration
        get_background_concentration(background_concentration_file, \
                                         current_date, street_list_eff)
    # From Polair3d output
    elif (config.option_background == 2):
        set_bkgd_bin(street_list_eff, bkgd_data, current_date, config.date_min_bkgd,
                     config.delta_t_bkgd, config.Nt_bkgd, config.x_min_bkgd,
                     config.y_min_bkgd, config.delta_x_bkgd, config.delta_y_bkgd,
                     config.Nx_bkgd, config.Ny_bkgd)
        
        
    # From Chimere output
    elif (config.option_background == 3):
        get_chimere_background_concentration(current_date, street_list_eff, \
                                             melchior_spec_list, \
                                             molar_mass_melchior2, \
                                             chimere_dir,chimout_lab)

    
    # Compute emissions in grid cells.
    emission.fill(0.0)
    for st in street_list_eff:
            indx, indy = get_polair_ind_v2(st.lon_cen, st.lat_cen, x_min, Delta_x, y_min, Delta_y, Nx, Ny)
            lat = polair_lat[indy]
            lon = polair_lon[indx]

            surface_polair = earth_radius_2 * np.cos(lat * pi / 180.) * Delta_x * (pi / 180.) * Delta_y * (pi / 180.);

            for i, species in enumerate(emis_species_list):
                    emission_polair = st.eff_emission[i] / surface_polair # ug/s to ug/m2/s
                    emission[i, 0, indy, indx] = emission[i, 0, indy, indx] + emission_polair

                    
# ------------
# Write output
# ------------
    wind_dir, wind_speed, pblh, ust, lmo, psfc, t2, sh, attenuation, background, wind_dir_inter_, wind_speed_inter_, pblh_inter_, ust_inter_, lmo_inter_, emission_array, lwc, rain = write_output(node_list, street_list, node_list_eff, street_list_eff, current_date, textfile_dir, emis_species_list)

    append_binary(wind_dir, file_wind_dir)
    append_binary(wind_speed, file_wind_speed)
    append_binary(pblh, file_pblh)
    append_binary(ust, file_ust)
    append_binary(lmo, file_lmo)
    append_binary(psfc, file_psfc)
    append_binary(t2, file_t2)
    append_binary(sh, file_sh)
    append_binary(attenuation, file_attenuation)
    append_binary(lwc, file_lwc)
    append_binary(rain, file_rain)

    append_binary(wind_dir_inter_, file_wind_dir_inter)
    append_binary(wind_speed_inter_, file_wind_speed_inter)
    append_binary(pblh_inter_, file_pblh_inter)
    append_binary(ust_inter_, file_ust_inter)
    append_binary(lmo_inter_, file_lmo_inter)

    for i, species in enumerate(emis_species_list):
            file_emission = output_dir + "/emission/" + species + ".bin"
            append_binary(emission_array[:,i], file_emission)

    for spec in list(background.keys()):
        filename = output_dir + "/background/" + spec + ".bin"
        append_binary(background[spec], filename)

        
    # Write grid-averaged emissions.
    for i, species in enumerate(emis_species_list):
            file_emission = output_dir + "/grid_emission/" + species + ".bin"
            append_binary(emission[i], file_emission)

# VOC speciated
if (config.is_voc_speciated):
        import speciation_aggregation
        speciation_aggregation.speciation_voc(sys.argv[1])
        
# NOx speciated
if (config.is_nox_speciated):
        import speciation_aggregation
        speciation_aggregation.speciation_nox(sys.argv[1])

#LL----------
# PM10 speciated
# if (config.is_pm10_speciated):
#         import speciation_aggregation
#         speciation_aggregation.speciation_pm10_emis()
#         speciation_aggregation.speciation_pm10_bg()
        
# # ISVOC speciated
# if (config.is_isvoc_speciated):
#         import speciation_aggregation
#         speciation_aggregation.speciation_isvoc()
#------------        
