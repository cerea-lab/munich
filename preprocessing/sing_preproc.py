#! /usr/bin/env python

import numpy as np
import sys, datetime, math, os
from atmopy import *
from street_network import *
from optparse import OptionParser

parser = OptionParser(usage = "%prog configuration_file")
(options, args) = parser.parse_args()

if not args:
	parser.error("A configuration file is required.")

content = [("emission_dir_weekday", "[input]", "String"), \
           ("emission_dir_weekend", "[input]", "String"), \
           ("epsg_code", "[input]", "String"), \
           ("weekday_file_prefix", "[input]", "String"), \
           ("weekend_file_prefix", "[input]", "String"), \
           ("emission_species", "[input]", "StringList"), \
           ("geog_info", "[input]", "String"), \
           ("background_concentration", "[input]", "String"), \
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
           ("chimere_bg","[background]","Bool"), \
           ("chimout_dir", "[background]", "String"), \
           ("chimout_lab", "[background]", "String"), \
           ("Species","[background]","StringList"), \
           ("melch2molmass_file","[background]","String"), \
]
config = talos.Config(sys.argv[1], content)

################
# Main program #
################

output_dir = config.Output_dir
if not os.path.exists(output_dir):
        os.makedirs(output_dir)
emis_dir = output_dir + "/emission/"
if not os.path.exists(emis_dir):
        os.makedirs(emis_dir)
meteo_dir = output_dir + "/meteo/"
if not os.path.exists(meteo_dir):
        os.makedirs(meteo_dir)
background_dir = output_dir + "/background/"
if not os.path.exists(background_dir):
        os.makedirs(background_dir)
textfile_dir = output_dir + "/textfile/"
if not os.path.exists(textfile_dir):
        os.makedirs(textfile_dir)
emis_dir = output_dir + "/grid_emission/"
if not os.path.exists(emis_dir):
        os.makedirs(emis_dir)
        
emis_species_list = config.emission_species
print emis_species_list

# Number of emitted species
ns_emis = len(emis_species_list)

# emis_no_all = []
# emis_no2_all = []
emission_array_all = []

wind_dir_all = []
wind_speed_all = []
pblh_all = []
ust_all = []
lmo_all = []
psfc_all = []
t2_all = []
sh_all = []
attenuation_all = []


is_chimere=config.chimere_bg

#bg_o3_all = []
#bg_no2_all = []
#bg_no_all = []
background_all={}


# For the intersections
wind_dir_inter = []
wind_speed_inter = []
pblh_inter = []
ust_inter = []
lmo_inter = []

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

emission_no = np.zeros([config.Nt_polair, 1, Ny, Nx], 'float')
emission_no2 = np.zeros([config.Nt_polair, 1, Ny, Nx], 'float')
emission = np.zeros([config.Nt_polair, ns_emis, 1, Ny, Nx], 'float')

earth_radius_2 = 6371229. * 6371229.
pi = 3.14159265358979323846264





#################################
#                               #
#    IF BACKGROUND FROM CHIMERE #
#                               #
#################################


if is_chimere :


    melchior_spec_list=config.Species
    chimout_dir=config.chimout_dir
    chimout_lab=config.chimout_lab
    
    molar_mass_melchior2={}

    with open(config.melch2molmass_file,'r') as f:
        for ln in f:
            spec=ln.split()[0]
            val=np.float(ln.split()[1])
            molar_mass_melchior2[spec]=val
    f.close()


# Get date info
begin_date = config.t_min
delta_t = config.Delta_t # in hr
nt = config.Nt
read_lut = False
date_list = []
for t in range(nt):
    current_date = begin_date + datetime.timedelta(hours = (delta_t * t))
    print "\n====================================================="
    print "Current date (UTC): ", current_date
    date_list.append(current_date)

    # Set time index for polair
    print config.Date_min_polair
    time_diff = current_date - config.Date_min_polair
    time_diff_seconds = time_diff.days * 24 * 60 * 60 + time_diff.seconds
    indt = int(time_diff_seconds / 3600)
    print "Time index: ", indt 

    current_date_local = utc_to_local(current_date, 'Europe/Paris')

    print "Current date (local hour): ", current_date_local
    str_date = current_date_local.strftime("%Y%m%d%H")
    date = str_date[0:8]
    hour = str_date[8:10]
    # Should take into account the French holidays during the period.
    if (current_date.weekday() >= 0) and (current_date.weekday() <= 4):
            print "The current day is a weekday."
            if is_holiday(current_date):
                    print "The current day is a holiday."
                    input_file = config.emission_dir_weekend + "/" + \
                    config.weekend_file_prefix + hour
            else:
                    print "The current day is a weekday."
                    input_file = config.emission_dir_weekday + "/" + \
                    config.weekday_file_prefix + hour
    else:
            print "The current day is a weekend."
            input_file = config.emission_dir_weekend + "/" + \
            config.weekend_file_prefix + hour

    print "Read the input data (segment coordinates and emission rates) from the file --- ",input_file
    
    street_list, node_list = read_traffic_data(input_file,
                                               emis_species_list, config.epsg_code)

    # Get street geographical informations: length, width, builiding height.
    geog_info = config.geog_info
    get_street_geog(geog_info, street_list)


    # Merging streets if they have same (or very near) nodes.
    if (config.is_street_merged):
            n_street = len(street_list)
            print " - Initial number of streets: ", n_street

            # Automatic merging for the separated roads
            lut_file = "street-merging-lookup-table.txt"
            read_lut = os.path.isfile(lut_file)
            if (read_lut):
                    ntemp = lut_merging_street(lut_file, street_list)
            else:
                    ntemp = merging_street(lut_file, node_list, street_list)
                    print " - Number of streets after merging the same streets: %d" % (n_street - ntemp)

            # Manual merging for the separated roads.
            if (config.is_street_manually_merged):
                    ntemp2 = manual_merging_street(street_list)
                    print " - Number of streets after manual merging of the streets: %d" % (n_street - ntemp - ntemp2)

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
    print " === Merging nodes === "
    n_node = len(node_list_temp)
    print "Initial number of nodes: ", n_node
    n_node1 = merging_node(node_list_temp)
    print "Number of nodes after removing the same nodes: ", (n_node - n_node1)
    if config.is_near_node_merged:
            n_node2 = merging_near_node(node_list_temp)
            print "Number of nodes after removing the nearest nodes: ", (n_node - n_node1 - n_node2)
    if config.is_node_manually_merged:
            n_node3 = manual_merging_node(node_list_temp)
            print "Number of nodes after manually removing the nearest nodes: ", (n_node - n_node1 - n_node2 - n_node3)

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

    get_meteo_data(config.meteo_dir, current_date, \
                           street_list_eff, node_list_eff, config.wrfout_prefix)

    background_concentration_file = config.background_concentration


    if is_chimere :
        get_chimere_background_concentration(current_date, street_list_eff, \
                             melchior_spec_list,molar_mass_melchior2, \
                             chimout_dir,chimout_lab)


    else :
        get_background_concentration(background_concentration_file, \
                                         current_date, street_list_eff)
    # Compute emissions in grid cells.
    for st in street_list_eff:
            indx, indy = get_polair_ind_v2(st.lon_cen, st.lat_cen, x_min, Delta_x, y_min, Delta_y, Nx, Ny)
            lat = polair_lat[indy]
            lon = polair_lon[indx]

            surface_polair = earth_radius_2 * np.cos(lat * pi / 180.) * Delta_x * (pi / 180.) * Delta_y * (pi / 180.);

            for i, species in enumerate(emis_species_list):
                    emission_polair = st.eff_emission[i] / surface_polair # ug/s to ug/m2/s
                    emission[indt, i, 0, indy, indx] = emission[indt, i, 0, indy, indx] + emission_polair
                
# ------------
# Write output
# ------------
    wind_dir, wind_speed, pblh, ust, lmo, psfc, t2, sh, attenuation, background, wind_dir_inter_, wind_speed_inter_, pblh_inter_, ust_inter_, lmo_inter_, emission_array = write_output(node_list, street_list, node_list_eff, street_list_eff, current_date, textfile_dir, emis_species_list)

    emission_array_all.append(emission_array)

    wind_dir_all.append(wind_dir)
    wind_speed_all.append(wind_speed)
    pblh_all.append(pblh)
    ust_all.append(ust)
    lmo_all.append(lmo)
    psfc_all.append(psfc)
    t2_all.append(t2)
    sh_all.append(sh)
    attenuation_all.append(attenuation)

    wind_dir_inter.append(wind_dir_inter_)
    wind_speed_inter.append(wind_speed_inter_)
    pblh_inter.append(pblh_inter_)
    ust_inter.append(ust_inter_)
    lmo_inter.append(lmo_inter_)

    for spec in background.keys():
        val=background[spec]
        if spec in background_all.keys():
           background_all[spec].append(val)
        else:
           background_all[spec]=[val]


emission_array_all = np.array(emission_array_all)
wind_dir_all = np.array(wind_dir_all)
wind_speed_all = np.array(wind_speed_all)
pblh_all = np.array(pblh_all)
ust_all = np.array(ust_all)
lmo_all = np.array(lmo_all)
psfc_all = np.array(psfc_all)
t2_all = np.array(t2_all)
sh_all = np.array(sh_all)
attenuation_all = np.array(attenuation_all)

wind_dir_inter = np.array(wind_dir_inter)
wind_speed_inter = np.array(wind_speed_inter)
pblh_inter = np.array(pblh_inter)
ust_inter = np.array(ust_inter)
lmo_inter = np.array(lmo_inter)

for spec in background_all.keys():
    background_all[spec]=np.array(background_all[spec])




for i, species in enumerate(emis_species_list):
    file_emission = output_dir + "/emission/" + species + ".bin"
    io.save_binary(emission_array_all[:,:,i], file_emission)

file_wind_dir = output_dir + "/meteo/WindDirection.bin"
file_wind_speed = output_dir + "/meteo/WindSpeed.bin"
file_pblh = output_dir + "/meteo/PBLH.bin"
file_ust = output_dir + "/meteo/UST.bin"
file_lmo = output_dir + "/meteo/LMO.bin"
file_psfc = output_dir + "/meteo/SurfacePressure.bin"
file_t2 = output_dir + "/meteo/SurfaceTemperature.bin"
file_sh = output_dir + "/meteo/SpecificHumidity.bin"
file_attenuation = output_dir + "/meteo/Attenuation.bin"

file_wind_dir_inter = output_dir + "/meteo/WindDirectionInter.bin"
file_wind_speed_inter = output_dir + "/meteo/WindSpeedInter.bin"
file_pblh_inter = output_dir + "/meteo/PBLHInter.bin"
file_ust_inter = output_dir + "/meteo/USTInter.bin"
file_lmo_inter = output_dir + "/meteo/LMOInter.bin"

for spec in background_all.keys():
    filename=output_dir + "/background/"+spec+".bin"
    io.save_binary(background_all[spec], filename)


io.save_binary(wind_dir_all, file_wind_dir)
io.save_binary(wind_speed_all, file_wind_speed)
io.save_binary(pblh_all, file_pblh)
io.save_binary(ust_all, file_ust)
io.save_binary(lmo_all, file_lmo)
io.save_binary(psfc_all, file_psfc)
io.save_binary(t2_all, file_t2)
io.save_binary(sh_all, file_sh)
io.save_binary(attenuation_all, file_attenuation)

io.save_binary(wind_dir_inter, file_wind_dir_inter)
io.save_binary(wind_speed_inter, file_wind_speed_inter)
io.save_binary(pblh_inter, file_pblh_inter)
io.save_binary(ust_inter, file_ust_inter)
io.save_binary(lmo_inter, file_lmo_inter)

for i, species in enumerate(emis_species_list):
        file_emission = output_dir + "/grid_emission/" + species + ".bin"
        io.save_binary(emission[:, i], file_emission)

# VOC speciated
if (config.is_voc_speciated):
        import speciation_aggregation

# NOx speciated
if (config.is_nox_speciated):
        import speciation_aggregation
        speciation_aggregation.speciation_nox()
        
