#!/usr/bin/env python3

import os, sys
import shutil
import argparse
import numpy as np
import datetime as dt
from atmopy import talos
from utils_sing_preproc import *

parser = argparse.ArgumentParser(usage='%(prog)s configuration_file')
parser.add_argument('config_file')
config_file = parser.parse_args().config_file

##### Read configuration file

content = [('indir', '[input]', 'String'),
           ('outdir', '[output]', 'String'),
           ('date_min', '[option]', 'DateTime'),
           ('delta_t', '[option]', 'Float'),
           ('Nt', '[option]', 'Int'),
           # [network]
           ('generate_network', '[network]', 'Bool'),
           ('network_outdir', '[network]', 'String'),
           ('epsg_code', '[network]', 'Int'),
           ('create_network', '[network]', 'Bool'),
           ('street_file', '[network]', 'String'),
           ('min_distance', '[network]', 'Float'),
           ('manual_street_merging', '[network]', 'Bool'),
           ('manual_street_merging_file', '[network]', 'String'),
           ('manual_node_merging', '[network]', 'Bool'),
           ('manual_node_merging_file', '[network]', 'String'),
           ('write_luts', '[network]', 'Bool'),
           ('street_lut_file', '[network]', 'String'),
           ('node_lut_file', '[network]', 'String'),
           # [emission]
           ('generate_emission', '[emission]', 'Bool'),
           ('emission_outdir', '[emission]', 'String'),
           ('emission_species', '[emission]', 'StringList'),
           ('is_local_time', '[emission]', 'Bool'),
           ('time_zone', '[emission]', 'String'),
           ('emission_input', '[emission]', 'String'),
           ('emission_indir_weekday', '[emission]', 'String'),
           ('emission_indir_weekend', '[emission]', 'String'),
           ('country_code', '[emission]', 'String'),
           ('emission_indir', '[emission]', 'String'),
           ('date_min_emission', '[emission]', 'DateTime'),
           ('delta_t_emission', '[emission]', 'Float'),
           ('Nt_emission', '[emission]', 'Int'),
           ('emission_type', '[emission]', 'String'),
           # [meteo]
           ('generate_meteo', '[meteo]', 'Bool'),
           ('meteo_indir', '[meteo]', 'String'),
           ('meteo_outdir', '[meteo]', 'String'),
           ('meteo_type', '[meteo]', 'String'),
           ('date_min_meteo', '[meteo]', 'DateTime'),
           ('delta_t_meteo', '[meteo]', 'Float'),
           ('Nt_meteo', '[meteo]', 'Int'),
           ('x_min_meteo', '[meteo]', 'Float'),
           ('y_min_meteo', '[meteo]', 'Float'),
           ('delta_x_meteo', '[meteo]', 'Float'),
           ('delta_y_meteo', '[meteo]', 'Float'),
           ('Nx_meteo', '[meteo]', 'Int'),
           ('Ny_meteo', '[meteo]', 'Int'),
           ('Nz_meteo', '[meteo]', 'Int'),
           # [background]
           ('generate_bkgd', '[background]', 'Bool'),
           ('bkgd_indir', '[background]', 'String'),
           ('bkgd_outdir', '[background]', 'String'),
           ('bkgd_species', '[background]', 'StringList'),
           ('bkgd_type', '[background]', 'String'),
           ('date_min_bkgd', '[background]', 'DateTime'),
           ('delta_t_bkgd', '[background]', 'Float'),
           ('Nt_bkgd', '[background]', 'Int'),
           ('x_min_bkgd', '[background]', 'Float'),
           ('y_min_bkgd', '[background]', 'Float'),
           ('delta_x_bkgd', '[background]', 'Float'),
           ('delta_y_bkgd', '[background]', 'Float'),
           ('Nx_bkgd', '[background]', 'Int'),
           ('Ny_bkgd', '[background]', 'Int'),
           ('Nz_bkgd', '[background]', 'Int')
           # [gridded_emission]
           ]

config = talos.Config(config_file, content)

# Clear and create output directories
#shutil.rmtree(config.outdir, ignore_errors=True)
os.makedirs(config.outdir, exist_ok=True)
if config.generate_emission:
    shutil.rmtree(config.emission_outdir, ignore_errors=True)
    os.makedirs(config.emission_outdir)
if config.generate_meteo:
    shutil.rmtree(config.meteo_outdir, ignore_errors=True)
    os.makedirs(config.meteo_outdir)
if config.generate_bkgd:
    shutil.rmtree(config.bkgd_outdir, ignore_errors=True)
    os.makedirs(config.bkgd_outdir)

##### Street network
print('============================== Street network')

if config.create_network:
    print('-> Create network from scratch')

    # Read data from file
    print('--------------------')
    print('Read data')
    node_list, street_list = read_street_data(config.street_file, config.epsg_code)
    n_street = len(street_list)
    print('-- Number of nodes: ', len(node_list))
    print('-- Number of streets: ', n_street)

    # Merge street
    print('--------------------')
    print('Merge streets')
    print('-- Initial number of streets: ', n_street)
    print('-- Number of streets after automatic merging: ', end='')
    n_street1 = auto_street_merging(node_list, street_list, config.min_distance)
    print(n_street - n_street1)
    print('Manual street merging: ', end='')
    if config.manual_street_merging:
        print('yes')
        print('-- Number of streets after manual merging: ', end='')
        n_street2 = manual_street_merging(street_list, config.manual_street_merging_file)
        print(n_street - n_street1 - n_street2)
    else:
        print('no')

    # Make final list of streets and set eff_id of corresponding nodes
    street_list_eff = []
    for i in range(n_street):
        if street_list[i].removed == False:
            street_list_eff.append(street_list[i])
            node_list[2 * i].eff_id = node_list[2 * i].id # begin
            node_list[2 * i + 1].eff_id = node_list[2 * i + 1].id # end
    # Make temporary list of nodes
    node_list_tmp = []
    for i in range(len(node_list)):
        if node_list[i].id == node_list[i].eff_id:
            node_list_tmp.append(node_list[i])
    n_node = len(node_list_tmp)

    # Merge node
    print('--------------------')
    print('Merge nodes')
    print('-- Initial number of nodes: ', n_node)
    print('-- Number of nodes after automatic merging: ', end='')
    n_node1 = auto_node_merging(node_list_tmp, config.min_distance)
    print(n_node - n_node1)
    print('Manual node merging: ', end='')
    if config.manual_node_merging:
        print('yes')
        print('-- Number of nodes after manual merging: ', end='')
        n_node2 = manual_node_merging(node_list_tmp, config.manual_node_merging_file)
        print(n_node - n_node1 - n_node2)
    else:
        print('no')

    # Make final list of nodes
    node_list_eff = []
    for i in range(n_node):
        if node_list_tmp[i].removed == False:
            node_list_eff.append(node_list_tmp[i])
    # Correct effective begin and end of the streets
    for i in range(len(node_list)):
        for j in range(len(street_list_eff)):
            if street_list_eff[j].begin == node_list[i].id:
                street_list_eff[j].eff_begin = node_list[i].eff_id
            elif street_list_eff[j].end == node_list[i].id:
                street_list_eff[j].eff_end = node_list[i].eff_id

    print('--------------------')
    print('-- Number of effective nodes: ', len(node_list_eff))
    print('-- Number of effective streets: ', len(street_list_eff))

    # Write LUTs
    print('--------------------')
    print('Write LUTs: ', end='')
    if config.write_luts:
        write_luts(node_list_eff, street_list, config.indir)
        print('yes')
    else:
        print('no')

else:
    print('-> Create network from LUTs')

    # Read LUTs
    print('--------------------')
    print('Read LUTs')
    node_list_eff, street_list, street_list_eff = read_luts(config.node_lut_file,
                                                            config.street_lut_file)
    print('-- Number of streets: ', len(street_list))

    print('--------------------')
    print('-- Number of effective nodes: ', len(node_list_eff))
    print('-- Number of effective streets: ', len(street_list_eff))

# Generate network files for MUNICH
print('--------------------')
print('Generate network files for MUNICH: ', end='')
if config.generate_network:
    generate_network(node_list_eff, street_list_eff, config.network_outdir)
    print('yes\n')
else:
    print('no\n')

##### Generate files for MUNICH
print('============================== MUNICH input')

# Initialize emission information
if config.generate_emission:
    if config.emission_type == 'bin':
        if config.emission_input == 'day':
            # Load files using memmap and return dicts
            emission_data_weekday, emission_data_weekend = read_emission_bin_day(config.emission_indir_weekday,
                                                                                 config.emission_indir_weekend,
                                                                                 len(street_list),
                                                                                 config.emission_species)
        elif config.emission_input == 'period':
            # Load file using memmap and return a dict
            emission_data = read_emission_bin_period(config.emission_indir, config.Nt_emission,
                                                     len(street_list), config.emission_species)

    elif config.emission_type == 'csv':
        sys.exit('WARNING: emission_type csv not implemented yet.')

    else:
        sys.exit('ERROR: emission_type {} not valid.'.format(config.emission_type))

# Initialize meteo information
if config.generate_meteo:
    # Meteo variables
    # WARNING: order is important
    meteo_var_street = ['Attenuation', 'LiquiWaterContent', 'LMO', 'PARdb', 'PARdiff',
                        'PBLH', 'Rain', 'SoilWater', 'SolarRadiation', 'SpecificHumidity',
                        'SurfacePressure', 'SurfaceRichardson', 'SurfaceTemperature',
                        'UST', 'WindDirection', 'WindSpeed']
    meteo_var_node = ['LMOInter', 'PBLHInter', 'USTInter', 'WindDirectionInter', 'WindSpeedInter']
    # Output binary files to append to for streets and nodes
    meteo_files_street = [config.meteo_outdir + x + '.bin' for x in meteo_var_street]
    meteo_files_node = [config.meteo_outdir + x + '.bin' for x in meteo_var_node]

    if config.meteo_type == 'bin':
        # Load files using memmap and return a dict
        meteo_data = read_meteo_bin(config.meteo_indir, config.Nt_meteo, config.Nx_meteo,
                                    config.Ny_meteo, config.Nz_meteo)

    elif config.meteo_type == 'wrf':
        sys.exit('WARNING: meteo_type wrf not implemented yet.')

    else:
        sys.exit('ERROR: meteo_type {} not valid.'.format(config.meteo_type))

# Initialize background information
if config.generate_bkgd:
    # Output binary files to append to for streets
    bkgd_files = {x: config.bkgd_outdir + x + '.bin' for x in config.bkgd_species}

    if config.bkgd_type == 'bin':
        # Load files using memmap and return a dict
        bkgd_data = read_bkgd_bin(config.bkgd_species, config.bkgd_indir, config.Nt_bkgd,
                                  config.Nx_bkgd, config.Ny_bkgd, config.Nz_bkgd)

    elif config.bkgd_type == 'csv':
        sys.exit('WARNING: bkgd_type csv not implemented yet.')

    elif config.bkgd_type == 'chimere':
        sys.exit('WARNING: bkgd_type chimere not implemented yet.')

    else:
        sys.exit('ERROR: bkgd_type {} not valid.'.format(config.bkgd_type))

start_date = config.date_min
delta_t = config.delta_t
nt = config.Nt
for t in range(nt):
    current_date = start_date + dt.timedelta(seconds=(t * delta_t))
    print('====================')
    print('Current date (UTC): ', current_date.strftime('%Y%m%d_%H'))

    ##### Emission
    print('---------- Emission')
    if config.generate_emission:
        if config.is_local_time:
            print('Conversion from {} to UTC.'.format(config.time_zone))
            current_date_local = utc_to_local(current_date, config.time_zone)
        else:
            current_date_local = current_date
        print('Current date (local): ', current_date_local.strftime('%Y%m%d_%H'))

        print('Emission_type: ', config.emission_type)
        if config.emission_type == 'bin':
            if config.emission_input == 'day':
                if current_date.weekday() >= 0 and current_date.weekday() <= 4:
                    if is_holiday(current_date, config.country_code):
                        print('Current day is a holiday.')
                        emission_data = emission_data_weekend
                    else:
                        print('Current day is a weekday.')
                        emission_data = emission_data_weekday
                else:
                    print('Current day is a weekend.')
                    emission_data = emission_data_weekend
                # Hour index in the emission dataset
                hour_id = current_date_local.hour
                set_emission_bin_day(street_list, street_list_eff, emission_data, hour_id,
                                     config.emission_species)

            else:
                set_emission_bin_period(street_list, street_list_eff, emission_data,
                                        current_date_local, config.date_min_emission,
                                        config.delta_t_emission, config.Nt_emission,
                                        config.emission_species)

        elif config.emission_type == 'csv':
            pass

        # Append data to output files
        append_emission_data(street_list_eff, config.emission_outdir, config.emission_species)
        print('-> Emission files generated.')

    else:
        print('-> Emission files not generated.')

    ##### Meteo
    print('---------- Meteo')
    if config.generate_meteo:
        print('Meteo_type: ', config.meteo_type)
        if config.meteo_type == 'bin':
            set_meteo_bin(node_list_eff, street_list_eff, meteo_data, current_date,
                          config.date_min_meteo, config.delta_t_meteo, config.Nt_meteo,
                          config.x_min_meteo, config.y_min_meteo, config.delta_x_meteo,
                          config.delta_y_meteo, config.Nx_meteo, config.Ny_meteo)

        elif config.meteo_type == 'wrf':
            pass

        # Append data to output files
        append_meteo_data(node_list_eff, street_list_eff, meteo_files_street, meteo_files_node)
        print('-> Meteo files generated.')

    else:
        print('-> Meteo files not generated.')

    ##### Background
    print('---------- Background')
    if config.generate_bkgd:
        print('Bkgd_type: ', config.bkgd_type)
        if config.bkgd_type == 'bin':
            set_bkgd_bin(street_list_eff, bkgd_data, current_date, config.date_min_bkgd,
                         config.delta_t_bkgd, config.Nt_bkgd, config.x_min_bkgd,
                         config.y_min_bkgd, config.delta_x_bkgd, config.delta_y_bkgd,
                         config.Nx_bkgd, config.Ny_bkgd)

        elif config.bkgd_type == 'csv':
            pass

        elif config.bkgd_type == 'chimere':
            pass

        # Append data to output files
        append_bkgd_data(street_list_eff, bkgd_files)
        print('-> Background files generated.')

    else:
        print('-> Background files not generated.')

    ##### Gridded emission
    #print('---------- Gridded emission')

    print('\n')