#! /usr/bin/env python

import math
import numpy as np
from io_tools import *
from street_network import *

# Read emission binary files of shape (24, Nstreet)
def read_emission_bin_day(indir_weekday, indir_weekend, n_street, species_list):
    # Put data in dicts, one for each typical day
    emission_data_weekday = {}
    emission_data_weekend = {}
    for species in species_list:
        infile = indir_weekday + species + '.bin'
        emission_data_weekday[species] = np.memmap(infile, dtype='float32', mode='r',
                                                   shape=(24, n_street))
        infile = indir_weekend + species + '.bin'
        emission_data_weekend[species] = np.memmap(infile, dtype='float32', mode='r',
                                                   shape=(24, n_street))

    return emission_data_weekday, emission_data_weekend

# Read emission binary files of shape (Nt_emission, Nstreet)
def read_emission_bin_period(indir, Nt_emission, n_street, species_list):
    # Put data in a dict
    emission_data = {}
    for species in species_list:
        infile = indir + species + '.bin'
        emission_data[species] = np.memmap(infile, dtype='float32', mode='r',
                                           shape=(Nt_emission, n_street))

    return emission_data

# Set emission for streets using typical days emission data
def set_emission_bin_day(street_list, street_list_eff, emission_data, hour_id, species_list):
    # Initialize emission for effective streets
    for i in range(len(street_list_eff)):
        for species in species_list:
            street_list_eff[i].emission[species] = 0.
    # Add emission to effective streets
    for i in range(len(street_list)):
        for j in range(len(street_list_eff)):
            if street_list[i].eff_id == street_list_eff[j].id:
                for species in species_list:
                    # Convert from g/km/h to µg/s
                    emission = emission_data[species][hour_id, i] * 1e6
                    emission *= (street_list[i].length / 1000.) / 3600.
                    street_list_eff[j].emission[species] += emission

# Set emission for streets using period emission data
def set_emission_bin_period(street_list, street_list_eff, emission_data, current_date,
                            date_min, delta_t, Nt, species_list):
    # Get index of current date
    c_id = int((current_date - date_min).total_seconds() / delta_t)
    if c_id < 0 or c_id >= Nt:
        sys.exit('ERROR: emission data not available for this date.')

    # Initialize emission for effective streets
    for i in range(len(street_list_eff)):
        for species in species_list:
            street_list_eff[i].emission[species] = 0.
    # Add emission to effective streets
    for i in range(len(street_list)):
        for j in range(len(street_list_eff)):
            if street_list[i].eff_id == street_list_eff[j].id:
                for species in species_list:
                    # Convert from g/km/h to µg/s
                    emission = emission_data[species][c_id, i] * 1e6
                    emission *= (street_list[i].length / 1000.) / 3600.
                    street_list_eff[j].emission[species] += emission

# Append emission to binary files
def append_emission_data(street_list, outdir, species_list):
    for species in species_list:
        data = np.zeros((len(street_list)), 'float')
        for j in range(len(street_list)):
            data[j] = street_list[j].emission[species]
        outfile = outdir + species + '.bin'
        append_binary(data, outfile)

##### Gridded emission

# Append background to binary files
def append_grid_emission(emission,
                         emis_species_list,
                         grid_emission_files):
    
    for i, var in enumerate(emis_species_list):
        data = emission[i]
        append_binary(data, grid_emission_files[var])



def read_street_emission(input_file,
                         emis_species_list,
                         epsg_code,
                         is_street_manually_merged,
                         manual_street_merging_file,
                         is_node_manually_merged,
                         manual_node_merging_file,
                         output_dir):


    street_list, node_list = read_traffic_data(input_file,
                                               emis_species_list,
                                               epsg_code)

    # # Get street geographical informations: length, width, builiding height.
    # geog_info = config.geog_info
    # get_street_geog(geog_info, street_list)

    # Merging streets if they have same (or very near) nodes.
    is_street_merged = True
    if (is_street_merged):
            n_street = len(street_list)
            print(" - Initial number of streets: ", n_street)

            # Automatic merging for the separated roads
            lut_file = "street-merging-lookup-table.txt"
            read_lut = os.path.isfile(lut_file)
            if (read_lut):
                    ntemp = lut_merging_street(lut_file,
                                               street_list,
                                               emis_species_list)
            else:
                    ntemp = merging_street(lut_file,
                                           node_list,
                                           street_list,
                                           emis_species_list)
            print(" - Number of streets after merging the same streets: %d" % (n_street - ntemp))

            # Manual merging for the separated roads.
            if (is_street_manually_merged):
                    ntemp2 = manual_merging_street(street_list, emis_species_list, manual_street_merging_file)
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
    is_near_node_merged = True
    if is_near_node_merged:
            n_node2 = merging_near_node(node_list_temp)
            print("Number of nodes after removing the nearest nodes: ", (n_node - n_node1 - n_node2))

    if is_node_manually_merged:
            n_node3 = manual_merging_node(node_list_temp,
                                          manual_node_merging_file)
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

    outfile = output_dir + 'textfile/street_all.csv'
    shutil.rmtree(output_dir + 'textfile')
    os.makedirs(output_dir + 'textfile')
    with open(outfile, 'w') as f:
        f.write('id' + ',' + 'eff_id' + ',' + 'node_begin'\
            + ',' + 'node_eff_begin' + ',' + 'node_end' + ','\
            + 'node_eff_end' + '\n')

        for street in street_list:
            f.write(str(street.id) + ',' + str(street.eff_id) + ',' + str(street.begin)\
                    + ',' + str(street.eff_begin) + ',' + str(street.end) + ','\
                    + str(street.eff_end) + '\n')

    return street_list_eff, node_list_eff
        
