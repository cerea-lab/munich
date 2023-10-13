#! /usr/bin/env python

import re, sys, os
from atmopy import *
from optparse import OptionParser
import numpy as np


################
# Main program #
################

# parser = OptionParser(usage = "%prog configuration_file")
# (options, args) = parser.parse_args()

content = [("outdir", "[output]", "String"), \
           ("speciation_dir", "[emission]", "String"), \
           ("meca","[emission]","String"), \
           ('generate_grid_emission', '[gridded_emission]', 'Bool'),
           ("Nt_polair", "[gridded_emission]", "Int")
]
#config = talos.Config("sing_preproc.cfg", content)


content_spec = [
    # [speciation_pm10]
    ('speciate_pm10', '[speciation_pm10]', 'Bool'),
    ('infile_pm25', '[speciation_pm10]', 'String'),
    ('infile_pm10', '[speciation_pm10]', 'String'),
    ('N_bin', '[speciation_pm10]', 'Int'),
    ('fraction_pm10_exhaust', '[speciation_pm10]', 'Float'),
    ('fraction_pm10_brake', '[speciation_pm10]', 'Float'),
    ('fraction_pm10_tyre', '[speciation_pm10]', 'Float'),
    ('fraction_pm10_road', '[speciation_pm10]', 'Float'),
    ('spec_pm10_exhaust_pbc', '[speciation_pm10]', 'Float'),
    ('spec_pm10_exhaust_pmd', '[speciation_pm10]', 'Float'),
    ('spec_pm10_exhaust_ppoalp', '[speciation_pm10]', 'Float'),
    ('spec_pm10_brake', '[speciation_pm10]', 'String'),
    ('spec_pm10_tyre', '[speciation_pm10]', 'String'),
    ('spec_pm10_road', '[speciation_pm10]', 'String'),
]



def speciation_voc(configuration_file = "sing_preproc.cfg"):

    print(configuration_file)
    config = talos.Config(configuration_file, content)
    meca = config.meca
    home_dir = config.outdir

    speciation_file = config.speciation_dir + "/COVNM_"+meca+".dat"
    speciation = open(speciation_file)
    header = speciation.readline()

    # Get snap code for traffic source.
    snap_line = speciation.readline()
    line_info = [x for x in re.split('\t| ', snap_line) if x.strip() != '']
    for n, line in enumerate(line_info):
        if line[0] == "7":
                index = n
                break
    print("traffic snap index: ", index)

    sp_real_name = []
    speciation_coeff = []
    mw_real = []
    total = 0.0
    intoh=1.  # for melchior2 aggreagation
    for line in speciation.readlines():
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        sp_real_name.append(line_info[0])
        speciation_coeff.append(float(line_info[index]))
        total += float(line_info[index])
        mw_real.append(float(line_info[1]))
    speciation.close()

    print("total coeff.: ", total, "%")
    ns_real = len(sp_real_name)
    print("Number of the real species: ", ns_real)

    # Get VOC aggregation to model species.
    if meca == "cb05":
        aggregation_file = config.speciation_dir + "/aggregation_cb05-siream.dat"
    elif meca == "melchior2": 
        aggregation_file = config.speciation_dir + "/aggregation_melchior2.dat"
    else:
        print('For now the "meca" flag at the [option] section of the sing_preproc.cfg has to be set to either cb05 or melchior2')
        sys.exit()

    aggregation = open(aggregation_file)
    # Read the line for the list of the model species.
    line = aggregation.readline()
    line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
    model_species = []
    for n, line in enumerate(line_info):
        if n > 2:
                model_species.append(line.strip())
    ns_model = len(model_species)

    line = aggregation.readline()
    line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
    mw_model = []
    for n, line in enumerate(line_info):
        if n > 2:
                mw_model.append(float(line))

    header = aggregation.readline()

    if meca == 'melchior2' :
        line_info = [x for x in re.split('\t| ', header) if x.strip() != '']
        reac_model = []
        for n, line in enumerate(line_info):
            if n > 2:
                reac_model.append(float(line))




    species_factor = np.zeros([ns_model], 'float')
    for s_real in range(ns_real):
        line = aggregation.readline()
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        if meca == 'melchior2' :
           reac_real=float(line_info[2])
        for s_model in range(ns_model):
                aggregation_coeff = float(line_info[s_model + 3])
                if meca=="cb05":
                    temp = speciation_coeff[s_real] * 0.01 * aggregation_coeff * mw_model[s_model] / mw_real[s_real]
                elif meca == "melchior2":
                    temp=speciation_coeff[s_real] * 0.01 * aggregation_coeff * mw_model[s_model] / mw_real[s_real] * \
                                               (1.-np.exp(-intoh*reac_real))/(1.-np.exp(-intoh*reac_model[s_model]))
                else:
                    print('For now the "meca" flag at the [option] section of the sing_preproc.cfg has to be set to either cb05 or melchior2')
                    sys.exit()

                species_factor[s_model] += temp
    aggregation.close()



    ### Compute speciated VOC emissions for street segments.

    input_dir = home_dir + "/emission/"
    input_file = input_dir + "NMHC.bin"

    # Get array shape of input file and load to an array.
    inputfile_size = io.get_filesize(input_file)
    array_shape = (config.Nt, int(inputfile_size / config.Nt / 4.0))
    input_array = io.load_binary(input_file, array_shape)

    total = 0.0
    # Compute and write to binary files;
    for s_model in range(ns_model):
        # print("species factor: ",model_species[s_model], species_factor[s_model])
        # command = "mult_nb_float " + input_file + " " + str(species_factor[s_model]) + " " + model_species[s_model] + ".bin"
        
        output_array = input_array * species_factor[s_model]
        output_filename = input_dir + model_species[s_model] + ".bin"

        if (model_species[s_model] == "CH4"):
            print(model_species[s_model], " is ignored.")
        else:
            io.save_binary(output_array, output_filename)
            total += species_factor[s_model]

    print("total species factor:", total)


    ### Compute speciated VOC emissions for grided data.

    if (config.generate_grid_emission == True):
        # Get array shape of input file and load to an array.
        input_dir = home_dir + "/grid_emission/"
        input_file = input_dir + "NMHC.bin"

        inputfile_size = io.get_filesize(input_file)
        array_shape = (config.Nt_polair, config.Ny, config.Nx)
        input_array = io.load_binary(input_file, array_shape)

        # Compute and write to binary files;
        for s_model in range(ns_model):
#            print("species factor: ",model_species[s_model], species_factor[s_model])
            # command = "mult_nb_float " + input_file + " " + str(species_factor[s_model]) + " " + model_species[s_model] + ".bin"
            # print command
            # os.system(command)

            output_array = input_array * species_factor[s_model]
            output_filename = input_dir + model_species[s_model] + ".bin"

            if (model_species[s_model] == "CH4"):
                print(model_species[s_model], " is ignored.")
            else:
                io.save_binary(output_array, output_filename)



def speciation_nox(configuration_file = "sing_preproc.cfg"):
    
        config = talos.Config(configuration_file, content)
    
        home_dir = config.outdir
        input_dir = home_dir + "/emission/"
        input_file = input_dir + "NOx.bin"

        ### Compute speciated VOC emissions for street segments.

        # Get array shape of input file and load to an array.
        inputfile_size = io.get_filesize(input_file)
        array_shape = (config.Nt, int(inputfile_size / config.Nt / 4.0))
        nox = io.load_binary(input_file, array_shape)

        # Compute and write to binary files for NO and NO2 species;
        # NO2 is assumed to be 20% of NOx.
        no2 = nox * 0.2;
        # NOx is given as NO2 equivalent.
        no = (nox - no2) *30.0 / 46.0
        
        output_filename = input_dir + "NO2.bin"
        io.save_binary(no2, output_filename)
        output_filename = input_dir + "NO.bin"
        io.save_binary(no, output_filename)

        ### Compute speciated VOC emissions for grided data.

        if (config.generate_grid_emission == True):
        
            # Get array shape of input file and load to an array.
            input_dir = home_dir + "/grid_emission/"
            input_file = input_dir + "NOx.bin"

            inputfile_size = io.get_filesize(input_file)
            array_shape = (config.Nt_polair, config.Ny, config.Nx)
            nox = io.load_binary(input_file, array_shape)

            # Compute and write to binary files  for NO and NO2 species;
            # NO2 is assumed to be 20% of NOx.
            no2 = nox * 0.2;
            # NOx is given as NO2 equivalent.
            no = (nox - no2) *30.0 / 46.0

            output_filename = input_dir + "NO2.bin"
            io.save_binary(no2, output_filename)
            output_filename = input_dir + "NO.bin"
            io.save_binary(no, output_filename);

"""
Author: Thibaud Sarica

"""            
def speciation_pm10(configuration_file = "sing_preproc.cfg"):

    
    config_main = talos.Config(configuration_file, content)
    config = talos.Config("speciation_aggregation.cfg",
                          content_spec)

    home_dir = config_main.outdir
    input_dir = home_dir + "/emission/"

    # Get array shape of input file and load to an array.
    # Read PM10
    input_file = input_dir + config.infile_pm10
    inputfile_size = io.get_filesize(input_file)
    array_shape = (config_main.Nt, int(inputfile_size / config_main.Nt / 4.0))
    pm10 = io.load_binary(input_file, array_shape)

    # Read PM2.5
    input_file = input_dir + config.infile_pm25
    inputfile_size = io.get_filesize(input_file)
    array_shape = (config_main.Nt, int(inputfile_size / config_main.Nt / 4.0))
    pm25 = io.load_binary(input_file, array_shape)

        
    # Split exhaust vs. non-exhaust
    exhaust_pm10 = config.fraction_pm10_exhaust * pm10
    brake_pm10 = config.fraction_pm10_brake * pm10
    tyre_pm10 = config.fraction_pm10_tyre * pm10
    road_pm10 = config.fraction_pm10_road * pm10

    # Final speciated PM species
    pm_species = {x: [np.zeros(array_shape) for i in range(config.N_bin)]
                  for x in ['PBC', 'PMD', 'PPOAlP']}
        
    # Exhaust
    pm_species['PBC'][1] += exhaust_pm10 * config.spec_pm10_exhaust_pbc
    pm_species['PMD'][1] += exhaust_pm10 * config.spec_pm10_exhaust_pmd
    pm_species['PPOAlP'][1] += exhaust_pm10 * config.spec_pm10_exhaust_ppoalp
    

    # Brake wear
    with open(config.spec_pm10_brake, 'r') as f:
        header = f.readline().replace('\n', '').split(',')
        if len(header)-1 != config.N_bin:
            sys.exit('ERROR: N_bin different in brake_wear.csv.')
        for line in f.readlines():
            line_info = line.replace('\n', '').split(',')
            for i in range(config.N_bin):
                pm_species[line_info[0]][i] += brake_pm10 * float(line_info[i+1])

    # Tyre wear
    with open(config.spec_pm10_tyre, 'r') as f:
        header = f.readline().replace('\n', '').split(',')
        if len(header)-1 != config.N_bin:
            sys.exit('ERROR: N_bin different in tyre_wear.csv.')
        for line in f.readlines():
            line_info = line.replace('\n', '').split(',')
            for i in range(config.N_bin):
                pm_species[line_info[0]][i] += tyre_pm10 * float(line_info[i+1])

    # Road wear
    with open(config.spec_pm10_road, 'r') as f:
        header = f.readline().replace('\n', '').split(',')
        if len(header)-1 != config.N_bin:
            sys.exit('ERROR: N_bin different in road_wear.csv.')
        for line in f.readlines():
            line_info = line.replace('\n', '').split(',')
            for i in range(config.N_bin):
                pm_species[line_info[0]][i] += road_pm10 * float(line_info[i+1])
    
    # Write files
    for species, value in pm_species.items():
        for i in range(config.N_bin):
            output_filename = input_dir + "/" + species + '_{}.bin'.format(i)
            io.save_binary(value[i], output_filename)
        

    ### Compute speciated PM emissions for grided data.

    if (config_main.generate_grid_emission == True):
        
        # Get array shape of input file and load to an array.
        input_dir = home_dir + "/grid_emission/"

        # Get array shape of input file and load to an array.
        # Read PM10
        input_file = input_dir + config.infile_pm10
        inputfile_size = io.get_filesize(input_file)
        array_shape = (config_main.Nt_polair, config_main.Ny, config_main.Nx)
        pm10 = io.load_binary(input_file, array_shape)

        # Read PM2.5
        input_file = input_dir + config.infile_pm25
        inputfile_size = io.get_filesize(input_file)
        array_shape = (config_main.Nt_polair, config_main.Ny, config_main.Nx)
        pm25 = io.load_binary(input_file, array_shape)


        # Split exhaust vs. non-exhaust
        exhaust_pm10 = config.fraction_pm10_exhaust * pm10
        brake_pm10 = config.fraction_pm10_brake * pm10
        tyre_pm10 = config.fraction_pm10_tyre * pm10
        road_pm10 = config.fraction_pm10_road * pm10

        # Final speciated PM species
        pm_species = {x: [np.zeros(array_shape) for i in range(config.N_bin)]
                      for x in ['PBC', 'PMD', 'PPOAlP']}

        # Exhaust
        pm_species['PBC'][1] += exhaust_pm10 * config.spec_pm10_exhaust_pbc
        pm_species['PMD'][1] += exhaust_pm10 * config.spec_pm10_exhaust_pmd
        pm_species['PPOAlP'][1] += exhaust_pm10 * config.spec_pm10_exhaust_ppoalp


        # Brake wear
        with open(config.spec_pm10_brake, 'r') as f:
            header = f.readline().replace('\n', '').split(',')
            if len(header)-1 != config.N_bin:
                sys.exit('ERROR: N_bin different in brake_wear.csv.')
            for line in f.readlines():
                line_info = line.replace('\n', '').split(',')
                for i in range(config.N_bin):
                    pm_species[line_info[0]][i] += brake_pm10 * float(line_info[i+1])

        # Tyre wear
        with open(config.spec_pm10_tyre, 'r') as f:
            header = f.readline().replace('\n', '').split(',')
            if len(header)-1 != config.N_bin:
                sys.exit('ERROR: N_bin different in tyre_wear.csv.')
            for line in f.readlines():
                line_info = line.replace('\n', '').split(',')
                for i in range(config.N_bin):
                    pm_species[line_info[0]][i] += tyre_pm10 * float(line_info[i+1])

        # Road wear
        with open(config.spec_pm10_road, 'r') as f:
            header = f.readline().replace('\n', '').split(',')
            if len(header)-1 != config.N_bin:
                sys.exit('ERROR: N_bin different in road_wear.csv.')
            for line in f.readlines():
                line_info = line.replace('\n', '').split(',')
                for i in range(config.N_bin):
                    pm_species[line_info[0]][i] += road_pm10 * float(line_info[i+1])

        # Write files
        for species, value in pm_species.items():
            for i in range(config.N_bin):
                output_filename = input_dir + "/" + species + '_{}.bin'.format(i)
                io.save_binary(value[i], output_filename)
        
def main():

    speciation_voc("sing_preproc.cfg")

    speciation_nox("sing_preproc.cfg")

if __name__ == '__main__':

    main()
