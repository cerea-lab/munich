#! /usr/bin/env python

import re, sys, os
from atmopy import *
from optparse import OptionParser
import numpy as np


################
# Main program #
################


home_dir = "/profils_cerea/kimy/work/StreetInGrid/munich/preprocessing/output/"
input_dir = home_dir + "/emission/"
input_file = "NMHC.bin"

speciation_file = "COVNM.dat"
speciation = open(speciation_file)
header = speciation.readline()

# Get snap code for traffic source.
snap_line = speciation.readline()
line_info = [x for x in re.split('\t| ', snap_line) if x.strip() != '']
for n, line in enumerate(line_info):
        if line[0] == "7":
                index = n
                break
print "traffic snap index: ", index

sp_real_name = []
speciation_coeff = []
mw_real = []
total = 0.0
for line in speciation.readlines():
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        sp_real_name.append(line_info[0])
        speciation_coeff.append(float(line_info[index]))
        total += float(line_info[index])
        mw_real.append(float(line_info[1]))
speciation.close()

print "total coeff.: ", total, "%"
ns_real = len(sp_real_name)
print "Number of the real species: ", ns_real

# Get VOC aggregation to model species.
aggregation_file = "aggregation_cb05-siream.dat"
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
species_factor = np.zeros([ns_model], 'float')
for s_real in range(ns_real):
        line = aggregation.readline()
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        for s_model in range(ns_model):
                aggregation_coeff = float(line_info[s_model + 3])
                temp = speciation_coeff[s_real] * 0.01 * aggregation_coeff * mw_model[s_model] / mw_real[s_real]
                species_factor[s_model] += temp
aggregation.close()


os.chdir(input_dir)
total = 0.0
for s_model in range(ns_model):
        total += species_factor[s_model]
        print "species factor: ",model_species[s_model], species_factor[s_model]
        command = "mult_nb_float " + input_file + " " + str(species_factor[s_model]) + " " + model_species[s_model] + ".bin"
        print command
        os.system(command)
print "total species factor:", total

# Grided data
input_dir = home_dir
os.chdir(input_dir)
for s_model in range(ns_model):
        total += species_factor[s_model]
        print "species factor: ",model_species[s_model], species_factor[s_model]
        command = "mult_nb_float " + input_file + " " + str(species_factor[s_model]) + " " + model_species[s_model] + ".bin"
        print command
        os.system(command)
