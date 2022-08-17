# coding=utf-8
#!/usr/bin/python3

# import matplotlib as mpl
# mpl.use('agg')
import matplotlib.pyplot as plt

# This script calculates the hot emissions of gasoline passenger car
# for each arc during after-work rush hours
import numpy as np
import os
import sys
import re
from atmopy import *
from atmopy.display import *

content = [("Species", "[input]", "String"),
	   ("Directory", "[input]", "String"),
           ("Date_min", "[input]", "DateTime"),
           ("Delta_t", "[input]", "Float")]

config = talos.Config(sys.argv[1], content)

config.Nx = 1
config.Ny = 1
config.Nt = 0

date_min = config.Date_min
delta_t = config.Delta_t

species = config.Species
directory = config.Directory

# Input data 
# -----------

input_file = open("intersection.dat")
node = {}
input_file.readline()
for line in input_file.readlines():
    line_info = [x for x in re.split(';', line) if x.strip() != '']
    # Converting from text to number.
    node[int(line_info[0])] = (float(line_info[1]), float(line_info[2]))
input_file.close()

node_begin = []
node_end = []
arc_id = []
input_file = open("street.dat")
input_file.readline()
for line in input_file.readlines():
    line_info = [x for x in re.split(';', line) if x.strip() != '']
    node_begin.append(int(line_info[1]))
    node_end.append(int(line_info[2]))
    arc_id.append(int(line_info[0]))
input_file.close()
   
nstreet = len(arc_id)

data = getd(config, directory + species + ".bin")

if (nstreet != data.shape[1]):
    raise Exception("Incompatible data. Please set Nz to ", nstreet)

# Dimension: nt * nstreet
data_reshape = np.reshape(data, (data.shape[0], data.shape[1]))
nt = data.shape[0]

print("Nt: ", nt)
print("Number of streets: ", data.shape[1])
print("Begining date: ", date_min)

## ==========
import datetime
sim_date = []
disp_conc = []

for t in range(nt):
    current_date = date_min + datetime.timedelta(seconds = t * delta_t)
    sim_date.append(current_date)
    street_data = data_reshape[t]
        
# Street fragments for Boulevard Alsace-Lorraine 
    AL_list = [631, 598, 685, 668, 687, 688, 689, 595, 692, 693, 826]
    v_AL = []
    for li in AL_list:
        for i in range(nstreet):
            if li == arc_id[i]:
                v_AL.append(street_data[i])

    disp_conc.append(sum(v_AL) / float(len(v_AL)))

fig, ax = plt.subplots(figsize = (16, 6))

ax.plot(sim_date, disp_conc, 'k-', linewidth = 2, label = "MUNICH")
ax.set_title(species)
print("******* average : ", np.mean(disp_conc))

# Formats the date axis.
locator = HourLocator(interval = 48)
ax.xaxis.set_major_locator(locator)

import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%d/%m %H')
ax.xaxis.set_major_formatter(myFmt)

ltick = ax.get_xticks()
ltick_center = []
ltick_center_lab = []
for i in range(len(ltick)-1):
     ltick_center.append(ltick[i] + 0.5)
     ltick_center_lab.append(ltick[i] + 0.5)

loc = 0.45,0.7
leg = legend(loc=1, shadow=True)
frame = leg.get_frame()
frame.set_facecolor('0.80')
for t in leg.get_texts():
	t.set_fontsize('12')
for l in leg.get_lines():
	l.set_linewidth(1.5)

plt.grid()
plt.show()
