# coding=utf-8
#!/usr/bin/python

# import matplotlib as mpl
# mpl.use('agg')
import matplotlib.pyplot as plt
from atmopy.display import *

# This script calculates the hot emissions of gasoline passenger car
# for each arc during after-work rush hours
import numpy as np
import os
import sys
import re
from atmopy import *
from atmopy.display import *

content = [("Species", "[input]", "StringList"),
	   ("Directory", "[input]", "String")]

config = talos.Config(sys.argv[1], content)

shape = (config.Ny, config.Nx)

input_node = open("node.txt")
node = {}

for line in input_node.readlines():
    line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
    # Converting from text to number.
    node[int(line_info[0])] = (float(line_info[1]), float(line_info[2]))

input_node.close()

# Input data 
# -----------

plot_data = np.loadtxt("emission.txt")
plot_arc_size = np.size(plot_data[:,0])
node_begin = []
node_end = []
arc_id = []

for i in range(0,plot_arc_size):
    node_begin.append(int(plot_data[i,1]))
    node_end.append(int(plot_data[i,2]))
    arc_id.append(int(plot_data[i,0]))

species = config.Species[0]

NOxratio = False

# Simulation 1
directory = "/net/libre/coba/kimy/munich-testcase/dynamic/results/"
if species == "NOx":
    data = getd(config, directory + "NO.bin") * 46. / 30.
    data2 = getd(config, directory + "NO2.bin")
    
else:
    data = getd(config, directory + species + ".bin")
    data2 = 0.0 
if (NOxratio):
    data_sum = data2 / (data + data2)
else:
    data_sum = data + data2
data_reshape = np.reshape(data_sum, (data_sum.shape[0], data_sum.shape[1]))
print data_reshape.shape

# Simulation 2
directory = "/net/libre/coba/kimy/munich-testcase/stationary/results/"
if species == "NOx":
    data = getd(config, directory + "NO.bin") * 46. / 30.
    data2 = getd(config, directory + "NO2.bin")
    
else:
    data = getd(config, directory + species + ".bin")
    data2 = 0.0 
if (NOxratio):
    data_sum = data2 / (data + data2)
else:
    data_sum = data + data2
data_reshape2 = np.reshape(data_sum, (data_sum.shape[0], data_sum.shape[1]))
print data_reshape2.shape


## ==========
import datetime
date_begin = datetime.datetime(2014,3,16,1)
sim_date = []
disp_conc = []
disp_conc2 = []

for t in range(data_sum.shape[0]):
    current_date = date_begin + datetime.timedelta(hours = t)
    sim_date.append(current_date)
    street_data = []
    street_data2 = [] 

    # ID, O3, NO, NO2
    for i in range(0,plot_arc_size):
        street_data.append(data_reshape[t, i])
        street_data2.append(data_reshape2[t, i])

# Street fragments for Boulevard Alsace-Lorraine 
    AL_list = [631, 598, 685, 668, 687, 688, 689, 595, 692, 693, 826]
    v_AL = []
    v_AL2 = []
    for li in AL_list:
        for i in range(0,plot_arc_size):
            if li == arc_id[i]:
                v_AL.append(street_data[i])
                v_AL2.append(street_data2[i])

    disp_conc.append(sum(v_AL) / float(len(v_AL)))
    disp_conc2.append(sum(v_AL2) / float(len(v_AL)))

fig, ax = plt.subplots(figsize = (16, 6))

ax.plot(sim_date, disp_conc, 'k-', linewidth = 2, label = "Dynamic")
ax.plot(sim_date, disp_conc2, 'r-', linewidth = 2, label = "Stationary")
ax.set_title(species)
print "******* average : ", np.mean(disp_conc), np.mean(disp_conc2)

# Formats the date axis.
locator = HourLocator(interval = 24)
ax.xaxis.set_major_locator(locator)

import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%d/%m')
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
