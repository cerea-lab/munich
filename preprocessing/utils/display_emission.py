# coding=utf-8
#!/usr/bin/python
import matplotlib
#matplotlib.use('agg')


# This script calculates the hot emissions of gasoline passenger car
# for each arc during after-work rush hours
import csv
import numpy as np
import pylab as plt
import scipy as sy
import os
import sys
import matplotlib
import re
import matplotlib.lines as mlines

#plt.ioff()

# pre-setting figures for articles

width = 25 / 2.54

matplotlib.rcParams["figure.subplot.left"] = 0.08

WithLegend = True
if WithLegend:
    matplotlib.rcParams["figure.figsize"] = (0.9 * width, .5 * width)
    matplotlib.rcParams["figure.figsize"] = (0.9 * width, width)
    matplotlib.rcParams["figure.subplot.right"] = 0.75
else:
    matplotlib.rcParams["figure.figsize"] = (0.7 * width, width)
    matplotlib.rcParams["figure.subplot.right"] = 0.95

matplotlib.rcParams["figure.subplot.bottom"] = 0.08
matplotlib.rcParams["figure.subplot.top"] = 0.95
font_size = 11
matplotlib.rcParams["font.size"] = font_size
matplotlib.rcParams["axes.titlesize"] = font_size + 2
matplotlib.rcParams["axes.labelsize"] = font_size +2
matplotlib.rcParams["xtick.labelsize"] = font_size
matplotlib.rcParams["ytick.labelsize"] = font_size
matplotlib.rcParams["legend.fontsize"] = font_size
matplotlib.rcParams["xtick.major.size"] = 2
matplotlib.rcParams["xtick.minor.size"] = 1
matplotlib.rcParams["ytick.major.size"] = 2
matplotlib.rcParams["ytick.minor.size"] = 1

# display the emissions map

input_node = open("../output/textfile/intersection.dat")

#
node = {}

for line in input_node.readlines():
    line_info = [x for x in re.split(';| ', line) if x.strip() != '']
    # Converting from text to number.
    node[int(line_info[0])] = (float(line_info[1]), float(line_info[2]))

input_node.close()

#!!!!!!!!!!!!!!!!!#
# User input data # 
#!!!!!!!!!!!!!!!!!#

# date
# see file names in ../output/textfile/emission.DATE.txt
date = "20140316-00"

# default species_ind is 0
# For the species list, please see emission_species in ../sing_preproc.cfg
# For example,
# emission_species: NO NO2
# species_ind = 1 to display NO2 species
species_ind = 0
    
plot_data_visum = np.loadtxt("../output/textfile/emission." + date + ".txt")

plot_arc_size_visum = np.size(plot_data_visum[:,0])
node_begin_visum = []
node_end_visum = []
emission_visum = []
arc_id_visum = []

for i in range(0,plot_arc_size_visum):
    node_begin_visum.append(int(plot_data_visum[i,1]))
    node_end_visum.append(int(plot_data_visum[i,2]))
    emission_visum.append(float(plot_data_visum[i,species_ind]))
    arc_id_visum.append(int(plot_data_visum[i,0]))


# Display intervals
size = 15
vmax = max(emission_visum)
interv = range(0, int(vmax), int(vmax / size))
alpha = 1.0
beta = 1.
interval = []
interval_l = []

# Interval for plotting
for i in range(0,size):
    interval.append(alpha*interv[i])
    interval_l.append(alpha*interv[i]/beta)
s = 1000

    

# PLOTTING
# ---------

w = 2
w_0 = 1


# Plot Visum
# ----------

fig = plt.figure()
visum = plt.subplot(111)

# Display the node numbers.
# ------------------------
disp_node_number = False
disp_street_number = False
if (disp_node_number):
  for inode in node:
      visum.text(node[inode][0], node[inode][1], str(inode), size=10, color='b')

lw_max = 2
for i in range(len(node_begin_visum)):
    xy_begin = node[node_begin_visum[i]]
    xy_end = node[node_end_visum[i]]

    # Display the street numbers.
    if (disp_street_number):
        x_ = (xy_begin[0] + xy_end[0]) / 2.0
        y_ = (xy_begin[1] + xy_end[1]) / 2.0

        visum.text(x_, y_, str(arc_id_visum[i]), size=10, color='r')

    if emission_visum[i] != 1.0 and emission_visum[i] >= interval[0] and emission_visum[i] < interval[1]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'indigo',linestyle = '-',
                     lw = 0.1)
#                     lw = emission_visum[i]/s)
    elif emission_visum[i] >= interval[1] and emission_visum[i] < interval[2]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'midnightblue',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[2] and emission_visum[i] < interval[3]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'navy',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[3] and emission_visum[i] <interval[4] :
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'blue',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[4] and emission_visum[i] < interval[5]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'royalblue',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[5] and emission_visum[i] < interval[6]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'deepskyblue',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[6] and emission_visum[i] < interval[7]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'cyan',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[7] and emission_visum[i] < interval[8]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'lime',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[8] and emission_visum[i] < interval[9]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'yellow',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[9] and emission_visum[i] < interval[10]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'gold',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[10] and emission_visum[i] <interval[11]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'orange',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[11] and emission_visum[i] < interval[12]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'salmon',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[12] and emission_visum[i] < interval[13]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'red',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    elif emission_visum[i] >= interval[13] and emission_visum[i] < interval[14]:
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                 color = 'firebrick',linestyle = '-',
                     lw = lw_max) #emission_visum[i]/s)
    else :
            visum.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color='darkred',linestyle='-',
                     lw = lw_max) #emission_visum[i]/s)

# plot legend
# -----------

intro = mlines.Line2D([],[],color = 'w',linestyle='-',lw = 0.5 )
in_0 = mlines.Line2D([],[],color = 'k',linestyle='-',lw = 0.5 )
in_1 = mlines.Line2D([],[],color = 'indigo',linestyle='-',lw = 0.8 )
in_2 = mlines.Line2D([],[],color = 'midnightblue',linestyle='-',lw = 1 )
in_3 = mlines.Line2D([],[],color = 'navy',linestyle='-',lw = 1.2 )
in_4 = mlines.Line2D([],[],color = 'blue',linestyle='-',lw = 1.4 )
in_5 = mlines.Line2D([],[],color = 'royalblue',linestyle='-',lw = 1.6 )
in_6 = mlines.Line2D([],[],color = 'deepskyblue',linestyle='-',lw = 1.8 )
in_7 = mlines.Line2D([],[],color = 'cyan',linestyle='-',lw = 2 )
in_8 = mlines.Line2D([],[],color = 'lime',linestyle='-',lw = 2.2 )
in_9 = mlines.Line2D([],[],color = 'gold',linestyle='-',lw = 2.4 )
in_10 = mlines.Line2D([],[],color = 'orange',linestyle='-',lw = 2.6 )
in_11 = mlines.Line2D([],[],color = 'salmon',linestyle='-',lw = 2.8 )
in_12 = mlines.Line2D([],[],color = 'tomato',linestyle='-',lw = 3 )
in_13 = mlines.Line2D([],[],color = 'red',linestyle='-',lw = 3.2 )
in_14 = mlines.Line2D([],[],color = 'firebrick',linestyle='-',lw = 3.4 )
in_15 = mlines.Line2D([],[],color = 'darkred',linestyle='-',lw = 3.6 )

handles = [intro, in_0, in_1, in_2, in_3, in_4, in_5, in_6, in_7, in_8, in_9, in_10, in_11, in_12, in_13, in_14, in_15]

labels = ['in $\mu$g/s','0.0', '] '+str(interval_l[0])+', '+str(interval_l[1])+' [','[ '+str(interval_l[1])+', '+str(interval_l[2])+' [','[ '+str(interval_l[2])+', '+str(interval_l[3])+' [', '[ '+str(interval_l[3])+', '+str(interval_l[4])+' [',  '[ '+str(interval_l[4])+', '+str(interval_l[5])+' [','[ '+str(interval_l[5])+', '+str(interval_l[6])+' [','[ '+str(interval_l[6])+', '+str(interval_l[7])+' [','[ '+str(interval_l[7])+', '+str(interval_l[8])+' [','[ '+str(interval_l[8])+', '+str(interval_l[9])+' [','[ '+str(interval_l[9])+', '+str(interval_l[10])+' [', '[ '+str(interval_l[10])+', '+str(interval_l[11])+' [','[ '+str(interval_l[11])+', '+str(interval_l[12])+' [','[ '+str(interval_l[12])+', '+str(interval_l[13])+' [','[ '+str(interval_l[13])+', '+str(interval_l[14])+' [','[ '+str(interval_l[14])+', ~ [' ]

plt.legend(handles,labels,bbox_to_anchor=(1.02,0.02),loc = 'lower left')


plt.subplots_adjust(bottom = 0.1)

a = visum.get_yticks().tolist()
visum.set_yticklabels(a)

plt.show()

