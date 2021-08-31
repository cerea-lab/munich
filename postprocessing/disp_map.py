# coding=utf-8
#!/usr/bin/python
import matplotlib
#matplotlib.use('agg')

import csv
import numpy as np
import pylab as plt
import scipy as sy
import os
import sys
import matplotlib
import re
import matplotlib.lines as mlines
import datetime
from atmopy import *
from atmopy.display import *

content = [("Species", "[input]", "StringList"),\
	       ("Directory", "[input]", "String"),\
           ("Date_begin", "[input]", "DateTime"), \
           ("Delta_t", "[input]", "Float"), \
	       ("file_street", "[input]", "String"), \
	       ("file_intersection", "[input]", "String"), \
           ("disp_node_number", "[option]", "Bool"), \
           ("disp_street_number", "[option]", "Bool")
]

config = talos.Config(sys.argv[1], content)
# shape = (config.Ny, config.Nx)

# pre-setting figures for articles

width = 25 / 2.54

matplotlib.rcParams["figure.figsize"] = (0.7 * width, width)
matplotlib.rcParams["figure.subplot.left"] = 0.10
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

# Read intersection location data.
# ------------------

input_file = open(config.file_intersection)
node = {}
input_file.readline()
for line in input_file.readlines():
    line_info = [x for x in re.split(';', line) if x.strip() != '']
    # Converting from text to number.
    node[int(line_info[0])] = (float(line_info[1]), float(line_info[2]))
input_file.close()

# Input street location data
# ------------------

node_begin = []
node_end = []
street_id = []
input_file = open(config.file_street)
input_file.readline()
for line in input_file.readlines():
    line_info = [x for x in re.split(';', line) if x.strip() != '']
    node_begin.append(int(line_info[1]))
    node_end.append(int(line_info[2]))
    street_id.append(int(line_info[0]))
input_file.close()

nstreet = len(street_id)

# Read binary data
# ------------------

# species = config.Species
# directory = config.Directory
# data = getd(config, directory + species[0] + ".bin")	
# data2 = np.reshape(data, (data.shape[0], data.shape[1]))


# date_begin = config.Date_begin
# t_ind = 1
# current_date = date_begin + datetime.timedelta(hours = t_ind)
street_concentration = []
# hour = current_date.strftime("%Y%m%d_%H")

for i in range(0,nstreet):
    street_concentration.append(10)


# Display options
vmin = 0
vmax = int(max(street_concentration))
ninterv = 15
step = (vmax - vmin) / ninterv
interv = np.arange(vmin,vmax,step)

alpha = 1. # 0.5
beta = 1.
size = np.size(interv)

interval = []
interval_l = []
# Interval for plotting
for i in range(0,size):
    interval.append(alpha*interv[i])
    interval_l.append(alpha*interv[i]/beta)

# Set the line thickness
s = 800

# PLOTTING
# ---------

w = 2
w_0 = 1


fig = plt.figure()
ax = plt.subplot(111)


# Options to display the node/street id.
# ------------------------
disp_node_number = False
disp_street_number = False

# Display the node id
if (disp_node_number):
  for inode in node:
      ax.text(node[inode][0], node[inode][1], str(inode), size=10, color='b')


lw_max = 1.5
for i in range(len(node_begin)):
    xy_begin = node[node_begin[i]]
    xy_end = node[node_end[i]]


    # Display the street id
    if (disp_street_number):
        x_ = (xy_begin[0] + xy_end[0]) / 2.0
        y_ = (xy_begin[1] + xy_end[1]) / 2.0
        ax.text(x_, y_, str(street_id[i]), size=10, color='r')        
    
    if street_concentration[i] >= interval[0] and street_concentration[i] < interval[1]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'indigo',linestyle = '-',
                     lw = 0.8)
#                     lw = street_concentration[i]/s)
    elif street_concentration[i] >= interval[1] and street_concentration[i] < interval[2]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'midnightblue',linestyle = '-',
                     lw = 1) #street_concentration[i]/s)
    elif street_concentration[i] >= interval[2] and street_concentration[i] < interval[3]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'navy',linestyle = '-',
                     lw = lw_max) #street_concentration[i]/s)
    elif street_concentration[i] >= interval[3] and street_concentration[i] <interval[4] :
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'blue',linestyle = '-',
                     lw = lw_max) # street_concentration[i]/s)
    elif street_concentration[i] >= interval[4] and street_concentration[i] < interval[5]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'royalblue',linestyle = '-',
                     lw = lw_max) # street_concentration[i]/s)
    elif street_concentration[i] >= interval[5] and street_concentration[i] < interval[6]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'deepskyblue',linestyle = '-',
                     lw = lw_max) # street_concentration[i]/s)
    elif street_concentration[i] >= interval[6] and street_concentration[i] < interval[7]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'cyan',linestyle = '-',
                     lw = lw_max) #street_concentration[i]/s)
    elif street_concentration[i] >= interval[7] and street_concentration[i] < interval[8]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'lime',linestyle = '-',
                     lw = lw_max) #street_concentration[i]/s)
    elif street_concentration[i] >= interval[8] and street_concentration[i] < interval[9]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'yellow',linestyle = '-',
                     lw = lw_max) # street_concentration[i]/s)
    elif street_concentration[i] >= interval[9] and street_concentration[i] < interval[10]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'gold',linestyle = '-',
                     lw = lw_max) # street_concentration[i]/s)
    elif street_concentration[i] >= interval[10] and street_concentration[i] <interval[11]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'orange',linestyle = '-',
                     lw = lw_max) # street_concentration[i]/s)
    elif street_concentration[i] >= interval[11] and street_concentration[i] < interval[12]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'salmon',linestyle = '-',
                     lw = lw_max) #street_concentration[i]/s)
    elif street_concentration[i] >= interval[12] and street_concentration[i] < interval[13]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color = 'red',linestyle = '-',
                     lw = lw_max) # street_concentration[i]/s)
    elif street_concentration[i] >= interval[13] and street_concentration[i] < interval[14]:
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                 color = 'firebrick',linestyle = '-',
                     lw = lw_max) #street_concentration[i]/s)
    else :
            ax.plot([xy_begin[0], xy_end[0]],
                     [xy_begin[1], xy_end[1]],
                     color='darkred',linestyle='-',
                     lw = lw_max) #street_concentration[i]/s)

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

#handles = [intro, in_0, in_1, in_2, in_3, in_4, in_5, in_6, in_7, in_8, in_9, in_10, in_11, in_12, in_13, in_14, in_15]
handles = [intro, in_1, in_2, in_3, in_4, in_5, in_6, in_7, in_8, in_9, in_10, in_11, in_12, in_13, in_14, in_15]

labels = ['in $\mu$g/m$^3$',\
#          '0.0',\
          '] '+str(interval_l[0])+', '+str(interval_l[1])+' [',\
          '[ '+str(interval_l[1])+', '+str(interval_l[2])+' [',\
          '[ '+str(interval_l[2])+', '+str(interval_l[3])+' [',\
          '[ '+str(interval_l[3])+', '+str(interval_l[4])+' [',\
          '[ '+str(interval_l[4])+', '+str(interval_l[5])+' [',\
          '[ '+str(interval_l[5])+', '+str(interval_l[6])+' [',\
          '[ '+str(interval_l[6])+', '+str(interval_l[7])+' [',\
          '[ '+str(interval_l[7])+', '+str(interval_l[8])+' [',\
          '[ '+str(interval_l[8])+', '+str(interval_l[9])+' [',\
          '[ '+str(interval_l[9])+', '+str(interval_l[10])+' [',\
          '[ '+str(interval_l[10])+', '+str(interval_l[11])+' [',\
          '[ '+str(interval_l[11])+', '+str(interval_l[12])+' [',\
          '[ '+str(interval_l[12])+', '+str(interval_l[13])+' [',\
          '[ '+str(interval_l[13])+', '+str(interval_l[14])+' [',\
          '[ '+str(interval_l[14])+', ~ [' ]

# bbox_to_anchor(0.5,-0.1)
#plt.legend(handles,labels,bbox_to_anchor=(1.05,0.0),loc = 'lower left')
plt.subplots_adjust(bottom = 0.1)

# a = ax.get_yticks().tolist()
# ax.set_yticklabels(a)
# a = [48.835, 48.84, 48.845, 48.85, 48.855, 48.86, 48.865, 48.87]
# ax.set_yticklabels(a)

# coor_y = np.arange(48.840, 48.870, 0.01)
# for y_ in coor_y:
#     l = plt.axhline(y = y_, color = 'k')
# #coor_x = np.arange(2.4763, 2.5263, 0.01)
# coor_x = np.arange(2.48, 2.52, 0.01)
# for x_ in coor_x:
#     l = plt.axvline(x = x_, color = 'k')


# ax.set_title(species[0])

plt.show()

