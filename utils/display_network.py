#!/usr/bin/env python3

# Author: Thibaud SARICA
# Date: 2021-04-02

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

"""
Plot street network with options for node id and/or street id.
"""

##### CONFIGURATION

# These are the files used by MUNICH
indir = '../processing/ssh-aerosol/'
node_file = indir + 'node.dat'
street_file = indir + 'street.dat'

# Display or save plot
save_plot = False # True/False
outfile = 'plot.png'

# Display node id and/or street id
display_node_id = False # True/False
display_street_id = False # True/False

# Size of node and street id text
size = 10

# Colors of node id, street id and street
color_node_id = 'b'
color_street_id = 'r'
color_street = 'g'

##### MAIN

# Read node file
nodes = {} # dict of tuples id: (lon, lat)
with open(node_file, 'r') as f:
    for line in f.readlines():
        line_info = line.replace('\n', '').split(';')
        nodes[int(line_info[0])] = (float(line_info[1]), float(line_info[2]))

# Read street file
streets = [] # list of tuples (id, begin, end)
with open(street_file, 'r') as f:
    for line in f.readlines():
        line_info = line.replace('\n', '').split(';')
        streets.append((int(line_info[0]), int(line_info[1]), int(line_info[2])))

# Plot
fig = plt.figure()
ax = plt.subplot(111)

# Display node id
if display_node_id:
    for i, node in nodes.items():
        ax.text(node[0], node[1], str(i), size=size, color=color_node_id)

# Display street id
if display_street_id:
    for i in range(len(streets)):
        lon1 = nodes[streets[i][1]][0]
        lat1 = nodes[streets[i][1]][1]
        lon2 = nodes[streets[i][2]][0]
        lat2 = nodes[streets[i][2]][1]
        lon_cen = (lon1 + lon2) * 0.5
        lat_cen = (lat1 + lat2) * 0.5
        ax.text(lon_cen, lat_cen, str(i), size=size, color=color_street_id)

# Display streets
for i in range(len(streets)):
    ax.plot([nodes[streets[i][1]][0], nodes[streets[i][2]][0]],
            [nodes[streets[i][1]][1], nodes[streets[i][2]][1]],
            lw=0.5, color=color_street)

if save_plot:
    plt.savefig(outfile)
else:
    plt.show()
