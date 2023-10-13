import re, sys, os, shutil
import datetime
import numpy as np
from misc import *

##### Class for a node (=intersection)
class Node:
    def __init__(self, node_id, lon, lat):
        self.id = node_id
        self.eff_id = 0 # Starts automatically at 1
        self.removed = False
        self.lon = lon
        self.lat = lat
        self.connected_street = []
        # meteo
        self.lmo = 0.
        self.pblh = 0.
        self.ust = 0.
        self.wdir = 0.
        self.wspd = 0.

##### Class for a street
class Street:
    def __init__(self, street_id, begin_node, end_node,
                 length, width, height, lon_cen, lat_cen,
                 emission):
        self.id = street_id
        self.begin = begin_node
        self.end = end_node
        self.eff_id = street_id
        self.eff_begin = begin_node
        self.eff_end = end_node
        self.removed = False
        self.length = length
        self.width = width
        self.height = height
        self.lon_cen = lon_cen
        self.lat_cen = lat_cen
        self.typo = 0
        # emission
        self.emission = emission
        # meteo
        self.attenuation = 0.
        self.liquidwaterc = 0.
        self.lmo = 0.
        self.pardb = 0.
        self.pardiff = 0.
        self.pblh = 0.
        self.rain = 0.
        self.soilwater = 0.
        self.solarrad = 0.
        self.spechumid = 0.
        self.surfpressure = 0.
        self.surfrichard = 0.
        self.surftemp = 0.
        self.ust = 0.
        self.wdir = 0.
        self.wspd = 0.
        # background
        self.bkgd = {} # in µg/m3

# ---------------------------
# Check if two nodes are the same or near
def are_nodes_same(node1, node2, min_distance):
    if node1.lon == node2.lon and node1.lat == node2.lat:
        return True
    else:
        geod = pyproj.Geod(ellps='WGS84')
        _, _, length = geod.inv(node1.lon, node1.lat, node2.lon, node2.lat) # in meter
        if length < min_distance:
            return True
        else:
            return False

# Check if two streets are the same or near based on their nodes
def are_streets_same(node_list, i, j, min_distance):
    begin1 = 2 * i
    end1 = 2 * i + 1
    begin2 = 2 * j
    end2 = 2 * j + 1
    # begin1 == begin2 and end1 == end2
    if are_nodes_same(node_list[begin1], node_list[begin2], min_distance) and\
       are_nodes_same(node_list[end1], node_list[end2], min_distance):
        return True
    # begin1 == end2 and end1 == begin2
    elif are_nodes_same(node_list[begin1], node_list[end2], min_distance) and\
         are_nodes_same(node_list[end1], node_list[begin2], min_distance):
        return True
    else:
        return False

    
# Merge streets in the case where there are two emission data sets for both ways. 
# --------------------

def merging_street(output_file, node_list, street_list,
                   species_list):

    ntemp = 0
    n_street = len(street_list)
    for i in range(n_street - 1):
      if (street_list[i].removed == False):
        is_street_found = False
        j = i + 1
        while (is_street_found == False and j < n_street):
            # if are_streets_same(node_list, street_list[i], street_list[j]):
            if are_streets_same(node_list, i, j):                
                street_list[j].eff_begin = street_list[i].begin
                street_list[j].eff_end = street_list[i].end
                street_list[j].eff_id = street_list[i].id
              
                for spe in species_list:
                    street_list[i].emission[spe] = street_list[i].emission[spe] + street_list[j].emission[spe]
                    street_list[j].emission[spe] = 0.0
                street_list[j].removed = True
                ntemp = ntemp + 1
                is_street_found = True
            j = j + 1

    f = open(output_file, 'w')
    f.write("# Street id \t Effective street id \n")
    for i in range(n_street):
        street = street_list[i]
        f.write(str(street.id) + "\t" + str(street.eff_id) + "\n")
    f.close()

    return ntemp                

# Manual merging of intersections.
# -------------------------------
def manual_merging_street(street_list, species_list,
                          input_file):
#    input_file = "street-merging.txt"
    input_merging = open(input_file)
    print("Manual merging using the street list in " + input_file)
    header = input_merging.readline()
    ntemp = 0
    for line in input_merging.readlines():
        line = line.replace('\n','')
        line_info = [x for x in re.split(',', line) if x.strip() != '']
        if len(line_info) != 2:
            break
        else:
            removed_street_id, remained_street_id = int(line_info[0]), int(line_info[1])
            for i in range(len(street_list)):
                if (street_list[i].id == remained_street_id):
                    for j in range(len(street_list)):
                        if (street_list[j].id == removed_street_id):
                            street_list[j].eff_id = street_list[i].id
                            street_list[j].eff_begin = street_list[i].begin
                            street_list[j].eff_end = street_list[i].end
                            for spe in species_list:
                                street_list[i].emission[spe] = street_list[i].emission[spe] + street_list[j].emission[spe]
                                street_list[j].emission[spe] = 0.0
                            street_list[j].removed = True
                            ntemp = ntemp + 1
    input_merging.close()
    return ntemp

# Merging using a look-up table.
# -------------------------------
def lut_merging_street(lut_file, street_list,
                       species_list):
    n_street = len(street_list)
    ntemp = 0
    print("Read the lookup-table: ", lut_file)    
    input_merging = open(lut_file)
    header = input_merging.readline()
    for line in input_merging.readlines():
        line = line.replace('\n','')
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        if len(line_info) != 2:
            break
        else:
            street_id, effective_street_id = int(line_info[0]), int(line_info[1])
            if street_id != effective_street_id:
                for nst in range(n_street):
                    if (street_list[nst].id == effective_street_id):
                        i = nst 
                    elif (street_list[nst].id == street_id):
                        # the emissions in this street are merged into
                        # the street having the index i.
                        j = nst 
                street_list[j].eff_begin = street_list[i].begin
                street_list[j].eff_end = street_list[i].end
                street_list[j].eff_id = street_list[i].id

                for spe in species_list:
                    street_list[i].emission[spe] = street_list[i].emission[spe] + street_list[j].emission[spe]
                    street_list[j].emission[spe] = 0.0
                street_list[j].removed = True
                ntemp = ntemp + 1
    input_merging.close()
    return ntemp

# Remove the same nodes
# --------------------
def merging_node(node_list):
    n_node = 0
    for i in range(len(node_list) - 1):
      if (node_list[i].removed == False):  
        for j in range(i + 1, len(node_list)):
          if (node_list[j].removed == False):  
            if ((node_list[i].lon == node_list[j].lon) and 
                (node_list[i].lat == node_list[j].lat)) :
                node_list[j].eff_id = node_list[i].eff_id
                node_list[j].removed = True
                n_node = n_node + 1
                for st_ in node_list[j].connected_street:
                    node_list[i].connected_street.append(st_)
    return n_node


# Remove the near nodes: 
# the distance between the nodes 
# is smaller than 10 m
# --------------------
def merging_near_node(node_list):
    n_node = 0
    for i in range(len(node_list) - 1):
      if (node_list[i].removed == False):  
        lon1, lat1 = node_list[i].lon, node_list[i].lat
        for j in range(i + 1, len(node_list)):
          if (node_list[j].removed == False):  
            lon2, lat2 = node_list[j].lon, node_list[j].lat
            if ((lon1 != lon2) and (lat1 != lat2)):
                length = distance_on_unit_sphere(lat1, lon1, lat2, lon2) / 1000. # in km
                if ((length < 0.01)):
                    id_old = node_list[j].eff_id
                    node_list[j].eff_id = node_list[i].eff_id
                    id_new = node_list[j].eff_id
                    node_list[j].removed = True
                    n_node = n_node + 1
                    for st_ in node_list[j].connected_street:
                        # check if st_ is an element of node_list[i]
                        if (st_ in node_list[i].connected_street):
                            sys.exit("Error: street " + str(st_) + " is " \
                                     "already in node_list.")
                        else:
                            node_list[i].connected_street.append(st_)

                    # Node A was merged to Node B.
                    # Node B was merged to Node C.
                    # Node A should be merged to Node C.
                    for node_ in node_list:
                        if (node_.eff_id == id_old):
                            node_.eff_id = id_new
    return n_node

# Manual merging of intersections.
# -------------------------------
def manual_merging_node(node_list, input_file):
    n_node = 0
#    input_file = "intersection-merging.txt"
    print("Manual merging of the nodes using the node list in " + input_file)
    input_merging = open(input_file)
    header = input_merging.readline()
    for line in input_merging.readlines():
        line = line.replace('\n','')
        line_info = [x for x in re.split(',', line) if x.strip() != '']
        if len(line_info) != 2:
            break
        else:
            removed_node_id, remained_node_id = int(line_info[0]), int(line_info[1])
            for i in range(len(node_list)):
                if (node_list[i].id == removed_node_id): 
                    for j in range(len(node_list)):
                        if (node_list[j].id == remained_node_id): 
                            id_old = node_list[i].eff_id
                            node_list[i].eff_id = node_list[j].eff_id
                            id_new = node_list[i].eff_id
                            node_list[i].removed = True
                            n_node = n_node + 1
                            for st_ in node_list[i].connected_street:
                                node_list[j].connected_street.append(st_)
                            #     
                            for node_ in node_list:
                                if (node_.eff_id == id_old):
                                    node_.eff_id = id_new
    input_merging.close()
    return n_node

# Get the street width and the builiding height from the input file.
# ------------------------------------------------------------------
def get_street_geog(input_file, street_list):
    # Read street width, length and builiding height.
    print("Read the geographical informations: ", input_file)
    input_street_geog = open(input_file)
    header = input_street_geog.readline()
    nstreet = 0
    for line in input_street_geog.readlines():
        line = line.replace('\n','')        
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        street_id = int(line_info[0])    
        for i in range(len(street_list)):
            street = street_list[i]
            if street.id == street_id:
                street.width = float(line_info[2])
                street.height = float(line_info[3])
                nstreet = nstreet + 1
    input_street_geog.close()
    if nstreet != len(street_list):
        print("Warning: Missing input data: data given for %d streets, but data are needed for %d streets." % (nstreet, len(street_list)))
    return 0



def read_traffic_data(input_file, emis_species_list, epsg_code):

    #
    # Define the projections input/output
    # -----------------------------------

    # For convertion to WGS84
    transformer = pyproj.Transformer.from_crs(epsg_code, 4326, always_xy=True)

    # Define Geod to compute distance on Earth and midpoint
    geod = pyproj.Geod(ellps='WGS84')
    
    print("=================", input_file)
    node_id = 0
    street_id = 0
    input_street = open(input_file)
    header = input_street.readline()
    header = header.strip("\n")
    header_info = [x for x in re.split('\t| ', header) if x.strip() != '']

    species_ind = []
    species_list = []
    for emis_species in emis_species_list:
        species_check = False
        for inds, species_name in enumerate(header_info):
            if emis_species == species_name:
                print(emis_species + " found with the index:", inds) 
                species_ind.append(inds)
                species_list.append(species_name)
                species_check = True
        if (species_check == False):
            sys.exit(emis_species + " is not found in " + input_file)
                
    node_list = []
    street_list = []
    for line in input_street.readlines():
        line = line.replace('\n','')
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        street_id = int(line_info[0])
        node_id = node_id + 1
        id_begin = node_id
        x = float(line_info[3])
        y = float(line_info[4])
        lon1, lat1 = transformer.transform(x, y)
        node = Node(node_id, lon1, lat1)
        node.connected_street.append(street_id)
        node_list.append(node)
        node_id = node_id + 1
        id_end = node_id
        x = float(line_info[5])
        y = float(line_info[6])
        lon2, lat2 = transformer.transform(x, y)
        node = Node(node_id, lon2, lat2)
        node.connected_street.append(street_id)
        node_list.append(node)

        # Check input data.
        if ((lon1 == lon2) and (lat1 == lat2)):
            sys.exit("Error: a street has two same intersection coordinate " + \
                     "for the node " + str(node_id) + ", lon: " + str(lon1) + \
                     ", lat: " + str(lat1))
        
        # Street length
        length = distance_on_unit_sphere(lat1, lon1, lat2, lon2) # in meter
        width = 0.0
        height = 0.0
        lon_cen = (lon1 + lon2) * 0.5
        lat_cen = (lat1 + lat2) * 0.5

        # Conversion of the unit of the emission input data
        # from ug/km/h to ug/s
#        emission = np.zeros([len(species_ind)], 'float')
        emission = {}
        for i, ind in enumerate(species_ind):
            emission[species_list[i]] = float(line_info[species_ind[i]]) * (length / 1000.) / 3600.0

        street = Street(street_id, id_begin, id_end, length, width, height, lon_cen, lat_cen, emission)        

        street_list.append(street)
    print(" --------------------------")
    print(" --- Number of the nodes: ", len(node_list))
    print(" --- Number of the streets: ", len(street_list))
    print(" --------------------------")
    input_street.close()
    return street_list, node_list    




def is_holiday(date, country_code):
    # Requirement: holidays python library (https://pypi.org/project/holidays/)
    # Install: pip install holidays
    try:
        import holidays
        isCountryFound = False
        for country in holidays.list_supported_countries():
            if country_code == country:
                print(('Found country code "{}"'.format(country_code)))
                country_holidays = holidays.CountryHoliday(country_code)
                isCountryFound = True
                break
        
        if (isCountryFound == False):
            print(('Error: given country code "{}" is not found in supported country codes'.format(country_code)))
            print(holidays.list_supported_countries())
            return False
        
        return date.date() in country_holidays
        
    except ImportError as err:
        print('*** Could not import "holidays" Python module ***')
        print('Please install "holidays", or Holidays are not taken into account.')
        print('To install, type pip install holidays')
        return False
  
    


def utc_to_local(utc, zone = 'Europe/Paris'):
    # Need to install dateutil
    # pip install python-dateutil
    from dateutil import tz
    
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz(zone)

    # Tell the datetime object that it's in UTC time zone since 
    # datetime objects are 'naive' by default
    utc = utc.replace(tzinfo=from_zone)

    # Convert time zone
    return utc.astimezone(to_zone)


# def projection_conversion(epsg_code):
# #
# # Define the projections input/output
# # -----------------------------------
#     import pyproj
#     wgs84 = pyproj.Proj('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
#     lambert93 = pyproj.Proj('+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')

def read_street_data(street_file, epsg_code):
    # For convertion to WGS84
    transformer = pyproj.Transformer.from_crs(epsg_code, 4326, always_xy=True)

    # Define Geod to compute distance on Earth and midpoint
    geod = pyproj.Geod(ellps='WGS84')

    # Read input file and create street/node lists
    node_id = 1
    node_list = []
    street_list = []
    with open(street_file, 'r') as f:
        header = f.readline()
        for line in f.readlines():
            line_info = line.replace('\n', '').split(',')
            street_id = int(line_info[0])
            id_begin = node_id
            node_id += 1
            x1 = float(line_info[1])
            y1 = float(line_info[2])
            lon1, lat1 = transformer.transform(x1, y1)
            node = Node(id_begin, lon1, lat1)
            node.connected_street.append(street_id)
            node_list.append(node)
            id_end = node_id
            node_id += 1
            x2 = float(line_info[3])
            y2 = float(line_info[4])
            lon2, lat2 = transformer.transform(x2, y2)
            node = Node(id_end, lon2, lat2)
            node.connected_street.append(street_id)
            node_list.append(node)
            # Check if nodes are superimposed
            if lon1 == lon2 and lat1 == lat2:
                sys.exit('ERROR: nodes of street {} are superimposed.'.format(street_id))
            # Street morphology
            _, _, length = geod.inv(lon1, lat1, lon2, lat2) # in meter
            width = float(line_info[5])
            height = float(line_info[6])
            mid_point = geod.npts(lon1, lat1, lon2, lat2, 1)
            street = Street(street_id, id_begin, id_end, length, width, height,
                            mid_point[0][0], mid_point[0][1], 0)
            street_list.append(street)

    return node_list, street_list


def create_network(street_file, epsg_code):

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


# Automatically merge streets assuming there are two
# emissions sets for both ways of a street.
def auto_street_merging(node_list, street_list, min_distance):
    n_node = len(node_list)
    n_street = len(street_list)
    n_removed_street = 0
    for i in range(n_street - 1):
        if street_list[i].removed == False:
            street_found = False
            j = i + 1
            while street_found == False and j < n_street:
                if are_streets_same(node_list, i, j, min_distance):
                    street_found = True
                    street_list[j].removed = True
                    street_list[j].eff_id = street_list[i].id
                    street_list[j].eff_begin = street_list[i].begin
                    street_list[j].eff_end = street_list[i].end
                    n_removed_street += 1
                j += 1

    return n_removed_street

# Manually merge streets
def manual_street_merging(street_list, street_file):
    n_street = len(street_list)
    n_removed_street = 0
    with open(street_file, 'r') as f:
        header = f.readline()
        for line in f.readlines():
            line_info = line.replace('\n', '').split(',')
            if len(line_info) != 2:
                sys.exit('ERROR: more than 2 columns in manual street merging file.')
            removed_id = int(line_info[0])
            remained_id = int(line_info[1])
            for i in range(n_street):
                if street_list[i].id == remained_id:
                    if street_list[i].removed:
                        sys.exit('ERROR: street {} has been automatically removed.'\
                                     .format(street_list[i].id))
                    for j in range(n_street):
                        if street_list[j].id == removed_id:
                            if street_list[j].removed:
                                sys.exit('ERROR: street {} has been automatically removed.'\
                                         .format(street_list[j].id))
                            street_list[j].removed = True
                            street_list[j].eff_id = street_list[i].id
                            street_list[j].eff_begin = street_list[i].begin
                            street_list[j].eff_end = street_list[i].end
                            n_removed_street += 1

    return n_removed_street


# Automatically merge same and near nodes
def auto_node_merging(node_list, min_distance):
    n_node = len(node_list)
    n_removed_node = 0
    for i in range(n_node - 1):
        if node_list[i].removed == False:
            for j in range(i + 1, n_node):
                if node_list[j].removed == False:
                    if are_nodes_same(node_list[i], node_list[j], min_distance):
                        node_list[j].removed = True
                        node_list[j].eff_id = node_list[i].id
                        for c_street in node_list[j].connected_street:
                            node_list[i].connected_street.append(c_street)
                        for k in range(n_node):
                            if node_list[k].eff_id == node_list[j].id:
                                node_list[k].eff_id = node_list[i].id
                        n_removed_node += 1

    return n_removed_node


# Manually merge nodes
def manual_node_merging(node_list, node_file):
    n_node = len(node_list)
    n_removed_node = 0
    with open(node_file, 'r') as f:
        header = f.readline()
        for line in f.readlines():
            line_info = line.replace('\n', '').split(',')
            if len(line_info) != 2:
                sys.exit('ERROR: more than 2 columns in manual node merging file.')
            removed_id = int(line_info[0])
            remained_id = int(line_info[1])
            for i in range(n_node):
                if node_list[i].id == remained_id:
                    if node_list[i].removed:
                        sys.exit('ERROR: node {} has been automatically removed.'\
                                     .format(node_list[i].id))
                    for j in range(n_node):
                        if node_list[j].id == removed_id:
                            if node_list[j].removed:
                                sys.exit('ERROR: node {} has been automatically removed.'\
                                         .format(node_list[j].id))
                            node_list[j].removed = True
                            node_list[j].eff_id = node_list[i].id
                            for c_street in node_list[j].connected_street:
                                node_list[i].connected_street.append(c_street)
                            for k in range(n_node):
                                if node_list[k].eff_id == node_list[j].id:
                                    node_list[k].eff_id = node_list[i].id
                            n_removed_node += 1

    return n_removed_node


# Write LUTs
def write_luts(node_list, street_list, outdir):
    # Write node LUT
    node_file = outdir + 'node_lut.csv'
    with open(node_file, 'w') as f:
        f.write('id,eff_id,lon,lat,N_street,[list_street]\n')
        for i in range(len(node_list)):
            node_id = node_list[i].id
            eff_id = node_list[i].eff_id
            lon = node_list[i].lon
            lat = node_list[i].lat
            n_c_street = len(node_list[i].connected_street)
            c_streets = ''
            for c_street in node_list[i].connected_street:
                c_streets += ',' + str(c_street)
            f.write(str(node_id) + ',' + str(eff_id) + ',' + str(lon) + ',' + str(lat)\
                    + ',' + str(n_c_street) + c_streets + '\n')

    # Write street LUT
    street_file = outdir + 'street_lut.csv'
    with open(street_file, 'w') as f:
        f.write('id,eff_id,begin,eff_begin,end,eff_end,length,width,height,lon_cen,lat_cen\n')
        for i in range(len(street_list)):
            street_id = str(street_list[i].id)
            eff_id = str(street_list[i].eff_id)
            begin = str(street_list[i].begin)
            eff_begin = str(street_list[i].eff_begin)
            end = str(street_list[i].end)
            eff_end = str(street_list[i].eff_end)
            length = str(street_list[i].length)
            width = str(street_list[i].width)
            height = str(street_list[i].height)
            lon_cen = str(street_list[i].lon_cen)
            lat_cen = str(street_list[i].lat_cen)
            f.write(str(street_id) + ',' + str(eff_id) + ',' + str(begin) + ',' + str(eff_begin)\
                    + ',' + str(end) + ',' + str(eff_end) + ',' + str(length) + ',' + str(width)\
                    + ',' + str(height) + ',' + str(lon_cen) + ',' + str(lat_cen) + '\n')

# Read LUTs
def read_luts(node_file, street_file):
    # Read node LUT
    node_list = []
    with open(node_file, 'r') as f:
        header = f.readline()
        for line in f.readlines():
            line_info = line.replace('\n', '').split(',')
            node_id = int(line_info[0])
            eff_id = int(line_info[1])
            lon = float(line_info[2])
            lat = float(line_info[3])
            n_c_street = int(line_info[4])
            node = Node(node_id, lon, lat)
            node.eff_id = eff_id
            if node_id != eff_id:
                node.removed = True
            for i in range(n_c_street):
                node.connected_street.append(int(line_info[5 + i]))
            node_list.append(node)

    # Read street LUT
    street_list = []
    street_list_eff = []
    with open(street_file, 'r') as f:
        header = f.readline()
        for line in f.readlines():
            line_info = line.replace('\n', '').split(',')
            street_id = int(line_info[0])
            eff_id = int(line_info[1])
            begin = int(line_info[2])
            eff_begin = int(line_info[3])
            end = int(line_info[4])
            eff_end = int(line_info[5])
            length = float(line_info[6])
            width = float(line_info[7])
            height = float(line_info[8])
            lon_cen = float(line_info[9])
            lat_cen = float(line_info[10])
            street = Street(street_id, begin, end, length, width, height,
                            lon_cen, lat_cen, {})
            street.eff_id = eff_id
            street.eff_begin = eff_begin
            street.eff_end = eff_end
            if street_id == eff_id:
                street_list.append(street)
                street_list_eff.append(street)
            else:
                street.removed = True
                street_list.append(street)

    return node_list, street_list, street_list_eff

# Generate network files for MUNICH
def generate_network(node_list, street_list, outdir):
    # Write node file
    node_file = outdir + 'node.dat'
    with open(node_file, 'w') as f:
        header = "#id;lon;lat;number_of_streets;1st_street_id;2nd_street_id;...\n"
        f.write(header)
        for i in range(len(node_list)):
            node_id = node_list[i].id
            lon = node_list[i].lon
            lat = node_list[i].lat
            n_c_street = len(node_list[i].connected_street)
            c_streets = ''
            for c_street in node_list[i].connected_street:
                c_streets += ';' + str(c_street)
            f.write(str(node_id) + ';' + str(lon) + ';' + str(lat) +\
                    ';' + str(n_c_street) + c_streets + '\n')

    # Write street file
    street_file = outdir + 'street.dat'
    with open(street_file, 'w') as f:
        header = "#id;begin_inter;end_inter;length;width;height;typo\n"
        f.write(header)
        for i in range(len(street_list)):
            street_id = street_list[i].id
            begin = street_list[i].eff_begin
            end = street_list[i].eff_end
            length = street_list[i].length
            width = street_list[i].width
            height = street_list[i].height
            typo = street_list[i].typo
            f.write(str(street_id) + ';' + str(begin) + ';' +  str(end) + ';' + \
                    str(length) + ';' + str(width) + ';' + str(height) + ';' + \
                    str(typo) + '\n')
