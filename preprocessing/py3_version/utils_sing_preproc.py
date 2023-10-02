import sys
import pyproj
import numpy as np
import datetime as dt
from dateutil import tz

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
                 length, width, height, lon_cen, lat_cen):
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
        # emission
        self.emission = {}
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

##### Misc.

def get_polair_id(lon, lat, x_min, y_min, dx, dy, Nx, Ny):
    Xid = max(int((lon - x_min + dx / 2.) / dx), 0)
    Xid = min(Xid, Nx - 1)
    Yid = max(int((lat - y_min + dy / 2.) / dy), 0)
    Yid = min(Yid, Ny - 1)

    return Xid, Yid # YK, Majorlaine

def append_binary(arrayToSave, filename, type = 'f'):
    with open(filename, 'ab') as f:
        np.array(arrayToSave, dtype=type).tofile(f)

##### Street network

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
                            mid_point[0][0], mid_point[0][1])
            street_list.append(street)

    return node_list, street_list

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
                            lon_cen, lat_cen)
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
        for i in range(len(street_list)):
            street_id = street_list[i].id
            begin = street_list[i].eff_begin
            end = street_list[i].eff_end
            length = street_list[i].length
            width = street_list[i].width
            height = street_list[i].height
            f.write(str(street_id) + ';' + str(begin) + ';' +  str(end) + ';' +\
                    str(length) + ';' + str(width) + ';' + str(height) + '\n')

##### Emission

def utc_to_local(utc, zone='Europe/Paris'):
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz(zone)

    # Tell the datetime object that it is in UTC since
    # datetime object are 'naive' by default.
    utc = utc.replace(tzinfo=from_zone)

    # Convert time zone
    return utc.astimezone(to_zone)

# Requires holidays python library (https://pypi.org/project/holidays/)
def is_holiday(date, country_code):
    # Try import holidays
    try:
        import holidays
    except ImportError:
        print('WARNING: holidays library could not be imported.')
        return False
    # Return True if it is a holiday
    country_found = False
    for country in holidays.list_supported_countries():
        if country_code == country:
            country_found = True
            country_holidays = holidays.CountryHoliday(country_code)
            break
    if country_found == False:
        print('WARNING: country_code {} not supported.'.format(country_code))
        return False

    return date.date() in country_holidays

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

##### Meteo

# Compute wind direction from zonal (u) and meridional (v) data
def compute_wdir(u, v):
    wdir = np.arctan2(v, u)
    hpi = np.pi / 2. # Half of Pi
    wdir = np.where(wdir<=hpi, hpi-wdir, (hpi-wdir)+2.*np.pi)
    return wdir

# Read meteo binary files from Polair3D input
def read_meteo_bin(indir, Nt, Nx, Ny, Nz):
    # Put data in a dict
    meteo_data = {}
    infile = indir + 'Attenuation.bin'
    meteo_data['attenuation'] = np.memmap(infile, dtype='float32', mode='r',
                                          shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]
    infile = indir + 'LiquidWaterContent.bin'
    meteo_data['liquidwaterc'] = np.memmap(infile, dtype='float32', mode='r',
                                           shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]
    infile = indir + 'LMO.bin'
    meteo_data['lmo'] = np.memmap(infile, dtype='float32', mode='r',
                                  shape=(Nt, Ny, Nx))
    infile = indir + 'PARdb.bin'
    meteo_data['pardb'] = np.memmap(infile, dtype='float32', mode='r',
                                    shape=(Nt, Ny, Nx))
    infile = indir + 'PARdiff.bin'
    meteo_data['pardiff'] = np.memmap(infile, dtype='float32', mode='r',
                                      shape=(Nt, Ny, Nx))
    infile = indir + 'BoundaryHeight.bin'
    meteo_data['pblh'] = np.memmap(infile, dtype='float32', mode='r',
                                   shape=(Nt, Ny, Nx))
    infile = indir + 'Rain.bin'
    meteo_data['rain'] = np.memmap(infile, dtype='float32', mode='r',
                                   shape=(Nt, Ny, Nx))
    infile = indir + 'SoilWater.bin'
    meteo_data['soilwater'] = np.memmap(infile, dtype='float32', mode='r',
                                        shape=(Nt, Ny, Nx))
    infile = indir + 'SolarRadiation.bin'
    meteo_data['solarrad'] = np.memmap(infile, dtype='float32', mode='r',
                                       shape=(Nt, Ny, Nx))
    infile = indir + 'SpecificHumidity.bin'
    meteo_data['spechumid'] = np.memmap(infile, dtype='float32', mode='r',
                                        shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]
    infile = indir + 'SurfacePressure.bin'
    meteo_data['surfpressure'] = np.memmap(infile, dtype='float32', mode='r',
                                           shape=(Nt, Ny, Nx))
    infile = indir + 'SurfaceRichardson.bin'
    meteo_data['surfrichard'] = np.memmap(infile, dtype='float32', mode='r',
                                          shape=(Nt, Ny, Nx))
    infile = indir + 'SurfaceTemperature.bin'
    meteo_data['surftemp'] = np.memmap(infile, dtype='float32', mode='r',
                                       shape=(Nt, Ny, Nx))
    infile = indir + 'FrictionModule.bin'
    meteo_data['ust'] = np.memmap(infile, dtype='float32', mode='r',
                                  shape=(Nt, Ny, Nx))
    infile = indir + 'WindModule.bin'
    meteo_data['wspd'] = np.memmap(infile, dtype='float32', mode='r',
                                   shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]
    infile = indir + 'ZonalWind.bin'
    meteo_data['zonalwind'] = np.memmap(infile, dtype='float32', mode='r',
                                        shape=(Nt, Nz, Ny, Nx+1))[:, 0, :, :]
    infile = indir + 'MeridionalWind.bin'
    meteo_data['meridiowind'] = np.memmap(infile, dtype='float32', mode='r',
                                          shape=(Nt, Nz, Ny+1, Nx))[:, 0, :, :]

    return meteo_data

# Set meteo for streets and nodes
def set_meteo_bin(node_list, street_list, meteo_data, current_date, date_min, delta_t,
                  Nt, x_min, y_min, delta_x, delta_y, Nx, Ny):
    # Get index of current date
    c_id = int((current_date - date_min).total_seconds() / delta_t)
    if c_id < 0 or c_id >= Nt:
        sys.exit('ERROR: meteo data not available for this date.')

    # Set meteo for the streets
    for i in range(len(street_list)):
        # Get cell indices (X, Y) in the Polair3D grid
        # Xid = (street_list[i].lon_cen - x_min) / delta_x
        # Yid = (street_list[i].lat_cen - y_min) / delta_y
        # Xid = int(Xid)-1 if Xid%1 < 0.5 else int(Xid)
        # Yid = int(Yid)-1 if Yid%1 < 0.5 else int(Yid)
        Xid, Yid = get_polair_id(street_list[i].lon_cen, street_list[i].lat_cen,
                                 x_min, y_min, delta_x, delta_y, Nx, Ny)
        street_list[i].attenuation = meteo_data['attenuation'][c_id, Yid, Xid]
        street_list[i].liquidwaterc = meteo_data['liquidwaterc'][c_id, Yid, Xid]
        street_list[i].lmo = meteo_data['lmo'][c_id, Yid, Xid]
        street_list[i].pardb = meteo_data['pardb'][c_id, Yid, Xid]
        street_list[i].pardiff = meteo_data['pardiff'][c_id, Yid, Xid]
        street_list[i].pblh = meteo_data['pblh'][c_id, Yid, Xid]
        street_list[i].rain = meteo_data['rain'][c_id, Yid, Xid]
        street_list[i].soilwater = meteo_data['soilwater'][c_id, Yid, Xid]
        street_list[i].solarrad = meteo_data['solarrad'][c_id, Yid, Xid]
        street_list[i].spechumid = meteo_data['spechumid'][c_id, Yid, Xid]
        street_list[i].surfpressure = meteo_data['surfpressure'][c_id, Yid, Xid]
        street_list[i].surfrichard = meteo_data['surfrichard'][c_id, Yid, Xid]
        street_list[i].surftemp = meteo_data['surftemp'][c_id, Yid, Xid]
        street_list[i].ust = meteo_data['ust'][c_id, Yid, Xid]
        u = meteo_data['zonalwind'][c_id, Yid, Xid]
        v = meteo_data['meridiowind'][c_id, Yid, Xid]
        street_list[i].wdir = compute_wdir(u, v)
        street_list[i].wspd = meteo_data['wspd'][c_id, Yid, Xid]

    # Set meteo for the nodes
    for i in range(len(node_list)):
        # Get cell indices (X, Y) in the Polair3D grid
        # Xid = (node_list[i].lon - x_min) / delta_x
        # Yid = (node_list[i].lat - y_min) / delta_y
        # Xid = int(Xid)-1 if Xid%1 < 0.5 else int(Xid)
        # Yid = int(Yid)-1 if Yid%1 < 0.5 else int(Yid)
        Xid, Yid = get_polair_id(node_list[i].lon, node_list[i].lat,
                                 x_min, y_min, delta_x, delta_y, Nx, Ny)
        node_list[i].lmo = meteo_data['lmo'][c_id, Yid, Xid]
        node_list[i].pblh = meteo_data['pblh'][c_id, Yid, Xid]
        node_list[i].ust = meteo_data['ust'][c_id, Yid, Xid]
        u = meteo_data['zonalwind'][c_id, Yid, Xid]
        v = meteo_data['meridiowind'][c_id, Yid, Xid]
        node_list[i].wdir = compute_wdir(u, v)
        node_list[i].wspd = meteo_data['wspd'][c_id, Yid, Xid]

# Append meteo to binary files
def append_meteo_data(node_list, street_list, street_files, node_files):
    # Street
    street_data = [np.zeros((len(street_list)), 'float') for i in range(len(street_files))]
    for i in range(len(street_list)):
        street_data[0][i] = street_list[i].attenuation
        street_data[1][i] = street_list[i].liquidwaterc
        street_data[2][i] = street_list[i].lmo
        street_data[3][i] = street_list[i].pardb
        street_data[4][i] = street_list[i].pardiff
        street_data[5][i] = street_list[i].pblh
        street_data[6][i] = street_list[i].rain
        street_data[7][i] = street_list[i].soilwater
        street_data[8][i] = street_list[i].solarrad
        street_data[9][i] = street_list[i].spechumid
        street_data[10][i] = street_list[i].surfpressure
        street_data[11][i] = street_list[i].surfrichard
        street_data[12][i] = street_list[i].surftemp
        street_data[13][i] = street_list[i].ust
        street_data[14][i] = street_list[i].wdir
        street_data[15][i] = street_list[i].wspd
    for i in range(len(street_files)):
        append_binary(street_data[i], street_files[i])

    # Node
    node_data = [np.zeros((len(node_list)), 'float') for i in range(len(node_files))]
    for i in range(len(node_list)):
        node_data[0][i] = node_list[i].lmo
        node_data[1][i] = node_list[i].pblh
        node_data[2][i] = node_list[i].ust
        node_data[3][i] = node_list[i].wdir
        node_data[4][i] = node_list[i].wspd
    for i in range(len(node_files)):
        append_binary(node_data[i], node_files[i])

##### Background

# Read background binary files from Polair3D output
def read_bkgd_bin(bkgd_species, indir, Nt, Nx, Ny, Nz):
    # Put data in a dict
    bkgd_data = {}
    for species in bkgd_species:
        infile = indir + species + '.bin'
        bkgd_data[species] = np.memmap(infile, dtype='float32', mode='r',
                                       shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]

    return bkgd_data

# Set background for streets
def set_bkgd_bin(street_list, bkgd_data, current_date, date_min, delta_t,
                  Nt, x_min, y_min, delta_x, delta_y, Nx, Ny):
    # Get index of current date
    c_id = int((current_date - date_min).total_seconds() / delta_t)
    if c_id < 0 or c_id >= Nt:
        sys.exit('ERROR: background data not available for this date.')

    for i in range(len(street_list)):
        # Get cell indices (X, Y) in the Polair3D grid
        # Xid = (street_list[i].lon_cen - x_min) / delta_x
        # Yid = (street_list[i].lat_cen - y_min) / delta_y
        # Xid = int(Xid)-1 if Xid%1 < 0.5 else int(Xid)
        # Yid = int(Yid)-1 if Yid%1 < 0.5 else int(Yid)
        Xid, Yid = get_polair_id(street_list[i].lon_cen, street_list[i].lat_cen,
                                 x_min, y_min, delta_x, delta_y, Nx, Ny)
        for key, value in bkgd_data.items():
            street_list[i].bkgd[key] = value[c_id, Yid, Xid]

# Append background to binary files
def append_bkgd_data(street_list, bkgd_files):
    for i, var in enumerate(street_list[0].bkgd.keys()):
        data = np.zeros((len(street_list)), 'float')
        for j in range(len(street_list)):
            data[j] = street_list[j].bkgd[var]
        append_binary(data, bkgd_files[var])

##### Gridded emission
