import numpy

def append_binary(arrayToSave, filename, type = 'f'):
    """
    Saves a numpy in a binary file using specified type.

    @type arrayToSave: numpy.array
    @param arrayToSave: The array to save.

    @type filename: string or python file object
    @param filename: The name of the file to save the array into.

    @type type: string
    @param type: Format of data to save the array in file.
    """
    f = open(filename, 'ab')
    numpy.array(arrayToSave, dtype = type).tofile(f)
    f.close()
    
def write_binary(arrayToSave, filename, type = 'f'):
    """
    Saves a numpy in a binary file using specified type.

    @type arrayToSave: numpy.array
    @param arrayToSave: The array to save.

    @type filename: string or python file object
    @param filename: The name of the file to save the array into.

    @type type: string
    @param type: Format of data to save the array in file.
    """
    f = open(filename, 'wb')
    numpy.array(arrayToSave, dtype = type).tofile(f)
    f.close()
    
def write_output(node_list, street_list, node_list_eff, street_list_eff, current_date, output_dir, emis_species_list):

    str_date = current_date.strftime("%Y%m%d%H")
    date = str_date[0:8]
    hour = str_date[8:10]

    # Write 
    # street_ID, begin_node, end_node, -----
    file_emission = output_dir + "emission." + date + "-" + hour + ".txt"
    f = open(file_emission, 'w')

    for i in range(len(street_list_eff)):
        street = street_list_eff[i]
        street_id = street.id
        begin_node = street.eff_begin
        end_node = street.eff_end

        nemis = len(street.eff_emission)

        f.write(str(street_id) + "\t" + str(begin_node) + "\t" + 
                str(end_node) + "\t")        
        for s in range(nemis):
            f.write(str(float(street.eff_emission[s])) + "\t")
        f.write("\n")
    f.close()

    # Write
    # node_ID, longitude, latitude, number of connected segments, segment indices
    file_node = output_dir + "intersection.dat"
    f = open(file_node, 'w')
    header = "#id;lon;lat;number_of_streets;1st_street_id;2nd_street_id;...\n"
    f.write(header)
    for i in range(len(node_list_eff)):
        node = node_list_eff[i]
        if node.removed == False:
            lon = node.lon
            lat = node.lat
            node_id = node.id
            ns = len(node.connected_street)
            st_ = ""
            for st in node.connected_street:
                st_ += str(st) + ";"
            f.write(str(node_id) + ";" + str(lon) + ";" + str(lat) + 
                    ";" + str(ns) + ";" + st_ +  "\n")
    f.close()

    # Write
    file_node = output_dir + "street.dat"
    f = open(file_node, 'w')
    header = "#id;begin_inter;end_inter;length;width;height;typo\n"
    f.write(header)
    for i in range(len(street_list_eff)):
        street = street_list_eff[i]
        street_id = street.id
        begin_node = street.eff_begin
        end_node = street.eff_end
        length = street.length
        width = street.width
        height = street.height
        typo = street.typo
        f.write(str(street_id) + ";" + str(begin_node) + ";" + str(end_node) + 
                ";" + str(length) + ";" + str(width) + ";" + str(height) +
                ";" + str(typo) + "\n")
    f.close()
    
    # Write
    import os.path
    
    # Emission data
    emission_array = np.zeros([len(street_list_eff), len(emis_species_list)], 'float')

    wind_dir = np.zeros((len(street_list_eff)), 'float')
    wind_speed = np.zeros((len(street_list_eff)), 'float')
    pblh = np.zeros((len(street_list_eff)), 'float')
    ust = np.zeros((len(street_list_eff)), 'float')
    lmo = np.zeros((len(street_list_eff)), 'float')
    psfc = np.zeros((len(street_list_eff)), 'float')
    t2 = np.zeros((len(street_list_eff)), 'float')
    attenuation = np.zeros((len(street_list_eff)), 'float')
    sh = np.zeros((len(street_list_eff)), 'float')
    lwc = np.zeros((len(street_list_eff)), 'float')
    rain = np.zeros((len(street_list_eff)), 'float')
    background={}
    street0=street_list_eff[0]
    for spec in list(street0.background.keys()):
        background[spec]=np.zeros((len(street_list_eff)), 'float')

    for i in range(len(street_list_eff)):
        street = street_list_eff[i]

        emission_array[i] = street.eff_emission

        wind_dir[i] = street.wind_dir
        wind_speed[i] = street.wind_speed
        pblh[i] = street.pblh
        ust[i] = street.ust
        lmo[i] = street.lmo
        psfc[i] = street.psfc
        t2[i] = street.t2
        sh[i] = street.sh
        attenuation[i] = street.attenuation
        lwc[i] = street.lwc
        rain[i] = street.rain

        for spec in list(background.keys()):
            background[spec][i]=street.background[spec]

    wind_dir_inter = np.zeros((len(node_list_eff)), 'float')
    wind_speed_inter = np.zeros((len(node_list_eff)), 'float')
    pblh_inter = np.zeros((len(node_list_eff)), 'float')
    ust_inter = np.zeros((len(node_list_eff)), 'float')
    lmo_inter = np.zeros((len(node_list_eff)), 'float')
    for i in range(len(node_list_eff)):
        node = node_list_eff[i]
        wind_dir_inter[i] = node.wind_dir
        wind_speed_inter[i] = node.wind_speed
        pblh_inter[i] = node.pblh
        ust_inter[i] = node.ust
        lmo_inter[i] = node.lmo

    return wind_dir, wind_speed, pblh, ust, lmo, psfc, t2, sh, attenuation, background, wind_dir_inter, wind_speed_inter, pblh_inter, ust_inter, lmo_inter, emission_array, lwc, rain
