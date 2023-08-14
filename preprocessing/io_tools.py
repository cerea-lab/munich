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
    
