

# Compile munich-ssh

compile

# Clean munich-ssh

clean

# Edit the path to the downloaded input data. 

Open the file munich-data.cfg
Modify the first line to put the path to the input data.

# Run the simulation

munich-ssh munich.cfg  (single processor)
or
mpirun -n 8 munich-ssh munich.cfg (multi processors)
