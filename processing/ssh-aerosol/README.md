

# Compile munich-ssh in two steps

## Step 1. Compile SSH-Aerosol

cd munich/include/ssh-aerosol
compile -s=yes

## Step 2. Compile MUNICH

cd ../../processing/ssh-aerosol
scons -j6 mpi=yes

# Edit the path to the downloaded input data. 

Open the file munich-data.cfg
Modify the first line to put the path to the input data.

# Run the simulation

munich-ssh munich.cfg  (single processor)
or
mpirun -n 8 munich-ssh munich.cfg (multi processors)

