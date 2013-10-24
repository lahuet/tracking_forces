This code calculates dynamic reaction forces from tracked whisker data. 

To install the module, cd into the top directory and run:
$ sudo python setup.py install

The four basic steps of the process are:
    1. Converting point data from images to angle configurations.
    2. Filtering the trajectories.
    3. Calculating reaction forces (using trep).
    4. Post-processing (plotting and animation).

The entire process can be run from the command line using:
$ tracking_forces point_data_file.mat
where point_data_file contains data for the shape of the whisker. (If the setup
script is run, this can be run from any directory.)

There are many options to choose when running this code. All are read from a
.cfg file. If no .cfg file exists in the current directory, default values are
used from a file installed with the rest of the module. For an example
configuration file, see examples.

The output is a set of files saved in a folder with the same name as the input
file. All converted and filtered data are stored in file.p and file_filtered.p.
The output forces and moments are stored in file_forces.mat. Optional plots and
movies are stored in the same location.

This file structure is used because Step 1 can take a long time for large input
files. This allows one to recalculate forces after doing the conversion only
once. (See the .cfg file.)
