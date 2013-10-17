#!/usr/bin/env python
import os
import numpy as np
import argparse
import ConfigParser

from util import load_converted_data, get_reference_frame
from whisker import make_whisker
from convert import convert_frames
from qfilter import filter_data
from forces import calc_forces
from plot import plot_and_save
from animate import animate_whisker

def main(data_file):
    """ 
    Runs the entire process of computing reaction forces from raw input to 
    final output. The file name containing the raw data should be given 
    as the only argument. The various options are read from the configuration
    file.
    """
    print 20*'+'+'TRACKING FORCES'+21*'+'

    if os.path.splitext(data_file)[1]:
        data_file = os.path.splitext(data_file)[0]  

    # Read the options from the configuration file.
    config = ConfigParser.SafeConfigParser()
    config.read(['tracking_forces.cfg'])

    dim = config.getint('general', 'dimension')
    converted_file_path = config.get('general', 'converted_dir')
    output_file_path = config.get('general', 'output_dir')

    # Load the raw data, sample the images based on the number of links and 
    # convert to configurations in terms of angles.
    variable_names = {'x': config.get('convert', 'mat_xname'),
                      'y': config.get('convert', 'mat_yname'),
                      'z': config.get('convert', 'mat_zname'),
                      'cp': config.get('convert', 'mat_cpname')}
    conversion_args = (data_file+'.mat', dim,
                       config.getint('convert', 'N'),
                       config.getint('convert','start'),
                       config.getint('convert','stop'),
                       config.getint('convert','convert_k'),
                       config.getfloat('convert','convert_s'),
                       config.get('convert', 'rx0_method'),
                       True, converted_file_path, variable_names)

    # Perform a new conversion.
    if config.getboolean('convert', 'new_conversion'):
        converted_data = convert_frames(*conversion_args)

    # Load previously converted data.
    else:
        try:
            converted_data = load_converted_data(data_file+'.p')
        except:
            print "WARNING: Converted data not found.\
                    Performing new conversion."
            converted_data = convert_frames(*conversion_args)

    # Get the reference shape (index of frame where whisker is undeformed).
    ref = config.getint('filter', 'ref_index')
    if  ref < 0:
        ref = get_reference_frame(converted_data['CP'])   

    # Filter the trajectories.
    filter_args = (data_file, dim, ref,
                   config.get('filter', 'filter_type'),
                   config._sections['filter'], converted_file_path)
    filtered_data = filter_data(*filter_args)

    # Build the whisker system in trep.
    whisker = make_whisker(dim, filtered_data['ref'], 
                           converted_data['link_length'],
                           config.getfloat('whisker', 'rbase'),
                           config.getfloat('whisker', 'taper'),
                           config.getfloat('whisker', 'damping_ratio'),
                           config.getfloat('whisker', 'rho'),
                           config.getfloat('whisker', 'E'))

    # Compute the reaction forces and moments.
    forces_args = (whisker, filtered_data['q'], 
                            filtered_data['v'],
                            filtered_data['a'],
                            filtered_data['cp'],
                            data_file, output_file_path,
                            config.getboolean('forces','overwrite'))
    force_data = calc_forces(dim, *forces_args)

    # Plot the results.
    if config.getboolean('plot', 'run_plot'):
        plot_args = ('./%s/%s_forces.p' %(data_file, data_file),
                     './output_data/%s' %data_file, 
                     '/%s.mat' %data_file, 
                     config.getboolean('plot', 'show_plots'))
        plot_and_save(dim, *plot_args)

    # Animate the whisker's motion.
    if config.getboolean('animate', 'run_animation'):
        animate_args = (data_file, 
                        config.getboolean('animate', 'show_animation'),
                        config.getboolean('animate', 'save_animation'),
                        config.getboolean('animate', 'debug_conversion'))
        animate_whisker(dim, *animate_args)

if __name__ == "__main__":
    # Get the input file(s) with the raw data. 
    parser = argparse.ArgumentParser(description='Compute reaction forces from\
            whisker tracking images')
    parser.add_argument('data_files', help=".mat file with tracked image data",
            nargs='+')

    args = parser.parse_args()

    for file_name in args.data_files:
        main(file_name)
