#!/usr/bin/env python
import os
import numpy as np
import argparse
import ConfigParser
from pkg_resources import Requirement, resource_filename
import scipy.io as sio

from util import load_file, get_reference_frame
from whisker import make_whisker
from convert import convert_frames
from qfilter import filter_data
from forces import calc_forces
from plot import plot_and_save
from animate import animate_whisker

def save_output_file(path_name, file_name):
    """Reads all of the data from a single run and saves to a .mat file."""

    converted_data = load_file('%s%s.p'%(path_name, file_name))
    filtered_data = load_file('%s%s_filtered.p'%(path_name, file_name))
    force_data = load_file('%s%s_forces.p'%(path_name, file_name))

    # Save all the data in a useful way.
    if force_data['dim'] == 2:
        converted_data['z'] = []
    output_data = {'config': filtered_data['q'],
                   'velocity': filtered_data['v'],
                   'acceleration': filtered_data['a'],
                   'reference_config': filtered_data['ref'],
                   'dimension': force_data['dim'],
                   'time': force_data['t'],
                   'link_length': filtered_data['link_length'],
                   'contact_point': filtered_data['cp'],
                   'reaction_forces': {'static': force_data['static'],
                                       'dynamic': force_data['dyn'], 
                                       'discrete': force_data['discrete']},
                   'sampled_points': {'x': converted_data['x'],
                                      'y': converted_data['y'],
                                      'z': converted_data['z']}}
    sio.savemat('%s%s_output.mat' %(path_name, file_name), output_data)

def main():
    """ 
    Runs the entire process of computing reaction forces from raw input to 
    final output. The file name containing the raw data should be given 
    as the only argument. The various options are read from the configuration
    file. Since the dimension might be changed often, you can override the
    dimension on the command line.
    """
    print 20*'+'+'TRACKING FORCES'+21*'+'

    # Get the input file(s) with the raw data. 
    parser = argparse.ArgumentParser(description='Compute reaction forces from\
            whisker tracking images.')
    parser.add_argument('data_file', help=".mat file with tracked image data")
    parser.add_argument('--dim', type=int, default=-1)

    args = parser.parse_args()

    if os.path.splitext(args.data_file)[1]:
        data_file = os.path.splitext(args.data_file)[0]  
    else:
        data_file = args.data_file

    # Read the options from the configuration file. If no configuration file,
    # use default.
    file_name = None
    for files in os.listdir('.'):
        if files.endswith('.cfg'):
            file_name = files
    if file_name is None: 
        file_name = resource_filename(Requirement.parse("tracking_forces"),
                                     'tracking_forces/tracking_forces.cfg')
    config = ConfigParser.SafeConfigParser()
    config.read([file_name])
    
    # Create a folder to save the outputs.
    output_dir = './%s/' %data_file
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    dim = args.dim
    if dim < 0:
        dim = config.getint('general', 'dimension')
    mm_per_pixel = config.getfloat('general', 'mm_per_pixel')
    scale = mm_per_pixel/1000 
    dt = config.getfloat('general', 'dt')     

    # Load the raw data, sample the images based on the number of links and 
    # convert to configurations in terms of angles.
    if dim == 2:
        try:
            variable_names = {'x': config.get('convert', 'mat_xname_2d'),
                              'y': config.get('convert', 'mat_yname_2d'),
                              'z': config.get('convert', 'mat_zname_2d'),
                              'cp': config.get('convert', 'mat_cpname_2d')}
        except:
            variable_names = {'x': 'xw', 'y': 'yw', 'z': None, 'cp': 'CP'}

    else:
        try:
            variable_names = {'x': config.get('convert', 'mat_xname_3d'),
                              'y': config.get('convert', 'mat_yname_3d'),
                              'z': config.get('convert', 'mat_zname_3d'),
                              'cp': config.get('convert', 'mat_cpname_3d')}
        except:
            variable_names = {'x': 'xw3d', 'y': 'yw3d', 'z': 'zw3d', 'cp': 'CP'}

    conversion_args = (data_file+'.mat', scale, dim,
                       config.getint('convert', 'N'),
                       config.getint('convert','start'),
                       config.getint('convert','stop'),
                       config.getint('convert','convert_k'),
                       config.getfloat('convert','convert_s'),
                       config.get('convert', 'rx0_method'),
                       True, output_dir, variable_names)
    
    # Perform a new conversion.
    if config.getboolean('convert', 'new_conversion'):
        converted_data = convert_frames(*conversion_args)

    # Load previously converted data.
    else:
        try:
            converted_data = load_file('%s%s.p'%(output_dir,data_file))
        except:
            print "WARNING: Converted data not found. Performing new conversion."
            converted_data = convert_frames(*conversion_args)

    # Get the reference shape (index of frame where whisker is undeformed).
    ref = config.getint('filter', 'ref_index')
    if  ref < 0:
        ref = get_reference_frame(converted_data['CP'])   

    # Filter the trajectories.
    filter_args = (data_file+'.p', output_dir, dt, dim, ref,
                   config.get('filter', 'filter_type'),
                   config._sections['filter'])
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
                            filtered_data['cp'], dt,
                            data_file, output_dir,
                            config.getboolean('forces','overwrite'))
    force_data = calc_forces(dim, *forces_args)

    # Plot the results.
    if config.getboolean('plot', 'run_plot'):
        plot_args = ('%s_forces.p' %(data_file), output_dir, './%s.mat' %data_file, 
                     config.getboolean('plot', 'show_plots'))
        plot_and_save(dim, *plot_args)

    # Animate the whisker's motion.
    if config.getboolean('animate', 'run_animation'):
        animate_args = (data_file, output_dir, dt, 
                        config.getboolean('animate', 'show_animation'),
                        config.getboolean('animate', 'save_animation'),
                        config.getboolean('animate', 'debug_conversion'))
        animate_whisker(dim, *animate_args)

    save_output_file(output_dir, data_file)        

if __name__ == "__main__":
    main()
