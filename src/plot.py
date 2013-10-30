"""
This file loads saved reaction force data and plots them. It also saves image
files of the results to the /output_data/ directory.
"""
import sys
import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})    # Adjust margins for labels.
import argparse

from util import *


def draw_spike_times(spike_times):  
    """Draws vertical lines denoting spike times on the current figure."""
    for line in spike_times:
        plt.axvline(x=line, color='y')

def plot_and_save_2d(file_name, path_name, raw_data_file, show=False):
    """ 
    Given a data file with output forces and moments, this method creates plots
    of the results and saves them to image files. The figures can also
    optionally be displayed.
    """
    print '-'*23+'PLOT (2d)'+'-'*24
    
    print 'Loading data...',
    data = load_file(path_name+file_name)
    t = data['t']
    
    pic_path = path_name+'pics/'
    if not os.path.exists(pic_path):
        os.makedirs(pic_path)
    print 'done'
    print 'Creating and saving plots...',    

    # Moment.
    plt.figure(1)
    plt.plot(t, data['dyn']['M'], t, data['static']['M'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('M')
    plt.title('Moment')
    plt.grid()
    plt.savefig('%sM.png' %pic_path)

    # Axial force.
    plt.figure(2)
    plt.plot(t, data['dyn']['FY'], t, data['static']['FY'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fa')
    plt.title('Fa')
    plt.grid()
    plt.savefig('%sFa.png' %pic_path)

    # Transverse force.
    plt.figure(3)
    plt.plot(t, data['dyn']['FZ'], t, data['static']['FZ'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Ft')
    plt.title('Ft')
    plt.grid()
    plt.savefig('%sFt.png' %pic_path)

    # Resultant force.
    plt.figure(4)
    plt.plot(t, np.sqrt(data['dyn']['FY']**2+data['dyn']['FZ']**2),
             t, np.sqrt(data['static']['FY']**2+data['static']['FZ']**2))
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fr')
    plt.title('Fr')
    plt.grid()
    plt.savefig('%sFr.png' %pic_path)
    print 'done'

    if show:
        plt.show()

def plot_and_save_3d(file_name, path_name, raw_data_file, show=False):
    """
    Given a data file with output forces and moments, this method creates plots
    of the results and saves them to files.
    """
    print '-'*23+'PLOT (3d)'+'-'*24
    
    print 'Loading force data...',   
    data = load_file(path_name+file_name)
    t = data['t']
    dyn = 1.0
    
    pic_path = path_name+'pics/'
    if not os.path.exists(pic_path):
        os.makedirs(pic_path)
    print 'done'
    print 'Creating and saving plots...', 

    # x-moment
    plt.figure(1)
    plt.plot(t, dyn*data['dyn']['MX'], t, data['static']['MX'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Mx')
    plt.title('Moment (x)')
    plt.grid()
    plt.savefig('%sMx.png' %pic_path)

    # y-moment
    plt.figure(2)
    plt.plot(t, dyn*data['dyn']['MY'], t, data['static']['MY'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('M')
    plt.title('Moment (y)')
    plt.grid()
    plt.savefig('%sMy.png' %pic_path)

    # z-moment
    plt.figure(3)
    plt.plot(t, dyn*data['dyn']['MZ'], t, data['static']['MZ'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Mz')
    plt.title('Moment (z)')
    plt.grid()
    plt.savefig('%sMz.png' %pic_path)
   
    # x-force
    plt.figure(4)
    plt.plot(t, dyn*data['dyn']['FX'], t, data['static']['FX'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fx')
    plt.title('Fx')
    plt.grid()
    plt.savefig('%sFx.png' %pic_path)

    # y-force
    plt.figure(5)
    plt.plot(t, dyn*data['dyn']['FY'], t, data['static']['FY'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fy')
    plt.title('Fy')
    plt.grid()
    plt.savefig('%sFy.png' %pic_path)

    # z-force
    plt.figure(6)
    plt.plot(t, dyn*data['dyn']['FZ'], t, data['static']['FZ'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fz')
    plt.title('Fz')
    plt.grid()
    plt.savefig('%sFz.png' %pic_path)
    print 'done'

    #nice_looking_plots(t, data['dyn'], data['static'])

    if show:
        plt.show()

def nice_looking_plots(t, dyn_data, static_data):
    plt.rc('text', usetex=True)
    fig = plt.figure(facecolor='white')
    ax = plt.axes()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.plot(t, static_data['FY']/(1.0e-6), t, dyn_data['FY']/(1.0e-6))
    plt.xlabel('$t\;\; (\mathrm{sec})$')
    plt.ylabel('$F_z\;\; (\mathrm{\mu N})$')
    l = ax.legend(['$\mathrm{static}$', 
                   '$\mathrm{dynamic}$'])
    l.draw_frame(False)
    plt.show()

def plot_and_save(dim, *args, **kwargs):
    """Wrapper for plotting in the correct dimension."""
    if dim == 2:
        plot_and_save_2d(*args, **kwargs)
    elif dim == 3:
        plot_and_save_3d(*args, **kwargs)


