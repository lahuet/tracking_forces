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

def plot_and_save_2d(data_file, data_path, raw_data_file, show=False):
    """ 
    Given a data file with output forces and moments, this method creates plots
    of the results and saves them to image files. The figures can also
    optionally be displayed.
    """
    print '-'*23+'PLOT (2d)'+'-'*24
    
    print 'Loading data...',
    data = load_output_data(data_file)
    t = data['t']

    spike_times = get_spike_times(raw_data_file, t[-1], t[1]-t[0])
    print 'done'

    # The first 3 figures are the raw force/moment signals (both static and 
    # dynamic forces are plotted on the same axes). More plots are also 
    # generated overlaying the spike times.
    print 'Creating and saving plots...',    

    # Moment.
    plt.figure(1)
    plt.plot(t, data['dyn']['M'], t, data['static']['M'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('M')
    plt.title('Moment')
    plt.grid()
    plt.savefig('%s/M.png' %data_path)
    draw_spike_times(spike_times)
    plt.savefig('%s/M_spikes.png' %data_path)

    # Axial force.
    plt.figure(2)
    plt.plot(t, data['dyn']['FA'], t, data['static']['FA'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fa')
    plt.title('Fa')
    plt.grid()
    plt.savefig('%s/Fa.png' %data_path)
    draw_spike_times(spike_times)
    plt.savefig('%s/Fa_spikes.png' %data_path)

    # Transverse force.
    plt.figure(3)
    plt.plot(t, data['dyn']['FT'], t, data['static']['FT'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Ft')
    plt.title('Ft')
    plt.grid()
    plt.savefig('%s/Ft.png' %data_path)
    draw_spike_times(spike_times)
    plt.savefig('%s/Ft_spikes.png' %data_path)

    # Resultant force.
    plt.figure(4)
    plt.plot(t, np.sqrt(data['dyn']['FT']**2+data['dyn']['FA']**2),
             t, np.sqrt(data['static']['FT']**2+data['static']['FA']**2))
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fr')
    plt.title('Fr')
    plt.grid()
    plt.savefig('%s/Fr.png' %data_path)
    draw_spike_times(spike_times)
    plt.savefig('%s/Fr_spikes.png' %data_path)

    print 'done'

    if show:
        plt.show()


def plot_and_save_3d(data_file, data_path, raw_data_file, show=False):
    """
    Given a data file with output forces and moments, this method creates plots
    of the results and saves them to files.
    """
    print '-'*23+'PLOT (3d)'+'-'*24
    
    dyn=1.0
    print 'Loading force data...',   
    data = load_output_data(data_file)
    t = data['t']
    spike_times = get_spike_times(raw_data_file, t[-1], t[1]-t[0])
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
    plt.savefig('%s/Mx.png' %data_path)
    draw_spike_times(spike_times)
    plt.savefig('%s/Mx_spikes.png' %data_path)

    # y-moment
    plt.figure(2)
    plt.plot(t, dyn*data['dyn']['MY'], t, data['static']['MY'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('M')
    plt.title('Moment (y)')
    plt.grid()
    plt.savefig('%s/My.png' %data_path)
    draw_spike_times(spike_times)
    plt.savefig('%s/My_spikes.png' %data_path)

    # z-moment
    plt.figure(3)
    plt.plot(t, dyn*data['dyn']['MZ'], t, data['static']['MZ'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Mz')
    plt.title('Moment (z)')
    plt.grid()
    plt.savefig('%s/Mz.png' %data_path)
    draw_spike_times(spike_times)
    plt.savefig('%s/Mz_spikes.png' %data_path)
   
    # x-force
    plt.figure(4)
    plt.plot(t, dyn*data['dyn']['FX'], t, data['static']['FX'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fx')
    plt.title('Fx')
    plt.grid()
    plt.savefig('%s/Fx.png' %data_path)
    draw_spike_times(spike_times)
    plt.savefig('%s/Fx_spikes.png' %data_path)

    # y-force
    plt.figure(5)
    plt.plot(t, dyn*data['dyn']['FY'], t, data['static']['FY'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fy')
    plt.title('Fy')
    plt.grid()
    plt.savefig('%s/Fy.png' %data_path)
    draw_spike_times(spike_times)
    plt.savefig('%s/Fy_spikes.png' %data_path)

    # z-force
    plt.figure(6)
    plt.plot(t, dyn*data['dyn']['FZ'], t, data['static']['FZ'])
    #plt.plot(t, data['static']['FZ'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fz')
    plt.title('Fz')
    plt.grid()
    plt.savefig('%s/Fz.png' %data_path)
    draw_spike_times(spike_times)
    plt.savefig('%s/Fz_spikes.png' %data_path)

    print 'done'

    if show:
        plt.show()

def plot_and_save(dim, *args, **kwargs):
    if dim == 2:
        plot_and_save_2d(*args, **kwargs)
    elif dim == 3:
        plot_and_save_3d(*args, **kwargs)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='whisker data plotting \
            options')
    parser.add_argument('data_file', help=".mat file with tracked image data")
    parser.add_argument('--dim', help="dimension (2 or 3)", type=int)
    parser.add_argument('--show', help="show figures", action='store_true')
    args = parser.parse_args()

    if os.path.splitext(args.data_file)[1]:
        args.data_file = os.path.splitext(args.data_file)[0]  

    if args.dim == 2:
        plot_and_save_2d('./%s/%s_forces.p' %(args.data_file, args.data_file),
                         './output_data/%s' %args.data_file, 
                         '/%s.mat' %args.data_file, show=args.show)
    
    else:
        plot_and_save_3d('./%s/%s_forces.p' %(args.data_file, args.data_file),
                         './output_data/%s' %args.data_file,
                         '/%s.mat' %args.data_file, show=args.show)

