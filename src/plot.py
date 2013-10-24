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
    print 'path', data_path
    
    print 'Loading data...',
    data = load_file('%s/%s' %(data_path, data_file))
    t = data['t']

    #spike_times = get_spike_times(raw_data_file, t[-1], t[1]-t[0])
    print 'done'
    if not os.path.exists(data_path+'/pics'):
        os.makedirs(data_path+'/pics')

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
    plt.savefig('%s/pics/M.png' %data_path)
    #draw_spike_times(spike_times)
    #plt.savefig('%s/M_spikes.png' %data_path)

    # Axial force.
    plt.figure(2)
    plt.plot(t, data['dyn']['FY'], t, data['static']['FY'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fa')
    plt.title('Fa')
    plt.grid()
    plt.savefig('%s/pics/Fa.png' %data_path)
    #draw_spike_times(spike_times)
    #plt.savefig('%s/Fa_spikes.png' %data_path)

    # Transverse force.
    plt.figure(3)
    plt.plot(t, data['dyn']['FZ'], t, data['static']['FZ'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Ft')
    plt.title('Ft')
    plt.grid()
    plt.savefig('%s/pics/Ft.png' %data_path)
    #draw_spike_times(spike_times)
    #plt.savefig('%s/Ft_spikes.png' %data_path)

    # Resultant force.
    plt.figure(4)
    plt.plot(t, np.sqrt(data['dyn']['FT']**2+data['dyn']['FA']**2),
             t, np.sqrt(data['static']['FT']**2+data['static']['FA']**2))
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fr')
    plt.title('Fr')
    plt.grid()
    plt.savefig('%s/pics/Fr.png' %data_path)
    #draw_spike_times(spike_times)
    #plt.savefig('%s/Fr_spikes.png' %data_path)

    print 'done'

    if show:
        plt.show()


def plot_and_save_3d(file_name, path_name, raw_data_file, show=False):
    """
    Given a data file with output forces and moments, this method creates plots
    of the results and saves them to files.
    """
    print '-'*23+'PLOT (3d)'+'-'*24
    
    dyn=1.0
    print 'Loading force data...',   
    data = load_file(path_name+file_name)
    t = data['t']
    #spike_times = get_spike_times(raw_data_file, t[-1], t[1]-t[0])
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
    #draw_spike_times(spike_times)
    #plt.savefig('%s/Mx_spikes.png' %data_path)

    # y-moment
    plt.figure(2)
    plt.plot(t, dyn*data['dyn']['MY'], t, data['static']['MY'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('M')
    plt.title('Moment (y)')
    plt.grid()
    plt.savefig('%sMy.png' %pic_path)
    #draw_spike_times(spike_times)
    #plt.savefig('%s/My_spikes.png' %data_path)

    # z-moment
    plt.figure(3)
    plt.plot(t, dyn*data['dyn']['MZ'], t, data['static']['MZ'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Mz')
    plt.title('Moment (z)')
    plt.grid()
    plt.savefig('%sMz.png' %pic_path)
    #draw_spike_times(spike_times)
    #plt.savefig('%s/Mz_spikes.png' %data_path)
   
    # x-force
    plt.figure(4)
    plt.plot(t, dyn*data['dyn']['FX'], t, data['static']['FX'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fx')
    plt.title('Fx')
    plt.grid()
    plt.savefig('%sFx.png' %pic_path)
    #draw_spike_times(spike_times)
    #plt.savefig('%s/Fx_spikes.png' %data_path)

    # y-force
    plt.figure(5)
    plt.plot(t, dyn*data['dyn']['FY'], t, data['static']['FY'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fy')
    plt.title('Fy')
    plt.grid()
    plt.savefig('%sFy.png' %pic_path)
    #draw_spike_times(spike_times)
    #plt.savefig('%s/Fy_spikes.png' %data_path)

    # z-force
    plt.figure(6)
    plt.plot(t, dyn*data['dyn']['FZ'], t, data['static']['FZ'])
    #plt.plot(t, data['static']['FZ'])
    plt.legend(["Dynamic", "Static"])
    plt.xlabel('t')
    plt.ylabel('Fz')
    plt.title('Fz')
    plt.grid()
    plt.savefig('%sFz.png' %pic_path)
    #draw_spike_times(spike_times)
    #plt.savefig('%s/Fz_spikes.png' %data_path)

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

    data = load_file('./output_data/%s/%s' %(args.data_file,
        args.data_file+'_forces.p'))
    t = data['t']
    nice_looking_plots(t, data['dyn'], data['static'])
    '''
    if args.dim == 2:
        plot_and_save_2d('/%s_forces.p' %(args.data_file),
                         './output_data/%s' %args.data_file, 
                         '/%s.mat' %args.data_file, show=args.show)
    
    else:
        plot_and_save_3d('/%s_forces.p' %(args.data_file),
                         './output_data/%s' %args.data_file,
                         '/%s.mat' %args.data_file, show=args.show)
    '''
