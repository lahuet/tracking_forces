"""
This file loads tracked whisker data from a .mat file and samples and converts
point data to get configurations in terms of joint angles. It then saves these
data to a .p file in the specified directory.
"""
import pickle
import os
import numpy as np
from numpy.linalg import norm
import scipy.io as sio
import scipy.optimize
from scipy.interpolate import splprep, splev, interp1d, splrep, interp2d
from scipy.spatial.distance import cdist
import argparse
from multiprocessing import Pool
import time

from util import *

_multiprocessing = True
#_multiprocessing = False

PIC_COUNTER = 0

def sample_points(n, d, x, y, z=None, k=1, s=10):
    """
    Samples a 2d or 3d curve to get n points that are each a distance d from
    one another. The first point is always at t=0. (Fixed initial point).
    Points on the curve are found by spline interpolation parameterized by (k,
    s) where k is the order and s is a smoothing factor. Returns the x, y (and
    z) coordinates of the sampled points. If no z-coordinates are provided, the
    conversion is carried out in 2d.
    """
    N = len(x)
    if z is None:
        dim = 2
        tck, u = splprep([x, y], k=k, s=s)
    else:
        dim = 3
        tck, u = splprep([x, y, z], k=k, s=s)

    if scipy.isscalar(d):
        d = d*np.ones(n-1)
    else:
        d = interp1d(np.linspace(0,1,N-1),d)(np.linspace(0,1,n-1))
            
    def F(t_vec):
        if dim == 2:
            xt, yt = splev(np.append([0.0], t_vec), tck)
            p = np.reshape(zip(xt,yt), (n,2))
        else:
            xt, yt, zt = splev(np.append([0.0], t_vec), tck)
            p = np.reshape(zip(xt,yt,zt), (n,3))

        return norm(cdist(p[1:], p[:-1]).diagonal()-d)

    t_guess = interp1d(np.linspace(0,1,N),u)(np.linspace(0,1,n))[1:]
    sol = scipy.optimize.fmin(F, t_guess, disp=0)

    if dim == 2:
        x, y = splev(np.append([0.0], sol), tck)
        z = None
    else:
        x, y, z = splev(np.append([0.0], sol), tck)

    return (x, y, z)   

def convert_single_frame(n, d, scale, x_raw, y_raw, z_raw=None, cp_raw=None,
                         k=1, s=10, rx0_method='zero'):
    """
    Converts a single video frame into an angle configuration and a scaled
    contact point. Returns sampled (x,y,z) coordinates, contact point
    coordinates, and the (dynamic) system configuration.
    """
    # Convert the input points and initialize output arrays.
    x = x_raw.astype(np.float).reshape((1,-1))[0]
    y = y_raw.astype(np.float).reshape((1,-1))[0]
    if z_raw is None:
        dim = 2
        Q = np.empty(n-1+2)
        z = None
    else:
        dim = 3
        Q = np.empty(4+2*(n-1))
        z = -z_raw.astype(np.float).reshape((1,-1))[0]
    CP = np.empty(dim)

    # Sample the curves.
    X, Y, Z = sample_points(n, d, x, y, z, k, s)

    # Subtract an offset so that the base is at (0,0) at time 0. Scale to 
    # units of meters.
    if (not np.isnan(cp_raw[0])) and (cp_raw[0] != 0):
        CP[0] = scale*(cp_raw[0]-X[0])
        CP[1] = scale*(cp_raw[1]-Y[0])  
        if dim == 3: CP[2] = scale*(cp_raw[2]-Z[0])
    else:
        CP.fill(np.nan)

    X = scale*(X-X[0])
    Y = scale*(Y-Y[0])
    if dim == 3: Z = scale*(Z-Z[0])

    # Compute the angle configuration.
    if dim == 2:
        #plot_sampled_points(x,y,X,Y,CP[0],CP[1],scale)
        ang = points_to_angles(zip(X, Y), dim=2)
        #Q = np.hstack((Y[0], X[0], ang))
        Q = np.hstack((ang[0], Y[0], X[0], ang[1:]))

    else: 
        theta_x, theta_y, theta_z =\
                points_to_angles(zip(X, Y, Z), dim=3, rx0_method=rx0_method)
        Q[3:6] = (0.0, 0.0, 0.0)
        Q[2] =  theta_x[0]
        Q[6::2] = theta_z[1:]; Q[0] = theta_z[0]
        Q[7::2] = theta_y[1:]; Q[1] = theta_y[0]

    return (X, Y, Z, CP, Q) 

def convert_single_frame_wrapper(X):
    """Function callable from the multiprocessing methods."""
    n, d, scale, x_raw, y_raw, z_raw, cp_raw, k, s, rx0_method = X
    return convert_single_frame(n, d, scale, x_raw, y_raw, z_raw, 
                                cp_raw, k, s, rx0_method)
 
def convert_frames(file_name, scale, dim, N, start_index=0, stop_index=-1, k=1, s=10,
                   rx0_method='zero', save=True, path_name='./',
                   variable_names=None):
    """ 
    Loads raw whisker data, samples the given points and saves to a file. This
    method computes configurations both in terms of points and angles. 
    Multiprocessing is used to speed up the computations.
    """
    # Note: N is the number of points, including the base. There are N-1 links.
    assert rx0_method in ['zero', 'min_angles'], "invalid rx0 method"
    assert dim in [2,3], "dimension must be 2 or 3"
    print '-'*22+'CONVERT (%dd)' %(dim)+'-'*22

    # Get the data from the .mat file.
    print 'Loading MATLAB data...',
    if variable_names is None:
        # Try to use the defaults.
        if dim == 2: names = {'x': 'xw', 'y': 'yw', 'cp': 'CP'}
        else: names = {'x': 'xw3d', 'y': 'yw3d', 'z': 'zw3d', 'cp': 'CP'}
    else:
        names = variable_names

    try:
        data_raw = sio.loadmat(file_name)
    except IOError:
        raise Exception('File %s not found. Point data could not be loaded.' %file_name)

    x_raw = data_raw[names['x']]
    y_raw = data_raw[names['y']]
    if dim == 3: z_raw = data_raw[names['z']]
    else: z_raw = [[None for i in range(len(x_raw[0]))]]
    try:
        cp_raw = data_raw[names['cp']]
    except KeyError:
        print 'WARNING: No contact point found. Assuming no contacts.'
        cp_raw = np.empty((dim,len(range(start_index,stop_index))))
        cp_raw.fill(np.nan)

    print 'done (Loaded from %s)' %file_name

    # Get the indices over which to convert. If no stopping index, convert 
    # all of the data points.
    if stop_index < 0:
        stop_index = len(x_raw[0])
    indices = range(start_index, stop_index)    

    # Next we will sample the whisker to get equally-spaced points along 
    # the curve. These will be the trajectories that are fed into the whisker
    # system.
    print 'Sampling whisker images...',
    start_time = time.time()

    # Get a link length by averaging several configurations over the entire
    # time interval.
    d = np.mean([arc_length(x_raw[0][i], y_raw[0][i], z_raw[0][i]) 
                for i in np.linspace(start_index, stop_index-1, 
                                     min(len(indices),20)).astype(int)])/(N-1)

    # Set and run the multiprocessing of the conversions.
    args = ((N, d, scale, x_raw[0][i], y_raw[0][i], z_raw[0][i],
             cp_raw[:,i], k, s, rx0_method) for i in indices)
    if _multiprocessing:
        p = Pool()
        res = p.map(convert_single_frame_wrapper, args)
    else:
        res = map(convert_single_frame_wrapper, args)  
    var = zip(*res)
    X = np.asarray(var[0]) 
    Y = np.asarray(var[1]) 
    Z = np.asarray(var[2]) 
    CP = np.asarray(var[3]) 
    Q = np.asarray(var[4]) 

    CONVERTED_DATA = {'x': X, 'y': Y, 'z': Z, 'CP': CP, 'q': Q,
                      'link_length': scale*d}
    print 'done (Converted %d frames in %.1f sec)' %(len(indices),
                                                     time.time()-start_time)

    # Finally, save the data to a file.
    if save:
        c_dir = path_name
        if not os.path.exists(c_dir):
            os.makedirs(c_dir)
        print 'Saving to file...',
        f = open('%s%s.p' %(path_name, os.path.splitext(file_name)[0]), 'w')
        pickle.dump(CONVERTED_DATA, f)
        f.close()
        print 'done (Saved to %s%s.p)' %(path_name, os.path.splitext(file_name)[0])

    return CONVERTED_DATA 

def plot_sampled_points(x,y,xs,ys,cpx,cpy,scale):
    """Get nice-looking pictures of the conversion results. (debug)"""
    import matplotlib.pyplot as plt
    plt.rc('text', usetex=True)
    fig = plt.figure(facecolor='white')
    ax = plt.axes()
    plt.ylim((-2.5, 1.0))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.plot(scale*(x-x[0])*1000,scale*(y-y[0])*1000,xs*1000,ys*1000,'-o',
            cpx*1000,cpy*1000, 'o', )
    plt.xlabel('$x\;\; (\mathrm{mm})$')
    plt.ylabel('$y\;\; (\mathrm{mm})$')
    #l = ax.legend(['$\mathrm{image\;data}$', 
    #               '$\mathrm{sampled\;points}$',
    #               '$\mathrm{contact\;point}$'],numpoints=1)
    #l.draw_frame(False)
    #plt.show()
    global PIC_COUNTER
    plt.savefig('/home/elliot/wpics/im%d.pdf' %PIC_COUNTER)
    PIC_COUNTER+=1


