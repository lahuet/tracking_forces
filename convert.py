"""
This file loads tracked whisker data from a .mat file and samples and converts
point data to get trajectories in terms of joint angles. It then saves these
data to a .p file.
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

import matplotlib.pyplot as plt

_multiprocessing = True
#_multiprocessing = False

def sample_points(n, d, x, y, z=None, k=1, s=10):
    """
    Samples a 2d or 3d curve to get n points that are each a distance d 
    from one another. The first point is always at t=0. (Fixed initial point).
    Points on the curve are found by spline interpolation parameterized by
    (s, k). Returns the x, y (and z) coordinates of the sampled points.
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
       
def convert_single_frame(n, d, x_raw, y_raw, z_raw=None, cp_raw=None,
                         k=1, s=10, rx0_method='zero'):
    """
    Converts a single video frame into an angle configuration and a scaled
    contact point.
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
        CP[0] = SCALE*(cp_raw[0]-X[0])
        CP[1] = SCALE*(cp_raw[1]-Y[0])  
        if dim == 3: CP[2] = SCALE*(cp_raw[2]-Z[0])
    else:
        CP.fill(np.nan)

    X = SCALE*(X-X[0])
    Y = SCALE*(Y-Y[0])
    if dim == 3: Z = SCALE*(Z-Z[0])

    # Compute the angle configuration.
    if dim == 2:
        #plot_sampled_points(x,y,X,Y,CP[0],CP[1])
        ang = points_to_angles(zip(X, Y), dim=2)
        Q = np.hstack((Y[0], X[0], ang))

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
    n, d, x_raw, y_raw, z_raw, cp_raw, k, s, rx0_method = X
    return convert_single_frame(n, d, x_raw, y_raw, z_raw, 
                                cp_raw, k, s, rx0_method)
 
def load_and_convert(file_name, dim, N, start_index=0, stop_index=-1, 
                     k=1, s=10, rx0_method='zero', save=True):
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
    if dim == 2: names = {'x': 'xw', 'y': 'yw', 'cp': 'CP'}
    else: names = {'x': 'xw3d', 'y': 'yw3d', 'z': 'zw3d', 'cp': 'CP'}

    data_raw = sio.loadmat('./raw_data/%s' %file_name)
    x_raw = data_raw[names['x']]
    y_raw = data_raw[names['y']]
    if dim == 3: z_raw = data_raw[names['z']]
    else: z_raw = [[None for i in range(len(x_raw[0]))]]
    cp_raw = data_raw[names['cp']]
    print 'done (Loaded from /raw_data/%s)' %file_name

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
    args = ((N, d, x_raw[0][i], y_raw[0][i], z_raw[0][i],
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
                      'link_length': SCALE*d}
    print 'done (Converted %d frames in %.1f sec)' %(len(indices),
                                                     time.time()-start_time)

    # Finally, save the data to a file.
    if save:
        print 'Saving to file...',
        f = open('./converted_data/%s.p' %os.path.splitext(file_name)[0], 'w')
        pickle.dump(CONVERTED_DATA, f)
        f.close()
        print 'done (Saved to /converted_data/%s.p)' \
                %os.path.splitext(file_name)[0]

    return CONVERTED_DATA 

def plot_sampled_points(x,y,xs,ys,cpx,cpy):
    import matplotlib.pyplot as plt
    plt.rc('text', usetex=True)
    fig = plt.figure(facecolor='white')
    ax = plt.axes()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.plot(SCALE*(x-x[0])*1000,SCALE*(y-y[0])*1000,xs*1000,ys*1000,'-o',
            cpx*1000,cpy*1000, 'o', )
    plt.xlabel('$x\;\; (\mathrm{mm})$')
    plt.ylabel('$y\;\; (\mathrm{mm})$')
    l = ax.legend(['$\mathrm{image\;data}$', 
                   '$\mathrm{sampled\;points}$',
                   '$\mathrm{contact\;point}$'],numpoints=1)
    l.draw_frame(False)
    plt.show()




if __name__ == "__main__":
      
    parser = argparse.ArgumentParser(description='whisker data conversion\
                                                  options')
    parser.add_argument('data_file', help=".mat file with tracked image data")
    parser.add_argument('--start', help="start index", type=int, default=0)
    parser.add_argument('--stop', help="stop index", type=int, default=-1)
    parser.add_argument('--N', help="number of links in the whisker",
                               type=int, default=14)
    parser.add_argument('--k', help="k for interp", type=int, default=1)
    parser.add_argument('--s', help="s for interp", type=float, default=10.0)
    parser.add_argument('--dim', help="dimension (2 or 3)",
                                 type=int, default=2)
    parser.add_argument('--rx0_method', help="method for computing x-axis \
                         rotation at base", default='zero')
    parser.add_argument('--save', help="save (overwrite) converted data",
            action='store_true')
    args = parser.parse_args()

    if not os.path.splitext(args.data_file)[1]:
        args.data_file += '.mat' 

    load_and_convert(args.data_file, dim=args.dim, N=args.N,
            start_index=args.start, stop_index=args.stop,
            k=args.k, s=args.s, rx0_method=args.rx0_method, save=args.save)    

