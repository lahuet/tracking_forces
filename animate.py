import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.animation as animation
import numpy as np
import argparse
import sys
from StringIO import StringIO

from util import *
from whisker import make_whisker

VIDEO_TIME_SCALING = 0.1/20     # (real time)/(video time)

def load_points(file_name, dim):
    """Returns x, y, (and z) coordinates of whisker."""
    data = load_converted_data(file_name)
    X = data['x']
    Y = data['y']
    if dim == 2:
        return X, Y
    else:
        Z = data['z']
        return X, Y, Z

def load_contact_points(file_name, dim):
    """Returns the x, y, (and z) locations of the contact points."""
    data = load_converted_data(file_name)
    CP = data['CP']
    x = CP[:,0]
    y = CP[:,1]
    if dim == 2:
        return x, y
    else:
        z = CP[:,3]
        return x, y, z

def animate_whisker_2d(file_name, show=True, save_movie=False, debug=False):
    """Draws 2d animation of whisker motion."""
    print '-'*22+'ANIMATE (2d)'+'-'*22

    X, Y = load_points(file_name+'.p', 2)
    if debug: Xf, Yf = get_filtered_points(file_name, 2)
    cp_x, cp_y = load_contact_points(file_name+'.p', 2)

    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(111)

    pts = ax.plot([], [], '-o')
    if debug: fpts = ax.plot([], [], 'r-o')
    cp_pt = ax.plot([], [], '*', ms=14)
    time_text = ax.text(0.01, 0.01, '', transform=ax.transAxes)

    ax.set_xlim((np.min(X), np.max(X)))
    ax.set_ylim((np.min(Y), np.max(Y)))

    def init():
        if debug: fpts[0].set_data([], [])
        pts[0].set_data([], [])
        cp_pt[0].set_data([], [])
        time_text.set_text('')

    def animate(i):
        if debug: fpts[0].set_data(Xf[i], Yf[i])
        pts[0].set_data(X[i],Y[i])
        if not np.isnan(cp_x[i]):    
            cp_pt[0].set_data(cp_x[i], cp_y[i])
        else:
            cp_pt[0].set_data([], [])

        time_text.set_text('t = %.3f s' %(DT*i))
        fig.canvas.draw()

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.grid()
    ax.set_title(file_name)
    
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=len(X), interval=10, blit=False)

    if save_movie:
        anim.save('./output_data/%s/%s.mp4' %(file_name, file_name), fps=15)
                 # writer=animation.FFMpegFileWriter(),
                 # extra_args=['-vcodec', 'libx264'])

    if show:
        sys.exit(plt.show())


def animate_whisker_3d(file_name, show=True, save_movie=False, debug=False):
    """ Draws 3d animation of whisker motion. """
    print '-'*22+'ANIMATE (3d)'+'-'*22

    X, Y, Z = load_points(file_name+'.p', 3)
    if debug: Xf, Yf, Zf = get_filtered_points(file_name, 3)

    plt.close('all')
    fig = plt.figure()
    ax = Axes3D(fig)

    pts = ax.plot([], [], [], '-o')
    if debug: fpts = ax.plot([], [], [], 'r-o')
    time_text = ax.text(0.01, 0.01, 0.0, '', transform=ax.transAxes)

    ax.set_xlim((np.min(X), np.max(X)))
    ax.set_ylim((np.min(Y), np.max(Y)))
    ax.set_zlim((np.min(Z), np.max(Z)))

    def init():
        pts[0].set_data([], [])
        pts[0].set_3d_properties([])
        if debug:
            fpts[0].set_data([], [])
            fpts[0].set_3d_properties([])
        time_text.set_text('')

    def animate(i):
        pts[0].set_data(X[i], Y[i])
        pts[0].set_3d_properties(Z[i])
        if debug:
            fpts[0].set_data(Xf[i], Yf[i])
            fpts[0].set_3d_properties(Zf[i])
        time_text.set_text('t = %.3f s' %(DT*i))
        fig.canvas.draw()

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(file_name)
    
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(X), interval=10, blit=False)

    if save_movie:
        anim.save('./output_data/%s/%s.mp4' %(file_name, file_name), fps=15)
                  #writer=animation.FFMpegFileWriter(),
                  #extra_args=['-vcodec', 'libx264'])

    if show:
        print 'done'
        sys.exit(plt.show())

def animate_whisker(dim, *args, **kwargs):
    if dim == 2:
        animate_whisker_2d(*args, **kwargs)
    elif dim == 3:
        animate_whisker_3d(*args, **kwargs)

def get_filtered_points(file_name, dim):
    converted_data = load_converted_data(file_name+'.p')
    filtered_data = load_converted_data(file_name+'_filtered.p')
    ref = filtered_data['ref']
    realstdout = sys.stdout
    sys.stdout = StringIO()
    whisker = make_whisker(dim, filtered_data['ref'], converted_data['link_length'])
    sys.stdout = realstdout
    N = whisker.num_links  
    n = len(converted_data['x'])

    X = np.empty((n, N+1))
    Y = np.empty((n, N+1))
    if dim > 2: Z = np.empty((n,N+1))
    for i in range(n):
        whisker.qd = filtered_data['q'][i]
        if dim == 2:
            X[i], Y[i] = zip(*whisker.config_points)  
        else:
            X[i], Y[i], Z[i] = zip(*whisker.config_points)        
       
    if dim == 2:
        return X, Y 

    return X, Y, Z

   
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='whisker animation options')
    parser.add_argument('data_file', help=".mat file with tracked image data")
    parser.add_argument('--dim', help="dimension (2 or 3)", default=2, type=int)
    parser.add_argument('--save', help="save movie", action='store_true')
    parser.add_argument('--show', help="show animation", action='store_true')
    parser.add_argument('--debug', help="debug config conversion", 
                        action='store_true')
    args = parser.parse_args()

    assert args.dim in [2,3], "dimension must be 2 or 3"

    if args.dim == 2:
        animate_whisker_2d(args.data_file, save_movie=args.save,
                           show=args.show, debug=args.debug)

    else:
        animate_whisker_3d(args.data_file, save_movie=args.save, 
                           show=args.show, debug=args.debug)

