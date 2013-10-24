import numpy as np
import scipy
import scipy.optimize
import scipy.io as sio
from scipy.interpolate import splprep, splev, interp1d, splrep
from scipy.spatial.distance import cdist
import pickle
import ConfigParser


def dist_from_points(x, y):
    """Array of distances between a series of points."""
    p = np.reshape(zip(x,y), (len(x),2))
    return cdist(p[1:], p[:-1]).diagonal()

def arc_length(x, y, z=None):
    """Sum of euclidean distances between points."""
    if z is None:
        return sum(dist_from_points(x,y))
    return sum(dist_from_points_3d(x,y,z))

def closest_point_spline(x, y, CP):
    """ 
    Returns the point on a spline constructed from the points (x,y) that is 
    closest to CP. 
    """
    tck, u = splprep([x, y], k=5)
    X, Y = CP[0], CP[1]
    def F(t):
        xt, yt = splev(t, tck)
        return np.linalg.norm((xt-X, yt-Y))
    t_guess = 0.5
    sol = scipy.optimize.fmin(F, t_guess, disp=0)
    return splev(sol, tck)

def rotate(x, y, theta):
    """Rotates the points (x,y) through an angle theta about the origin."""
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]])
    p = np.row_stack((x,y))
    rp = np.dot(R, p)
    return (rp[0,:], rp[1,:])  

def closest_point_2d(x, y, p):
    """ 
    Given a list of points (x,y), returns the index of the point that is
    closest to p.
    """
    XY = np.column_stack((x,y))
    p = np.reshape(p, (1,2))
    d = cdist(XY, p)
    return np.argmin(d)

#def wrap_to_pi(num):
#    """Converts angle to between -pi and pi."""
#    return (num + np.pi) % (2*np.pi)-np.pi

def closest_point_3d(x, y, z, p):
    """ 
    Given a list of points (x,y), returns the index of the point that is 
    closest to p.
    """
    XYZ = np.column_stack((x,y,z))
    p = np.reshape(p, (1,3))
    d = cdist(XYZ, p)
    return np.argmin(d)

def points_to_angles_2d(points):
    """ 
    Takes a set of points in R^2 and converts to a unique angle configuration.
    """
    x, y = zip(*points)
    angles = np.arctan2(np.diff(y),np.diff(x)) 
    for i in range(1,len(points)-1):
        angles[i] = angles[i] - sum(angles[:i])
    return angles

def get_spike_times(raw_data_file, t_max, dt):
    """Reads raw data file to find all spikes that occur before time t_max."""
    try:
        data = sio.loadmat('./raw_data/%s' %raw_data_file)
        spikes = data['AllSpikesOut'].astype(int).flatten()
    except:
        return []
    spike_times = dt*(np.where(spikes)[0])
    return spike_times[np.where(spike_times<t_max)[0]]

def load_file(file_name):
    """Loads a pickled data file."""
    try:
        f = open(file_name)
        data = pickle.load(f)
        f.close()
    except IOError:
        raise Exception("Could not load file %s" %file_name)
    return data

def dist_from_points_3d(x, y, z):
    """Gives a list of the distances between a series of points."""
    p = np.reshape(zip(x,y,z), (len(x),3))
    return cdist(p[1:], p[:-1]).diagonal()

def arc_length_3d(x, y, z):
    """Sum of euclidean distances between points (x,y)."""
    return sum(dist_from_points_3d(x,y,z))

def closest_point_spline_3d(x, y, z, CP):
    """ 
    Returns the point on a spline constructed from the points 
    (x,y,z) that is closest to CP.
    """
    tck, u = splprep([x, y, z], k=5)
    def F(t):
        xt, yt, zt = splev(t, tck)
        return np.linalg.norm((xt-CP[0], yt-CP[1], zt-CP[2]))
    t_guess = 0.5
    sol = scipy.optimize.fmin(F, t_guess, disp=0)
    return splev(sol, tck)

# Homogeneous transformations.
def g_rz(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0, 0],
                     [np.sin(theta),  np.cos(theta), 0, 0],
                     [            0,              0, 1, 0],
                     [            0,              0, 0, 1]])

def g_ry(theta):
    return np.array([[ np.cos(theta), 0, np.sin(theta), 0],
                     [             0, 1,             0, 0],
                     [-np.sin(theta), 0, np.cos(theta), 0],
                     [             0, 0,             0, 1]])
def g_rx(theta):
    return np.array([[1,             0,              0, 0],
                     [0, np.cos(theta), -np.sin(theta), 0],
                     [0, np.sin(theta),  np.cos(theta), 0],
                     [0,             0,              0, 1]])

def g_tx(x):
    return np.array([[1,0,0,x], [0,1,0,0], [0,0,1,0], [0,0,0,1]])

def g_frame(i, theta_x, theta_y, theta_z, L):
    """Returns the global rotation for the i-th frame."""
    grx = g_rx(theta_x[i])
    gry = g_ry(theta_y[i])
    grz = g_rz(theta_z[i])
    gtx = g_tx(L[i])
    g = grz.dot(gry).dot(grx).dot(gtx)
    if i == 0:
        return g
    return np.dot(g_frame(i-1, theta_x, theta_y, theta_z, L), g)
    
def points_to_angles_3d(points, rx0=0.0):
    """ 
    Converts a set of 3d (x,y,z) points to an angle configuration. The angles
    are unique up to a rotation about the x-axis at the base (rx0). The angles
    follow the convention R = Rz(theta_z)*Ry(theta_y)*Rx(theta_x).
    """
    n = len(points)
    p = np.asarray(points)

    x, y, z = zip(*points)
    L = dist_from_points_3d(x,y,z)

    theta_x = np.empty(n-1)
    theta_y = np.empty(n-1)
    theta_z = np.empty(n-1)

    r = p[1]
    theta_x[0] = rx0
    theta_y[0] = -np.arctan2(r[2], np.sqrt(r[0]**2+r[1]**2))
    theta_z[0] = np.arctan2(r[1],r[0])
 
    for i in range(1,n-1):
        g = g_frame(i-1, theta_x, theta_y, theta_z, L)
        g_inv = np.linalg.inv(g)
        r = np.dot(g_inv, np.append(p[i+1],1))
        theta_x[i] = 0.0
        theta_y[i] = -np.arctan2(r[2], np.sqrt(r[0]**2+r[1]**2))
        theta_z[i] = np.arctan2(r[1],r[0])

    return theta_x, theta_y, theta_z, L   

def points_to_angles(points, dim=2, rx0_method='zero'):
    """ 
    Wrapper for converting point data to angle trajectories. 
    Returned as three lists with theta_x, theta_y, theta_z.
    """
    assert dim in [2,3], "dimension must be 2 or 3"

    if dim == 2:
        return points_to_angles_2d(points)

    assert rx0_method in ['zero', 'min_angles']

    if rx0_method == 'zero':
        return points_to_angles_3d(points)[:-1]
    else:
        def func(rx0):
            theta_x, theta_y, theta_z, L = points_to_angles_3d(points, rx0=rx0)
            x = np.hstack((theta_x, theta_y, theta_z))
            return np.linalg.norm(x)

        #res = scipy.optimize.minimize_scalar(func, bracket=[-np.pi, np.pi])
        (xmin, fval, niter, funcalls) = scipy.optimize.brent(func,
                brack=[-np.pi, np.pi], full_output=True)
        return points_to_angles_3d(points, rx0=xmin)[:-1] # rx0=res.x

def angles_to_points_2d(theta, L):
    """Takes a set of angles and link lengths and returns a set of points."""
    n = len(theta)+1
    x = np.empty(n)
    y = np.empty(n)
    x[0], y[0] = (0,0)
    for i in range(1,n):
        g = g_frame(i-1, np.zeros(n), np.zeros(n), theta, L)
        p = np.dot(g, np.array((0,0,0,1)))
        x[i], y[i] = p[:-2]
    return x, y 

def angles_to_points_3d(theta_x, theta_y, theta_z, L):
    """Takes a set of angles and link lengths and returns a set of points."""
    n = len(theta_x)+1
    x = np.empty(n)
    y = np.empty(n)
    z = np.empty(n)
    x[0], y[0], z[0] = (0,0,0)
    for i in range(1,n):
        g = g_frame(i-1, theta_x, theta_y, theta_z, L)
        p = np.dot(g, np.array((0,0,0,1)))
        x[i], y[i], z[i] = p[:-1]
    return x, y, z   

def rotation_to_euler(R, zy_only=True):
    """ 
    Takes a rotation matrix and returns the euler angles. If zy_only is True,
    the matrix is assumed to be R=Rz*Ry. If zy_only is false, R=Rz*Ry*Rx. 
    """
    # Take theta_x = 0, only find rotations about z- and y-axes. 
    if zy_only:
        theta_z = np.arcsin(-R[0,1])
        theta_y = np.arcsin(-R[2,0])
        theta_x = 0.0
        return theta_x, theta_y, theta_z

    # Otherwise assume all rotations are nonzero. Find rotations about all axes.
    if (R[2,0] != 1) and (R[2,0] != -1):
        theta_y = -np.arcsin(R[2,0]); c1 = np.cos(theta_y)
        theta_x = np.arctan2(R[2,1]/c1, R[2,2]/c1)
        theta_z = np.arctan2(R[1,0]/c1, R[0,0]/c1)
    else:
        theta_z = 0.0
        if R[2,0] == -1:
            theta_y = np.pi/2
            theta_x = theta_z + np.arctan2(R[0,1], R[0,2])
        else:
            theta_y = -np.pi/2
            theta_x = -theta_z + np.arctan2(-R[0,1], -R[0,2])
    return theta_x, theta_y, theta_z        

def offset_euler_angles(theta_x, theta_y, theta_z, 
                        theta_x0, theta_y0, theta_z0):
    """Adjusts euler angles to be relative to some known offest."""
    Rb = np.dot(np.dot(g_rz(theta_z0), g_ry(theta_y0)), g_rx(theta_x0))[:-1,:-1]
    Ri = np.dot(np.dot(g_rz(theta_z), g_ry(theta_y)), g_rx(theta_x))[:-1,:-1]
    R = np.dot(Rb.T, Ri)
    if theta_x == 0:
        return rotation_to_euler(R)
    else:
        return rotation_to_euler(R, zy_only=False)

def relative_angle_config_2d(q, q_ref):
    """
    Offsets the angles in a 2D whisker relative to a reference configuration.
    """
    ang = q[2:]
    ref_ang = q_ref[2:]

    rel_config = np.empty(len(q))
    rel_config[0] = q[0]
    rel_config[1] = q[1]

    for i in range(len(q)-2):
        Ri = g_rz(ang[i])[:-1,:-1]
        Rb = g_rz(ref_ang[i])[:-1,:-1]
        R = np.dot(Rb.T, Ri)
        rel_config[i+2] = np.arcsin(-R[0,1])
   
    return rel_config

def relative_angle_config(q, q_ref):
    """
    Converts a full configuration to a configuration relative to a reference.
    """
    x, y, z, thx, thy, thz = extract_angles(q)
    x0, y0, z0, thx0, thy0, thz0 = extract_angles(q_ref)

    n = (len(q)-4)/2
    theta_x = np.empty(n)
    theta_y = np.empty(n)
    theta_z = np.empty(n)

    for i in range(n):
        theta_x[i], theta_y[i], theta_z[i]\
            = offset_euler_angles(thx[i],thy[i],thz[i],thx0[i],thy0[i],thz0[i])

    q_rel = np.empty(len(q))
    q_rel[3:6] = q[3:6]
    q_rel[0] = theta_z[0]
    q_rel[1] = theta_y[0]
    q_rel[2] =  theta_x[0]
    q_rel[6::2] = theta_z[1:]  
    q_rel[7::2] = theta_y[1:]
    return q_rel

def get_relative_configs(q, q_ref, dim):
    """Converts an entire trajectory to relative angles."""
    n = len(q)
    q_rel = np.empty((n, len(q_ref)))
    for i in range(n):
        if dim == 3:
            q_rel[i] = relative_angle_config(q[i], q_ref)
        elif dim == 2:
            q_rel[i] = relative_angle_config_2d(q[i], q_ref)
    return q_rel    

def extract_angles(config):
    """
    Given a system configuration, returns the base displacement and the x, y,
    and z rotations.
    """
    x, y, z = config[3:6]
    th_z = np.append(config[0], config[6::2])
    th_y = np.append(config[1], config[7::2])
    th_x = np.append(config[2], np.zeros(len(th_z)-1))
    return x, y, z, th_x, th_y, th_z

def get_reference_frame(cp):
    """
    Returns the index of the frame just before first contact. If there are 
    no contacts, returns 0. Return i-1 where i is the index of the first
    non-nan element. 
    """
    index = np.where(~np.isnan(cp[:,1]))[0]
    if index.any():
        return index[0]
    else:
        return 0


def get_base_points_2d(xbw, ybw, theta_0):
    """Gets the base point of the whisker in the rotated base coordinates."""
    pass

def get_base_points_3d():
    pass

