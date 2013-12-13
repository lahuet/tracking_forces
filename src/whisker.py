import trep
from trep import tx, ty, tz, rx, ry, rz
import numpy as np
from scipy.interpolate import splev, splprep, interp1d
from scipy.optimize import fmin

from util import extract_angles

class FixedConfig(trep.Constraint):
    """A constraint that fixes a configuration variable to a constant value."""
    def __init__(self, system, config_name, reference):
        trep.Constraint.__init__(self, system, name=config_name+'_constraint')
        self.config = system.get_config(config_name)
        self.reference = reference
        
    def h(self):
        return (self.config.q-self.reference)

    def h_dq(self, q1):
        if q1 == self.config:
            return 1.0
        return 0.0

    def h_dqdq(self, q1, q2):
        return 0.0

    def h_dt(self):
        return 0.0

class Whisker2D(trep.System):
    """A two-dimensional whisker."""
    def __init__(self, parameters):
        trep.System.__init__(self)

        self.dim = 2
        self.num_links = parameters['N']
        self.lengths = parameters['L'] 
        self.import_frames(self.whisker_frames(parameters['L']))
        self.add_base_constraints()

        self.add_node_springs(parameters['k'])
        self.add_node_damping(parameters['c'])
        self.set_masses(parameters['m'])
        self._q0 = None
        
    def whisker_frames(self, L):
        """Creates a list of frames that define the whisker."""
        frames = []
        for j in reversed(range(1, self.num_links)):
            frames = [ty(L, name='Link-%d' %j), frames]
            frames = [rx('theta-%d' %j, name='f%d' %j), frames]
            frames = [rx('curvature-%d' %j, kinematic=True), frames]
        frames = [ty(L, name='Link-0'), frames]
        frames = [tz('z'), [ty('y', name='Base_Point'), frames]]
        frames = [rx('theta-0', name='f0'), frames]
        frames = [tz('z0', kinematic=True), [ty('y0', kinematic=True), frames]]
        frames = [rx('curvature-0', kinematic=True), frames]
        return frames
            
    def set_masses(self, m):
        """Sets masses at each of the whisker's nodes."""
        for j in range(self.num_links):
            self.get_frame('Link-%d' %j).set_mass((m[j],0,0,0))
        self.get_frame('Base_Point').set_mass((m[0],0,0,0))
    
    def add_node_springs(self, k):
        """Adds configuration (torsional) springs at each node."""
        for j in range(1,self.num_links):
            trep.potentials.ConfigSpring(self, 'theta-%d' %j, k=k[j], 
                    name='spring-%d' %j)

    def add_node_damping(self, c, default=0.0):
        """Adds damping at each node."""
        damp_dict = {}
        for j in range(1,self.num_links):
            damp_dict.update({'theta-%d' %j: c[j]})   
        trep.forces.Damping(self, default, damp_dict)

    def add_base_constraints(self):
        """ 
        Adds constraints on the rotation and translation at whisker's base.
        """
        FixedConfig(self, 'theta-0', 0.0)
        trep.constraints.PointOnPlane(self, 'Base_Point', (0,0,1), 
                'World', name='FZ') 
        trep.constraints.PointOnPlane(self, 'Base_Point', (0,-1,0),
                'World', name='FY')
     
    @property
    def arc_length(self):
        """ Returns arc length of the whisker by summing link lengths. """
        return arc_length(self.config_points)  

    @property
    def config_points(self):
        """ 
        Returns list of points of the whisker's nodes in the current 
        configuration. 
        """
        return np.array([f.p()[1:3] for f in self.masses])

    @config_points.setter
    def config_points(self, points):
        """ 
        Sets the positions of the whisker's nodes to match a given list 
        of points.
        """
        ang = points_to_angles(points)
        for i in range(self.num_links):
            self.get_config('theta-%d' %i).q = ang[i]

    @property
    def config_angles(self):
        return np.array([cnfg.q for cnfg in self.dyn_configs
                                           if 'theta' in cnfg.name])

    @config_angles.setter
    def config_angles(self, angles):
        for i in range(self.num_links):
            self.get_config('theta-%d' %i).q =\
                    (angles[i]-self.get_config('curvature-%d' %i).q)

    @property
    def reference_shape(self):
        return self._q0

    @reference_shape.setter
    def reference_shape(self, ref_shape):
        self.qk = ref_shape
        self._q0 = ref_shape

    @property
    def in_contact(self):
        """
        Returns true if whisker is in contact with the peg, 
        false otherwise. 
        """
        for con in self.constraints:
            if con.name == 'peg':
                return True
        return False        
       
    @property
    def base_angle(self):
        """Angle of the first link in the world frame."""
        p0 = self.masses[0].p()[1:3]
        p1 = self.masses[1].p()[1:3]
        return np.arctan2(p1[1]-p0[1], p1[0]-p0[0])

    def get_angle(self, num):
        return self.get_config('theta-%d' %num)

class Whisker3D(trep.System):
    """A three-dimensional whisker."""
    def __init__(self, parameters, ref):
        trep.System.__init__(self)

        self.dim = 3
        self.num_links = parameters['N']
        self._link_length = parameters['L'] 
        self._ref = ref
        self.import_frames(self.whisker_frames())
        self.add_base_constraints()

        self.add_node_springs(parameters['k'])
        self.add_node_damping(parameters['c'])
        self.set_masses(parameters['m'])
        
    def whisker_frames(self):
        """Creates a list of frames that define the whisker."""
        x0, y0, z0, th_x0, th_y0, th_z0 = extract_angles(self._ref)
        L = self._link_length
        frames = []
        for j in reversed(range(1, self.num_links)):
            frames = [tx(L, name='Link-%d' %j), frames]
            frames = [rz('theta_z-%d' %j), [ry('theta_y-%d' %j), frames]]
            frames = [rz(th_z0[j]),[ry(th_y0[j]),frames]]
        frames = [tx(L, name='Link-0'), frames]
        frames = [tx('xb'), [ty('yb'), [tz('zb', name='Rotated_Base_Point'), 
                   frames]]]
        frames = [rz('theta_z-0'), [ry('theta_y-0'), [rx('theta_x-0'),
                   frames]]]
        frames = [rz(th_z0[0]),[ry(th_y0[0]), [rx(th_x0[0]), frames]]]    
        frames = [tx(x0), [ty(y0), [tz(z0, name='Base_Point'), frames]]]
        return frames

    def set_masses(self, m):
        """Sets masses at each of the whisker's nodes."""
        for j in range(self.num_links):
            self.get_frame('Link-%d' %j).set_mass((m[j],0,0,0))
        self.get_frame('Rotated_Base_Point').set_mass((m[0],0,0,0))
        #self.get_frame('Base_Point').set_mass((m[0],0,0,0))
    
    def add_node_springs(self, k):
        """Adds configuration (torsional) springs at each node."""
        for j in range(1,self.num_links):
            trep.potentials.ConfigSpring(self, 'theta_y-%d' %j, k=k[j],
                    name='spring_y-%d' %j)
            trep.potentials.ConfigSpring(self, 'theta_z-%d' %j, k=k[j],
                    name='spring_z-%d' %j)

    def add_node_damping(self, c, default=0.0):
        """Adds damping at each node."""
        damp_dict = {}
        for j in range(1,self.num_links):
            damp_dict.update({'theta_y-%d' %j: c[j], 'theta_z-%d' %j: c[j]})   
        trep.forces.Damping(self, default, damp_dict)

    def add_base_constraints(self):
        """
        Adds constraints on the rotation and translation at whisker's base.
        """
        FixedConfig(self, 'theta_x-0', 0.0)
        FixedConfig(self, 'theta_y-0', 0.0)
        FixedConfig(self, 'theta_z-0', 0.0)
        trep.constraints.PointOnPlane(self, 'Rotated_Base_Point',
                (0,0,1), 'Base_Point', name='FZ') 
        trep.constraints.PointOnPlane(self, 'Rotated_Base_Point',
                (0,1,0), 'Base_Point', name='FY')
        trep.constraints.PointOnPlane(self, 'Rotated_Base_Point',
                (1,0,0), 'Base_Point', name='FX')
     
    @property
    def arc_length(self):
        """Returns arc length of the whisker by summing link lengths."""
        return arc_length_3d(self.config_points)  

    @property
    def config_points(self):
        """
        Returns list of points of the whisker's nodes in the current 
        configuration. 
        """
        return np.array([f.p()[:-1] for f in self.masses])

    @property
    def config_angles(self):
        return np.array([cnfg.q for cnfg in self.dyn_configs 
                                             if 'theta' in cnfg.name])

    @property
    def reference_shape(self):
        return self._ref

    @property
    def in_contact(self):
        """ 
        Returns true if whisker is in contact with the peg, false otherwise.
        """
        for con in self.constraints:
            if con.name == 'peg':
                return True
        return False        
       
    @property
    def base_angle(self):
        """Angle of the first link in the world frame."""
        return self.get_config('theta-0').q + self.get_config('curvature-0').q

    def get_angle(self, num, axis):
        return self.get_config('theta_%s-%d' %(axis, num))


def make_whisker(dim, q0, L, rbase=100e-6, taper=1.0/15, damping_ratio=None,
        rho=1.0, E=3.3e9):
    """ 
    Assembles the whisker based on material properties and initial 
    configuration. The default material and geometric properties from 
    elastica2d are used by default. 
    """
    assert dim in [2,3], "dimension must be 2 or 3"    
    print 19*'-'+'BUILD WHISKER (%dd)'%dim+19*'-'

    print 'Getting parameters...',
    if dim==2: N = len(q0)-2
    else: N = (len(q0)-4)/2
    I = calc_inertia(N, rbase, taper)
    K = E*I/L
    M = calc_mass(L, N, rho, rbase, taper)
    if damping_ratio<0:
        damping_ratio = np.append([1.0], get_interp_damping(L,N))
    C = calc_damping(N, K, M, L, damping_ratio)
    parameters = {'L': L, 'k': K[:-1], 'c': C, 'm': M, 'N':N}
    print 'done'

    print 'Building system...',
    if dim==2:
        whisker = Whisker2D(parameters)
        whisker.reference_shape = q0
        peg_frames = [ty('py', kinematic=True),
                      [tz('pz', kinematic=True, name='peg_frame')]]
    else:
        whisker = Whisker3D(parameters, q0)
        peg_frames = [tx('px', kinematic=True), [ty('py', kinematic=True),
                      [tz('pz', kinematic=True, name='peg_frame')]]]
    whisker.world_frame.import_frames(peg_frames)
    print 'done'
    return whisker

def calc_inertia(N, rbase, taper):
    """The moment of inertia for a tapered cylinder with N points."""
    R = np.linspace(rbase, rbase*taper, N+1)
    return 0.25*np.pi*R**4

def calc_mass(L, N, rho, rbase, taper):
    """The mass at each node."""
    R = np.linspace(rbase, rbase*taper, N+1)
    return rho*np.pi*(R**2)*L

def calc_damping(N, k, m, L, zeta):
    """
    Calculates damping coefficients at each node from a given damping ratio.
    """
    return 2*L*np.sqrt(k*m)*zeta

def get_interp_mass(L, N):
    """Calculates mass by interpolation of SWD values."""
    f = interp1d(np.linspace(0, sum(SWD['L']), SWD['N']), SWD['m'], 
            fill_value=0.0, bounds_error=False)
    return f(np.linspace(L, L*N, N))

def get_interp_damping(L, N):
    """Calculates damping by interpolaton of SWD values. """
    f = interp1d(np.linspace(0, sum(SWD['L']), SWD['N']), SWD['zeta'], 
            fill_value=0.0, bounds_error=False)
    return f(np.linspace(L, L*N, N))

# Default parameters for a single whisker, used for interpolation of values for 
# all whiskers. These are the latest values from CND.
SWD_old = {'N': 13,
       'L': [4.08e-3]*13,
       'k': [3.26e-5, 3.09e-5, 2.68e-5, 2.16e-5, 1.64e-5, 1.15e-5, 7.47e-6,
             4.40e-5, 2.28e-6, 9.80e-7, 3.17e-7, 5.97e-8, 2.73e-9],
       'c': [1.19e-7, 7.85e-8, 4.98e-8, 3.02e-8, 1.72e-8, 9.18e-9, 4.46e-9,
             1.91e-9, 6.81e-10, 1.81e-10, 2.59e-11, 1.03e-12, 4.50e-12],
       'm': [3.00e-7, 2.35e-7, 1.80e-7, 1.34e-7, 9.72e-8, 6.76e-8, 4.47e-8,
             2.77e-8, 1.56e-8, 7.71e-9, 3.06e-9, 8.01e-10, 6.83e-11]}        

SWD = {'N': 13, 
       'L': [3.78e-3]*13,
       'k': [65.3e-6, 46.5e-6, 32.2e-6, 21.5e-6, 13.9e-6, 8.53e-6, 4.94e-6,
             2.64e-6, 1.27e-6, 0.529e-6, 0.175e-6, 0.00358e-6],
       'c': [],
       'zeta': [0.988, 0.909, 0.829, 0.750, 0.670, 0.590, 0.511, 0.431, 0.352,
           0.272, 0.192, 0.113, 0.033],
       'm': [104, 88.7, 74.2, 61.1, 49.9, 38.6, 29.3, 21.3, 14.5, 9.09, 4.92,
           2.03, 0.42]}
