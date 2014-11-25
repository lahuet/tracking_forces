import numpy as np
import os
import argparse
import scipy.io as sio
import trep

from util import *

"""
20140221 - Lucie edits code to call self.nodes instead of self.masses -
a fallout of moving the masses to link center of masses instead of nodes.
"""


def update_contact(whisker, CP):
    """Adjusts the contact constraint to match current contact point."""
    
    # If the whisker is not in contact, remove constraint if necessary.
    if np.isnan(CP[0]):
        if whisker.in_contact:
            remove_peg_constraint(whisker)

    # Otherwise, the whisker is in contact and we need to update the point of
    # contact.
    else:
        if len(CP) == 2:
            whisker.get_config('py').q = CP[0]
            whisker.get_config('pz').q = CP[1]
            plane_frame = get_closest_frame(whisker, CP)
        else:
            whisker.get_config('px').q = CP[0]
            whisker.get_config('py').q = CP[1]
            whisker.get_config('pz').q = CP[2]
            plane_frame = get_closest_frame(whisker, CP)

        # If the peg constraint already exists in the system, simply move the
        # plane frame.
        if whisker.in_contact:
            whisker.get_constraint('peg')._plane_frame = plane_frame

        # Otherwise we need to add the constraint to the system.
        else:
            point_frame = whisker.get_frame('peg_frame')
            if len(CP) == 2:
                # In 2d the normal to the whisker is always in the z-direction
                # of the local frame.
                n = (0,0,-1)
            else:
                # The normal to the whisker at the point of contact is chosen 
                # to lie entirely in the xy plane, normal to the projection in 
                # the xy plane. To do this, we calculate the normal in the
                # world frame and then rotate back to the local frame.
                # This means there is no force in the z-direction due to the peg
                # constraint.
                x1, y1 = plane_frame.p()[:2]
                # next_frame = whisker.masses[whisker.masses.index(plane_frame)+1]
                next_frame = whisker.nodes[whisker.nodes.index(plane_frame)+1]
                x2, y2 = next_frame.p()[:2]
                nw = (-(y2-y1), x2-x1, 0.0)
                RT = np.transpose(plane_frame.g()[:-1,:-1])
                n = np.dot(RT, nw)

            trep.constraints.PointOnPlane(whisker, plane_frame, n, 
                                          point_frame, name='peg')

def get_closest_frame(whisker, CP):
    """Gets the whisker frame that is the closest to the point of contact."""
    if len(CP) == 2:
        x, y = zip(*whisker.config_points)
        index = closest_point_2d(x, y, CP)
    else:
        x, y, z = zip(*whisker.config_points)
        index = closest_point_3d(x, y, z, CP)
    # return whisker.masses[index]
    return whisker.nodes[index]

def remove_peg_constraint(whisker):
    """Removes the peg constraint from the whisker system."""
    new_cons = tuple()
    for c in whisker.constraints:
        if c.name != 'peg':
            new_cons += (c,)
    whisker._constraints = new_cons
    whisker._structure_changed()

def constraint_vec(whisker, length):
    """Returns two points on the normal of the point constraint. (debug)"""
    if not whisker.in_contact:
        return ([], [])
    plane_frame = whisker.contact_frame
    p1 = plane_frame.p()[1:3]
    p2 = np.dot(plane_frame.g(), np.array((0,0,length,1)))[1:3]
    return ([p1[0], p2[0]], [p1[1], p2[1]])

def static_constraint_forces(sys):
    """
    Static constraint forces are calculted assuming dq = ddq = 0. This solves
    the EL equations for lambda in the least squares sense.
    """
    dqd = sys.dqd
    sys.dqd = 0.0

    A = np.array([[con.h_dq(q1) for q1 in sys.dyn_configs] 
                                for con in sys.constraints])
    AT = np.transpose(A)
    F = np.array([sum([force.f(q1) for force in sys.forces])
                                   for q1 in sys.dyn_configs])
    dLdq = np.array([sys.L_dq(q1) for q1 in sys.dyn_configs])

    # B = np.dot(np.linalg.inv(np.dot(A, AT)), A)
    # Lucie: pseudoinverse insertion!!
    B = np.linalg.pinv(np.transpose(A))
    D = dLdq + F
    LAM = -np.dot(B, D)

    # Reset velocity so system state is unchanged by this method
    sys.dqd = dqd

    return LAM

def dynamic_constraint_forces(sys, ddq):
    """ 
    Dynamic forces for system with nonzero velocity/acceleration. Uses current
    dq of system and specified ddq. Solves the EL equations for lambda in the
    least squares sense.
    """
    # dh/dq
    A = np.array([[con.h_dq(q1) for q1 in sys.dyn_configs] 
                                for con in sys.constraints])
    AT = np.transpose(A)

    # array of external forcing for each dynamic config
    F = np.array([sum([force.f(q1) for force in sys.forces]) 
                                   for q1 in sys.dyn_configs])
    dLdq = np.array([sys.L_dq(q1) for q1 in sys.dyn_configs])
    dLddqdq = np.array([[sys.L_ddqdq(q2, q1) for q1 in sys.dyn_configs] 
                                             for q2 in sys.dyn_configs])
    dLddqddq = np.array([[sys.L_ddqddq(q1, q2) for q1 in sys.dyn_configs] 
                                               for q2 in sys.dyn_configs])

    # B = np.dot(np.linalg.inv(np.dot(A, AT)), A)
    # Lucie: pseudoinverse insertion!!
    B = np.linalg.pinv(np.transpose(A))
    D = np.dot(dLddqddq, ddq) + np.dot(dLddqdq, sys.dqd) - dLdq - F
    LAM = np.dot(B, D)

    return LAM

def discrete_constraint_forces(sys, qd0, qd1, qd2, dt, dynamic=True):
    """
    Calculates an approximation to the continuous constraint forces by
    solving the DEL equation for lambda.
    """
    # To calculate the static forces, we set the configuration to q1 for the
    # previous and following steps.
    q, dq = sys.q, sys.dq
    if not dynamic:
        qd0 = qd1.copy()
        qd2 = qd1.copy()
    
    sys.qd = 0.5*(qd1+qd2)
    sys.dqd = (qd2-qd1)/dt
    D1Ld = np.array([0.5*dt*sys.L_dq(q1)-sys.L_ddq(q1) for q1 in sys.dyn_configs])
    Fm = np.array([sum([dt*F.f(q1) for F in sys.forces]) for q1 in sys.dyn_configs])
    
    sys.qd = 0.5*(qd0+qd1)
    sys.dqd = (qd1-qd0)/dt
    D2Ld = np.array([0.5*dt*sys.L_dq(q1)+sys.L_ddq(q1) for q1 in sys.dyn_configs])
    
    sys.qd = qd1
    sys.dqd = 0.0
    Dh = np.array([[c.h_dq(q1) for q1 in sys.dyn_configs] for c in
                    sys.constraints])
    
    # dh_inv = np.linalg.inv(np.dot(Dh, Dh.T))
    # Lucie: pseudoinverse insertion!!
    dh_inv = np.linalg.pinv(np.dot(Dh, Dh.T))
    lam_bar = np.dot(dh_inv, np.dot(Dh,(D2Ld + D1Ld + Fm)))
    sys.set_state(q, dq) # Reset the initial system state.

    LAM = lam_bar/dt
    return -LAM

def save_forces_to_file(file_name, path_name, data, overwrite=True):
    """Saves force data to .p and .mat files."""
    if not os.path.exists(path_name):
        os.makedirs(path_name)
        
    f = open('%s%s_forces.p' %(path_name, file_name), 'w')
    pickle.dump(data, f)
    f.close()
    sio.savemat('%s%s_forces.mat' %(path_name, file_name), data)

def calc_forces_2d(whisker, qd, dqd, ddqd, CP, dt, file_name, path_name, ow):
    """Calculates the constraint forces given the dynamics."""
    print '-'*20+'CALC FORCES (2d)'+'-'*20
    n_steps = len(qd) 

    RXN_DYN = {'M': np.empty(n_steps), 'FY': np.empty(n_steps),
              'FZ': np.empty(n_steps)}
    RXN_STC = {'M': np.empty(n_steps), 'FY': np.empty(n_steps),
              'FZ': np.empty(n_steps)}
    RXN_DSC = {'M': np.empty(n_steps), 'FY': np.empty(n_steps),
              'FZ': np.empty(n_steps)}
    T = np.empty(n_steps)
    m_index = whisker.get_constraint('theta-0_constraint').index
    fy_index = whisker.get_constraint('FY').index
    fz_index = whisker.get_constraint('FZ').index
    t = 0.0

    # Iterate over the time steps.
    print 'Calculating constraint forces...',
    for i in range(n_steps):

        # Set system configuration and velocity.
        whisker.qd  = qd[i]
        whisker.dqd = dqd[i]

        # Update the contact point/constraint.
        update_contact(whisker, CP[i])
        
        # Calculate constraint forces.
        lam_dyn = dynamic_constraint_forces(whisker, ddqd[i])
        lam_static = static_constraint_forces(whisker)
        if (i > 0) and (i < (n_steps-1)):
            lam_discrete = discrete_constraint_forces(whisker, qd[i-1], qd[i],
                    qd[i+1], dt, dynamic=True)
        else:
            lam_discrete = np.zeros(3)

        # Save the reaction forces and time. 
        RXN_DYN['M'][i]  = lam_dyn[m_index]
        RXN_DYN['FY'][i] = lam_dyn[fy_index]
        RXN_DYN['FZ'][i] = lam_dyn[fz_index]
        RXN_STC['M'][i]  = lam_static[m_index]
        RXN_STC['FY'][i] = lam_static[fy_index]
        RXN_STC['FZ'][i] = lam_static[fz_index]
        RXN_DSC['M'][i]  = lam_discrete[m_index]
        RXN_DSC['FY'][i] = lam_discrete[fy_index]
        RXN_DSC['FZ'][i] = lam_discrete[fz_index]

        T[i] = t
        t += dt
    print 'done'    

    FORCES = {'dyn': RXN_DYN, 'static': RXN_STC, 'discrete': RXN_DSC, 't': T, 'dim': 2}

    # Save the results to a file.
    print 'Saving to file...',
    save_forces_to_file(file_name, path_name, FORCES, ow)    
    print 'done (Saved to %s/%s/%s_forces.mat(.p))' %(path_name, file_name, file_name)
    return FORCES

def calc_forces_3d(whisker, qd, dqd, ddqd, CP, dt, file_name, path_name, ow):
    """ Calculates the constraint forces given the dynamics. """
    print '-'*20+'CALC FORCES (3d)'+'-'*20
    n_steps = len(qd) 

    RXN_DYN = {'MX': np.empty(n_steps), 'MY': np.empty(n_steps), 
               'MZ': np.empty(n_steps),
               'FX': np.empty(n_steps), 'FY': np.empty(n_steps), 
               'FZ': np.empty(n_steps)}
    RXN_STC = {'MX': np.empty(n_steps), 'MY': np.empty(n_steps), 
               'MZ': np.empty(n_steps),
               'FX': np.empty(n_steps), 'FY': np.empty(n_steps), 
               'FZ': np.empty(n_steps)}
    RXN_DSC = {'MX': np.empty(n_steps), 'MY': np.empty(n_steps), 
               'MZ': np.empty(n_steps),
               'FX': np.empty(n_steps), 'FY': np.empty(n_steps), 
               'FZ': np.empty(n_steps)}
    T = np.empty(n_steps)
    mx_index = whisker.get_constraint('theta_x-0_constraint').index
    my_index = whisker.get_constraint('theta_y-0_constraint').index
    mz_index = whisker.get_constraint('theta_z-0_constraint').index
    fx_index = whisker.get_constraint('FX').index
    fy_index = whisker.get_constraint('FY').index
    fz_index = whisker.get_constraint('FZ').index
    t = 0.0

    # Iterate over the time steps.
    print 'Calculating constraint forces...',
    for i in range(n_steps):

        # Set system configuration and velocity.
        whisker.qd  = qd[i]
        whisker.dqd = dqd[i]

        # Update the contact point/constraint.
        update_contact(whisker, CP[i])
        
        # Calculate constraint forces.
        lam_dyn = dynamic_constraint_forces(whisker, ddqd[i])
        lam_static = static_constraint_forces(whisker)
        if (i > 0) and (i < (n_steps-1)):
            lam_discrete = discrete_constraint_forces(whisker, qd[i-1], qd[i],
                    qd[i+1], dt, dynamic=True)
        else:
            lam_discrete = np.zeros(6)
                
        # Save the reaction forces and time.     
        RXN_DYN['MX'][i] = lam_dyn[mx_index]
        RXN_DYN['MY'][i] = lam_dyn[my_index]
        RXN_DYN['MZ'][i] = lam_dyn[mz_index]
        RXN_DYN['FX'][i] = lam_dyn[fx_index]
        RXN_DYN['FY'][i] = lam_dyn[fy_index]
        RXN_DYN['FZ'][i] = lam_dyn[fz_index]

        RXN_STC['MX'][i] = lam_static[mx_index]
        RXN_STC['MY'][i] = lam_static[my_index]
        RXN_STC['MZ'][i] = lam_static[mz_index]
        RXN_STC['FX'][i] = lam_static[fx_index]
        RXN_STC['FY'][i] = lam_static[fy_index]
        RXN_STC['FZ'][i] = lam_static[fz_index]

        RXN_DSC['MX'][i] = lam_discrete[mx_index]
        RXN_DSC['MY'][i] = lam_discrete[my_index]
        RXN_DSC['MZ'][i] = lam_discrete[mz_index]
        RXN_DSC['FX'][i] = lam_discrete[fx_index]
        RXN_DSC['FY'][i] = lam_discrete[fy_index]
        RXN_DSC['FZ'][i] = lam_discrete[fz_index]
        
        T[i] = t
        t += dt
    print 'done'    

    FORCES = {'dyn': RXN_DYN, 'static': RXN_STC, 'discrete': RXN_DSC, 't': T, 'dim': 3}

    # Save the results to a file.
    print 'Saving to file...',
    save_forces_to_file(file_name, path_name, FORCES, ow)
    print 'done (Saved to %s%s_forces.mat(.p))' %(path_name, file_name)
    return FORCES

def calc_forces(dim, *args, **kwargs):
    """Wrapper for calculating the forces in the correct dimension."""
    if dim == 2:
        calc_forces_2d(*args, **kwargs)
    elif dim == 3:
        calc_forces_3d(*args, **kwargs)

