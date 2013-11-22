import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np

def plot_stuff1():
    f = open('./Rat1206_N04_Tr09/Rat1206_N04_Tr09.p')
    data = pickle.load(f)
    f.close()

    X = data['x']*1000
    Y = data['y']*1000

    plt.rc('text', usetex=True)
    fig = plt.figure(facecolor='white')
    ax = plt.axes()

    for i in xrange(0,200,1):
        ax.plot(X[i],Y[i])

    x1, y1 = 2.5, -1.75
    x2, y2 = 1.0, 1.25

    #ax.plot([x1, x2], [y1, y2], ".")
    ax.annotate("",
            xy=(x2,y2), xycoords='data',
            xytext=(x1,y1), textcoords='data',
            arrowprops=dict(arrowstyle="simple",
                            color="0.5", 
                            shrinkA=5, 
                            shrinkB=5, 
                            patchA=None, 
                            patchB=None,
                            connectionstyle="angle3,angleA=90,angleB=0"),
            )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.xlabel('$x\;\;(\mathrm{mm})$')
    plt.ylabel('$y\;\;(\mathrm{mm})$')

    plt.show()

def plot_stuff2():
    f = open('./Rat1206_N04_Tr09/Rat1206_N04_Tr09_forces.p')
    data = pickle.load(f)
    f.close()

    # get data.
    N = 600
    N = min(N,len(data['t']))
    scale = 1.0e6
    T = data['t']
    FZ_static = data['static']['FZ']*scale
    FY_static = data['static']['FY']*scale
    M_static = data['static']['M']*scale
    FZ_dyn = data['dyn']['FZ']*scale
    FY_dyn = data['dyn']['FY']*scale
    M_dyn = data['dyn']['M']*scale
    M_dyn_dis = data['discrete']['M']*scale

    # matplotlib settings.
    plt.rc('text', usetex=True)
    fig = plt.figure(num=1, facecolor='white',figsize=(12,5))
    ax = plt.axes()

    # Axial force.
    ax.plot(T[:N],FY_static[:N],T[:N],FY_dyn[:N])
    ax.legend(['$\mathrm{quasistatic}$','$\mathrm{dynamic}$'],4)
    plt.title('$\mathrm{Axial\;Force},\;F_x$')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.xlabel('$t\;\;(\mathrm{sec})$')
    plt.ylabel('$F\;\;(\mu\mathrm{N})$')
    
    # transverse force
    fig = plt.figure(num=2, facecolor='white',figsize=(12,5))
    ax = plt.axes()
    ax.plot(T[:N],-1.0*FZ_static[:N],T[:N],-1.0*FZ_dyn[:N])
    ax.legend(['$\mathrm{quasistatic}$','$\mathrm{dynamic}$'],4)
    plt.title('$\mathrm{Transverse\;Force},\;F_y$') 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.xlabel('$t\;\;(\mathrm{sec})$')
    plt.ylabel('$F\;\;(\mu\mathrm{N})$')


    # moment
    fig = plt.figure(num=3, facecolor='white',figsize=(12,5))
    ax = plt.axes()
    ax.plot(T[:N],-1.0*M_static[:N],T[:N],-1.0*M_dyn[:N])
    ax.legend(['$\mathrm{quasistatic}$','$\mathrm{dynamic}$'],4)
    plt.title('$\mathrm{Moment},\;M_z$')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.xlabel('$t\;\;(\mathrm{sec})$')
    plt.ylabel('$M\;\;(\mathrm{mN\mbox{-}mm})$')

    # discrete vs continuous
    fig = plt.figure(num=4, facecolor='white',figsize=(12,5))
    ax = plt.axes()
    ax.plot(T[:N],-1.0*M_dyn[:N],T[:N],-1.0*M_dyn_dis[:N])
    ax.legend(['$\mathrm{continuous}$','$\mathrm{discrete}$'],4)
    plt.title('$\mathrm{Moment\;(Continuous\;vs.\;Discrete)}$')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.xlabel('$t\;\;(\mathrm{sec})$')
    plt.ylabel('$M\;\;(\mathrm{mN\mbox{-}mm})$')
    
    plt.show()

def plot_stuff3():
    f = open('./Rat1206_N04_Tr09_3d/Rat1206_N04_Tr09_3d_forces.p')
    data = pickle.load(f)
    f.close()

    # get data.
    N = 600
    N = min(N,len(data['t']))
    scale = 1.0e6
    T = data['t']
    FZ_static = data['static']['FZ']*scale
    FY_static = data['static']['FY']*scale
    FX_static = data['static']['FX']*scale
    MZ_static = data['static']['MZ']*scale
    MY_static = data['static']['MY']*scale
    MX_static = data['static']['MX']*scale
    FZ_dyn = data['dyn']['FZ']*scale
    FY_dyn = data['dyn']['FY']*scale
    FX_dyn = data['dyn']['FX']*scale
    MZ_dyn = data['dyn']['MZ']*scale
    MY_dyn = data['dyn']['MY']*scale
    MX_dyn = data['dyn']['MX']*scale

    M_dyn = np.sqrt(MX_dyn**2+MY_dyn**2+MZ_dyn**2)
    M_stc = np.sqrt(MX_static**2+MY_static**2+MZ_static**2)

    F_dyn = np.sqrt(FX_dyn**2+FY_dyn**2+FZ_dyn**2)
    F_stc = np.sqrt(FX_static**2+FY_static**2+FZ_static**2)
    
    # matplotlib settings.
    plt.rc('text', usetex=True)
    fig = plt.figure(num=1, facecolor='white',figsize=(12,5))
    ax = plt.axes()
    ax.plot(T[:N],F_dyn[:N])

    ax2 = ax.twinx()
    ax2.plot(T[:N],M_dyn[:N])
    #ax.legend(['$\mathrm{continuous}$','$\mathrm{discrete}$'],4)
    plt.title('$\mathrm{Moment\;(Continuous\;vs.\;Discrete)}$')
    ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax2.get_yaxis().tick_right()
    ax.set_xlabel('$t\;\;(\mathrm{sec})$')
    ax2.set_ylabel('$\|M\|\;\;(\mathrm{mN\mbox{-}mm})$',color='b')
    ax.set_ylabel('$\|F\|\;\;(\mu\mathrm{N})$')
    plt.show()

def plot_stuff4():
    f = open('./Rat1206_N04_Tr09_3d/Rat1206_N04_Tr09_3d.p')
    data = pickle.load(f)
    f.close()

    X = data['x']
    Y = data['y']
    Z = data['z']

    plt.rc('text', usetex=True)
    fig = plt.figure()
    ax = Axes3D(fig)

    for i in xrange(0,100,5):
        ax.plot(X[i], Y[i], Z[i])

    ax.set_xlabel('$x\;\;(\mathrm{mm})$')
    ax.set_ylabel('$y\;\;(\mathrm{mm})$')
    ax.set_zlabel('$z\;\;(\mathrm{mm})$')
    
    plt.show()

def plot_stuff5():
    f = open('./Rat1206_N04_Tr09/Rat1206_N04_Tr09.p')
    data=pickle.load(f)
    f.close()
    ba = data['q'][2]
    plt.plot(ba)
    plt.show()


def plot_stuff6():

    f = open('./Rat1206_N04_Tr09_3d/Rat1206_N04_Tr09_3d_forces.p')
    data=pickle.load(f)
    f.close()

    md = np.sqrt(data['dyn']['MX']**2+data['dyn']['MY']**2+data['dyn']['MZ']**2)
    mq =\
    np.sqrt(data['static']['MX']**2+data['static']['MY']**2+data['static']['MZ']**2)
    fd = data['dyn']['FX']
    fq = data['static']['FX']

    t = data['t']
    plt.plot(t,mq,t,md)
    plt.show()

if __name__ == "__main__":
    plot_stuff6()






