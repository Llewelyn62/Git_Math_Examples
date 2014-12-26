import numpy
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def ftcs(T, nt, alpha, dt, dx, dy):
    ''' T is the 2D array, which is updated through ftcs.
        nt: number of time steps.
        alpha: diffusivity value = k/c_p*rho: rho is density, c_p specific heat & k thermal constant.
        dt: time step.
        dx, dy: array steps, can be set equal. '''

    #Takes shape of matrix and finds midpoint.
    j_mid = (numpy.shape(T)[0])/2
    i_mid = (numpy.shape(T)[1])/2

    for n in range(nt):
        Tn = T.copy()
        #Forward time difference and centred spatial change across x and y(i,j).
        #Note that the y-value is the FIRST (row reference), x the col. to preserve orientation.
        T[1:-1,1:-1] = Tn[1:-1,1:-1] + alpha *\
            (dt/dy**2 * (Tn[2:,1:-1] - 2*Tn[1:-1,1:-1] + Tn[:-2,1:-1]) +\
             dt/dx**2 * (Tn[1:-1,2:] - 2*Tn[1:-1,1:-1] + Tn[1:-1,:-2]))

        # Enforce Neumann BCs: The last must be the gradient equal to a function q(x,y).
        T[-1,:] = T[-2,:] #-dy*5500 ## Gradient dyQ(x,y) must be negative. Stasis around flux 5500.
        T[:,-1] = T[:,-2] #-dx*5500

        # Check if we reached T=70C
        if T[j_mid, i_mid] >= 70:
            print ("Center of plate reached 70C at time {0:.2f}s.".format(dt*n))
            break

    if T[j_mid, i_mid]<70:
        print ("Center has not reached 70C yet, it is only {0:.2f}C.".format(T[j_mid, i_mid]))

    return T
#Initial values and parameters.
L = 1.0e-2
H = 1.0e-2
#Steps for each dimension.
nx = 21
ny = 21
nt = 500
#Size of each step for spatial co-ordinates.
dx = L/(nx-1)
dy = H/(ny-1)
#Set two arrays to later combine into a mesh.
x = numpy.linspace(0,L,nx)
y = numpy.linspace(0,H,ny)
#Diffusivity value.
alpha = 1e-4
#Create initial array. The 20 reverses calculation of dx, dy.
Ti = numpy.ones((ny, nx))*20
Ti[0,:]= 100
Ti[:,0] = 100
#Sigma is the value to maintain stability of the system across time.
sigma = 0.25
dt = sigma * min(dx, dy)**2 / alpha

#Set up the problem and create the array for the time taken to reach 70C.
T=Ti.copy()
T = ftcs(T, nt, alpha, dt, dx, dy)

mx, my = numpy.meshgrid(x,y)

##Now plot the figure.
#plt.figure(figsize=(8,5))
#plt.contourf(my,mx,T,20)
#plt.xlabel('$x$')
#plt.ylabel('$y$')
#plt.colorbar();
#plt.show()

from mayavi.mlab import *
imshow(T, colormap='spectral')
scalarbar()
show()
