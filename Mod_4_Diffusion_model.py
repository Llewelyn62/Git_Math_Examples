import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#Set parameters.
n = 192
Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.065 # Bacteria 1
#Du, Dv, F, k = 0.00016, 0.00008, 0.060, 0.062 # Coral
#Du, Dv, F, k = 0.00014, 0.00006, 0.035, 0.065 # Bacteria 2
#Du, Dv, F, k = 0.00016, 0.00008, 0.060, 0.062 # Coral
#Du, Dv, F, k = 0.00019, 0.00005, 0.060, 0.062 # Fingerprint
#Du, Dv, F, k = 0.00010, 0.00010, 0.018, 0.050 # Spirals
#Du, Dv, F, k = 0.00012, 0.00008, 0.020, 0.050 # Spirals Dense
#Du, Dv, F, k = 0.00010, 0.00016, 0.020, 0.050 # Spirals Fast
#Du, Dv, F, k = 0.00016, 0.00008, 0.020, 0.055 # Unstable
#Du, Dv, F, k = 0.00016, 0.00008, 0.050, 0.065 # Worms 1
#Du, Dv, F, k = 0.00016, 0.00008, 0.054, 0.063 # Worms 2
#Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.060 # Zebrafish
dh = 5./(n-1)
T = 8000
dt = dt = .9 * dh**2 / (4*max(Du,Dv))
nt = int(T/dt)
#Get the U and V data.
uvinitial = np.load('/Users/Llewelyn/Dropbox/Python code/Git_Math_Examples/uvinitial.npz')
U = uvinitial['U']
V = uvinitial['V']
#print np.shape((U,V))
#
#Create discretisation method.
def Forward_time_centre_space(Uin,Vin, nt, Du,Dv,F,k, dt, dh):
    ''' Un & Vn are updated through ftcs.
        nt: number of time steps.
        dt: time step.
        dh: array steps, dh=dx=dy. '''

    for j in xrange(nt):
        Un = Uin.copy()
        Vn = Vin.copy()
        #Forward time difference and centred spatial change across x and y(i,j).
        Uin[1:-1,1:-1] = Du*dt/dh**2*((Un[2:,1:-1]-2*Un[1:-1,1:-1]+Un[:-2,1:-1])+(Un[1:-1,2:]-2*Un[1:-1,1:-1]+Un[1:-1,:-2]))\
        -dt*Un[1:-1,1:-1]*Vn[1:-1,1:-1]**2 + dt*F*(1-Un[1:-1,1:-1]) + Un[1:-1,1:-1]
        # Enforce Neumann BCs: The border must be the gradient equal to a function q(x,y)=0.
        Uin[-1,:] = Uin[-2,:] #Last item of top = last-1
        Uin[0,:] = Uin[1,:] #Last item of bottom = last-1
        Uin[:,-1] = Uin[:,-2] #Last item of rhs = last-1
        Uin[:,0] = Uin[:,1] #Last item of lhs = last-1

        Vin[1:-1,1:-1] = Dv*dt/dh**2*((Vn[2:,1:-1]-2*Vn[1:-1,1:-1]+Vn[:-2,1:-1])+(Vn[1:-1,2:]-2*Vn[1:-1,1:-1]+Vn[1:-1,:-2]))\
        +dt*Un[1:-1,1:-1]*Vn[1:-1,1:-1]**2 - dt*(F+k)*Vn[1:-1,1:-1] + Vn[1:-1,1:-1]
        # Enforce Neumann BCs: The border must be the gradient equal to a function q(x,y)=0.
        Vin[-1,:] = Vin[-2,:] #Last item of top = last-1
        Vin[0,:] = Vin[1,:] #Last item of bottom = last-1
        Vin[:,-1] = Vin[:,-2] #Last item of rhs = last-1
        Vin[:,0] = Vin[:,1] #Last item of lhs = last-1
        #print j
    return Uin, Vin

#Create data
Unew, Vnew = Forward_time_centre_space(U, V,nt,Du,Dv,F,k,dt,dh)
print Unew[100,::40]
#Plot results
plt.close()
fig = plt.figure(figsize=(8,5))
plt.subplot(121)
plt.imshow(Unew, cmap = cm.RdBu)
plt.xticks([]), plt.yticks([]);
plt.subplot(122)
plt.imshow(Vnew, cmap = cm.RdBu)
plt.xticks([]), plt.yticks([]);
plt.show()
