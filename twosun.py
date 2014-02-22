from scipy.integrate import odeint  # for integrate.odeint
from pylab import *  # for plotting commands
#1 is earth 2 is jup
GM = (4*pi**2)
#same order for x input into planets
xi = 1.1
yi = 1
vxi = sqrt(GM/1)
vyi = sqrt(GM/1)
#
 
def planets(x,t):
    dxdt = zeros(size(x))
    r1 = sqrt((x[0] - 0)**2 + (x[1] - 0)**2)
    r2 = sqrt((x[0] - 2)**2 + (x[1] - 0)**2)
    dxdt[0] = x[2]
    dxdt[1] = x[3]
    dxdt[2] = -(GM/r1**3)*x[0] - (GM/r2**3)*(x[0] - 2)
    dxdt[3] = -(GM/r1**3)*x[1] - (GM/r2**3)*(x[1])
    return dxdt
 
#initial conditions
state1 = array([xi,yi,vxi,vyi])
state2 = array([xi,yi,sqrt(GM/1.1),sqrt(GM/1)])
state3 = array([xi,yi,sqrt(GM/1),sqrt(GM/1.1)])
state4 = array([xi,yi,sqrt(GM/0.9),sqrt(GM/1)])
 
times = linspace(0,100,10000)
 
#python ODE solver solution
states1 = odeint(planets,state1,times)
states2 = odeint(planets,state2,times)
states3 = odeint(planets,state3,times)
states4 = odeint(planets,state4,times)
 
E_x1 = states1[:,0]
E_y1 = states1[:,1]
E_x2 = states2[:,0]
E_y2 = states2[:,1]
E_x3 = states3[:,0]
E_y3 = states3[:,1]
E_x4 = states4[:,0]
E_y4 = states4[:,1]
 
#plotting
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax1.plot(E_x1,E_y1,'or')
plt.xlabel('Distances (A.U.)')
plt.ylabel('Distances (A.U.)')
legend(["Vxi = sqrt(GM), Vyi = sqrt(GM)"])
xlim(-7,7)
ylim(-7,7)
ax2 = fig.add_subplot(222)
ax2.plot(E_x2,E_y2,'or')
plt.xlabel('Distances (A.U.)')
plt.ylabel('Distances (A.U.)')
legend(["Vxi = sqrt(GM/1.1), Vyi = sqrt(GM)"])
xlim(-7,7)
ylim(-7,7)
ax3 = fig.add_subplot(223)
ax3.plot(E_x3,E_y3,'or')
plt.xlabel('Distances (A.U.)')
plt.ylabel('Distances (A.U.)')
legend(["Vxi = sqrt(GM), Vyi = sqrt(GM/1.1)"])
xlim(-7,7)
ylim(-7,7)
ax4 = fig.add_subplot(224)
ax4.plot(E_x4,E_y4,'or')
plt.xlabel('Distances (A.U.)')
plt.ylabel('Distances (A.U.)')
legend(["Vxi = sqrt(GM/0.9), Vyi = sqrt(GM)"])
xlim(-7,7)
ylim(-7,7)
plt.show()