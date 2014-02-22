from scipy.integrate import odeint  # for integrate.odeint
from pylab import *  # for plotting commands
#1 is earth 2 is jup
GM = (4*pi**2)
#same order for x input into planets
x1 = 2.52
y1 = 0
vx1 = 0
vy1 = sqrt(GM/2.52)
x2 = 5.24
y2 = 0
vx2 = 0
vy2 = sqrt(GM/5.24)
#
#for python ODE
def planets(x,t):
    dxdt = zeros(size(x))
    r1 = sqrt(x[0]**2 + x[1]**2)
    r2 = sqrt(x[4]**2 + x[5]**2)
    r21 = sqrt((x[4] - x[0])**2 + (x[5] - x[1])**2)
    dxdt[0] = x[2]
    dxdt[1] = x[3]
    dxdt[2] = -(GM/r1**3)*x[0] + (0.04*GM/r21**3)*(x[4] - x[0])
    dxdt[3] = -(GM/r1**3)*x[1] + (0.04*GM/r21**3)*(x[5] - x[1])
    dxdt[4] = x[6]
    dxdt[5] = x[7]
    dxdt[6] = -(GM/r2**3)*x[4] - (0.001*GM/r21**3)*(x[4] - x[0])
    dxdt[7] = -(GM/r2**3)*x[5] - (0.001*GM/r21**3)*(x[5] - x[1])
    return dxdt
#for my rungkutta because x and t are switched
def planet(t,x):
    dxdt = zeros(size(x))
    r1 = sqrt(x[0]**2 + x[1]**2)
    r2 = sqrt(x[4]**2 + x[5]**2)
    r21 = sqrt((x[4] - x[0])**2 + (x[5] - x[1])**2)
    dxdt[0] = x[2]
    dxdt[1] = x[3]
    dxdt[2] = -(GM/r1**3)*x[0] + (0.04*GM/r21**3)*(x[4] - x[0])
    dxdt[3] = -(GM/r1**3)*x[1] + (0.04*GM/r21**3)*(x[5] - x[1])
    dxdt[4] = x[6]
    dxdt[5] = x[7]
    dxdt[6] = -(GM/r2**3)*x[4] - (0.001*GM/r21**3)*(x[4] - x[0])
    dxdt[7] = -(GM/r2**3)*x[5] - (0.001*GM/r21**3)*(x[5] - x[1])
    return dxdt
 
#ODE int
def RungKutta(x,f,t,dt):
    RK1 = f(t,x) * dt
    RK2 = f(t + dt/2, x + RK1/2) * dt
    RK3 = f(t + dt/2, x + RK2/2) * dt
    RK4 = f(t+ dt, x + RK3) * dt    
    return x + (1.0/6.0)*(RK1 + 2*RK2 + 2*RK3 + RK4)
 
#initial conditions
state = array([x1,y1,vx1,vy1,x2,y2,vx2,vy2])
times = linspace(0,100,1000)
 
#Rungkutta solution
rk_sol = [state]
dt = 0.1
rk_time = [0]
while rk_time[-1] < 100:
    rk_sol.append(RungKutta(rk_sol[-1],planet,rk_time[-1],dt))
    rk_time.append(rk_time[-1] + dt)
 
#python ODE solver solution
states = odeint(planets,state,times)
 
E_x = states[:,0]
E_y = states[:,1]
J_x = states[:,4]
J_y = states[:,5]
 
#energy conservation
def Energy(x,t):
    r1 = sqrt(x[0]**2 + x[1]**2)
    r2 = sqrt(x[4]**2 + x[5]**2)
    r21 = sqrt((x[4] - x[0])**2 + (x[5] - x[1])**2)
    v1sq = (x[2]**2 + x[3]**2)
    v2sq = (x[6]**2 + x[7]**2)
    return (0.5*0.001*v1sq) + (0.5*0.04*v2sq) - (0.001*GM/r1) - (0.04*GM/r2) - ((0.001*0.04*GM)/r21)
#angular momentum
def Momentum(x,t):
    return 0.001*(x[0]*x[3] - x[1]*x[2]) + 0.04*(x[4]*x[7] - x[5]*x[6])
 
EoverM = []
for i in range(len(times)):
    EoverM.append(Energy(states[i],times[i]))
 
LoverM = []
for i in range(len(times)):
    LoverM.append(Momentum(states[i],times[i]))
 
#plotting    
fig = plt.figure()
ax1 = fig.add_subplot(311)
ax1.plot(E_x,E_y,'or')
ax1.plot(J_x,J_y,'ob')
plt.xlabel('Distances (A.U.)')
plt.ylabel('Distances (A.U.)')
legend(["Earth", "Jupiter"])
ax2 = fig.add_subplot(312)
ax2.plot(times,EoverM,'ob')
plt.xlabel('Time (yrs)')
plt.ylabel('Energy/M')
ax3 = fig.add_subplot(313)
ax3.plot(times,LoverM,'ob')
plt.xlabel('Time (yrs)')
plt.ylabel('Angular Momentum/M')
plt.show()