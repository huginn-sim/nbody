# -*- coding: utf-8 -*-
"""
.. module:: nbody
   :synopsis:   Simulates the motion of two bodies
                orbiting a fixed star.

.. moduleauthor:: Evin Özer
"""

#~ Modules
from pylab import *
from scipy.integrate import ode, odeint
import matplotlib.image as mpimg
from PIL import Image
#/~ Modules

#~ Custom Modules
from viz import configure
#/~ Custom Modules

#~ ODE Solvers
def rk4(f,t,dt,x):
    k1 = f(t        , x        )*dt
    k2 = f(t + dt/2., x + k1/2.)*dt
    k3 = f(t + dt/2., x + k2/2.)*dt
    k4 = f(t + dt   , x + k3   )*dt

    return x + (1./6.)*(k1 + 2*k2 + 2*k3 + k4)
#/~ ODE Solvers

def energy(x):
    r1 = np.linalg.norm(array(x[:2]))
    r2 = np.linalg.norm(array(x[4:6]))
    r21 = np.linalg.norm(array(x[4:6])-array(x[:2]))
    v1 = np.linalg.norm(array(x[2:4]))
    v2 = np.linalg.norm(array(x[6:]))
    
    return .5*m1*v1**2 + .5*m2*v2**2 - (GM*m1)/r1 - (GM*m2)/r2 - (GM*m2*m1)/r21

def angular_momentum(x):
    return m1*(x[0]*x[3] - x[1]*x[2]) + m2*(x[4]*x[7] - x[5]*x[6])

def _move(t, y, GM, Gm1, Gm2):
    r1 = array(y[:2]); r1_mag = np.linalg.norm(r1)**3
    rv1 = array(y[2:4])

    r2 = array(y[4:6]); r2_mag = np.linalg.norm(r2)**3
    rv2 = array(y[6:])

    r21 = r2 - r1
    r21_mag = sqrt((r2[0] - r1[0])**2 + (r2[1] - r1[1])**2)**3

    r1_next = rv1
    rv1_next = -((GM)/(r1_mag))*r1 + ((Gm2)/(r21_mag))*r21

    r2_next = rv2
    rv2_next = -((GM)/(r2_mag))*r2 - ((Gm1)/(r21_mag))*r21

    state = array([ r1_next[0], r1_next[1], rv1_next[0], rv1_next[1],
                    r2_next[0], r2_next[1], rv2_next[0], rv2_next[1]])
    return state

def move_(x, t):
    r1 = array(x[:2]); r1_mag = np.linalg.norm(r1)**3
    rv1 = array(x[2:4])

    r2 = array(x[4:6]); r2_mag = np.linalg.norm(r2)**3
    rv2 = array(x[6:])

    r21 = r2 - r1
    r21_mag = sqrt((r2[0] - r1[0])**2 + (r2[1] - r1[1])**2)**3

    r1_next = rv1
    rv1_next = -((GM)/(r1_mag))*r1 + ((Gm2)/(r21_mag))*r21

    r2_next = rv2
    rv2_next = -((GM)/(r2_mag))*r2 - ((Gm1)/(r21_mag))*r21

    state = array([ r1_next[0], r1_next[1], rv1_next[0], rv1_next[1],
                    r2_next[0], r2_next[1], rv2_next[0], rv2_next[1]])
    return state

def move(t, x):
    r1 = array(x[:2]); r1_mag = np.linalg.norm(r1)**3
    rv1 = array(x[2:4])

    r2 = array(x[4:6]); r2_mag = np.linalg.norm(r2)**3
    rv2 = array(x[6:])

    r21 = r2 - r1
    r21_mag = sqrt((r2[0] - r1[0])**2 + (r2[1] - r1[1])**2)**3

    r1_next = rv1
    rv1_next = -((GM)/(r1_mag))*r1 + ((Gm2)/(r21_mag))*r21

    r2_next = rv2
    rv2_next = -((GM)/(r2_mag))*r2 - ((Gm1)/(r21_mag))*r21

    state = array([ r1_next[0], r1_next[1], rv1_next[0], rv1_next[1],
                    r2_next[0], r2_next[1], rv2_next[0], rv2_next[1]])
    return state

def plot_move_rk4():
    states = [array([x1, y1, vx1, vy1, x2, y2, vx2, vy2])]
    for t in times[1:]:
        states.append(rk4(move, t, dt, states[-1]))

    states = array(states)

    earth = [states[:,0], states[:,1]]
    jupiter = [states[:,4], states[:,5]]
    sun = mpimg.imread('sun.png')
    
    fig, ax = subplots()

    ax.imshow(sun, aspect='auto', extent=(-.1, .1, -.1, .1), zorder=1)
    e_mark0, = ax.plot(earth[0], earth[1], 'b-', lw=2, alpha=.4)
    e_mark1, = ax.plot(earth[0], earth[1], 'bo', ms=5, alpha=.4)

    j_mark0, = ax.plot(jupiter[0], jupiter[1], 'r-', lw=2, alpha=.4)
    j_mark1, = ax.plot(jupiter[0], jupiter[1], 'ro', ms=7, alpha=.4)

    legend( [(e_mark0, e_mark1), (j_mark0, j_mark1)],
            ['Earth Orbit:\t$m_{\mathcal{E}}='+str(m1)+'$',
             'Jupiter Orbit:\t$m_{\mathcal{J}}='+str(m2)+'$'],
             numpoints=1)

    configure(  ax=ax,
                title=r'$'+str(tf-t0)+'$ Years; $\Delta t='+str(dt)+'$',
                xlabel=r'Displacement (AU)',
                ylabel=r'Displacement (AU)',
                xbounds=None, ybounds=None)

    fig.suptitle('RK4 - Orbit of Two Planets About a Fixed Star', size=30)
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.08)

    show()

def plot_move_ode():
    states = [array([x1, y1, vx1, vy1, x2, y2, vx2, vy2])]

    i = ode(_move)
    i.set_integrator('vode', method='adams')
    i.set_initial_value(states[-1], t0)
    i.set_f_params(GM, Gm1, Gm2)
    for t in times[1:]:
        i.integrate(i.t+dt)
        states.append(i.y)

    states = array(states)

    earth = [states[:,0], states[:,1]]
    jupiter = [states[:,4], states[:,5]]
    sun = mpimg.imread('sun.png')
    
    fig, ax = subplots()

    ax.imshow(sun, aspect='auto', extent=(-.1, .1, -.1, .1), zorder=1)
    e_mark0, = ax.plot(earth[0], earth[1], 'b-', lw=2, alpha=.4)
    e_mark1, = ax.plot(earth[0], earth[1], 'bo', ms=5, alpha=.4)

    j_mark0, = ax.plot(jupiter[0], jupiter[1], 'r-', lw=2, alpha=.4)
    j_mark1, = ax.plot(jupiter[0], jupiter[1], 'ro', ms=7, alpha=.4)

    legend( [(e_mark0, e_mark1), (j_mark0, j_mark1)],
            ['Earth Orbit:\t$m_{\mathcal{E}}='+str(m1)+'$',
             'Jupiter Orbit:\t$m_{\mathcal{J}}='+str(m2)+'$'],
             numpoints=1)

    configure(  ax=ax,
                title=r'$'+str(tf-t0)+'$ Years; $\Delta t='+str(dt)+'$',
                xlabel=r'Displacement (AU)',
                ylabel=r'Displacement (AU)',
                xbounds=None, ybounds=None)

    fig.suptitle('VODE - Orbit of Two Planets About a Fixed Star', size=30)
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.08)

    show()

def plot_move_odeint():
    states = [array([x1, y1, vx1, vy1, x2, y2, vx2, vy2])]

    states = odeint(move_, states[-1], times)
    states = array(states)

    earth = [states[:,0], states[:,1]]
    jupiter = [states[:,4], states[:,5]]
    sun = mpimg.imread('sun.png')
    
    fig, ax = subplots()

    ax.imshow(sun, aspect='auto', extent=(-.1, .1, -.1, .1), zorder=1)
    e_mark0, = ax.plot(earth[0], earth[1], 'b-', lw=2, alpha=.4)
    e_mark1, = ax.plot(earth[0], earth[1], 'bo', ms=5, alpha=.4)

    j_mark0, = ax.plot(jupiter[0], jupiter[1], 'r-', lw=2, alpha=.4)
    j_mark1, = ax.plot(jupiter[0], jupiter[1], 'ro', ms=7, alpha=.4)

    legend( [(e_mark0, e_mark1), (j_mark0, j_mark1)],
            ['Earth Orbit:\t$m_{\mathcal{E}}='+str(m1)+'$',
             'Jupiter Orbit:\t$m_{\mathcal{J}}='+str(m2)+'$'],
             numpoints=1)

    configure(  ax=ax,
                title=r'$'+str(tf-t0)+'$ Years; $\Delta t='+str(dt)+'$',
                xlabel=r'Displacement (AU)',
                ylabel=r'Displacement (AU)',
                xbounds=None, ybounds=None)

    fig.suptitle('VODE - Orbit of Two Planets About a Fixed Star', size=30)
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.08)

    show()

def plot_conservation():
    fig, ax = subplots()
    tax = ax.twinx()

    atols = [.00001, .000001, .0000001]
    y0 = array([x1, y1, vx1, vy1, x2, y2, vx2, vy2])
    for i,atol in enumerate(atols):
        states = odeint(move_, y0, times, atol=atol)

        enrgy = []
        for state in states:
            enrgy.append(energy(state))

        am = []
        for state in states:
            am.append(angular_momentum(state))

        alpha = (i+1.) / (len(atols)+1.)
        e_mark, = ax.plot(times, enrgy, 'b--', lw=3, alpha=alpha)
        am_mark, = tax.plot(times, am, 'r--', lw=3, alpha=alpha)

        legend( [e_mark, am_mark],
                [r'$E = \frac{1}{2}m_1v_1^2 + \frac{1}{2}m_2v_2^2 - \frac{GMm_1}{r_1} - \frac{GMm_2}{r_2} - \frac{Gm_2m1}{r_21}$',
                 r'$L = m_1 \left( x_1 v_{y1} - y_1 v_{x1} \right) + m_2 \left( x_2 v_{y2} - y_2 v_{x2} \right)$'],
                 numpoints=1, loc='upper right')

    configure(  ax=ax,
                    title=r'Absolute Tolerances: '+str(atols),
                    xlabel=r'Time (years)',
                    ylabel=r'Energy (Joules)',
                    colors=('k', 'b'))

    configure(  ax=tax,
                    ylabel=r'Angular Momentum $\left(kg\frac{m^2}{s}\right)$',
                    xbounds=(t0,tf+dt), colors=('k', 'r'), suppress=True)

    fig.suptitle('Conservation of Energy and Angular Momentum', size=30)
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.08)

    show()

#~ Entry point of the script.
if __name__ == "__main__":
    #~ State Variables
    GM = 4*math.pi**2

    x1 = 2.52; y1 = 0.0; vx1 = 0.0; vy1 = sqrt(GM/2.52)
    x2 = 5.24; y2 = 0.0; vx2 = 0.0; vy2 = sqrt(GM/5.24)

    m1 = .001
    m2 = .04

    Gm1 = m1 * GM
    Gm2 = m2 * GM

    t0 = 0; tf = 100; dt=.01;
    times = np.arange(t0, tf+dt, dt)
    #/~ State Variables
    
    plot_move_rk4()
    plot_move_ode()
    plot_move_odeint()
    plot_conservation()