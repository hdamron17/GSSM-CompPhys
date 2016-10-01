'''
Demonstrates numerical integration of the motion of linear damped driven 
    pendulums
Created on Sep 30, 2016
@author: Hunter Damron
'''

import numpy as np
import matplotlib.pyplot as plt
import math

def dw_dt(theta,ang_v, t, nat_freq, friction_coef, damping_freq, damping_torque):
    """
    Derivative of omega (angular velocity)
    :param theta: anglular position
    :param ang_v: angular velocity
    :param t: time from beginning of damping force oscillation 
    :param nat_freq: natural frequency of ideal pendulum
    :param friction_coeff: coefficient of friction force
    :param damping_freq: frequency of damping force
    :param damping_torque: torque due to damping force
    :param return: returns instantaneous acceleration
    """
    return -nat_freq**2 * theta - friction_coef * ang_v + \
                                    damping_torque * math.sin(damping_freq * t)

def shm_damped_driven(theta0, ang_v0, dt, end_t, nat_freq, friction_coef, 
                      damping_freq, damping_torque, m, plot_x_vs_time=False, 
                      plot_phase_space=False):
    """
    Calculates Leapfrog solution by Leapfrog algorithm
    :param theta0: initial position
    :param ang_v0: initial angular velocity
    :param dt: time delta (s)
    :param end_t: end time of calculation (s)
    :param nat_freq: natural frequency of ideal pendulum
    :param friction_coeff: coefficient of friction force
    :param damping_freq: frequency of damping force
    :param damping_torque: torque due to damping force
    :param m: pendulum mass
    :param plot_x_vs_t: if True, plots position vs time
    :param plot_phase_space: if True, plots momentum vs position
    :return: Returns tuple containing lists (time, position, angular velocity)
    """
    
    assert nat_freq > 0, "natural frequency must be positive\n" \
                                        + "(It's really boring at zero)"
    assert friction_coef >= 0, "coefficient of friction must be positive"
    assert damping_freq >= 0, "damping force frequency must be positive"
    assert damping_torque >= 0, "damping torque must be positive"
    assert dt > 0, "time delta must be positive"
    assert end_t > 0, "end time must be positive"
    
    t = np.arange(0, end_t+dt, dt)
    theta = [theta0]
    ang_v = [ang_v0]
    
    for time in t[1:]:
        ai = dw_dt(theta[-1], ang_v[-1], time, nat_freq, friction_coef, 
                                            damping_freq, damping_torque)
        theta.append(theta[-1] + ang_v[-1] * dt + 0.5 * ai * dt**2)
        af = dw_dt(theta[-1], ang_v[-1], time, nat_freq, friction_coef, 
                                            damping_freq, damping_torque)
        ang_v.append(ang_v[-1] + 0.5*(ai + af)*dt)
    
    if plot_x_vs_time:
        plt.figure()
        plt.plot(t, theta)
        plt.title("Damped Driven Pendulum Position vs Time")
        plt.ylabel("Position (radians)")
        plt.xlabel("Time (s)")
    
    if plot_phase_space:
        plt.figure()
        plt.plot(theta, np.multiply(ang_v, m))
        plt.title("Damped Driven Pendulum Phase Space Plot")
        plt.ylabel("Position (radians)")
        plt.xlabel("Time (s)")
        
    return t, theta, ang_v

if __name__ == '__main__':
    shm_damped_driven(
        theta0=1, ang_v0=0, dt=0.01, end_t=20, nat_freq=8, 
        friction_coef=0, damping_freq=0, damping_torque=0, m=3, 
        plot_x_vs_time=True, plot_phase_space=True)
    plt.show()
    pass