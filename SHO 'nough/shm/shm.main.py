"""
Project for Simple Harmonic Motion to show failure of Euler Method
Created September 19, 2016
Author: Hunter Damron
"""

import numpy as np
import matplotlib.pyplot as plt
import math

def dv_dt(k, m, x):
    """
    Calculates acceleration based on position x
    :param k: spring constant (kg/s^2)
    :param m: mass (kg)
    :param x: position (m)
    """
    return -(k / m)*x

def shm_exact(x0, k, m, dt, end_t, plot_x_vs_time=False, 
                                                plot_phase_space=False):
    """
    Calculates SHM exact solution (to be plotted over numerical solutions)
    :param x0: initial position
    :param k: spring constant
    :param m: mass (kg)
    :param dt: time delta (s)
    :param end_t: end time of calculation (s)
    :param plot_x_vs_t: if True, plots x vs time
    :param plot_phase_space: if True, plots momentum vs position
    :return: Returns tuple containing lists (time, position, velocity)
    """
    omega = math.sqrt(k / m)
    t = np.arange(0, end_t+dt, dt)
    x = x0 * np.cos(omega * t)
    v = -x0 * omega * np.sin(omega * t)
    
    if plot_x_vs_time:
        plt.figure()
        plt.plot(t, x)
        plt.title("SHM Exact Solution Position vs Time")
        plt.ylabel("Position (m)")
        plt.xlabel("Time (s)")
    
    if plot_phase_space:
        plt.figure()
        plt.plot(x, np.multiply(v, m))
        plt.title("SHM Exact Solution Phase Space Plot")
        plt.ylabel("Position (m)")
        plt.xlabel("Time (s)")
        
    return t, x, v

def shm(x0, k, m, dt, end_t=100, plot_x_vs_time=False, plot_phase_space=False,
                                           use_euler=False, plot_exact=False):
    """
    Calculates SHM using euler method (should fail)
    :param x0: initial position (m)
    :param k: spring constant (kg/s^2)
    :param m: mass (kg)
    :param dt: time delta between calculations (s)
    :param end_t: end time of calculation (s)
    :param plot_x_vs_t: if True, plots x vs time
    :param plot_phase_space: if True, plots momentum vs position
    :param use_euler: If True, use Euler Method instead of Euler Cromer Method
    :param plot_exact: If True, plots exact solution
    :return: Returns tuple containing lists (time, position, velocity)
    """
    assert k > 0, "k must be positive"
    assert m > 0, "mass must be positive"
    assert dt > 0, "time delta must be positive"
    assert end_t > 0, "End time must be positive"
    t = np.arange(0, end_t+dt, dt)
    x=[x0]
    v=[0]
    x_index = -2 if use_euler else -1
    cromer_string = "" if use_euler else "Cromer "
    for time in t[1:]:
        v.append(v[-1] + dv_dt(k, m, x[-1]) * dt)
        x.append(x[-1] + v[x_index] * dt)
    
    if plot_exact:
        t_exact, x_exact, v_exact = shm_exact(x0, k, m, dt, end_t)
    
    if plot_x_vs_time:
        plt.figure()
        plt.plot(t, x, 'k-', label="Euler %sMethod"% cromer_string)
        if plot_exact:
            plt.plot(t_exact, x_exact, 'r-', label="Exact Solution")
            plt.legend()
        plt.title("Position vs Time for SHM by Euler %sMethod"% cromer_string)
        plt.ylabel("Position (m)")
        plt.xlabel("Time (s)")
        
    if plot_phase_space:
        plt.figure()
        plt.plot(x, np.multiply(v, m), 'k-', label="Euler %sMethod"% cromer_string)
        if plot_exact:
            plt.plot(x_exact, np.multiply(v_exact, m), 'r-', 
                                                       label="Exact Solution")
            plt.legend()
        plt.title("Phase Space Plot for SHM by Euler %sMethod"% cromer_string)
        plt.ylabel("Momentum (kg*m/s)")
        plt.xlabel("Position (m)")
    
    return t, x, v

def shm_various_dt(x0, k, m, dt_log_range=(-5, -1), end_t=20, 
                        plot_x_vs_time=False, plot_phase_space=False, 
                        use_euler=False, plot_exact=False):
    """
    Calculates SHM using euler method (should fail)
    :param x0: initial position (m)
    :param k: spring constant (kg/s^2)
    :param m: mass (kg)
    :param dt_log_range: two value tuple containing logarithmic endpoints
                                        ( log(largest dt), log(smallest dt) )
    :param end_t: end time of calculation (s)
    :param plot_x_vs_t: if True, plots x vs time
    :param plot_phase_space: if True, plots momentum vs position
    :param use_euler: If True, use Euler Method instead of Euler Cromer Method
    :param plot_exact: If True, plots exact solution
    :return: Returns tuple containing lists (time, position, velocity)
    """
    assert k > 0, "k must be positive"
    assert m > 0, "mass must be positive"
    assert dt_log_range[1] > dt_log_range[0], "end of dt range must be \
                                            larger than end of dt range"
    assert dt_log_range[1] <= end_t / 5, "largest dt must be less than \
                                            one fifth of the end time"
    assert end_t > 0, "End time must be positive"
    
    t = []
    x = []
    v = []
    colors = ['k', 'b', 'c', 'g', 'y', 'r', 'm']
    dt = np.power(10, np.linspace(dt_log_range[0], dt_log_range[1], 
                                                len(colors), endpoint=True))
    for i in dt:
        t_new, x_new, v_new = shm(x0, k, m, i, end_t, False, False, use_euler)
        t.append(t_new)
        x.append(x_new)
        v.append(v_new)
    
    cromer_string = "" if use_euler else "Cromer "
    
    if plot_exact:
        t_exact, x_exact, v_exact = shm_exact(x0, k, m, dt[0], end_t)
    
    if plot_x_vs_time:
        plt.figure()
        for tp, xp, dtp, color in zip(t, x, dt, colors):
            plt.plot(tp, xp, color + ':', label="dt = %g" % dtp)
        if plot_exact:
            plt.plot(t_exact, x_exact, 'k-', label="Exact Solution")
        plt.title("Position vs Time for SHM by Euler %sMethod"% cromer_string)
        plt.ylabel("Position (m)")
        plt.xlabel("Time (s)")
        plt.legend()
        
    if plot_phase_space:
        plt.figure()
        for xp, vp, dtp, color in zip(x, v, dt, colors):
            plt.plot(xp, np.multiply(vp, m), color + ':', label="dt = %g"%dtp)
        if plot_exact:
            plt.plot(x_exact, np.multiply(v_exact, m), 'k-', 
                                                       label="Exact Solution")
        plt.title("Phase Space Plot for SHM by Euler %sMethod"% cromer_string)
        plt.ylabel("Momentum (kg*m/s)")
        plt.xlabel("Position (m)")
        plt.legend()

if __name__ == "__main__":
    shm(-10, 6, 0.7, 0.01, 30, plot_x_vs_time=True, plot_phase_space=True, 
                    use_euler=True, plot_exact=True)
    shm(-10, 6, 0.7, 0.01, 30, plot_x_vs_time=True, plot_phase_space=True, 
                    use_euler=False, plot_exact=True)
    
    shm_various_dt(-10, 6, 0.7, (-5, -1), 5, plot_x_vs_time=True, 
                    plot_phase_space=True, use_euler=True, plot_exact=True)    
    shm_various_dt(-10, 6, 0.7, (-5, -1), 5, plot_x_vs_time=True, 
                    plot_phase_space=True, use_euler=False, plot_exact=True)
    
    plt.show()