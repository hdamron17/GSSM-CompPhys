'''
Created on Sep 23, 2016

@author: hdamron1594
'''

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

def shm_exact(x0, v0, k, m, t):
    """
    Calculates SHM exact solution (to be plotted over numerical solutions)
    :param x0: initial position (m)
    :param v0: initial velocity (m/s)
    :param k: spring constant (kg/s^2)
    :param m: mass (kg)
    :param t: point in time (s)
    :return: Returns tuple containing (position, velocity) at time
    """
    A = math.sqrt((m * v0**2 + k * x0**2)/k)
    omega = math.sqrt(k / m)
    phi = math.atan2(-v0, (omega * x0))
    x = A * np.cos(omega * t + phi)
    v = -omega * A * np.sin(omega * t + phi)
    return x, v
    
def shm_exact_plot(x0, v0, k, m, dt, end_t, plot_x_vs_time=False, 
                                                plot_phase_space=False):
    """
    Calculates SHM exact solution (to be plotted over numerical solutions)
    :param x0: initial position (m)
    :param v0: initial velocity (m/s)
    :param k: spring constant (kg/s^2)
    :param m: mass (kg)
    :param dt: time delta (s)
    :param end_t: end time of calculation (s)
    :param plot_x_vs_t: if True, plots x vs time
    :param plot_phase_space: if True, plots momentum vs position
    :return: Returns tuple containing lists (time, position, velocity)
    """
    A = math.sqrt((m * v0**2 + k * x0**2)/k)
    omega = math.sqrt(k / m)
    t = np.arange(0, end_t+dt, dt)
    phi = math.atan2(-v0, (omega * x0))
    x = A * np.cos(omega * t + phi)
    v = -omega * A * np.sin(omega * t + phi)
    
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

def shm_euler_cromer(x0, v0, k, m, dt, end_t=100, plot_x_vs_time=False, 
                   plot_phase_space=False, use_euler=False, plot_exact=False):
    """
    Calculates SHM using Euler Cromer method (or Euler method)
    :param x0: initial position (m)
    :param v0: initial velocity (m/s)
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
    v=[v0]
    x_index = -2 if use_euler else -1
    cromer_string = "" if use_euler else "Cromer "
    for time in t[1:]:
        v.append(v[-1] + dv_dt(k, m, x[-1]) * dt)
        x.append(x[-1] + v[x_index] * dt)
    
    if plot_exact:
        t_exact, x_exact, v_exact = shm_exact_plot(x0,v0, k, m, dt, end_t)
    
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
        plt.plot(x, np.multiply(v, m), 'k-', label="Euler %sMethod"
                                                            % cromer_string)
        if plot_exact:
            plt.plot(x_exact, np.multiply(v_exact, m), 'r-', 
                                                       label="Exact Solution")
            plt.legend()
        plt.title("Phase Space Plot for SHM by Euler %sMethod"% cromer_string)
        plt.ylabel("Momentum (kg*m/s)")
        plt.xlabel("Position (m)")
    
    return t, x, v

def shm_leapfrog(x0, v0, k, m, dt, end_t=100, plot_x_vs_time=False,
                                    plot_phase_space=False, plot_exact=False):
    """
    Calculates SHM using leapfrog algorithm 
    :param x0: initial position (m)
    :param v0: initial velocity (m/s)
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
    v=[v0]
    for time in t[1:]:
        ai = dv_dt(k, m, x[-1])
        x.append(x[-1] + v[-1] * dt + 0.5 * ai * dt**2)
        v.append(v[-1] + 0.5*(ai + dv_dt(k, m, x[-1]))*dt)
    
    if plot_exact:
        t_exact, x_exact, v_exact = shm_exact_plot(x0, v0, k, m, dt, end_t)
    
    if plot_x_vs_time:
        plt.figure()
        plt.plot(t, x, 'k-', label="Leapfrog Algorithm")
        if plot_exact:
            plt.plot(t_exact, x_exact, 'r-', label="Exact Solution")
            plt.legend()
        plt.title("Position vs Time for SHM by Leapfrog Algorithm")
        plt.ylabel("Position (m)")
        plt.xlabel("Time (s)")
        
    if plot_phase_space:
        plt.figure()
        plt.plot(x, np.multiply(v, m), 'k-', label="Leapfrog Algorithm")
        if plot_exact:
            plt.plot(x_exact, np.multiply(v_exact, m), 'r-', 
                                                       label="Exact Solution")
            plt.legend()
        plt.title("Phase Space Plot for SHM by Leapfrog Algorithm")
        plt.ylabel("Momentum (kg*m/s)")
        plt.xlabel("Position (m)")
    
    return t, x, v

def shm_error(x0, v0, k, m, dt_power_range=(-5, -1), num_dt=100, end_t=100, 
                                    plot=False, plot_lin=False, loglog=True):
    """
    Calculates SHM using leapfrog algorithm 
    :param x0: initial position (m)
    :param v0: initial velocity (m/s)
    :param k: spring constant (kg/s^2)
    :param m: mass (kg)
    :param dt_power_range: two point tuple with powers to which to raise dt
    :param num_dt: number of dt points used in dt_power_range
    :param end_t: end time of calculation (s)
    :param plot: if True, plots a graph of error vs dt
    :param plot_lin: if True, plots a graph of linear regression lines
    :param loglog: if True, plots a loglog plot, else plots standard plot
    :return: Returns tuple with lists (dt, euler cromer error, leapfrog error)
    """
    dt = np.power(10, np.linspace(dt_power_range[0], dt_power_range[1],
                                                    num_dt, endpoint=True))
    euler_cromer_error = []
    leapfrog_error = []
    xf_exact = shm_exact(x0, v0, k, m, end_t)[0]
    for delta_t in dt:
        euler_cromer_xf =shm_euler_cromer(x0, v0, k, m, delta_t, end_t)[1][-1]
        euler_cromer_err = abs((euler_cromer_xf - xf_exact) / xf_exact)
        euler_cromer_error.append(euler_cromer_err)
        
        leapfrog_xf = shm_leapfrog(x0, v0, k, m, delta_t, end_t)[1][-1]
        leapfrog_err = abs((leapfrog_xf - xf_exact) / xf_exact)
        leapfrog_error.append(leapfrog_err)
    
    if plot:
        plt.figure()
        plot = plt.loglog if loglog else plt.plot
        if plot_lin:
            euler_cromer_lin_reg = lin_reg(np.log10(dt), 
                                                np.log10(euler_cromer_error))
            leapfrog_lin_reg = lin_reg(np.log10(dt), np.log10(leapfrog_error))
            plot(dt, np.power(dt, euler_cromer_lin_reg[0]) 
                                * np.power(10, euler_cromer_lin_reg[1]), 'r-')
            plot(dt, np.power(dt, leapfrog_lin_reg[0])
                                * np.power(10, leapfrog_lin_reg[1]), 'b-')
            print("Euler Cromer : y = %gx + %g on loglog plot" 
                                                    % (euler_cromer_lin_reg))
            print("Leapfrog : y = %gx + %g on loglog plot"
                                                         % (leapfrog_lin_reg))
        plot(dt, euler_cromer_error, 'r.', label="Euler Cromer Method")
        plot(dt, leapfrog_error, 'b.', label="Leapfrog Method")
        plt.title("Global Error for SHM after %g seconds" % end_t)
        plt.ylabel("Global Error")
        plt.xlabel("Time Step (s)")
        plt.legend(loc='lower right')
    
    
def lin_reg(x, y):
    """
    Calculates equation for linear regression using method of least squares
    :param x: x values of data set
    :param y: y values of data set
    :return: Returns tuple with (m, b) for line y=mx+b
    """
    assert (np.size(x) == np.size(y)), "x and y must have same dimension"
    assert (np.size(np.size(x)) == 1), "currently only supports 2D data"
    N = len(x)
    sum_x = np.sum(x)
    sum_y = np.sum(y)
    sum_x_y = np.sum(np.multiply(x, y))
    sum_x_squared = np.sum(np.square(x))
    m = (sum_x * sum_y - N * sum_x_y) / (sum_x**2 - N * sum_x_squared)
    b = (sum_x * sum_x_y - sum_x_squared * sum_y) \
                                         / (sum_x**2 - N * sum_x_squared)
    return (m , b)

if __name__ == '__main__':
    shm_euler_cromer(-10, 30, 1, 4, 0.1, 30, plot_x_vs_time=True, 
                                    plot_phase_space=False, plot_exact=True)
    shm_leapfrog(-10, 30, 1, 4, 0.1, 30, plot_x_vs_time=True, 
                                    plot_phase_space=False, plot_exact=True)
    shm_error(-10, 30, 1, 4, dt_power_range=(-3,-1), num_dt=100, end_t=30, 
                                    plot=True, plot_lin=True, loglog=True)
    plt.show()
    
    