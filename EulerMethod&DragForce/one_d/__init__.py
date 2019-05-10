'''
Created on Aug 31, 2016

@author: βrennan Cain (Hunter Damron)
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sci
import math

def drag_force_1d(v, b, a=sci.g):
    """
    Second order derivative dv/dt for a falling object with air resistance
    :param b: drag coefficient for object
    :param v: velocity of object
    :param g: acceleration due to gravity
    :return: returns the change in velocity
    """
    assert b >= 0, "b must be positive for a terminal velocity to be reached"
    return a - b * v * abs(v)

def without_drag(v, a=sci.g):
    """
    Second order derivative dv/dt for a falling object without air resistance
    (included only for calculation within with_and_without_drag() function)
    :param v: unused variable included for proper positioning
    :param a: acceleration due to gravity
    :return: returns the change in velocity
    """
    return a

def euler_method(fprime, fxi=0, dx=0.01, calc_range=(0, 10), plot=False,
                 save=False, title="Euler Method Function and Derivative", 
                 f_label="f(x)", fp_label="f\'(x)", x_label="x", y_label="y",
                 yp_label="y\'", *args, **kwargs):
    """
    Calculates f(x) and f'(x) using euler method within limits
    :param fprime: derivative f'(x) of f(x) for euler method calculation (must
        accept one numerical argument and return one numerical value
    :param fxi: initial value of f(x)
    :param dx: difference between x values
    :param calc_range: range to calculate in
    :param plot: if True, function is plotted in one_d import of pyplot
    :param save: if true, function is saved in rundata folder
    :param title: Title of graph
    :param f_label: label for f(x) in graph and save
    :param fp_label: label for f\'(x) in graph and save
    :param x_label: x value label on graph
    :param y_label: y value label on graph
    :param *args: unnamed arguments for fprime
    :param **kwargs: unnamed arguments for fprime
    :return: Returns a tuple containing ( x, f(x), f'(x) ) within limit
    """
    f_x = [fxi]
    fp_x = [fprime(calc_range[0], *args, **kwargs)]
    x = np.arange(calc_range[0], calc_range[1] + dx, dx)
    for x_var in x[1:]:
        f_x.append(f_x[-1] + fprime(f_x[-1], *args, **kwargs) * dx)
        fp_x.append(fprime(f_x[-2], *args, **kwargs))
    if plot:
        plt.figure()
        plt.plot(x, f_x, 'b-', label=f_label)
        plt.title(title + " f(x)")
        plt.ylabel(y_label)
        plt.xlabel(x_label)
        
        plt.figure()
        plt.plot(x, fp_x, 'r-', label=fp_label)
        plt.title(title + " f\'(x)")
        plt.ylabel(yp_label)
        plt.xlabel(x_label)
        
        
    if save:
        np.savetxt("../rundata/euler_method(%s, %s, %s).csv" % (fprime.__name__,
                dx, calc_range), np.array([x, f_x, fp_x]).T, delimiter=",", 
                header="%s,%s,%s" % (x_label, f_label, fp_label))
    return (x, f_x, fp_x)

def term_vel(fxi=0, dx=0.01, calc_range=(0, 10), *args, **kwargs):
    """
    Calculates terminal velocity from theoretical method
    """
    return euler_method(drag_force, fxi, dx, calc_range, *args, **kwargs)

def horizontal_asymptote(fx):
    """
    Determines the horizontal asymptote of the function array provided
    :param fx: array containing values from function
    :return: Returns a tuple containing the left and right horizontal 
        asymptotes or for either if the asymptote does not exist
    """
    left = fx[0] if math.isclose(fx[0], fx[1]) else None
    right = fx[-1] if math.isclose(fx[-1], fx[-2]) else None
    return (left, right)

def euler_method_various_fxi(fprime, fxi_range, dfxi=0.5, dx=0.01, 
                        calc_range=(0, 10), plot=False, save=False, 
                        title="Euler Method Function and Derivative", 
                        f_label="f(x)", x_label="x", y_label="y", 
                        *args, **kwargs):
    """
    Calculates f(x) and f'(x) using euler method for various fxi
    :param fprime: derivative f'(x) of f(x) for euler method calculation (must 
        accept one numerical argument and return one numerical value
    :param fxi_range: tuple containing start and end values for fxi
    :param dfxi: change between each initial value of f(x)
    :param dx: difference between x values
    :param calc_range: range to calculate in
    :param plot: if True, function is plotted in one_d import of pyplot
    :param save: if true, function is saved in rundata folder
    :param title: Title of graph
    :param f_label: label for f(x) in graph and save
    :param fp_label: label for f\'(x) in graph and save
    :param x_label: x value label on graph
    :param y_label: y value label on graph
    :param *args: unnamed arguments for fprime
    :param **kwargs: unnamed arguments for fprime
    :return: Returns a tuple containing ( x, f(x) ) within limit where f(x)
        is an array containing series for all values xfi
    """
    x = np.arange(calc_range[0], calc_range[1] + dx, dx)
    f_x = []
    xfi = np.arange(fxi_range[0], fxi_range[1] + dfxi, dfxi)
    
    for i in xfi:
        f_x.append(euler_method(fprime=fprime, fxi=i, dx=dx, 
                                calc_range=calc_range, *args, **kwargs)[1])
    
    if plot:
        plt.figure()
        for f_x_series in f_x:
            plt.plot(x, f_x_series)
        plt.title(title)
        plt.ylabel(y_label)
        plt.xlabel(x_label)
    if save:
        np.savetxt("../rundata/euler_method_various_fxi(%s, %s, %s, %s, %s)" \
                   ".csv" %(fprime.__name__, fxi_range, dfxi, dx, calc_range),
                np.append([x], f_x, axis=0).T, delimiter=",")
    return (x, f_x)
    
def with_and_without_drag(fxi=0, dt=0.01, calc_range=(0, 10), plot=False,
                 save=False, title="Euler Method Function and Derivative", 
                 f_label="f(x)", fp_label="f\'(x)", x_label="x", y_label="y",
                 yp_label="y\'", *args, **kwargs):
    """
    Calculates f(x) and f'(x) using euler method within limits
    :param fprime: derivative f'(x) of f(x) for euler method calculation (must
        accept one numerical argument and return one numerical value
    :param fxi: initial value of f(x)
    :param dx: difference between x values
    :param calc_range: range to calculate in
    :param plot: if True, function is plotted in one_d import of pyplot
    :param save: if true, function is saved in rundata folder
    :param title: Title of graph
    :param f_label: label for f(x) in graph and save
    :param fp_label: label for f\'(x) in graph and save
    :param x_label: x value label on graph
    :param y_label: y value label on graph
    :param *args: unnamed arguments for fprime
    :param **kwargs: unnamed arguments for fprime
    :return: Returns a tuple containing 
                ( (t_drag, x_drag, v_drag, a_drag), (t, x, v, a) within limit
    """
    
    t_drag, v_drag, a_drag = euler_method(drag_force, fxi, dt, calc_range, 
                                                            *args, **kwargs)
    t, v, a = euler_method(without_drag, fxi, dt, calc_range)
    
    x_drag = [0]
    x = [0]
    
    w_drag_increasing = True
    wo_drag_increasing = True
    
    w_drag_t_max = -1
    wo_drag_t_max = -1
    
    w_drag_t_final = -1
    wo_drag_t_final = -1
    
    for time, vel in zip(t_drag, v_drag):
        x_drag.append(x_drag[-1] + vel * dt)
        if w_drag_increasing:
            if vel > 0:
                w_drag_increasing = False
                w_drag_t_max = time
        elif x_drag[-1] < 0:
            w_drag_t_final = time
        
    for time, vel in zip(t, v):
        x.append(x[-1] + vel * dt)
        if wo_drag_increasing:
            if vel > 0:
                wo_drag_increasing = False
                wo_drag_t_max = time
        elif x[-1] < 0:
            wo_drag_t_final = time
            
    print(w_drag_t_max, w_drag_t_final - w_drag_t_max) # TODO remove print statements
    print(wo_drag_t_max, wo_drag_t_final - wo_drag_t_max)
    
    if plot:
        plt.figure()
        plt.plot(t_drag, np.negative(x_drag[1:]), label="With Drag") # negated for viewing pleasure
        plt.plot(t, np.negative(x[1:]), label="Without Drag") # negated for viewing pleasure
        plt.title("Euler Method with and without Drag Force")
        plt.ylabel("Position")
        plt.xlabel("Time")
        plt.legend()
        plt.text(1, -40, "Without Drag: up=%s, down=%s"
                        % (wo_drag_t_max, wo_drag_t_final - wo_drag_t_max))
        plt.text(1, -50, "With Drag: up=%s," 
              " down=%s" % (w_drag_t_max, w_drag_t_final - w_drag_t_max))
    
    return ((t_drag, x_drag, v_drag, a_drag), (t, x, v, a))
    
def drag_force_various_dx(b, dx_start=5, dx_precision=1e-2, 
                        calc_range=(0, 10), plot=False, save=False, 
                        title="Euler Method Function and Derivative", 
                        f_label="f(x)", x_label="x", y_label="y"):
    """
    Calculates f(x) and f'(x) using euler method for various dx
    :param fprime: derivative f'(x) of f(x) for euler method calculation (must 
        accept one numerical argument and return one numerical value
    :param fxi_range: tuple containing start and end values for fxi
    :param dfxi: change between each initial value of f(x)
    :param dx: difference between x values
    :param calc_range: range to calculate in
    :param plot: if True, function is plotted in one_d import of pyplot
    :param save: if true, function is saved in rundata folder
    :param title: Title of graph
    :param f_label: label for f(x) in graph and save
    :param fp_label: label for f\'(x) in graph and save
    :param x_label: x value label on graph
    :param y_label: y value label on graph
    :param *args: unnamed arguments for fprime
    :param **kwargs: unnamed arguments for fprime
    :return: Returns a tuple containing ( x, f(x) ) within limit where f(x)
        is an array containing series for all values xfi
    """
    assert dx_start <= (calc_range[1] - calc_range[0])/2, "Too large of a " \
                                                           " dx malfunctions"
    dx = [dx_start]
    t = []
    v = []
    a = []
    
    plt.figure()
    while dx[-1] >= dx_precision:
        dx.append(dx[-1] / 2)
        print(dx[-1])
        tpart, vpart, apart = euler_method(drag_force, fxi=0, dx=dx[-1], 
                               calc_range=calc_range, b=b, a=sci.g)
        t.append(tpart)
        v.append(vpart)
        a.append(apart)
        plt.plot(tpart, vpart, 'k:')
    t_theoretical = np.arange(calc_range[0], calc_range[1], 0.001)
    v_theoretical = math.sqrt(sci.g/b)*np.tanh(sci.g*t_theoretical*math.sqrt(b/sci.g))
    plt.plot(t_theoretical, v_theoretical, 'r-')
    
def max_dx_error(b, calc_range=(0,10), dt=0.1):
    """
    Calculates maximum error from theoretical solution
    :param b: drag coefficient
    :param calc_range: tuple containing range for calculation
    :param dt: delta t
    """
    v = euler_method(drag_force, fxi=0, dx=dt, calc_range=calc_range, b=b)[1]
    t = np.arange(calc_range[0], calc_range[1]+dt, dt)
    v_theoretical = math.sqrt(sci.g / b) * np.tanh(sci.g * t * math.sqrt(b / sci.g))
    error = np.divide(np.abs(np.subtract(v_theoretical[1:], v[1:])), v[1:])*100
    return np.max(error)
    
def dt_error_threshold(b, threshold=1, calc_range=(0,10), precision=1e-5):
    """
    Calculates maximum dt for error to be less than threshold
    :param b: drag coefficient
    :param threshold: percentage to calculate maximum delta t
    :param calc_range: pair containing range for calculation
    :param precision: decimal value below which for ddt solution is reported
    """
    ddt = 0.01
    dt = 0.1
    error = math.inf
    while ddt > precision:
        dt = dt - ddt
        error = max_dx_error(b, calc_range, dt)
        if error < threshold:
            dt = dt + ddt
            ddt = ddt / 10
    return dt

def euler_vs_theoretical_various_b(b_range, db, dt=0.01, calc_range=(0,10), plot=False):
    """
    Calculates maximum error for multiple values b
    :param b_range: range of values to calculate for b (first is not included)
    :param db: delta for drag coefficient
    :param dt: time delta
    :param calc_range: calculation range
    :param plot: creates a plot if true
    """
    b = np.arange(b_range[0]+db, b_range[1]+db, db)
    error = []
    for b_part in b:
        error.append(max_dx_error(b_part, calc_range=calc_range, dt=dt))
    
    if plot:
        plt.figure()
        plt.plot(b, error)
        plt.title("Error vs b for one dimensional drag force")
        plt.ylabel("Error (%)")
        plt.xlabel("b (kg/m)")
    return (b, error)
        
