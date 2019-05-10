'''
Created on Aug 31, 2016

@author: βrennan Cain (Hunter Damron)
'''

import numpy as np
import scipy.constants as sci
import matplotlib.pyplot as plt
import math

##### From one_d module:

'''
Created on Aug 31, 2016

@author: βrennan Cain (Hunter Damron)
'''

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
    
    t_drag, v_drag, a_drag = euler_method(drag_force_1d, fxi, dt, calc_range, 
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
        tpart, vpart, apart = euler_method(drag_force_1d, fxi=0, dx=dx[-1], 
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
    v = euler_method(drag_force_1d, fxi=0, dx=dt, calc_range=calc_range, b=b)[1]
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

##### From two_d module:

'''
Created on Aug 31, 2016

@author: βrennan Cain (Hunter Damron)
'''

def ax(v, vx, b, m):
    """
    Calculates acceleration in x direction
    :param v: total velocity (m/s)
    :param vx: x component of velocity (m/s)
    :param b: drag coefficient (kg/m)
    :param m: mass of object (kg)
    :return: Returns acceleration in x direction (m/s^2)
    """
    return - (b * v * vx) / m

def ay(v, vy, b, m):
    """
    Calculates acceleration in y direction
    :param v: total velocity (m/s)
    :param vy: y component of velocity (m/s)
    :param b: drag coefficient (kg/m)
    :param m: mass of object (kg)
    :return: Returns acceleration in y direction (m/s^2)
    """
    return - sci.g - (b * v * vy) / m

def drag_force_2d(vi, theta_i, C, A, m, rho=1.225, dt=0.01, object=None, 
                                                            plot=False):
    """
    Calculates trajectory of projectile
    :param vi: initial velocity (m/s)
    :param theta_i: initial angle (degrees)
    :param C: drag coefficient (dimensionless)
    :param a: cross-sectional area (m^2)
    :param m: mass (kg)
    :param dt: time delta (s)
    :param object: name of object being thrown
    :return: Returns a tuple containing arrays (x, y)
    """
    t = [0]
    x = [0]
    y = [0]
    vx = [vi * math.cos(theta_i * math.pi / 180)]
    vy = [vi * math.sin(theta_i * math.pi / 180)]
    
    b = 0.5 * C * rho * A
    
    while y[-1] >= 0:
        t.append(t[-1] + dt)
        v = math.sqrt((vx[-1]**2 + vy[-1]**2))
        vx.append(vx[-1] + ax(v, vx[-1], b, m) * dt)
        vy.append(vy[-1] + ay(v, vy[-1], b, m) * dt)
        x.append(x[-1] + vx[-1] * dt)
        y.append(y[-1] + vy[-1] * dt)
    
    if plot:
        plt.figure()
        plt.plot(x, y)
        plt.title("Projectile Motion%s" % 
                                 ("" if object is None else (" for %s" % object)))
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        
    return (x, y)

def projectile_motion(vi, theta_i, dt, object=None, plot=False):
    """
    Calculates trajectory with air resistance as negligible
    :param vi: initial velocity
    :param theta_i: initial angle
    :param dt: time delta
    :param object: name of object
    :param plot: plots if True
    """
    t = [0]
    x = [0]
    y = [0]
    vx = vi * math.cos(theta_i * math.pi / 180)
    vy = vi * math.sin(theta_i * math.pi / 180)
    
    while y[-1] >= 0:
        t.append(t[-1] + dt)
        x.append(vx * t[-1])
        y.append(vy * t[-1] - 0.5 * sci.g * t[-1]**2)
    
    if plot:
        plt.figure()
        plt.plot(x, y)
        plt.title("Projectile Motion without Drag for%s" % 
                                 ("" if object is None else (" for %s" % object)))
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        
    return x, y
        
def drag_force_various_theta(vi, dtheta_i, C, A, m, rho=1.225, dt=0.01, 
                                                        object=None, plot=False):
    """
    Calculates drag force for multiple angles
    :param vi: initial velocity
    :param dtheta_i: delta for initial angle
    :param C: drag coefficient
    :param A: cross sectional area
    :param m: mass
    :param rho: air mass density
    :param dt: time delta
    :param object: object title
    :param plot: plots if True
    """
    x = []
    y = []
    theta_i = np.arange(0, 90+dtheta_i, dtheta_i)
    for theta in theta_i:
        x_new, y_new = drag_force_2d(vi, theta, C, A, m, rho, dt, object)
        x.append(x_new)
        y.append(y_new)
    
    if plot:
        plt.figure()
        for xp, yp in zip(x, y):
            plt.plot(xp, yp)
        plt.title("Projectile Motion%s" % 
                            ("" if object is None else (" for %s" % object)))
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        
    return (theta_i, x, y)
        
def range_vs_theta(vi, dtheta_i, C, A, m, rho=1.225, dt=0.01, 
                            object=None, plot=False):
    """
    Calculates range vs initial angle using above parameters
    """
    theta_i, x, y = drag_force_various_theta(vi, dtheta_i, C, A, m, rho, dt)
    range = []
    for xp in x:
        range.append(xp[-1])
    
    if plot:
        plt.figure()
        plt.plot(theta_i, range)
        plt.title("Range vs Starting Angle for%s" % 
                            ("" if object is None else (" for %s" % object)))
        
    return theta_i, range
        
def range_vs_theta_wo_drag(vi, dtheta_i, dt, object=None, plot=False):
    """
    calculates the same range vs initial angle but without drag
    """
    theta_i = np.arange(0, 90+dtheta_i, dtheta_i)
    range = []
    for theta in theta_i:
        x_new = projectile_motion(vi, theta, dt)[0]
        range.append(x_new[-1])
        
    if plot:
        plt.figure()
        plt.plot(theta_i, range)
        plt.title("Range vs Starting Angle without Drag for%s" % 
                            ("" if object is None else (" for %s" % object)))
        
    return (theta_i, range)

def compare_trajectory(vi, theta_i, C, A, m, rho=1.225, dt=0.01, object=None, 
                                                            plot=False):
    """
    Calculates trajectory but without drag
    """
    x_drag, y_drag = drag_force_2d(vi, theta_i, C, A, m, rho, dt)
    x_wo_drag, y_wo_drag = projectile_motion(vi, theta_i, dt)
    
    if plot:
        plt.figure()
        plt.plot(x_drag, y_drag, 'r-', label="With Drag")
        plt.plot(x_wo_drag, y_wo_drag, 'k-', label="Without Drag")
        plt.title("Compare Trajectory with and without drag for%s" % 
                            ("" if object is None else (" for %s" % object)))
        plt.ylabel("y (m)")
        plt.xlabel("x (m)")
        plt.legend()
    
    return ((x_drag, y_drag), (x_wo_drag, y_wo_drag))

def compare_range_vs_theta(vi, dtheta_i, C, A, m, rho=1.225, dt=0.01, 
                                            object=None, plot=False):
    """
    Compares range vs theta for with and without drag
    """
    theta_drag, range_drag = range_vs_theta(vi, dtheta_i, C, A, m, rho, dt)
    theta_wo_drag, range_wo_drag = range_vs_theta_wo_drag(vi, dtheta_i, dt)
    
    if plot:
        plt.figure()
        plt.plot(theta_drag, range_drag, 'r-', label="With Drag")
        plt.plot(theta_wo_drag, range_wo_drag, 'k-', label="Without Drag")
        plt.title("Compare Range vs. Theta%s" % 
                            ("" if object is None else (" for %s" % object)))
        plt.ylabel("y (m)")
        plt.xlabel("x (m)")
        plt.legend()
    
    return ((theta_drag, range_drag), (theta_wo_drag, range_wo_drag))

if __name__ == '__main__':
    v = euler_method(drag_force_1d, calc_range=(0,3),
            plot=True, save=False, title="Euler Method Drag Force", 
            f_label="velocity", fp_label="acceleration", x_label="time", 
            y_label="velocity", yp_label="acceleration", b=1)[1]
    print("The function approaches %s from the left and %s from the right\n"
            % horizontal_asymptote(v))
      
    euler_method_various_fxi(drag_force_1d, 
        fxi_range=(-10, 10), calc_range=(0,3), plot=True, save=False, 
        title="Euler Method Drag Force for various vi", f_label="velocity",
        x_label="time", y_label="velocity", b=1)
  
    with_and_without_drag(fxi=-30, calc_range=(0,10), b=0.01, plot=True)
  
    drag_force_various_dx(0.4, dx_start=1, plot=True, calc_range=(0,3))
      
    print("Max time step under 1% error is", dt_error_threshold(0.2, 1, precision=1e-15))
  
    euler_vs_theoretical_various_b((0, 20), 0.1, plot=True)
  
    drag_force_2d(9000, 50, 0.1, 8.375*5.5*0.0254**2, 1.28, 
                                                    object="Syd", plot=True)
      
    drag_force_various_theta(9000, 0.1, 0.1, 8.375*5.5*0.0254**2, 1.28, 
                                            object="Lots of Syds", plot=True)
      
    range_vs_theta(9000, 1, 0.1, 8.375*5.5*0.0254**2, 1.28, 
                                             object="Lots of Syds", plot=True)
      
    projectile_motion(9000, 50, 0.1, "Syd",plot=True)
      
    range_vs_theta_wo_drag(9000, 1, 0.1, "Syd", plot=True)
  
    ### Plotting for copy of Anna Karenina using estimates ###
    vi = 20 # m/s
    theta_i = 30 # degrees
    C = 2.1
    A =  0.2127 * 0.1397 #m^2
    m = 0.6 # kg
    rho = 1.225 # kg/m^3
      
    compare_trajectory(vi, theta_i, C, A, m, rho=rho, object="Anna K", plot=True)
       
    compare_range_vs_theta(vi, 0.1, C, A, m, rho=rho, object="Anna K", plot=True)
    
    plt.show()