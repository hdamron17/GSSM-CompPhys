'''
Created on Aug 31, 2016

@author: βrennan Cain (Hunter Damron)
'''

import numpy as np
import scipy.constants as sci
import matplotlib.pyplot as plt
import math

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
        x_new, y_new = drag_force(vi, theta, C, A, m, rho, dt, object)
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
    x_drag, y_drag = drag_force(vi, theta_i, C, A, m, rho, dt)
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


