# SHO 'nough Documentation

## Running Instructions
This project is written in Python using Numpy and Matplotlib to render graphs of simple harmonic motion and the failure of the Euler Method to conserve energy. The main code is located in

shm.main.py

The script main.py executes a main which prints the following eight graphs of SHM:  
- Euler Method position vs time for x0=-10, k=6, m=0.7, dt=0.01 for 30 seconds  
- Euler Method phase space plot for x0=-10, k=6, m=0.7, dt=0.01 for 30 seconds  
- Euler Cromer Method position vs time for x0=-10, k=6, m=0.7, dt=0.01 for 30 seconds  
- Euler Cromer Method phase space plot for x0=-10, k=6, m=0.7, dt=0.01 for 30 seconds  
- Euler Method position vs time for x0=-10, k=6, m=0.7 with dt values in range (10^-5, 10^-1) for 5 seconds  
- Euler Method phase space plot for x0=-10, k=6, m=0.7 with dt values in range (10^-5, 10^-1) for 5 seconds  
- Euler Cromer Method position vs time for x0=-10, k=6, m=0.7 with dt values in range (10^-5, 10^-1) for 5 seconds  
- Euler Cromer Method phase space plot for x0=-10, k=6, m=0.7 with dt values in range (10^-5, 10^-1) for 5 seconds  

The script does not have a user interface because the main is intended only to demonstrate the functionality of the code. Feel free to change values in the main to view different graphs. The code also includes pydoc style function API documentation but the project does not include a generated HTML API.

## General Overview of Code
The project uses the Euler Method and the Euler Cromer Method to numerically integrate simple harmonic motion. The Euler Method is included only to show that it fails in this case. The project also allows the user to visualize each method using different time steps to show their convergence to the exact solution as the time step approaches zero. 

Thank You,  
The Damron Development Team
