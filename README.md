## Synopsis

This project uses the leapfrog algorithm to numerically integrate the motion of 
non-chaotic pendulums. It first solves the motion of a linear, damped, driven 
pendulum and plots graphs of position, phase space plots, and a plot of 
amplitude vs. frequency. It then plots the motion of a nonlinear, non-damped, 
non-driven pendulum, which also does not exhibit chaos. 

This project can be found on (Github)[

## Usage

To use, unpack tar.gz file to desired path
Ensure all dependencies are installed and part of linker path
In a terminal run 
```bash
$ cd path_to_project/"C++ Simple Ideal Pendulum"
$ make
$ ./dist/Debug/GNU-Linux/c___simple_ideal_pendulum
```
Twenty graphs should spam the screen with science

## Dependencies

This project uses:
* Boost libraries which may be installed by
    ```
    brew install boost
    ```
    or
    ```
    sudo apt-get install libboost-all-dev
    ```
    on Mac and Linux
* Gnuplot plotting software which can be installed by 
    ```
    brew install gnuplot
    ```
    or
    ```
    sudo apt-get install gnuplot
    ```
    on Mac and Linux. 
* Gnuplot iostream which is included in the project source. It can also be installed by
    ```
    sudo apt-get install libgnuplot-iostream-dev
    ```
    on Linux or by using [gnuplot-iostream.h](https://github.com/dstahlke/gnuplot-iostream/blob/master/gnuplot-iostream.h) on Github.

## Contributor

The one and only Hunter Damron