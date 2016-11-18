## Synopsis

This project uses the leapfrog algorithm to numerically integrate the motion of 
chaotic pendulums. It plots position versus time of a nonlinear, damped, driven 
pendulum which exhibits chaos. It then constructs the bifurcation diagram with 
angular position at long time as a function of the driving torque, showing the route to chaos
via period doubling. It lastly plots Lyapunov divergence as a function of time; future 
work will include calculating the Lyapunov coefficient.

This project can be found on [Github](https://github.com/hdamron17/Pendulum-Chaos).

More complete documentation can be found at (project_root)/doc/html/index.html.

## Usage

To use, unpack tar.gz file to desired path
Ensure all dependencies are installed and part of linker path
In a terminal run 
```bash
$ cd path_to_project/"Pendulum-Chaos"
$ make
$ ./dist/Debug/GNU-Linux/pendulum-chaos
```
Twelve graphs should spam the screen with science

## Dependencies

This project uses:
<ul>
<li>Boost libraries which may be installed by
```
brew install boost
```
or
```
sudo apt-get install libboost-all-dev
```
on Mac and Linux</li>
<li>Gnuplot plotting software which can be installed by 
```
brew install gnuplot
```
or
```
sudo apt-get install gnuplot
```
on Mac and Linux.</li>
<li>Gnuplot iostream which is included in the project source. It can also be installed by
```
sudo apt-get install libgnuplot-iostream-dev
```
on Linux or by using [gnuplot-iostream.h](https://github.com/dstahlke/gnuplot-iostream/blob/master/gnuplot-iostream.h) on Github.</li>
</ul>

## Contributor

The one and only Hunter Damron

