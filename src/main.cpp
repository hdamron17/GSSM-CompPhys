// Demo of sending data via temporary files.  The default is to send data to gnuplot directly
// through stdin.
//
// Compile it with:
//   g++ -o example-tmpfile example-tmpfile.cc -lboost_iostreams -lboost_system -lboost_filesystem

#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <typeinfo>

#include "gnuplot-iostream.h"

using namespace std;

bool is_close(double arg1, double arg2) {
    return abs(arg1 - arg2) < numeric_limits<double>::epsilon();
}

double dw_dt(double theta, double ang_v, double t, double nat_freq, 
        double friction_coef, double damping_freq, double damping_torque) {
    /**
     * Derivative of omega (angular velocity)
     * @param theta: anglular position
     * @param ang_v: angular velocity
     * @param t: time from beginning of damping force oscillation 
     * @param nat_freq: natural frequency of ideal pendulum
     * @param friction_coeff: coefficient of friction force
     * @param damping_freq: frequency of damping force
     * @param damping_torque: torque due to damping force
     * @return returns instantaneous acceleration
     */
    return -pow(nat_freq, 2) * theta - friction_coef * ang_v + \
                                    damping_torque * sin(damping_freq * t);
}

vector<tuple<double, double, double>> exact_damped_driven(double theta0, 
        double ang_v0, double dt, double end_t, double nat_freq, 
        double friction_coef, double damping_freq,  double damping_torque, 
        double m, bool plot_x_vs_time, bool plot_phase_space) {
    /**
     * Calculates exact damped driven pendulum solution
     * @param theta0: initial position
     * @param ang_v0: initial angular velocity
     * @param dt: time delta (s)
     * @param end_t: end time of calculation (s)
     * @param nat_freq: natural frequency of ideal pendulum
     * @param friction_coeff: coefficient of friction force
     * @param damping_freq: frequency of damping force
     * @param damping_torque: torque due to damping force
     * @param m: pendulum mass
     * @param plot_x_vs_t: if True, plots position vs time
     * @param plot_phase_space: if True, plots momentum vs position
     * @return: returns vector containing tuples (time, position, angular velocity)
     */
    vector<tuple<double,double,double>> data; //(t, theta, ang_v)
    data.push_back(make_tuple(0, theta0, ang_v0));
    while(get<0>(data.back()) < end_t) {
        tuple<double,double,double> last = data.back();
        double t_new = get<0>(last) + dt;
        if(is_close(friction_coef, 2 * nat_freq)) {
            //Critically damped
            double C = ang_v0 + theta0 * nat_freq; //TODO figure out what C needs to be
            double theta_new = (theta0 + C * t_new) \
                                            * exp(-friction_coef * t_new / 2);
            double ang_v_new = (C - friction_coef / 2 * (theta0 + C * t_new)) \
                                            * exp(-friction_coef * t_new / 2);
            data.push_back(make_tuple(t_new, theta_new, ang_v_new));
        } else if(friction_coef < 2 * nat_freq) {
            //Underdamped
            double phi = 0; //TODO figure out what phi needs to be
            double theta_new = theta0 * exp(-friction_coef * t_new / 2) \
                   * cos(sqrt(pow(nat_freq,2) - pow(friction_coef,2) / 4));
            double ang_v_new = (-friction_coef / 2 * cos(sqrt(pow(nat_freq,2) -\
                   pow(friction_coef,2) / 4) + phi) - sin(sqrt(pow(nat_freq,2)-\
                   pow(friction_coef,2) / 4) + phi));
            data.push_back(make_tuple(t_new, theta_new, ang_v_new));
        } else {
            //Overdamped
            double theta_new = theta0 * exp(-friction_coef * t_new / 2 + /*or - */ 
                      t_new * sqrt(pow(friction_coef,2) / 4 - pow(nat_freq,2)));
            double ang_v_new = theta0 * (-friction_coef / 2 + /*or - */ 
                      sqrt(pow(friction_coef,2) / 4 - pow(nat_freq,2)))
                      * exp(-friction_coef * t_new / 2 + /*or - */ 
                      sqrt(pow(friction_coef,2) / 4 - pow(nat_freq,2)) * t_new);
            data.push_back(make_tuple(t_new, theta_new, ang_v_new));
            //TODO figure out plus or minus
        }
        // TODO finish and check solutions
    }
    
    if(plot_x_vs_time) {
        Gnuplot gp;
        gp << "set autoscale xy\n"
           << "set title \'Position vs. Time for Damped Driven Pendulum\'\n"
           << "set ylabel \'Position (radians)\'\n"
           << "set xlabel \'Time (s)\'\n"
           << "unset key\n"
           << "plot" << gp.file1d(data) << "with lines lt rgb \"blue\"\n";
    }
    if(plot_phase_space) {
        vector<pair<double, double>> phase_space; //(position, momentum)
        for(tuple<double,double,double> point : data) {
            phase_space.push_back(make_pair(get<1>(point), m * get<2>(point)));
        }
        Gnuplot gp;
        gp << "set autoscale xy\n"
           << "set title \'Phase Space Diagram for Damped Driven Pendulum\'\n"
           << "set ylabel \'Angular Momentum (kg*radians/s)\'\n"
           << "set xlabel \'Position (radians)\'\n"
           << "unset key\n"
           << "plot" << gp.file1d(phase_space) << "with lines lt rgb \"red\"\n";
    }
    return data;
}

vector<tuple<double, double, double>> shm_damped_driven(double theta0, 
        double ang_v0, double dt, double end_t, double nat_freq, 
        double friction_coef, double damping_freq,  double damping_torque, 
        double m, bool plot_x_vs_time, bool plot_phase_space) {
    /**
     * Calculates damped driven pendulum solution by Leapfrog algorithm
     * @param theta0: initial position
     * @param ang_v0: initial angular velocity
     * @param dt: time delta (s)
     * @param end_t: end time of calculation (s)
     * @param nat_freq: natural frequency of ideal pendulum
     * @param friction_coeff: coefficient of friction force
     * @param damping_freq: frequency of damping force
     * @param damping_torque: torque due to damping force
     * @param m: pendulum mass
     * @param plot_x_vs_t: if True, plots position vs time
     * @param plot_phase_space: if True, plots momentum vs position
     * @return: returns vector containing tuples (time, position, angular velocity)
     */
    vector<tuple<double,double,double>> data; //formatted as (t, theta, ang_v)
    data.push_back(make_tuple(0, theta0, ang_v0));
    while(get<0>(data.back()) < end_t) {
        tuple<double,double,double> last = data.back();
        double ai = dw_dt(get<1>(last), get<2>(last), get<0>(last), nat_freq, 
                                 friction_coef, damping_freq, damping_torque);
        double theta_new = get<1>(last) + get<2>(last) * dt + 0.5 * ai \
                                                                * pow(dt, 2);
        double af = dw_dt(theta_new, get<2>(last), get<0>(last), nat_freq, 
                                friction_coef, damping_freq, damping_torque);
        double ang_v_new = get<2>(last) + 0.5 * (ai + af) * dt;
        double time_new = get<0>(last) + dt;
        data.push_back(make_tuple(time_new, theta_new, ang_v_new));
    }
    
    if(plot_x_vs_time) {
        Gnuplot gp;
        gp << "set autoscale xy\n"
           << "set title \'Position vs. Time for Damped Driven Pendulum\'\n"
           << "set ylabel \'Position (radians)\'\n"
           << "set xlabel \'Time (s)\'\n"
           << "unset key\n"
           << "plot" << gp.file1d(data) << "with lines lt rgb \"blue\"\n";
    }
    if(plot_phase_space) {
        vector<pair<double, double>> phase_space; //(position, momentum)
        for(tuple<double,double,double> point : data) {
            phase_space.push_back(make_pair(get<1>(point), m * get<2>(point)));
        }
        Gnuplot gp;
        gp << "set autoscale xy\n"
           << "set title \'Phase Space Diagram for Damped Driven Pendulum\'\n"
           << "set ylabel \'Angular Momentum (kg*radians/s)\'\n"
           << "set xlabel \'Position (radians)\'\n"
           << "unset key\n"
           << "plot" << gp.file1d(phase_space) << "with lines lt rgb \"red\"\n";
    }
    return data;
}

int main() {
        shm_damped_driven(/*theta0*/-0.5, /*ang_v0*/0, /*dt*/0.01, /*end_t*/10, 
                /*nat_freq*/8, /*friction_coef*/4, /*damping_freq*/0, 
                /*damping_torque*/0, /*mass*/3, /*plot_x_vs_y*/true, 
                /*plot_phase_space*/true);
        exact_damped_driven(/*theta0*/-0.5, /*ang_v0*/0, /*dt*/0.01, /*end_t*/10, 
                /*nat_freq*/8, /*friction_coef*/4, /*damping_freq*/0, 
                /*damping_torque*/0, /*mass*/3, /*plot_x_vs_y*/true, 
                /*plot_phase_space*/true);

#ifdef _WIN32
	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	std::cout << "Press enter to exit." << std::endl;
	std::cin.get();
#endif
}