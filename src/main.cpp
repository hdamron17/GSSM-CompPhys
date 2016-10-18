/**
 * Demo of sending data via temporary files.  The default is to send data to gnuplot directly
 * through stdin.
 * 
 * Compile it with:
 *   g++ -o main main.cpp -lboost_iostreams -lboost_system -lboost_filesystem
 * Or use the makefile attached to this project
 */

#include <vector>
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

vector<tuple<double, double>> exact_damped_driven(double theta0,
        double dt, double end_t, double nat_freq, 
        double friction_coef, double damping_freq,  double damping_torque, 
        bool plot_x_vs_time) {
    /**
     * Calculates exact damped driven pendulum solution
     * Note: initial angular velocity must be zero
     * @param theta0: initial position
     * @param dt: time delta (s)
     * @param end_t: end time of calculation (s)
     * @param nat_freq: natural frequency of ideal pendulum
     * @param friction_coeff: coefficient of friction force
     * @param damping_freq: frequency of damping force
     * @param damping_torque: torque due to damping force
     * @param plot_x_vs_t: if True, plots position vs time
     * @return: returns vector containing tuples (time, angular position)
     */
    
    vector<tuple<double,double>> data; //(t, theta)
    data.push_back(make_tuple(0, theta0));
    while(get<0>(data.back()) < end_t) {
        tuple<double,double> last = data.back();
        double t_new = get<0>(last) + dt;
        if(is_close(friction_coef, 2 * nat_freq)) {
            //Critically damped
//            double theta_new = theta0 * (friction_coef * t_new / 2) \
//                                              * exp(-friction_coef * t_new / 2);
            double temp = friction_coef * t_new / 2;
            double theta_new = theta0 * (1 + temp) * exp(-temp);
            data.push_back(make_tuple(t_new, theta_new));
        } else if(friction_coef < 2 * nat_freq) {
            //Underdamped
            double temp = sqrt(pow(nat_freq,2) - pow(friction_coef,2) / 4);
            double theta_new = theta0 * exp(-friction_coef * t_new / 2) * \
                                (cos(temp * t_new) + friction_coef / 2 / temp \
                                * sin(temp * t_new));
            data.push_back(make_tuple(t_new, theta_new));
        } else {
            //Overdamped
            double temp = sqrt(pow(friction_coef,2) / 4 - pow(nat_freq,2));
            double theta_new = theta0 * exp(-friction_coef * t_new / 2) * \
                                (cosh(temp * t_new) + friction_coef / 2 / temp \
                                * sinh(temp * t_new));
            data.push_back(make_tuple(t_new, theta_new));
        }
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
    return data;
}

vector<tuple<double, double, double>> shm_damped_driven(double theta0, 
        double ang_v0, double dt, double end_t, double nat_freq, 
        double friction_coef, double damping_freq,  double damping_torque, 
        bool plot_x_vs_time=false, bool plot_phase_space=false, 
        bool plot_exact=false) {
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
           << "set xlabel \'Time (s)\'\n";
        if(plot_exact) {
            vector<tuple<double,double>> exact = exact_damped_driven(theta0, dt,
                    end_t, nat_freq, friction_coef, damping_freq, 
                    damping_torque, false);
            gp << "set key\n"
               << "plot" << gp.file1d(exact) << "with lines lt rgb \"black\""
                                                       " title \"Exact\", \\\n";
        } else {
            gp << "unset key\nplot";
        }
        gp << gp.file1d(data) << "with lines lt rgb \"blue\""
                                                       " title \"Numerical\"\n";
    }
    if(plot_phase_space) {
        vector<pair<double, double>> phase_space; //(position, velocity)
        for(tuple<double,double,double> point : data) {
            phase_space.push_back(make_pair(get<1>(point), get<2>(point)));
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
        shm_damped_driven(/*theta0*/0.2, /*ang_v0*/0, /*dt*/0.01, /*end_t*/10, 
                /*nat_freq*/3, /*friction_coef*/10, /*damping_freq*/20, 
                /*damping_torque*/2, /*plot_x_vs_y*/true, 
                /*plot_phase_space*/false, /*plot_exact*/true);
#ifdef _WIN32
	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	std::cout << "Press enter to exit." << std::endl;
	std::cin.get();
#endif
}