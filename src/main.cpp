/**
 * Calculates pendulum behavior numerically using Leapfrog algorithm
 * 
 * Compile with:
 *   g++ -o main main.cpp -lboost_iostreams -lboost_system -lboost_filesystem
 * Or use the makefile attached to this project
 */
/////TODO Actual name for project

#include <iostream>
#include <vector>
#include <tuple>
#include "gnuplot-iostream.h"

using namespace std;

/**
 * Calculates if the difference between two args is less than epsilon
 * @param arg1 First argument
 * @param arg2 Second argument
 * @param epsilon Small number to account for numerical error
 * @return Returns true if numbers are close, else false
 */
bool is_close(double arg1, double arg2, double epsilon=numeric_limits<double>::epsilon() ) {
    return abs(arg1 - arg2) < epsilon;
}

/**
 * Derivative of omega for linear estimation (angular velocity)
 * @param theta anglular position
 * @param ang_v angular velocity
 * @param t time from beginning of driving force oscillation 
 * @param nat_freq natural frequency of ideal pendulum
 * @param friction_coef coefficient of friction force
 * @param driving_freq frequency of driving force
 * @param driving_torque torque due to driving force
 * @param linear if true, plots linear estimation, else nonlinear
 * @return returns instantaneous acceleration
 */
double dw_dt(double theta, double ang_v, double t, double nat_freq, 
        double friction_coef, double driving_freq, double driving_torque, 
        bool linear) {
    return -pow(nat_freq, 2) * (linear ? theta : sin(theta)) - friction_coef * \
                                ang_v + driving_torque * sin(driving_freq * t);
}

/**
 * Calculates one step of damped driven pendulum solution by Leapfrog algorithm
 * @param prev previous step (time, ang_position, ang_velocity)
 * @param dt time delta (s)
 * @param nat_freq natural frequency of ideal pendulum
 * @param friction_coef coefficient of friction force
 * @param driving_freq frequency of driving force
 * @param driving_torque torque due to driving force
 * @param linear plots linear estimation if true, otherwise nonlinear
 * @return returns tuple (time, ang_position, ang_velocity)
 */
tuple<double,double,double> one_step(tuple<double,double,double> prev,double dt, 
        double nat_freq, double friction_coef, double driving_freq, 
        double driving_torque, bool linear=true) {
    double theta_new, ang_v_new;
    
    double time_new = get<0>(prev) + dt;
    double ai = dw_dt(get<1>(prev), get<2>(prev), get<0>(prev), nat_freq, 
                       friction_coef, driving_freq, driving_torque, linear);
    theta_new = get<1>(prev) + get<2>(prev) * dt + 0.5 * ai \
                                                            * pow(dt, 2);
    ang_v_new = (get<2>(prev) + 0.5 * (ai - pow(nat_freq, 2) *  \
            (linear ? theta_new : sin(theta_new)) + driving_torque *  \
            sin(driving_freq * time_new) ) * dt) /  \
            (1 + friction_coef / 2 * dt);
    
    while(theta_new < -M_PI) {
        theta_new += 2 * M_PI; //add rotations until positive
    }
    while(theta_new > M_PI) {
        theta_new -= 2 * M_PI; //subtract rotations until negative
    }
    return make_tuple(time_new, theta_new, ang_v_new);
}

//////TODO REMOVE EXACT SOLUTION
/**
 * Calculates exact damped driven pendulum solution
 * Note initial angular velocity must be zero
 * @param theta0 initial position
 * @param dt time delta (s)
 * @param end_t end time of calculation (s)
 * @param nat_freq natural frequency of ideal pendulum
 * @param friction_coef coefficient of friction force
 * @param driving_freq frequency of driving force
 * @param driving_torque torque due to driving force
 * @param plot_x_vs_t if True, plots position vs time
 * @return returns vector containing tuples (time, angular position)
 */
vector<tuple<double, double>> exact_damped_driven(double theta0,
        double dt, double end_t, double nat_freq, 
        double friction_coef, double driving_freq,  double driving_torque, 
        bool plot_x_vs_time) {
    vector<tuple<double,double>> data; //(t, theta)
    data.push_back(make_tuple(0, theta0));
    while(get<0>(data.back()) < end_t) {
        tuple<double,double> last = data.back();
        double t_new = get<0>(last) + dt;
        if(is_close(friction_coef, 2 * nat_freq)) {
            //Critically damped
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
           << "plot" << gp.file1d(data) << "with points lc rgb \"red\" pt 7 ps 0.1\n";
    }
    return data;
}

/**
 * Calculates damped driven pendulum solution by Leapfrog algorithm
 * @param theta0 initial position
 * @param ang_v0 initial angular velocity
 * @param dt time delta (s)
 * @param end_t end time of calculation (s)
 * @param nat_freq natural frequency of ideal pendulum
 * @param friction_coef coefficient of friction force
 * @param driving_freq frequency of driving force
 * @param driving_torque torque due to driving force
 * @param plot_x_vs_t if True, plots position vs time
 * @param plot_phase_space if True, plots angular velocity vs position
 * @param plot_exact if true, plots exact solution (only valid with not damping/driving)
 * @param linear plots linear estimation if true, otherwise nonlinear
 * @return returns vector containing tuples (time, position, angular velocity)
 */
vector<tuple<double, double, double>> shm_damped_driven(double theta0, 
        double ang_v0, double dt, double end_t, double nat_freq, 
        double friction_coef, double driving_freq,  double driving_torque,
        bool plot_x_vs_time=false, bool plot_phase_space=false, 
        bool plot_exact=false, bool linear=true) {
    vector<tuple<double,double,double>> data; //formatted as (t, theta, ang_v)
    data.push_back(make_tuple(0, theta0, ang_v0));
    while(get<0>(data.back()) < end_t) {
        data.push_back(one_step(data.back(), dt, nat_freq, friction_coef, \
                                         driving_freq, driving_torque, linear));
    }
    
    if(plot_x_vs_time) {
        Gnuplot gp;
        gp << "set autoscale xy\n"
           << "set title \'Position vs. Time for Damped Driven Pendulum\'\n"
           << "set ylabel \'Position (radians)\'\n"
           << "set xlabel \'Time (s)\'\n"
           << "unset key\n"
           << "plot" << gp.file1d(data) << "with lines lt rgb \"blue\""
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
           << "set ylabel \'Angular Velocity (radians/s)\'\n"
           << "set xlabel \'Position (radians)\'\n"
           << "unset key\n"
           << "plot" << gp.file1d(phase_space) << "with points lc rgb \"red\" pt 7 ps 0.1\n";
    }
    return data;
}

vector<tuple<double,vector<double>>> bifurcation(double theta0, 
        double ang_v0, double dt, double end_t, double nat_freq, 
        double friction_coef, double driving_freq,  double driving_torque,
        bool plot_x_vs_time=false, bool plot_phase_space=false, 
        bool plot_exact=false, bool linear=true) {
    //TODO Figure out parameters and implement; also doxygen comment
}

/*
 * 
 */
int main(int argc, char** argv) {
    shm_damped_driven(/*theta0*/ 0.2, /*ang_v0*/ 0, /*dt*/ 0.01, /*end_t*/ 600, 
            /*nat_freq*/ 1, /*friction_coef*/ 1/2.0, /*driving_freq*/ 2/3.0, 
            /*driving_torque*/ 1.2, /*plot_x_vs_y*/ false, 
            /*plot_phase_space*/ true, /*plot_exact*/ false, /*linear*/ false);
    return 0;
    ///////TODO Figure out why this seems to be exhibiting non-chaos
}

