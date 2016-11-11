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
#include <fstream>
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
 * Writes 2 column data to file
 * @param data data to write
 * @param fname string path to file (defaults to base of project)
 * @param delimiter string to place between data entries
 * @return returns true if write is successful
 */
template <typename num> bool write_csv(vector<tuple<num,num>> data, 
        tuple<string,string> header, string fname, string delimiter=",") {
    ofstream out(fname);
    if(out) {
        out << get<0>(header) << delimiter << get<1>(header) << endl;
        for(auto point : data) {
            out << get<0>(point) << delimiter << get<1>(point) << endl;
        }
        out.close();
        return true;
    } else {
        return false;
    }
}

/**
 * Writes 2 column data to file
 * @param data data to write
 * @param fname string path to file (defaults to base of project)
 * @param delimiter string to place between data entries
 * @return returns true if write is successful
 */
template <typename num> bool write_csv(vector<tuple<num,num,num>> data, 
        tuple<string,string,string> header, string fname, string delimiter=","){
    ofstream out(fname, ios::out);
    if(out) {
        out << get<0>(header) << delimiter << get<1>(header) << delimiter 
                                                    << get<2>(header) << endl;
        for(auto point : data) {
            out << get<0>(point) << delimiter << get<1>(point) << delimiter 
                                                    << get<2>(point) << endl;
        }
        out.close();
        return true;
    } else {
        return false;
    }
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
 * @param ofile name of file for output
 * @return returns vector containing tuples (time, position, angular velocity)
 */
vector<tuple<double, double, double>> shm_damped_driven(double theta0, 
        double ang_v0, double dt, double end_t, double nat_freq, 
        double friction_coef, double driving_freq,  double driving_torque,
        bool plot_x_vs_time=false, bool plot_phase_space=false, 
        bool linear=true, string ofile="") {
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
    
    if(ofile != "") {
        write_csv<double>(data, make_tuple("t","theta","ang_v"), ofile);
    }
    return data;
}

vector<tuple<double,vector<double>>> bifurcation(double dt, double end_t, 
        double nat_freq, double friction_coef, double driving_freq, 
        double d_driving_torque, double end_driving_torque, double theta0, 
        double ang_v0, bool linear=false, string ofile="", bool plot=false) {
    //TODO Figure out parameters and implement; also doxygen comment
    
}

/* 
 * 
 */
int main(int argc, char** argv) {
    shm_damped_driven(/*theta0*/ 0.2, /*ang_v0*/ 0, /*dt*/ 0.001, /*end_t*/ 600, 
            /*nat_freq*/ 1, /*friction_coef*/ 1/2.0, /*driving_freq*/ 2/3.0, 
            /*driving_torque*/ 1.2, /*plot_x_vs_y*/ true, 
            /*plot_phase_space*/ true, /*linear*/ false, /*output_file*/ "");
    return 0;
    ///////TODO Figure out why this seems to be exhibiting non-chaos
}

