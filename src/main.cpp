/**
 * Calculates chaotic pendulum behavior numerically using Leapfrog algorithm
 * 
 * Compile with:
 *   g++ -o main main.cpp -lboost_iostreams -lboost_system -lboost_filesystem
 * Or use the makefile attached to this project
 */

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
 * Writes 3 column data to file
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

//zoom = (x,y,w,h)
/**
 * Creates a bifurcation diagram of damped, driven, nonlinear pendulum
 * @param dt Time step
 * @param end_t End time for each run
 * @param nat_freq Natural frequency of pendulum ( sqrt(g/l) )
 * @param friction_coef Coefficient of friction
 * @param driving_freq Driving frequency
 * @param d_driving_torque Step between each driving frequency on plot
 * @param end_driving_torque Endpoint for driving frequency
 * @param start_driving_torque Starting point for driving torque
 *      (must be smaller than driving torque)
 * @param theta0 Starting angular position (should not affect plot)
 * @param ang_v0 Starting angular velocity (should not affect plot)
 * @param points Number of points included for each iteration
 * @param linear If true, plots with linear estimation, else nonlinear
 * @param plot If true, plots to Gnuplot
 * @param zoom Numerous zoom with dimensions (x, y, width, height)
 * @return Returns vector containing the data of the plot in the form 
 *      [ (driving_torque, [positions at long time...] ) ]
 */
vector<tuple<double,vector<double>>> bifurcation(double dt, double end_t, 
        double nat_freq, double friction_coef, double driving_freq, 
        double d_driving_torque, double end_driving_torque, double start_driving_torque, double theta0=0, 
        double ang_v0=0, int points=10, bool linear=false,  
        bool plot=false, vector<tuple<double,double,double,double>> zoom={}) {
    
    vector<tuple<double,vector<double>>> data;
    
    for(double driving_torque = start_driving_torque; 
            driving_torque <= end_driving_torque; 
            driving_torque += d_driving_torque) {
        vector<tuple<double,double,double>> series = shm_damped_driven(theta0, ang_v0, 
                dt, end_t, nat_freq, friction_coef, driving_freq, 
                driving_torque, false, false, linear, "");
        vector<double> points_of_interest;
        int counter = 0;
        bool direction = get<2>(series.back()) > 0; //Gets initial direction
        for(auto iter = series.end(); iter > series.begin() && counter < points; iter--) {
            if( (direction && get<2>(*iter) < 0) || (!direction && get<2>(*iter) > 0) ) {
                points_of_interest.push_back(get<1>(*iter));
                direction = !direction;
                counter++;
            }
        }
        data.push_back(make_tuple(driving_torque, points_of_interest));
    }
    if(plot) {
        vector<tuple<double,double>> modified; //expands multiple points to one driving_torque
        for(auto series : data) {
            for(double value : get<1>(series)) {
                modified.push_back(make_tuple(get<0>(series), value));
            }
        }
        
        Gnuplot gp;
        gp  << "set autoscale xy\n"
            << "set title \'Position vs Driving Torque\'\n"
            << "set ylabel \'Position at long time (radians)\'\n"
            << "set xlabel \'Driving Torque (N*m)\'\n"
            << "unset key\n"
            << "plot" << gp.file1d(modified) << "with points lc rgb \"red\" pt 7 ps 0.1\n";
        int i = 1;
        for(auto rect : zoom) {
            Gnuplot zoom_gp;
            zoom_gp << "set title \'Bifurcation Zoom" << i << "\'\n"
                << "set ylabel \'Position at long time (radians)\'\n"
                << "set xlabel \'Driving Torque (N*m)\'\n"
                << "set xrange [" << get<0>(rect) << ":" 
                     << get<0>(rect) + get<2>(rect) << "]\n"
                << "set yrange [" << get<1>(rect) << ":"
                     << get<1>(rect) + get<3>(rect) << "]\n"
                << "unset key\n"
                << "plot" << zoom_gp.file1d(modified) << "with points lc rgb \"red\" pt 7 ps 0.1\n";
            i++;
        }
    }
    return data;
}

vector<tuple<double,double>> lyapunov (double theta1, double theta2, 
        double ang_v0, double dt, double end_t, double nat_freq, 
        double friction_coef, double driving_freq,  double driving_torque,
        bool plot=false, bool linear=true, string ofile="") {
    
    vector<tuple<double,double,double>> s1 =  shm_damped_driven(theta1, 
        ang_v0, dt, end_t, nat_freq, friction_coef, driving_freq, 
        driving_torque, false, false, false, "");
    vector<tuple<double,double,double>> s2 =  shm_damped_driven(theta2, 
        ang_v0, dt, end_t, nat_freq, friction_coef, driving_freq, 
        driving_torque, false, false, false, "");
    vector<tuple<double,double>> lyapunov_series(s1.size());
    for(size_t i = 0; i < lyapunov_series.size(); i++) {
        lyapunov_series[i] = make_tuple(get<0>(s1[i]), \
                                            abs(get<2>(s1[i]) - get<2>(s2[i])));
    }
    //TODO plot and save
    if(plot) {
        Gnuplot gp;
        gp  << "set autoscale xy\n"
            << "set title \'Lyapunov series divergence\'\n"
            << "set ylabel \'Difference between angular positions (m)\'\n"
            << "set xlabel \'Time (s)\'\n"
            << "unset key\n"
            << "plot" << gp.file1d(lyapunov_series) << "with points lc rgb \"red\" pt 7 ps 0.05\n";
    }
    return lyapunov_series;
}

/* 
 * Constructs a lot of plot and spams the screen
 */
int main() {
//    // Problem 1
//    shm_damped_driven(/*theta0*/ 0.3, /*ang_v0*/ 0, /*dt*/ 0.001, /*end_t*/ 100, 
//            /*nat_freq*/ 1, /*friction_coef*/ 1/2.0, /*driving_freq*/ 2/3.0, 
//            /*driving_torque*/ 1.6, /*plot_x_vs_y*/ true, 
//            /*plot_phase_space*/ true, /*linear*/ false, /*output_file*/ "");
//    
//    // Problem 3 - Modified initial conditions position vs time plot
//    shm_damped_driven(/*theta0*/ 0.2, /*ang_v0*/ 0, /*dt*/ 0.001, /*end_t*/ 100, 
//            /*nat_freq*/ 1, /*friction_coef*/ 1/2.0, /*driving_freq*/ 2/3.0, 
//            /*driving_torque*/ 1.6, /*plot_x_vs_y*/ true, 
//            /*plot_phase_space*/ true, /*linear*/ false, /*output_file*/ "");
//    
//    // Problem 2 - Bifurcation Base Plot
//    bifurcation(/*dt*/ 0.01, /*end_t*/ 400, 
//        /*nat_freq*/ 1, /*friction_coef*/ 0.5, /*driving_freq*/ 2/3.0, 
//        /*d_driving_torque*/ 0.01, /*end_driving_torque*/ 4, 
//        /*start_driving_torque*/ 0, /*theta0*/ 0.3, /*ang_v0*/ 0, /*points*/ 10, 
//        /*linear*/ false, /*plot*/ true, /*zoom*/ {});
//    
//    // Problem 2 - Modified Friction Coefficient
//    bifurcation(/*dt*/ 0.01, /*end_t*/ 400, 
//        /*nat_freq*/ 1, /*friction_coef*/ 1, /*driving_freq*/ 2/3.0, 
//        /*d_driving_torque*/ 0.01, /*end_driving_torque*/ 4, 
//        /*start_driving_torque*/ 0, /*theta0*/ 0.3, /*ang_v0*/ 0, /*points*/ 10, 
//        /*linear*/ false, /*plot*/ true, /*zoom*/ {});
//    
//    // Problem 3 - Modified initial conditions bifurcation
//    bifurcation(/*dt*/ 0.01, /*end_t*/ 400, 
//        /*nat_freq*/ 1, /*friction_coef*/ 0.5, /*driving_freq*/ 2/3.0, 
//        /*d_driving_torque*/ 0.01, /*end_driving_torque*/ 4, 
//        /*start_driving_torque*/ 0, /*theta0*/ 0.2, /*ang_v0*/ 0, /*points*/ 10, 
//        /*linear*/ false, /*plot*/ true, /*zoom*/ {});
//    
    // Problem 4 Initial View
    bifurcation(/*dt*/ 0.01, /*end_t*/ 200, 
        /*nat_freq*/ 1, /*friction_coef*/ 1/2.0, /*driving_freq*/ 2/3.0, 
        /*d_driving_torque*/ 0.005, /*end_driving_torque*/ 4, 
        /*start_driving_torque*/ 0, /*theta0*/ 0.2, /*ang_v0*/ 0, /*points*/ 10, 
        /*linear*/ false, /*plot*/ true, /*zoom*/ {});
    // Problem 4 - ZOOM 1
    bifurcation(/*dt*/ 0.005, /*end_t*/ 200, 
        /*nat_freq*/ 1, /*friction_coef*/ 1/2.0, /*driving_freq*/ 2/3.0, 
        /*d_driving_torque*/ 0.001, /*end_driving_torque*/ 1.491, 
        /*start_driving_torque*/ 1.28, /*theta0*/ 0.2, /*ang_v0*/ 0, /*points*/ 10, 
        /*linear*/ false, /*plot*/ true, /*zoom*/ {});
    // Problem 4 - ZOOM 2
    bifurcation(/*dt*/ 0.001, /*end_t*/ 200, 
        /*nat_freq*/ 1, /*friction_coef*/ 1/2.0, /*driving_freq*/ 2/3.0, 
        /*d_driving_torque*/ 0.0005, /*end_driving_torque*/ 1.491, 
        /*start_driving_torque*/ 1.415, /*theta0*/ 0.2, /*ang_v0*/ 0, /*points*/ 10, 
        /*linear*/ false, /*plot*/ true, /*zoom*/ {});
//
//    // TODO remove lyapunov
//    lyapunov (/*theta1*/0.2, /*theta2*/ 0.20001, /*ang_v0*/ 0, /*dt*/ 0.001, 
//        /*end_t*/ 40, /*nat_freq*/ 1, /*friction_coef*/ 0.5, 
//        /*driving_freq*/ 2/3.0,  /*driving_torque*/ 1.6, /*plot*/ true, 
//        /*linear*/ true, /*ofile*/ "");
}

