/**
 * Demo of sending data via temporary files.  The default is to send data to gnuplot directly
 * through stdin.
 * 
 * Compile it with:
 *   g++ -o main main.cpp -lboost_iostreams -lboost_system -lboost_filesystem
 * Or use the makefile attached to this project
 */

#include <vector>
#include <cmath>
#include <deque>
#include <map>
#include <typeinfo>

#include "gnuplot-iostream.h"

using namespace std;

bool is_close(double arg1, double arg2, double epsilon=numeric_limits<double>::epsilon() ) {
    return abs(arg1 - arg2) < epsilon;
}

/**
 * Derivative of omega (angular velocity)
 * @param theta: anglular position
 * @param ang_v: angular velocity
 * @param t: time from beginning of driving force oscillation 
 * @param nat_freq: natural frequency of ideal pendulum
 * @param friction_coef: coefficient of friction force
 * @param driving_freq: frequency of driving force
 * @param driving_torque: torque due to driving force
 * @return returns instantaneous acceleration
 */
double dw_dt(double theta, double ang_v, double t, double nat_freq, 
        double friction_coef, double driving_freq, double driving_torque) {
    return -pow(nat_freq, 2) * theta - friction_coef * ang_v + \
                                    driving_torque * sin(driving_freq * t);
}

/**
 * Derivative of omega for nonlinear, nondamped, nondriven pendulum
 * @param theta Angular position
 * @param nat_freq Natural frequency of pendulum
 * @return Returns instantaneous acceleration
 */
double dw_dt_nonlinear(double theta, double nat_freq) {
    return -pow(nat_freq, 2) * sin(theta);
}

/**
 * Calculates exact damped driven pendulum solution
 * Note: initial angular velocity must be zero
 * @param theta0: initial position
 * @param dt: time delta (s)
 * @param end_t: end time of calculation (s)
 * @param nat_freq: natural frequency of ideal pendulum
 * @param friction_coef: coefficient of friction force
 * @param driving_freq: frequency of driving force
 * @param driving_torque: torque due to driving force
 * @param plot_x_vs_t: if True, plots position vs time
 * @return: returns vector containing tuples (time, angular position)
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
           << "plot" << gp.file1d(data) << "with lines lt rgb \"blue\"\n";
    }
    return data;
}

/**
 * Calculates damped driven pendulum solution by Leapfrog algorithm
 * @param previous: previous step (time, ang_position, ang_velocity)
 * @param dt: time delta (s)
 * @param nat_freq: natural frequency of ideal pendulum
 * @param friction_coef: coefficient of friction force
 * @param driving_freq: frequency of driving force
 * @param driving_torque: torque due to driving force
 * @param linear: plots linear estimation if true, otherwise nonlinear
 * @return: returns tuple (time, ang_position, ang_velocity)
 */
tuple<double,double,double> one_step(tuple<double,double,double> prev,double dt, 
        double nat_freq, double friction_coef, double driving_freq, 
        double driving_torque, bool linear=true) {
    double theta_new, ang_v_new;
    double time_new = get<0>(prev) + dt;
    if(linear) {
        double ai = dw_dt(get<1>(prev), get<2>(prev), get<0>(prev), nat_freq, 
                                 friction_coef, driving_freq, driving_torque);
        theta_new = get<1>(prev) + get<2>(prev) * dt + 0.5 * ai \
                                                                * pow(dt, 2);
        double af = dw_dt(theta_new, get<2>(prev), get<0>(prev), nat_freq, 
                                friction_coef, driving_freq, driving_torque);
        ang_v_new = (get<2>(prev) + 0.5 * (ai - pow(nat_freq, 2) * theta_new \
                + driving_torque * sin(driving_freq * time_new) ) * dt) \
                / (1 + friction_coef / 2 * dt);
    } else {
        if(friction_coef != 0 || driving_freq != 0 || driving_torque != 0)
            cerr << "Nonlinear pendulum cannot have damping or driving\n";
        double ai = dw_dt_nonlinear(get<1>(prev), nat_freq);
        theta_new = get<1>(prev) + get<2>(prev) * dt + 0.5 * ai \
                                                                * pow(dt, 2);
        double af = dw_dt_nonlinear(theta_new, nat_freq);
        ang_v_new = get<2>(prev) + 0.5 * (ai + af) * dt;
    }
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
 * @param theta0: initial position
 * @param ang_v0: initial angular velocity
 * @param dt: time delta (s)
 * @param end_t: end time of calculation (s)
 * @param nat_freq: natural frequency of ideal pendulum
 * @param friction_coef: coefficient of friction force
 * @param driving_freq: frequency of driving force
 * @param driving_torque: torque due to driving force
 * @param linear: plots linear estimation if true, otherwise nonlinear
 * @param plot_x_vs_t: if True, plots position vs time
 * @param plot_phase_space: if True, plots angular velocity vs position
 * @param plot_exact: if true, plots exact solution (only valid with not damping/driving)
 * @return: returns vector containing tuples (time, position, angular velocity)
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
           << "set xlabel \'Time (s)\'\n";
        if(plot_exact) {
            vector<tuple<double,double>> exact = exact_damped_driven(theta0, dt,
                    end_t, nat_freq, friction_coef, driving_freq, 
                    driving_torque, false);
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
           << "set ylabel \'Angular Velocity (radians/s)\'\n"
           << "set xlabel \'Position (radians)\'\n"
           << "unset key\n"
           << "plot" << gp.file1d(phase_space) << "with lines lt rgb \"red\"\n";
    }
    return data;
}

/**
 * Calculates if deque is from stable oscillation
 * @param turns Queue containing (position, time) for function turns
 * @return Returns true if even indexes are all same and odd indexes are same
 */
bool _stable(deque<tuple<double,double>> *turns, double epsilon) {
    if(turns->size() <= 2) {
        cerr << "Queue must have at least 2 points to be stable";
    }
    double peak1 = get<0>((*turns)[0]);
    double peak2 = get<0>((*turns)[1]);
    int i = 0;
    for(tuple<double,double> turn : *turns) {
        if( (i % 2 == 0 && !is_close(get<0>(turn), peak1, epsilon)) //case: i is even
            || (i % 2 == 1 && !is_close(get<0>(turn), peak2, epsilon)) ) { //case: i is odd
            return false;
        }
        i++;
    }
    return true;
}

tuple<double,double> _calc_ampl_freq(deque<tuple<double,double>> *turns) {
    //Note: deque should be even size but is not required
    int i = 0;
    double prev_time; 
    double peak_total1 = 0; //sums either mins or maxes
    double peak_total2 = 0; //sums the opposite of peak_total1
    double time_total = 0; //sums periods
    for(tuple<double,double> turn : *turns) {
        if(i % 2 == 0) {
            peak_total1 += get<0>(turn);
            prev_time = get<1>(turn);
        } else {
            peak_total2 += get<0>(turn);
            time_total += get<1>(turn) - prev_time;
        }
        i++;
    }
    int size = turns->size();
    double avg_peak1 = peak_total1 / (size / 2 + size % 2); //avg for first set 
    double avg_peak2 = peak_total2 / (size / 2); //avg for second set
    double avg_time = time_total / (size / 2); //avg period
    return make_tuple( abs(avg_peak1 - avg_peak2) / 2, 1 / avg_time );
}

/**
 * Calculates amplitude from frequency
 * @param theta0: initial position
 * @param ang_v0: initial angular velocity
 * @param dt: time delta (s)
 * @param nat_freq: natural frequency of ideal pendulum
 * @param friction_coef: coefficient of friction force
 * @param driving_freq: frequency of driving force
 * @param driving_torque: torque due to driving force
 * @param linear: if true, plots linear estimation, else nonlinear
 * @param linear_time: time after which velocity can be confirmed zero
 * @param loops_precision: number wavelengths before confirming amplitude
 * @return: returns tuple with (amplitude, frequency) at long time
 */
tuple<double,double> _ampl_freq(double theta0, 
        double ang_v0, double dt, double nat_freq, 
        double friction_coef, double driving_freq, double driving_torque, 
        bool linear=true, double linear_time=20, int loops_precision=10) {
    int lin_max_steps = int(linear_time / dt) + 1;
    deque<tuple<double,double>> turns; //positions and times of mins and maxes
    tuple<double,double> initializer = make_tuple(theta0 + 1, 0); //value to fill deque before looping
    for(int i = 0; i < 2*loops_precision; i++) {
        turns.push_back(initializer);
    }
    int same_count = 0; //counts the number of points which are the same in a row
    bool direction = true;
    tuple<double,double,double> step = make_tuple(0, theta0, ang_v0);
    
    while(true) {
        step = one_step(step, dt, nat_freq, friction_coef, driving_freq, \
                                                        driving_torque, linear);
        double new_ang_v = get<2>(step);
        if(is_close(new_ang_v, 0)) {
            //case: next_step has same position as previousd
            same_count++;
        } else if(
            //case: next step is positive and direction is negative
                (new_ang_v > 0 && !direction)
            //case: ang_v is negative and direction is positive
                || (new_ang_v < 0 && direction) ) {
            direction = !direction; //reverse direction
            turns.pop_front(); //remove oldest min/max
            turns.push_back( make_tuple(get<1>(step), get<0>(step)) ); //new pos/time
            if(_stable(&turns, dt / 2)) { //precision capabilities based on dt
                return _calc_ampl_freq(&turns); //calculates amplitude and frequency
            }
        }
        if(same_count > lin_max_steps) {
            return make_tuple(0, 0); //linear has no amplitude or position
        }
    }
}

/**
 * Calculates amplitude as a function of frequency ratio
 * @param theta0: initial position
 * @param ang_v0: initial angular velocity
 * @param dt: time delta (s)
 * @param nat_freq: natural frequency of ideal pendulum
 * @param friction_coef: coefficient of friction force
 * @param d_driving_freq: step for frequencies of driving force
 * @param end_driving_freq: largest driving frequency
 * @param driving_torque: torque due to driving force
 * @param plot: if true, plots amplitude vs frequency
 * @return: returns vector containing tuples (frequency, amplitude)
 */
vector<tuple<double,double>> ampl_vs_freq(double theta0, 
        double ang_v0, double dt, double nat_freq, 
        double friction_coef, double d_driving_freq, double end_driving_freq, 
        double driving_torque, bool plot=false) {
    vector<tuple<double,double>> ret;
    double driving_freq = d_driving_freq;
    while(driving_freq < end_driving_freq) {
        double ampl = get<0>(_ampl_freq(theta0, ang_v0, dt, nat_freq, \
                                friction_coef, driving_freq, driving_torque));
        ret.push_back(make_tuple(driving_freq / nat_freq, ampl));
        driving_freq += d_driving_freq;
    }
    
    if(plot) {
        Gnuplot gp;
        gp << "set autoscale xy\n"
           << "set title \'Amplitude vs Frequency Ratio for Linear Damped Driven Pendulum\'\n"
           << "set ylabel \'Amplitude (rad)\'\n"
           << "set xlabel \'Natural Frequency - Driving Frequency Ratio\'\n"
           << "unset key\n"
           << "plot" << gp.file1d(ret) << "with lines lt rgb \"blue\"\n";
    }
    return ret;
}

/**
 * Calculates amplitude as a function of frequency ratio for various dampings
 * @param theta0: initial position
 * @param ang_v0: initial angular velocity
 * @param dt: time delta (s)
 * @param nat_freq: natural frequency of ideal pendulum
 * @param friction_coef: array of damping damping coefficients to plot
 * @param d_driving_freq: step for frequencies of driving force
 * @param end_driving_freq: largest driving frequency
 * @param driving_torque: torque due to driving force
 * @param plot: if true, plots amplitude vs frequency
 * @return: returns vector containing tuples (friction coeff, ampl_vs_freq data)
 */
map<double,vector<tuple<double,double>>> various_ampl_vs_freq(double theta0, 
        double ang_v0, double dt, double nat_freq, 
        vector<double> friction_coef, double d_driving_freq, double end_driving_freq, 
        double driving_torque, bool plot=false, vector<string> colors = 
        {"rgb \"red\"", "rgb \"blue\"", "rgb \"green\"", "rgb \"violet\""}) {
    map<double,vector<tuple<double,double>>> ret;
    for(size_t i = 0; i < friction_coef.size(); i++) {
        ret[friction_coef[i]] = ampl_vs_freq(theta0, ang_v0, dt, 
                nat_freq, friction_coef[i], d_driving_freq, 
                end_driving_freq, driving_torque);
    }
    if(plot) {
        Gnuplot gp;
        gp << "set autoscale xy\n"
           << "set title \'Position vs. Time for Damped Driven Pendulum\'\n"
           << "set ylabel \'Position (radians)\'\n"
           << "set xlabel \'Time (s)\'\n"
           << "set key\n"
           << fixed << setprecision(2)
           << "plot ";
        int i = 0;
        for(auto series : ret) {
            gp << gp.file1d(series.second) << "with lines lt ";
            if(i < colors.size()) {
                gp << colors[i];
            } else {
                gp << "black";
            }
            gp << " title \"q = " << series.first << "\", \\\n";
            i++;
        }
        gp << ";\n"; //to end plotting mutliple funcions
    }
    return ret;
}

/**
 * Calculates amplitude as a function of frequency ratio for various dampings
 * @param d_theta0: step size for initial position
 * @param ang_v0: initial angular velocity
 * @param dt: time delta (s)
 * @param nat_freq: natural frequency of ideal pendulum
 * @param friction_coef: array of damping damping coefficients to plot
 * @param driving_freq: frequency for driving force
 * @param driving_torque: torque due to driving force
 * @param plot: if true, plots amplitude vs frequency
 * @param linear: if true, plots linear estimation, else nonlinear
 * @return: returns vector containing tuples (friction coeff, ampl_vs_freq data)
 */
vector<tuple<double,double>> period_vs_ampl(double d_theta0, 
        double ang_v0, double dt, double nat_freq,
        double friction_coef, double driving_freq, double driving_torque, 
        bool linear=true, bool plot=false) {
    vector<tuple<double,double>> ret;
    if(ang_v0 != 0)
        cerr << "Cannot calculate if an initial velocity is present" << endl;
    for(double theta0 = d_theta0; theta0 <= M_PI-d_theta0; theta0 += d_theta0) {
        tuple<double,double> ampl_freq = _ampl_freq(theta0, ang_v0, dt, \
                    nat_freq, friction_coef, driving_torque, driving_torque, \
                    linear);
        //tuple contains (amplitude, frequency) - convert to (period, amplitude)
        ret.push_back(make_tuple(1 / get<1>(ampl_freq), get<0>(ampl_freq)));
    }
    if(plot) {
        Gnuplot gp;
        gp << "set autoscale xy\n"
           << "set title \'Period vs Amplitude\'\n"
           << "set ylabel \'Period (s)\'\n"
           << "set xlabel \'Amplitude (rad)\'\n"
           << "unset key\n";
        if(is_close(get<0>(ret.back()), get<0>(ret.front()))) {
            gp << "set xrange [" << M_PI-0.1 << ":" << M_PI+0.1 << "]\n";
            //point out of view to give series a domain
            ret.push_back(make_tuple(get<0>(ret.back())+0.001, get<1>(ret.back())));
        }
        gp << "plot" << gp.file1d(ret) << "with lines lt rgb \"blue\"\n";
    }
}

/**
 * Plots linear and nonlinear together
 * @param theta0: initial position
 * @param ang_v0: initial angular velocity
 * @param dt: time delta (s)
 * @param end_t: end time of calculation (s)
 * @param nat_freq: natural frequency of ideal pendulum
 * @param friction_coef: coefficient of friction force
 * @param driving_freq: frequency of driving force
 * @param driving_torque: torque due to driving force
 * @param plot_x_vs_t: if True, plots position vs time
 * @param plot_phase_space: if True, plots angular velocity vs position
 * Note: does not return the data, only plots it
 */
void plot_lin_and_nonlin(double theta0, double ang_v0, double dt, 
       double end_t, double nat_freq, bool plot_phase_space=false, 
       bool plot_exact=false) {
    vector<tuple<double,double,double>> lin = shm_damped_driven(theta0, ang_v0,
            dt, end_t, nat_freq, 0, 0,  0, false, false, false, true);
    vector<tuple<double,double,double>> nonlin = shm_damped_driven(theta0, 
            ang_v0, dt, end_t, nat_freq, 0, 0,  0, false, false, false, false);
    
    Gnuplot gp;
    gp << "set autoscale xy\n"
       << "set title \'Position vs. Time for Damped Driven Pendulum\'\n"
       << "set ylabel \'Position (radians)\'\n"
       << "set xlabel \'Time (s)\'\n"
       << "set key\n";
    gp << "plot" << gp.file1d(lin) << "with lines lt rgb \"blue\" title "
       << "\"Linear\", \\\n";
    if(plot_exact) {
        vector<tuple<double,double>> exact = exact_damped_driven(theta0, dt,
                end_t, nat_freq, 0, 0, 0, false);
        gp << gp.file1d(exact) << "with lines lt rgb \"black\""
                                                   " title \"Exact\", \\\n";
    }
    gp << gp.file1d(nonlin) << "with lines lt rgb \"red\" title "
       << "\"Nonlinear\"\n";
   if(plot_phase_space) {
       vector<pair<double,double>> lin_phase_space; //(position, velocity)
       vector<pair<double,double>> nonlin_phase_space;
       for(tuple<double,double,double> point : lin) {
           lin_phase_space.push_back(make_pair(get<1>(point), get<2>(point)));
       }
       for(tuple<double,double,double> point : nonlin) {
           nonlin_phase_space.push_back(make_pair(get<1>(point), get<2>(point)));
       }
       Gnuplot gp;
       gp << "set autoscale xy\n"
          << "set title \'Phase Space Diagram for Damped Driven Pendulum\'\n"
          << "set ylabel \'Angular Velocity (radians/s)\'\n"
          << "set xlabel \'Position (radians)\'\n"
          << "unset key\n"
          << "plot" << gp.file1d(lin_phase_space) << "with lines lt rgb \"blue\""
          << " title \"Linear\", \\\n"
          << gp.file1d(nonlin_phase_space) << "with lines lt rgb "
          << "\"red\" title \"Nonlinear\"\n";
   }
}

int main() {
//    shm_damped_driven(/*theta0*/3.14159, /*ang_v0*/0, /*dt*/0.01, /*end_t*/10, 
//            /*nat_freq*/5, /*friction_coef*/0, /*driving_freq*/0, 
//            /*driving_torque*/0, /*plot_x_vs_y*/true, 
//            /*plot_phase_space*/true, /*plot_exact*/false, /*linear*/false);
    
//    for(double nat_freq = 1; nat_freq <= 10; nat_freq += 1) {
//        shm_damped_driven(/*theta0*/5, /*ang_v0*/0, /*dt*/0.01, /*end_t*/10, 
//                /*nat_freq*/ nat_freq, /*friction_coef*/0, /*driving_freq*/0, 
//                /*driving_torque*/0, /*plot_x_vs_y*/true, 
//                /*plot_phase_space*/false, /*plot_exact*/false, /*linear*/false);
//        tuple<double,double> ampl_freq = _ampl_freq(/*theta0*/ 5, /*ang_v0*/ 0, /*dt*/ 0.01, /*nat_freq*/ nat_freq, 
//            /*friction_coef*/ 0, /*driving_freq*/ 0, /*driving_torque*/ 0, 
//            /*linear*/ false, /*linear_time*/ 1, /*loops_precision*/ 4);
//        cout << get<0>(ampl_freq) << " -> " << get<1>(ampl_freq) << endl;
//    }
    
//    ampl_vs_freq(/*theta0*/0.2, /*ang_v0*/0, /*dt*/0.0001, /*nat_freq*/5, 
//        /*friction_coef*/0.5, /*d_driving_freq*/0.1, /*end_driving_freq*/20, 
//        /*driving_torque*/4, /*plot*/true);
    
//    various_ampl_vs_freq(/*theta0*/0.2, /*ang_v0*/0, /*dt*/0.001, /*nat_freq*/5, 
//        /*friction_coef*/{0.1,1,4}, /*d_driving_freq*/1, /*end_driving_freq*/20, 
//        /*driving_torque*/4, /*plot*/true);

        //## COOL SCIENCE ##//
//    shm_damped_driven(/*theta0*/0.2, /*ang_v0*/0, /*dt*/0.01, /*end_t*/200, 
//        /*nat_freq*/5, /*friction_coef*/0.02, /*driving_freq*/0.15, 
//        /*driving_torque*/9, /*plot_x_vs_y*/true, 
//        /*plot_phase_space*/true, /*plot_exact*/false);
    
/******************************************************************************/
//    //Problem 1
//    shm_damped_driven(/*theta0*/ 0.2, /*ang_v0*/ 0, /*dt*/ 0.01, /*end_t*/ 20, 
//            /*nat_freq*/ 5, /*friction_coef*/ 2, /*driving_freq*/ 4, 
//            /*driving_torque*/ 0.5, /*plot_x_vs_y*/ true, 
//            /*plot_phase_space*/ true, /*plot_exact*/ false, /*linear*/ true);
//    
//    //Problems 3, 6, and 7
//    shm_damped_driven(/*theta0*/ 0.2, /*ang_v0*/ 0, /*dt*/ 0.01, /*end_t*/ 20, 
//            /*nat_freq*/ 4, /*friction_coef*/ 1, /*driving_freq*/ 0, 
//            /*driving_torque*/ 0, /*plot_x_vs_y*/ true, 
//            /*plot_phase_space*/ true, /*plot_exact*/ true, /*linear*/ true);
//
//    //Problems 4, 6, and 7
//    shm_damped_driven(/*theta0*/ 0.2, /*ang_v0*/ 0, /*dt*/ 0.01, /*end_t*/ 20, 
//            /*nat_freq*/ 4, /*friction_coef*/ 8, /*driving_freq*/ 0, 
//            /*driving_torque*/ 0, /*plot_x_vs_y*/ true, 
//            /*plot_phase_space*/ true, /*plot_exact*/ true, /*linear*/ true);
//    
//    //Problems 5, 6, and 7
//    shm_damped_driven(/*theta0*/ 0.2, /*ang_v0*/ 0, /*dt*/ 0.01, /*end_t*/ 20, 
//            /*nat_freq*/ 4, /*friction_coef*/ 15, /*driving_freq*/ 0, 
//            /*driving_torque*/ 0, /*plot_x_vs_y*/ true, 
//            /*plot_phase_space*/ true, /*plot_exact*/ true, /*linear*/ true);
//    
//    //Problems 10 - underdamped
//    shm_damped_driven(/*theta0*/ 0.2, /*ang_v0*/ 0, /*dt*/ 0.01, /*end_t*/ 20, 
//            /*nat_freq*/ 4, /*friction_coef*/ 1, /*driving_freq*/ 3, 
//            /*driving_torque*/ 0.5, /*plot_x_vs_y*/ true, 
//            /*plot_phase_space*/ true, /*plot_exact*/ false, /*linear*/ true);
//    
//    //Problems 10 - critically damped
//    shm_damped_driven(/*theta0*/ 0.2, /*ang_v0*/ 0, /*dt*/ 0.01, /*end_t*/ 20, 
//            /*nat_freq*/ 4, /*friction_coef*/ 8, /*driving_freq*/ 3, 
//            /*driving_torque*/ 0.5, /*plot_x_vs_y*/ true, 
//            /*plot_phase_space*/ true, /*plot_exact*/ false, /*linear*/ true);
//    
//    //Problems 10 - overdamped
//    shm_damped_driven(/*theta0*/ 0.2, /*ang_v0*/ 0, /*dt*/ 0.01, /*end_t*/ 20, 
//            /*nat_freq*/ 4, /*friction_coef*/ 15, /*driving_freq*/ 3, 
//            /*driving_torque*/ 0.5, /*plot_x_vs_y*/ true, 
//            /*plot_phase_space*/ true, /*plot_exact*/ false, /*linear*/ true);
//    
//    //Problem 11
//    various_ampl_vs_freq(/*theta0*/0.2, /*ang_v0*/0, /*dt*/0.001, /*nat_freq*/6, 
//        /*friction_coef*/{3, 12, 18}, /*d_driving_freq*/0.1, /*end_driving_freq*/20, 
//        /*driving_torque*/4, /*plot*/true);
//    
//    //Problem 14
//    shm_damped_driven(/*theta0*/3, /*ang_v0*/0, /*dt*/0.01, /*end_t*/10, 
//        /*nat_freq*/5, /*friction_coef*/0, /*driving_freq*/0, 
//        /*driving_torque*/0, /*plot_x_vs_y*/true, 
//        /*plot_phase_space*/false, /*plot_exact*/false, /*linear*/false);
//    
//    //Problem 15
//    plot_lin_and_nonlin(/*theta0*/ 0.1, /*ang_v0*/ 0, /*dt*/ 0.01, 
//       /*end_t*/ 20, /*nat_freq*/ 1, /*plot_phase_space*/ false, 
//       /*plot_exact*/ false);
//    
//    //Problem 16
//    plot_lin_and_nonlin(/*theta0*/ 3, /*ang_v0*/ 0, /*dt*/ 0.01, 
//       /*end_t*/ 20, /*nat_freq*/ 1, /*plot_phase_space*/ false, 
//       /*plot_exact*/ false);
//    
    //Problem 17 nonlinear
    period_vs_ampl(/*d_theta0*/ 0.1, /*ang_v0*/ 0, /*dt*/ 0.001, /*nat_freq*/ 1, 
        /*friction_coef*/ 0, /*driving_freq*/ 0, /*driving_torque*/ 0, 
        /*linear*/ false, /*plot*/ true);   
    
    //Problem 17 linear
    period_vs_ampl(/*d_theta0*/ 0.1, /*ang_v0*/ 0, /*dt*/ 0.001, /*nat_freq*/ 1, 
        /*friction_coef*/ 0, /*driving_freq*/ 0, /*driving_torque*/ 0, 
        /*linear*/ true, /*plot*/ true);
        
#ifdef _WIN32
    // For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
    // the gnuplot window doesn't get closed.
    std::cout << "Press enter to exit." << std::endl
    std::cin.get();
#endif
}