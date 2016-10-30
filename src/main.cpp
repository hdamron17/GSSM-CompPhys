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

#include "gnuplot-iostream.h"

using namespace std;

bool is_close(double arg1, double arg2, double epsilon=1e-3) {
    return abs(arg1 - arg2) < epsilon;
}

double dw_dt(double theta, double ang_v, double t, double nat_freq, 
        double friction_coef, double driving_freq, double driving_torque) {
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
    return -pow(nat_freq, 2) * theta - friction_coef * ang_v + \
                                    driving_torque * sin(driving_freq * t);
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

tuple<double,double,double> one_step(tuple<double,double,double> prev,double dt, 
        double nat_freq, double friction_coef, double driving_freq, 
        double driving_torque) {
    /**
     * Calculates damped driven pendulum solution by Leapfrog algorithm
     * @param previous: previous step (time, ang_position, ang_velocity)
     * @param dt: time delta (s)
     * @param nat_freq: natural frequency of ideal pendulum
     * @param friction_coef: coefficient of friction force
     * @param driving_freq: frequency of driving force
     * @param driving_torque: torque due to driving force
     * @return: returns tuple (time, ang_position, ang_velocity)
     */
    double time_new = get<0>(prev) + dt;
    double ai = dw_dt(get<1>(prev), get<2>(prev), get<0>(prev), nat_freq, 
                             friction_coef, driving_freq, driving_torque);
    double theta_new = get<1>(prev) + get<2>(prev) * dt + 0.5 * ai \
                                                            * pow(dt, 2);
    double af = dw_dt(theta_new, get<2>(prev), get<0>(prev), nat_freq, 
                            friction_coef, driving_freq, driving_torque);
    double ang_v_new = (get<2>(prev) + 0.5* (ai - pow(nat_freq, 2) * theta_new \
                    + driving_freq * sin(driving_freq * time_new)) * dt) \
                    / (1 + friction_coef / 2 * dt);
    return make_tuple(time_new, theta_new, ang_v_new);
}

vector<tuple<double, double, double>> shm_damped_driven(double theta0, 
        double ang_v0, double dt, double end_t, double nat_freq, 
        double friction_coef, double driving_freq,  double driving_torque, 
        bool plot_x_vs_time=false, bool plot_phase_space=false, 
        bool plot_exact=false) {
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
     * @param plot_x_vs_t: if True, plots position vs time
     * @param plot_phase_space: if True, plots momentum vs position
     * @return: returns vector containing tuples (time, position, angular velocity)
     */
    vector<tuple<double,double,double>> data; //formatted as (t, theta, ang_v)
    data.push_back(make_tuple(0, theta0, ang_v0));
    while(get<0>(data.back()) < end_t) {
        data.push_back(one_step(data.back(), dt, nat_freq, friction_coef, \
                                            driving_freq, driving_torque));
        //TODO use one_step() method instead of re-implementing it
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
           << "set ylabel \'Angular Momentum (kg*radians/s)\'\n"
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
bool _stable(deque<tuple<double,double>> *turns) {
    if(turns->size() <= 2) {
        //TODO throw error
        throw new exception(); 
        //TODO msg = "Queue must have at least 2 points to be stable"
    }
    double peak1 = get<0>((*turns)[0]);
    double peak2 = get<0>((*turns)[1]);
    int i = 0;
    for(tuple<double,double> turn : *turns) {
        if( (i % 2 == 0 && !is_close(get<0>(turn), peak1)) //case: i is even
            || (i % 2 == 1 && !is_close(get<0>(turn), peak2)) ) { //case: i is odd
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
 * @param linear_time: time after which velocity can be confirmed zero
 * @param loops_precision: number wavelengths before confirming amplitude
 * @return: returns tuple with (amplitude, frequency) at long time
 */
tuple<double,double> _ampl_freq(double theta0, 
        double ang_v0, double dt, double nat_freq, 
        double friction_coef, double driving_freq, double driving_torque, 
        double linear_time=1, int loops_precision=4) {
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
                                                                driving_torque);
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
            if(_stable(&turns)) {
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
    double driving_freq = nat_freq - 1; //TODO start at a different value
    while(driving_freq < end_driving_freq) {
        double ampl = get<0>(_ampl_freq(theta0, ang_v0, dt, nat_freq, \
                                friction_coef, driving_freq, driving_torque));
        cout << (driving_freq/nat_freq) << " -> " << ampl << endl; //TODO remove
        ret.push_back(make_tuple(ampl, driving_freq / nat_freq));
        driving_freq += d_driving_freq;
    }
    
    if(plot) {
        Gnuplot gp;
        //TODO proper titles and formatting
        gp << "set autoscale xy\n"
           << "set title \'Position vs. Time for Damped Driven Pendulum\'\n"
           << "set ylabel \'Position (radians)\'\n"
           << "set xlabel \'Time (s)\'\n"
           << "unset key\n"
           << "plot" << gp.file1d(ret) << "with lines lt rgb \"blue\""
                                                       " title \"Numerical\"\n";
    }
    return ret;
}

int main() {
        shm_damped_driven(/*theta0*/0.2, /*ang_v0*/0, /*dt*/0.01, /*end_t*/200, 
                /*nat_freq*/3, /*friction_coef*/2, /*driving_freq*/3, 
                /*driving_torque*/4, /*plot_x_vs_y*/true, 
                /*plot_phase_space*/false, /*plot_exact*/false);
        
//        tuple<double,double> single_ampl_freq = _ampl_freq(/*theta0*/0.2, /*ang_v0*/0, /*dt*/0.1, 
//                /*nat_freq*/1, /*friction_coef*/0, /*driving_freq*/0, 
//                /*driving_torque*/0); //TODO remove
//        cout << "Function has amplitude " << get<0>(single_ampl_freq) << " and frequency " << get<1>(single_ampl_freq) << endl; //TODO remove
        
        ampl_vs_freq(/*theta0*/0.2, /*ang_v0*/0, /*dt*/0.01, /*nat_freq*/5, 
        /*friction_coef*/2, /*d_driving_freq*/0.01, /*end_driving_freq*/10, 
        /*driving_torque*/4, /*plot*/true);
        
        //## COOL SCIENCE ##//
//        shm_damped_driven(/*theta0*/0.2, /*ang_v0*/0, /*dt*/0.01, /*end_t*/200, 
//                /*nat_freq*/5, /*friction_coef*/0.02, /*driving_freq*/0.15, 
//                /*driving_torque*/9, /*plot_x_vs_y*/true, 
//                /*plot_phase_space*/true, /*plot_exact*/false);
        
#ifdef _WIN32
	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	std::cout << "Press enter to exit." << std::endl;
	std::cin.get();
#endif
}