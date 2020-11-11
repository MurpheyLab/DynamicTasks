//an example class definition for a particular system
//To define a class that can be used in the template SAC and cost functions,
//it must include variables Xcurr, Ucurr. It must also have functions for the
//system dynamics and relevant derivatives named f, dfdx, and hx respectively.
//Finally, Xcurr should be managed only by this class through the use of a 
//function called step that evolves the system forward in time by dt seconds.

// Xcurr by element 
// thetaX, thetaX_dot, thetaY, thetaY_dot, bowl:X, bowl:X_dot, bowl:Y, bowl:Y_dot


#ifndef BALLBOWL_HPP
#define BALLBOWL_HPP
#include <vector>
#include "rk4_int.hpp"
#include <cmath>
using namespace std;
typedef vector<double> Vec;
typedef vector<Vec> Mat;
//define any relevant constants such as pi here
const double PI = 3.1415926535987;

class BallBowl {
    double m, g; //system parameters: mass, gravity
    public:
        double B; // damping
		double h; // pendulum length
        double dt; //time between system updates
        double tcurr=0.0; //initilaize system time
        Vec Xcurr, Ucurr, Xpotential; //current system state arranged as q1,dq1,q2,dq2
        //class member function prototypes
        BallBowl (double, double, double,double,double);
        inline Vec f(const Vec& x, const Vec& u);
        void step(void);
		void simulate(void);
};

BallBowl::BallBowl (double _m, double _B, double _g, double _h, double _dt){
    m = _m; B = _B; g = _g; h = _h; //system parameters
    dt = _dt; //step size
}


inline Vec BallBowl::f(const Vec& x, const Vec& u){ //system dynamics, xdot = f(x)
	Vec xdot{
	x[1] * dt, // xthetadot
	(g / h * sin(x[0]) + B * x[1] / (m * h * h) - u[0] * cos(x[0])) * dt, //xthetadotdot
	x[3] * dt, // ythetadot
	(g / h * sin(x[2]) + B * x[3] / (m * h * h) - u[1] * cos(x[2])) * dt, // ythetadotdot
	x[5] * dt, // bowlxdot
	u[0] * dt, // bowlxdotdot
	x[7] * dt, // bowlydot
	u[1] * dt }; // bowlydotdot
    return xdot;
}; 


void BallBowl::step(){ //step the system forward by dt seconds
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    //cout<<Xcurr[0]<<" ";
    tcurr = tcurr+dt;
};

void BallBowl::simulate() { 
	Xpotential = RK4_step(this, Xcurr, Ucurr, dt);
};


#endif