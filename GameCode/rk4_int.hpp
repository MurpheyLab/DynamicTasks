//This is a Runge-Kutte integration scheme that takes a class with 
//dynamics f and evolves it forward in time.
//Backwards integration can be performed by pass a negative value as
//the argument for dt.

#ifndef RK4_INT_HPP
#define RK4_INT_HPP
#include<vector>

std::vector<double> prod(std::vector<double> input, double coeff) {
	std::vector<double> temp=input;
	for (int n = 0; n < input.size(); n++) { temp[n] = temp[n] * coeff; }
	return temp;
};

std::vector<double> add(std::vector<double> v1, std::vector<double> v2) {
	std::vector<double> temp=v1;
	for (int n = 0; n < v1.size(); n++) { temp[n] = v1[n] + v2[n]; }
	return temp;
};

template <class T, class input> std::vector<double> RK4_step(T *sys, const std::vector<double>& x, input u, double dt){
  std::vector<double> k1, k2, k3, k4;
  k1 = sys->f(x, u);
  k2 = sys->f(add(x,prod(k1,0.5)), u); 
  k3 = sys->f(add(x,prod(k2,0.5)), u);
  k4 = sys->f(add(x,k3), u);
  //std::cout << k4[4]<<" ";
  return add(x,add(prod(k1,1./6.),add(prod(k2,1./3.),add(prod(k3,1./3.),prod(k4,1./6.)))));
    
};

#endif