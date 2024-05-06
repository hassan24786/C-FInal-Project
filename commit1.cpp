/*Description: This code provides a catalogue of all particles in the Standard Model of Physics. 
It can be used to view physical properties of particles, as well as their associated quantum numbers. 
The particles are split into separate classes, and some instances of particle decays are shown. 
Author: Hassan Hashmi
Date: 11/05/2024*/

#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <memory>
#include <cmath>
#include <list>
#include <algorithm>

const double speed_of_light = 1;

class FourMomentum
{
private:
  std::unique_ptr<std::vector<double>> momentum_vector;   //vector of smart pointers - easier to handle 
public:
  FourMomentum(double input_energy, double input_momentum_x, double input_momentum_y, double input_momentum_z) : momentum_vector(std::make_unique<std::vector<double>>(4))
  {
    set_energy(input_energy);
    set_momentum_x(input_momentum_x);
    set_momentum_y(input_momentum_y);
    set_momentum_z(input_momentum_z);
  }
  void set_energy(double input_energy)   //energy checks 
  {
    if(input_energy < 0)
      std::cout<<"The energy value cannot be negative."<<std::endl;
    else
      (*momentum_vector)[0] = input_energy / speed_of_light;
  }
  //setters and getters 
  void set_momentum_x(double input_momentum_x) {(*momentum_vector)[1] = input_momentum_x;}
  void set_momentum_y(double input_momentum_y) {(*momentum_vector)[2] = input_momentum_y;}
  void set_momentum_z(double input_momentum_z) {(*momentum_vector)[3] = input_momentum_z;}
  double get_energy() const {return (*momentum_vector)[0] * speed_of_light;}
  double get_momentum_x() const {return (*momentum_vector)[1];}
  double get_momentum_y() const {return (*momentum_vector)[2];}
  double get_momentum_z() const {return (*momentum_vector)[3];}
  std::vector<double> operator+(const FourMomentum &particle_two) 
  {
    std::vector<double> total_vector(4);                                //overloaded operators 
    total_vector[0] = get_energy() + particle_two.get_energy();
    total_vector[1] = get_momentum_x() + particle_two.get_momentum_x();
    total_vector[2] = get_momentum_y() + particle_two.get_momentum_y();
    total_vector[3] = get_momentum_z() + particle_two.get_momentum_z();
    return total_vector;
  }
  std::vector<double> operator-(const FourMomentum &particle) 
  {
    std::vector<double> subtracted_vector(4);
    subtracted_vector[0] = get_energy() - particle.get_energy();
    subtracted_vector[1] = get_momentum_x() - particle.get_momentum_x();
    subtracted_vector[2] = get_momentum_y() - particle.get_momentum_y();
    subtracted_vector[3] = get_momentum_z() - particle.get_momentum_z();
    return subtracted_vector;
  }
  double dot_product(const FourMomentum &particle)                     //dot product
  {
    return (*momentum_vector)[0] * (*particle.momentum_vector)[0] 
                      + (*momentum_vector)[1] * (*particle.momentum_vector)[1] 
                      + (*momentum_vector)[2] * (*particle.momentum_vector)[2] 
                      + (*momentum_vector)[3] * (*particle.momentum_vector)[3];
  }
  ~FourMomentum() {} 
};

class Particle           //plan to have all classes inherit from here 
{
protected:
  double mass;
  std::unique_ptr<FourMomentum> momentum;
public:
  Particle(double mass, double energy, double px, double py, double pz) :
    mass{mass}, momentum{std::make_unique<FourMomentum>(energy, px, py, pz)} {}
  virtual void print_particle_data() const = 0;   //makes base class abstract class 
}; 