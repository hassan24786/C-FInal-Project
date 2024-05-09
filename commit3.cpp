
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <memory>
const double speed_of_light = 1;

class FourMomentum
{
private:
  std::unique_ptr<std::vector<double>> momentum_vector;   //vector of smart pointers easier to handle 
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
  double compute_invariant_mass() const
  {
    double energy_squared = this->get_energy() * this->get_energy(); 
    double momentum_squared = (this->get_momentum_x() * this->get_momentum_x()) 
    + (this->get_momentum_y() * this->get_momentum_y()) 
    + (this->get_momentum_z() * this->get_momentum_z());
    return sqrt(energy_squared - momentum_squared);
  }
  ~FourMomentum() {}
};


class Calorimeter  //electron will inherit this 
{
private:
  double layer_1, layer_2, layer_3, layer_4;

public:
  Calorimeter() = default;
  void set_calorimeter_energies(const std::vector<double> &electron_energy_vector)
  {
    layer_1 = electron_energy_vector[0];
    layer_2 = electron_energy_vector[1];
    layer_3 = electron_energy_vector[2];
    layer_4 = electron_energy_vector[3];
  }
  double total_calorimeter_energy() const {return layer_1 + layer_2 + layer_3 + layer_4;}
};

class Electron:public Lepton
{
private:
  Calorimeter calorimeter;
public:
  Electron(double mass, double energy, double px, double py, double pz, bool anti_particle): 
    Lepton(mass, 0.5, energy, px, py, pz, -1, anti_particle ? 1: -1){}
  int get_lepton_number() const override {return anti_particle ? -1: 1;}
  std::string get_particle_type() const override {return "Electron";}
  double get_charge() const override {return anti_particle ? 1: -1;}
  double get_spin() const {return 0.5;}
  void print_particle_data() const override
  {
    Lepton::print_particle_data();
    std::cout <<"Total deposited energy in the calorimeter: "<<calorimeter.total_calorimeter_energy()<<std::endl;
  }
  void validate_deposited_energies(const std::vector<double> &deposited_energies) 
  {
  try 
  {
    calorimeter.set_calorimeter_energies(deposited_energies);
    if(std::abs(momentum->get_energy() - calorimeter.total_calorimeter_energy()) > 1e-3) //use closest acceptable values 
      throw std::runtime_error("The 4-momentum electron energy is not equal to deposited energy in calorimeter."); 
  }
  catch (const std::runtime_error&) 
  {
    std::cout<<"Incorrect deposited electron energies.  "<<std::endl;
  }
  };
}; 
