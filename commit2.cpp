#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <memory>
#include <cmath>


class Lepton:public Particle
{
protected:
  double spin;
  bool anti_particle; //true if antiparticle present 
  int charge;   

public:
  Lepton(double mass, double input_spin, double energy, double px, double py, double pz, int input_charge, bool anti_particle = false): 
    Particle(mass, energy, px, py, pz), spin{input_spin}, anti_particle{anti_particle}, charge{anti_particle ? input_charge: input_charge *(-1)}   //easier than if statements
  {
    double invariant_mass = momentum->compute_invariant_mass();
    if(std::abs(invariant_mass - mass) > 1e-3) 
    {
      double adjusted_energy = std::sqrt(mass * mass + (px * px + py * py + pz * pz)); //rearranged formula for invariant mass 
      momentum = std::make_unique<FourMomentum>(adjusted_energy, px, py, pz); 
      std::cout<<"The invariant mass is not equal to the rest mass of the particle. Calculated mass: "
      <<invariant_mass<<" MeV/c^2."<<" The energy in the four-momentum vector has been corrected to: "
      <<adjusted_energy<<" MeV."<<std::endl;
    }
  }
  double get_spin() const {return spin;}
  int get_lepton_number() const override {return anti_particle ? -1: 1;} //override abstract class virtual function
  double get_charge() const override {return charge;}
  std::string get_particle_type() const override{return "Lepton";}
  double get_energy() const {return momentum->get_energy();}   //access FourMomentum class 
  double get_momentum_x() const {return momentum->get_momentum_x();}
  double get_momentum_y() const {return momentum->get_momentum_y();}
  double get_momentum_z() const {return momentum->get_momentum_z();}  
  void set_energy(double energy) {momentum->set_energy(energy);}
  void set_momentum_x(double px) {momentum->set_momentum_x(px);}
  void set_momentum_y(double py) {momentum->set_momentum_y(py);}
  void set_momentum_z(double pz) {momentum->set_momentum_z(pz);}
  void print_particle_data() const override
  {
    std::cout<<"Lepton - Mass: "<<mass<<" MeV/c^2, Spin: "<<spin
    <<" Charge: "<<get_charge()<<" Lepton Number: "<<get_lepton_number() 
    <<", 4-Momentum: ("<<momentum->get_energy() << ", "<<momentum->get_momentum_x()<<", " 
    <<momentum->get_momentum_y()<<momentum->get_momentum_z()<<std::endl;
  }
};

class Quark:public Particle
{
protected:
  double spin;
  bool anti_particle;  
  double charge; 
  double baryon_number;
  std::string quark_flavour;
  std::string colour_charge;

public:
  Quark(double mass, double input_spin, double energy, double px, double py, double pz, const std::string &flavour, bool anti_particle, const std::string &colour):
    Particle(mass, energy, px, py, pz), spin{input_spin}, anti_particle{anti_particle}, quark_flavour{flavour}, colour_charge{colour}
  {
    double invariant_mass = momentum->compute_invariant_mass(); //using pointers 
    if(std::abs(invariant_mass - mass) > 1e-3) 
    {
      double adjusted_energy = std::sqrt(mass * mass + (px * px + py * py + pz * pz)); 
      momentum = std::make_unique<FourMomentum>(adjusted_energy, px, py, pz); 
      std::cout<<"The invariant mass is not equal to the rest mass of the particle. Calculated mass: "
      <<invariant_mass<<" MeV/c^2."<<" The energy in the four-momentum vector has been corrected to: "
      <<adjusted_energy<<" MeV."<<std::endl;
    }
    std::string start_of_string = "anti";
    if((anti_particle==true && (colour_charge.compare(0, 4, start_of_string) != 0)) || //comparing expected and required strings
    (anti_particle == false && (colour_charge.compare(0, 4, start_of_string) == 0))) //strings are equal
    {
      if(anti_particle == true)  
        std::cout<<"Error. Antiquarks should carry anticolour charge.";
      else if(anti_particle == false)
        std::cout<<"Error. Quarks must carry colour charge.";
    }
    baryon_number = anti_particle ? -1.0/3: 1.0/3;
  }
  //setters and getters for main physical quantities 
  double get_spin() const {return spin;}
  double get_charge() const override {return charge;}
  double get_baryon_number() const override {return baryon_number;}
  std::string get_colour_charge() const {return colour_charge;} 
  std::string get_particle_type() const override {return "Quark";}   //useful for catalogue class (currently working on)
  std::string get_quark_flavour() const override {return quark_flavour;}   
  void set_energy(double energy) {momentum->set_energy(energy);}
  void set_momentum_x(double px) {momentum->set_momentum_x(px);}
  void set_momentum_y(double py) {momentum->set_momentum_y(py);}
  void set_momentum_z(double pz) {momentum->set_momentum_z(pz);}
  double get_energy() const {return momentum->get_energy();}
  double get_momentum_x() const {return momentum->get_momentum_x();}
  double get_momentum_y() const {return momentum->get_momentum_y();}
  double get_momentum_z() const {return momentum->get_momentum_z();}
  void print_particle_data() const override
  {
    std::cout<<"Quark - Mass: "<<mass << " MeV/c^2, Spin: "<<spin 
    <<", Charge: "<<get_charge()<<", Baryon Number: "<<get_baryon_number()
    <<", 4-Momentum: ("<<momentum->get_energy()<<", "<<momentum->get_momentum_x() 
    <<", "<<momentum->get_momentum_y()<<", "<<momentum->get_momentum_z()<<")" << std::endl;
    std::cout<<std::endl;
  }
};
