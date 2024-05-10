/*Description: This code provides a catalogue of all particles in the Standard Model of Physics. 
It can be used to view physical properties of particles, as well as their associated quantum numbers. 
The particles are split into separate classes, and some instances of particle decays are shown. 
Author: Hassan Hashmi
Date: 10/05/2024*/

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
  std::unique_ptr<std::vector<double>> momentum_vector;   
public:
  FourMomentum(double input_energy, double input_momentum_x, double input_momentum_y, double input_momentum_z) : momentum_vector(std::make_unique<std::vector<double>>(4))
  {
    set_energy(input_energy);
    set_momentum_x(input_momentum_x);
    set_momentum_y(input_momentum_y);
    set_momentum_z(input_momentum_z);
  }
  void set_energy(double input_energy)   
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
    std::vector<double> total_vector(4);                               
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
  double dot_product(const FourMomentum &particle)                
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

class Particle
{
protected:
  double mass;
  std::unique_ptr<FourMomentum> momentum;
public:
  Particle(double mass, double energy, double px, double py, double pz) :
    mass{mass}, momentum{std::make_unique<FourMomentum>(energy, px, py, pz)} {}
  virtual void print_particle_data() const = 0;
  virtual int get_lepton_number() const {return 0;}    
  virtual double get_baryon_number() const {return 0;} 
  virtual double get_charge() const {return 0;}
  virtual std::string get_particle_type() const = 0;
  virtual std::string get_quark_flavour() const {return 0;}
  virtual double get_energy() {return momentum->get_energy();}
  FourMomentum &get_momentum() {return *momentum;}
  virtual ~Particle() {}
};

class Lepton:public Particle
{
protected:
  double spin;
  bool anti_particle; //true if antiparticle present 
  int charge;   

public:
  Lepton(double mass, double input_spin, double energy, double px, double py, double pz, int input_charge, bool anti_particle = false): 
    Particle(mass, energy, px, py, pz), spin{input_spin}, anti_particle{anti_particle}, charge{anti_particle ? input_charge: input_charge *(-1)} 
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
  double get_energy() const {return momentum->get_energy();}
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
    double invariant_mass = momentum->compute_invariant_mass();
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
  double get_spin() const {return spin;}
  double get_charge() const override {return charge;}
  double get_baryon_number() const override {return baryon_number;}
  std::string get_colour_charge() const {return colour_charge;} 
  std::string get_particle_type() const override {return "Quark";}
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

class Boson:public Particle
{
protected:
  double spin;
  bool anti_particle;  
  double charge; 
  std::string gauge_boson; 

public:
  Boson(double mass, double input_spin, double energy, double px, double py, double pz, const std::string &input_boson_type, bool anti_particle = false):
    Particle(mass, energy, px, py, pz), spin{input_spin}, anti_particle{anti_particle}, gauge_boson{input_boson_type}
  {
    double invariant_mass = momentum->compute_invariant_mass();
    if(std::abs(invariant_mass - mass) > 1e-3) 
    {
      double adjusted_energy = std::sqrt(mass *mass + (px * px + py * py + pz * pz)); 
      momentum = std::make_unique<FourMomentum>(adjusted_energy, px, py, pz); 
      std::cout<<"The invariant mass is not equal to the rest mass of the particle. Calculated mass: "
      <<invariant_mass<<" MeV/c^2."<<" The energy in the four-momentum vector has been corrected to: "
      <<adjusted_energy<<" MeV."<<std::endl;
    }
    if(input_boson_type == "W+" || input_boson_type == "W-")
      charge = (anti_particle ? -1: 1);
    else
      charge = 0; 
  }
  double get_spin() const {return spin;}
  double get_charge() const override {return charge;}
  std::string get_particle_type() const override {return "Boson";}
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
    std::cout<<"Gauge Boson - Mass: "<<mass<<" MeV/c^2 Spin: "<<spin
    << " Charge: "<<get_charge()<<", 4-Momentum: "<<momentum->get_energy() 
    <<", "<<momentum->get_momentum_x() << ", " << momentum->get_momentum_y() 
    <<", "<< momentum->get_momentum_z()<< std::endl;
  }
};

class HiggsBoson:public Particle
{
private:
  std::vector<std::unique_ptr<Particle>> higgs_decay_products;
public:
  HiggsBoson(double mass, double energy, double px, double py, double pz): 
    Particle(mass, energy, px, py, pz)
  {
    double invariant_mass = momentum->compute_invariant_mass();
    if(std::abs(invariant_mass - mass) > 1e-3) //allow a small tolerance 
    {
      double adjusted_energy = std::sqrt(mass * mass + (px * px + py * py + pz * pz)); 
      momentum = std::make_unique<FourMomentum>(adjusted_energy, px, py, pz); 
      std::cout<<"The invariant mass is not equal to the rest mass of the particle. Calculated mass: "
      <<invariant_mass<<" MeV/c^2."<<" The energy in the four-momentum vector has been corrected to: "
      <<adjusted_energy<<" MeV."<<std::endl;
    }
  }
  void add_higgs_decay_particles(std::vector<std::unique_ptr<Particle>> &higgs_decay_particles)
  {
    for(auto &decay_particle: higgs_decay_particles)
    {
      higgs_decay_products.push_back(std::move(decay_particle));  
    }
    higgs_decay_particles.clear();
  }

  void check_decay_charges() const
  {
    double total_higgs_decay_charge = 0.0;
    for(const auto &particle: higgs_decay_products)
    {
      total_higgs_decay_charge += particle->get_charge();
    }
    if(total_higgs_decay_charge == 0) 
      std::cout<<"The decay product charges sum to 0 as required."<<std::endl;
    else
      std::cout<<"The decay charges do not sum to 0."<<std::endl;
  }
  std::string get_particle_type() const override {return "HiggsBoson";}
  double get_charge() const {return 0;}   //does not need overriding in this case 
  void set_energy(double energy) {momentum->set_energy(energy);}
  void set_momentum_x(double px) {momentum->set_momentum_x(px);}
  void set_momentum_y(double py) {momentum->set_momentum_y(py);}
  void set_momentum_z(double pz) {momentum->set_momentum_z(pz);}
  double get_energy()const {return momentum->get_energy();}
  double get_momentum_x() const {return momentum->get_momentum_x();}
  double get_momentum_y() const {return momentum->get_momentum_y();}
  double get_momentum_z() const {return momentum->get_momentum_z();}
  double get_spin() const {return 0;} 
  void print_particle_data() const override 
  {
    std::cout<<"Higgs Boson - Mass: "<<mass<<" MeV/c^2, Spin: 0, Charge: 0, "
    <<"4-Momentum: "<<momentum->get_energy()<<", " 
    <<momentum->get_momentum_x()<< ", " 
    <<momentum->get_momentum_y()<< ", " 
    <<momentum->get_momentum_z()<< std::endl;
    for(const auto &decay_product: higgs_decay_products)
    {
      decay_product->print_particle_data();
    }
  }
};

class Calorimeter
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
    if(std::abs(momentum->get_energy() - calorimeter.total_calorimeter_energy()) > 1e-3) 
      throw std::runtime_error("The 4-momentum electron energy is not equal to deposited energy in calorimeter."); 
  }
  catch (const std::runtime_error&) 
  {
    std::cout<<"Incorrect deposited electron energies.  "<<std::endl;
  }
  };
}; 

class Muon:public Lepton
{
private:
  bool muon_isolation;
public:
  Muon(double mass, double energy, double px, double py, double pz, bool anti_particle, bool input_muon_isolation):
    Lepton(mass, 105.7, energy, px, py, pz, -1, anti_particle), muon_isolation{input_muon_isolation} {}
  void set_isolation_status(bool isolation_status) {muon_isolation = isolation_status;}
  bool get_isolation_status() const {return muon_isolation;}
  std::string get_particle_type() const override {return "Muon";}
  double get_charge() const override {return anti_particle ? 1: -1;}
  double get_spin() const {return 0.5;}
  int get_lepton_number() const override {return anti_particle ? -1: 1;}
  void print_particle_data() const override
  {
    Lepton::print_particle_data();
    if(muon_isolation == true) 
      std::cout<<"The muon is isolated."<<std::endl;
    else 
      std::cout<<"The muon is not isolated."<<std::endl;
  }
};

class ElectronNeutrino:public Lepton
{
private:
  bool has_interacted;
public:
  ElectronNeutrino(double mass, double energy, double px, double py, double pz, bool anti_particle, bool interaction_status):
    Lepton(mass, 0.5, energy, px, py, pz, 0, anti_particle), has_interacted{interaction_status} {}
  bool get_interaction_status() const {return has_interacted;}
  std::string get_particle_type() const override {return "ElectronNeutrino";}
  double get_charge() const {return 0;}
  double get_spin() const {return 0.5;}
  int get_lepton_number() const override {return anti_particle ? -1: 1;}
  void print_particle_data() const override
  {
    Lepton::print_particle_data();
    if(has_interacted == true) 
      std::cout<<"Interaction status: true"<<"\n"<<std::endl;
    else 
      std::cout<<"Interaction status: false"<<"\n"<<std::endl;
  }
};

class MuonNeutrino:public Lepton
{
private:
  bool has_interacted;
public:
  MuonNeutrino(double mass, double energy, double px, double py, double pz, bool anti_particle, bool interaction_status):
    Lepton(mass, 0.5, energy, px, py, pz, 0, anti_particle), has_interacted{interaction_status} {}
  bool get_interaction_status() const {return has_interacted;}
  std::string get_particle_type() const override {return "MuonNeutrino";}
  double get_charge() const {return 0;}
  double get_spin() const {return 0.5;}
  int get_lepton_number() const override {return anti_particle ? -1: 1;}
  void print_particle_data() const override
  {
    Lepton::print_particle_data();
    if(has_interacted == true) 
      std::cout<<"Interaction status: true"<<"\n"<<std::endl;
    else 
      std::cout<<"Interaction status: false"<<"\n"<<std::endl;
  }
};

class TauNeutrino:public Lepton
{
private:
  bool has_interacted;
public:
  TauNeutrino(double mass, double energy, double px, double py, double pz, bool anti_particle, bool interaction_status):
    Lepton(mass, 0.5, energy, px, py, pz, 0, anti_particle), has_interacted{interaction_status} {}
  bool get_interaction_status() const {return has_interacted;}
  std::string get_particle_type() const override {return "TauNeutrino";}
  double get_charge() const {return 0;}
  double get_spin() const {return 0.5;}
  int get_lepton_number() const override {return anti_particle ? -1: 1;}
  void print_particle_data() const override
  {
    Lepton::print_particle_data();
    if(has_interacted == true) 
      std::cout<<"Interaction status: true"<<"\n"<<std::endl;
    else 
      std::cout<<"Interaction status: false"<<"\n"<<std::endl;
  }
};

class Photon:public Boson
{
public:
  Photon(double mass, double energy, double px, double py, double pz, bool anti_particle = false):
    Boson(mass, 1, energy, px, py, pz, "Photon", anti_particle) {} 
  std::string get_particle_type() const override {return "Photon";}
  double get_charge() const {return 0;}
  double get_spin() const {return 1;}
};

class Gluon:public Boson
{
private:
  std::string colour_charge;
  std::string anticolour_charge;
public:
  Gluon(double mass, double energy, double px, double py, double pz, bool anti_particle, const std::string &input_colour, const std::string &input_anticolour):
    Boson(mass, 1, energy, px, py, pz, "Gluon", anti_particle), colour_charge{input_colour}, anticolour_charge{input_anticolour}
  {                                                                  
    std::string allowed_colour_charges[3] = {"red", "blue", "green"};                 
    std::string allowed_anticolour_charges[3] = {"antired", "antiblue", "antigreen"};
    bool valid_colour_inputted = false;  
    bool valid_anticolour_inputted = false;
    for(int i = 0; i < 3; i++)                 //iterating through lists above to check colour charges
    {
      if(input_colour == allowed_colour_charges[i])
        valid_colour_inputted = true;
      if(input_anticolour == allowed_anticolour_charges[i])
        valid_anticolour_inputted = true;
    }
    if (valid_colour_inputted == false || valid_anticolour_inputted == false)
      std::cout<<"Incorrect colour charges inputted. Gluons must carry a colour charge and an anticolour charge."<<std::endl;
  }
  std::string get_particle_type() const override {return "Gluon";}
  double get_spin() const {return 1;}
  double get_charge() const {return 0;}
  void print_particle_data() const override
  {
    Boson::print_particle_data();
    std::cout<<"Colour Charge: "<<colour_charge<< std::endl;
    std::cout<<"Anticolour Charge: "<<anticolour_charge<<std::endl; 
  };
};

class ZBoson:public Boson
{
public:
  ZBoson(double mass, double energy, double px, double py, double pz, bool anti_particle = false):
    Boson(mass, 1, energy, px, py, pz, "Z", anti_particle) {}
  std::string get_particle_type() const override {return "ZBoson";}
  double get_spin() const {return 1;}
  double get_charge() const {return 0;}
};

class WPlusBoson:public Boson
{
public:
  WPlusBoson(double mass, double energy, double px, double py, double pz, bool anti_particle = false): 
    Boson(mass, 1, energy, px, py, pz, "W+", anti_particle) {} 
  std::string get_particle_type() const override {return "WPlusBoson";}
  double get_spin() const {return 1;}
  double get_charge() const override {return 1;}
};

class WMinusBoson:public Boson
{
public:
  WMinusBoson(double mass, double energy, double px, double py, double pz, bool anti_particle = false): 
    Boson(mass, 1, energy, px, py, pz, "W-", anti_particle) {}  
  std::string get_particle_type() const override {return "WMinusBoson";}
  double get_charge() const override {return -1;}
  double get_spin() const {return 1;}
};

class UpQuark:public Quark
{
public:
  UpQuark(double mass, double energy, double px, double py, double pz, bool anti_particle, const std::string &colour): 
    Quark(mass, 0.5, energy, px, py, pz, "up", anti_particle, colour) {}
  std::string get_particle_type() const override {return "UpQuark";}
  double get_charge() const override {return anti_particle ? -2/3: 2/3;}
  double get_spin() const {return 0.5;}
  double get_baryon_number() const override {return anti_particle ? -1/3: 1/3;}
};

class DownQuark:public Quark
{
public:
  DownQuark(double mass, double energy, double px, double py, double pz, bool anti_particle, const std::string &colour): 
    Quark(mass, 0.5, energy, px, py, pz, "down", anti_particle, colour) {}
  std::string get_particle_type() const override {return "DownQuark";}
  double get_charge() const override {return anti_particle ? 1/3: -1/3;}
  double get_spin() const {return 0.5;}
  double get_baryon_number() const override { return anti_particle ? -1/3: 1/3; }
};

class Charm:public Quark
{
public:
  Charm(double mass, double energy, double px, double py, double pz, bool anti_particle, const std::string &colour):
    Quark(mass, 0.5, energy, px, py, pz, "charm", anti_particle, colour) {}
  std::string get_particle_type() const override {return "CharmQuark";}
  double get_charge() const override {return anti_particle ? -2/3: 2/3;}
  double get_baryon_number() const override {return anti_particle ? -1/3: 1/3;}
  double get_spin() const {return 0.5;}
};

class Strange:public Quark
{
public:
  Strange(double mass, double energy, double px, double py, double pz, bool anti_particle, const std::string &colour):
    Quark(mass, 0.5, energy, px, py, pz, "strange", anti_particle, colour) {}
  std::string get_particle_type() const override {return "StrangeQuark";}
  double get_charge() const override {return anti_particle ? 1/3: -1/3;}
  double get_baryon_number() const override {return anti_particle ? -1/3: 1/3;}
  double get_spin() const {return 0.5;}
};

class Top:public Quark
{
public:
  Top(double mass, double energy, double px, double py, double pz, bool anti_particle, const std::string &colour): 
    Quark(mass, 0.5, energy, px, py, pz, "top", anti_particle, colour) {}
  std::string get_particle_type() const override {return "TopQuark";}
  double get_charge() const override {return anti_particle ? -2/3: 2/3;}
  double get_baryon_number() const override {return anti_particle ? -1/3: 1/3;}
  double get_spin() const {return 0.5;}
};

class Bottom:public Quark
{
public:
  Bottom(double mass, double energy, double px, double py, double pz, bool anti_particle, const std::string &colour):
    Quark(mass, 0.5, energy, px, py, pz, "bottom", anti_particle, colour) {}
  std::string get_particle_type() const override {return "BottomQuark";}
  double get_charge() const override {return anti_particle ? 1/3: -1/3;}
  double get_baryon_number() const override {return anti_particle ? -1/3: 1/3;}
  double get_spin() const {return 0.5;}
};

class Tau:public Lepton
{
private:
  std::vector<std::unique_ptr<Particle>> tau_decay_products;
public:
  Tau(double mass, double energy, double px, double py, double pz, bool anti_particle):
    Lepton(mass, 0.5, energy, px, py, pz, anti_particle ? 1: -1, anti_particle)
  {add_tau_decay_products();}

  void add_tau_decay_products()
  {
  tau_decay_products.clear(); 
  if(anti_particle==false)
  {
    tau_decay_products.push_back(std::make_unique<UpQuark>(2.3, 0, 0, 0, 0, false, "red"));
    tau_decay_products.push_back(std::make_unique<DownQuark>(4.8, 0, 0, 0, 0, true, "antiblue"));
    tau_decay_products.push_back(std::make_unique<TauNeutrino>(0, 0, 0, 0, 0, false, false));
  }
  else
  {
    tau_decay_products.push_back(std::make_unique<Electron>(0.511, 0, 0, 0, 0, true));
    tau_decay_products.push_back(std::make_unique<ElectronNeutrino>(0, 0, 0, 0, 0, false, false));
    tau_decay_products.push_back(std::make_unique<TauNeutrino>(0, 0, 0, 0, 0, true, false));
  }
  check_tau_decay_total();
  }
  void check_tau_decay_total()
  {
    double total_charge = 0;
    for(const auto &decay_particle: tau_decay_products)
    {
      total_charge  += decay_particle->get_charge();
    }
    if(std::abs(total_charge)!= 1)
      std::cout<<"Tau decay products do not sum to expected value."<<std::endl; 
  }
  std::string get_particle_type() const override {return "Tau";}
  double get_charge() const override {return anti_particle ? 1 :-1;}
  int get_lepton_number() const override {return anti_particle ? -1 : 1;}
  double get_spin() const {return 0.5;}
  void print_particle_data() const override
  {
    Lepton::print_particle_data();
    std::cout<<"The Tau particle decayed into the following:"<<std::endl;
    std::cout<<std::endl; 
    for(const auto &decay_particle: tau_decay_products)
    {
      decay_particle->print_particle_data();
    }
  }
};

class Catalogue
{
private:
  std::list<std::shared_ptr<Particle>> catalogue_particles; //wrap particle catalogue in an STL list 
  int particle_count = 0, higgs_boson_count = 0, electron_count = 0, muon_count = 0, tau_count = 0,
  electron_neutrino_count = 0, muon_neutrino_count = 0, tau_neutrino_count = 0, up_quark_count = 0,
  down_quark_count = 0, charm_quark_count = 0, strange_quark_count = 0, top_quark_count = 0,
  bottom_quark_count = 0, gluon_count = 0, photon_count = 0, w_minus_count = 0, w_plus_count = 0,
  z_boson_count = 0;
public:
  void add_particle_to_catalogue(std::unique_ptr<Particle> particle) 
  {
    const std::string &type_of_particle = particle->get_particle_type();
    if(type_of_particle == "HiggsBoson")
      higgs_boson_count++;
    else if(type_of_particle == "Electron")
      electron_count++;
    else if(type_of_particle == "Muon")
      muon_count++;
    else if(type_of_particle == "Tau")
      tau_count++;
    else if(type_of_particle == "ElectronNeutrino")
      electron_neutrino_count++;
    else if(type_of_particle == "MuonNeutrino")
      muon_neutrino_count++;
    else if(type_of_particle == "TauNeutrino")
      tau_neutrino_count++;
    else if(type_of_particle == "UpQuark")
      up_quark_count++;
    else if(type_of_particle == "DownQuark")
      down_quark_count++;
    else if(type_of_particle == "CharmQuark")
      charm_quark_count++;
    else if(type_of_particle == "StrangeQuark")
      strange_quark_count++;
    else if(type_of_particle == "TopQuark")
      top_quark_count++;
    else if(type_of_particle == "BottomQuark")
      bottom_quark_count++;
    else if(type_of_particle == "Photon")
      photon_count++;
    else if(type_of_particle == "Gluon")
      gluon_count++;
    else if(type_of_particle == "WMinusBoson")
      w_minus_count++;
    else if(type_of_particle == "WPlusBoson")
      w_plus_count++;
    else if(type_of_particle == "ZBoson")
      z_boson_count++;
    catalogue_particles.push_back(std::move(particle));
    particle_count++;
    }
  int total_number_catalogue_particles() const {return particle_count;}
  void particles_sorted_by_descending_energy() //lambda function to sort particles by energy
  {
    catalogue_particles.sort([](const std::shared_ptr<Particle> &particle_one, const std::shared_ptr<Particle> &particle_two) 
    {
      return particle_one->get_energy() > particle_two->get_energy(); 
    });
  }
  void print_catalogue_by_energy() const 
  {
    std::cout<<"Particles sorted by energy:"<<std::endl;
    for(const auto &catalogue_particle: catalogue_particles) 
    {
      std::cout<<"Particle type, Energy: "<<catalogue_particle->get_particle_type()<<", "
      <<catalogue_particle->get_energy()<<" MeV"<<std::endl;
    }
  }
  void print_all_catalogue_particles() const 
  {
    for(const auto &catalogue_particle: catalogue_particles)
    {
      catalogue_particle->print_particle_data();
      std::cout<<std::endl; 
    }
  }
  std::list<std::shared_ptr<Particle>> get_particles_of_same_type(const std::string &type_of_particle) const 
  {
    std::list<std::shared_ptr<Particle>> sub_container;   
    for(const auto &catalogue_particle:catalogue_particles)
    {
      if(catalogue_particle->get_particle_type() == type_of_particle)
        sub_container.push_back(catalogue_particle);
    }
    return sub_container;
  }
  std::vector<int> get_each_particle_count() const    
  {
    return 
    {
      higgs_boson_count, electron_count, muon_count, tau_count,
      electron_neutrino_count, muon_neutrino_count, tau_neutrino_count,
      up_quark_count, down_quark_count, charm_quark_count, strange_quark_count, 
      top_quark_count, bottom_quark_count, gluon_count, photon_count, 
      w_minus_count, w_plus_count, z_boson_count
    };
  }
  std::vector<double> combined_catalogue_momentum() const  
  {
    std::vector<double> combined_momentum = {0, 0, 0, 0}; 
    for(const auto &catalogue_particle: catalogue_particles)
    {
      const FourMomentum &momentum = catalogue_particle->get_momentum();
      combined_momentum[0] += momentum.get_energy();
      combined_momentum[1] += momentum.get_momentum_x();
      combined_momentum[2] += momentum.get_momentum_y();
      combined_momentum[3] += momentum.get_momentum_z();
    }
    return combined_momentum;
  }
};

template <typename particle_class_one, typename particle_class_two>
bool particle_energy_comparison(const particle_class_one &particle1, const particle_class_two &particle2)
{
  return particle1.get_energy() < particle2.get_energy();
}

int main()
{
  //Testing calorimeter 
  Electron electron(0.511, 1, 0, 0, 0, false); 
  electron.set_energy(2.0); 
  std::vector<double> calorimeter_energies = {0.3, 0.4, 0.5, 0.3};  
  electron.validate_deposited_energies(calorimeter_energies);
  //test for antiquark having a regular colour charge
  UpQuark anti_up(2.3, 0, 0, 0, 0, true, "red");  
  //W plus decay products
  std::vector<std::unique_ptr<Particle>> wplus_decay_products;
  wplus_decay_products.push_back(std::make_unique<Electron>(0.511, 1, 0, 0, 0, false));
  wplus_decay_products.push_back(std::make_unique<Electron>(0.511, 1, 0, 0, 0, true));
  wplus_decay_products.push_back(std::make_unique<ElectronNeutrino>(0, 0, 0, 0, 0, true, true));
  wplus_decay_products.push_back(std::make_unique<ElectronNeutrino>(0, 0, 0, 0, 0, false, false));
  wplus_decay_products.push_back(std::make_unique<Muon>(105.7, 1, 0, 0, 0, true, true));
  wplus_decay_products.push_back(std::make_unique<Muon>(105.7, 1, 0, 0, 0, false, false));
  wplus_decay_products.push_back(std::make_unique<MuonNeutrino>(0, 0, 0, 0, 0, true, true));
  wplus_decay_products.push_back(std::make_unique<MuonNeutrino>(0, 0, 0, 0, 0, false, false));
  wplus_decay_products.push_back(std::make_unique<Tau>(1777, 1, 0, 0, 0, false));
  wplus_decay_products.push_back(std::make_unique<Tau>(1777, 1, 0, 0, 0, true));
  wplus_decay_products.push_back(std::make_unique<TauNeutrino>(0, 0, 0, 0, 0, true, true));
  wplus_decay_products.push_back(std::make_unique<TauNeutrino>(0, 0, 0, 0, 0, false, false));
  wplus_decay_products.push_back(std::make_unique<UpQuark>(2.3, 0, 0, 0, 0, false, "red"));
  wplus_decay_products.push_back(std::make_unique<DownQuark>(4.8, 0, 0, 0, 0, false, "blue"));
  wplus_decay_products.push_back(std::make_unique<Charm>(1275, 0, 0, 0, 0, false, "green"));
  wplus_decay_products.push_back(std::make_unique<Strange>(95, 0, 0, 0, 0, false, "red"));
  wplus_decay_products.push_back(std::make_unique<Top>(173000, 0, 0, 0, 0, false, "red"));
  wplus_decay_products.push_back(std::make_unique<Bottom>(4180, 0, 0, 0, 0, false, "blue"));
  //W minus decay products
  std::vector<std::unique_ptr<Particle>> wminus_decay_products;
  wminus_decay_products.push_back(std::make_unique<Electron>(0.511, 1, 0, 0, 0, true));
  wminus_decay_products.push_back(std::make_unique<Electron>(0.511, 1, 0, 0, 0, false));
  wminus_decay_products.push_back(std::make_unique<ElectronNeutrino>(0, 0, 0, 0, 0, false, false));
  wminus_decay_products.push_back(std::make_unique<ElectronNeutrino>(0, 0, 0, 0, 0, true, true));
  wminus_decay_products.push_back(std::make_unique<Muon>(105.7, 1, 0, 0, 0, false, false));
  wminus_decay_products.push_back(std::make_unique<Muon>(105.7, 1, 0, 0, 0, true, true));
  wminus_decay_products.push_back(std::make_unique<MuonNeutrino>(0, 0, 0, 0, 0, false, false));
  wminus_decay_products.push_back(std::make_unique<MuonNeutrino>(0, 0, 0, 0, 0, true, true));
  wminus_decay_products.push_back(std::make_unique<Tau>(1777, 1, 0, 0, 0, true));
  wminus_decay_products.push_back(std::make_unique<Tau>(1777, 1, 0, 0, 0, false));
  wminus_decay_products.push_back(std::make_unique<TauNeutrino>(0, 0, 0, 0, 0, false, false));
  wminus_decay_products.push_back(std::make_unique<TauNeutrino>(0, 0, 0, 0, 0, true, true));
  wminus_decay_products.push_back(std::make_unique<UpQuark>(2.3, 0, 0, 0, 0, true, "antired"));
  wminus_decay_products.push_back(std::make_unique<DownQuark>(4.8, 0, 0, 0, 0, true, "antiblue"));
  wminus_decay_products.push_back(std::make_unique<Charm>(1275, 0, 0, 0, 0, true, "antigreen"));
  wminus_decay_products.push_back(std::make_unique<Strange>(95, 0, 0, 0, 0, true, "antired"));
  wminus_decay_products.push_back(std::make_unique<Top>(173000, 0, 0, 0, 0, true, "antired"));
  wminus_decay_products.push_back(std::make_unique<Bottom>(4180, 0, 0, 0, 0, true, "antiblue"));
  //Z boson decay products
  std::vector<std::unique_ptr<Particle>> zboson_decay_products;
  zboson_decay_products.push_back(std::make_unique<Electron>(0.511, 1, 0, 0, 0, false));
  zboson_decay_products.push_back(std::make_unique<Electron>(0.511, 1, 0, 0, 0, true));
  zboson_decay_products.push_back(std::make_unique<ElectronNeutrino>(0, 0, 0, 0, 0, true, true));
  zboson_decay_products.push_back(std::make_unique<ElectronNeutrino>(0, 0, 0, 0, 0, false, false));
  zboson_decay_products.push_back(std::make_unique<Muon>(105.7, 1, 0, 0, 0, true, true));
  zboson_decay_products.push_back(std::make_unique<Muon>(105.7, 1, 0, 0, 0, false, false));
  zboson_decay_products.push_back(std::make_unique<MuonNeutrino>(0, 0, 0, 0, 0, true, true));
  zboson_decay_products.push_back(std::make_unique<MuonNeutrino>(0, 0, 0, 0, 0, false, false));
  zboson_decay_products.push_back(std::make_unique<Tau>(1777, 1, 0, 0, 0, false));
  zboson_decay_products.push_back(std::make_unique<Tau>(1777, 1, 0, 0, 0, true));
  zboson_decay_products.push_back(std::make_unique<TauNeutrino>(0, 0, 0, 0, 0, true, true));
  zboson_decay_products.push_back(std::make_unique<TauNeutrino>(0, 0, 0, 0, 0, false, false));
  zboson_decay_products.push_back(std::make_unique<UpQuark>(2.3, 0, 0, 0, 0, true, "antired"));
  zboson_decay_products.push_back(std::make_unique<UpQuark>(2.3, 0, 0, 0, 0, false, "red"));
  zboson_decay_products.push_back(std::make_unique<DownQuark>(4.8, 0, 0, 0, 0, true, "antiblue"));
  zboson_decay_products.push_back(std::make_unique<DownQuark>(4.8, 0, 0, 0, 0, false, "blue"));
  zboson_decay_products.push_back(std::make_unique<Charm>(1275, 0, 0, 0, 0, true, "antigreen"));
  zboson_decay_products.push_back(std::make_unique<Charm>(1275, 0, 0, 0, 0, false, "green"));
  zboson_decay_products.push_back(std::make_unique<Strange>(95, 0, 0, 0, 0, true, "antired"));
  zboson_decay_products.push_back(std::make_unique<Strange>(95, 0, 0, 0, 0, false, "red"));
  zboson_decay_products.push_back(std::make_unique<Top>(173000, 0, 0, 0, 0, true, "antired"));
  zboson_decay_products.push_back(std::make_unique<Top>(173000, 0, 0, 0, 0, false, "red"));
  zboson_decay_products.push_back(std::make_unique<Bottom>(4180, 0, 0, 0, 0, true, "antiblue"));
  zboson_decay_products.push_back(std::make_unique<Bottom>(4180, 0, 0, 0, 0, false, "blue"));
  //higgsboson decay testing 
  HiggsBoson higgs(125.0, 100, 0, 0, 0);
  std::vector<std::unique_ptr<Particle>> higgs_decay_particles;
  higgs_decay_particles.push_back(std::make_unique<ZBoson>(91.2, 0.5, 150, 0, 0, false));
  higgs_decay_particles.push_back(std::make_unique<ZBoson>(91.2, 0.5, 50, 0, 0, false));
  higgs_decay_particles.push_back(std::make_unique<WPlusBoson>(80.4, 0.5, 50, 0, 0, false));
  higgs_decay_particles.push_back(std::make_unique<WMinusBoson>(80.4, 0.5, 100, 0, 0, false));
  higgs_decay_particles.push_back(std::make_unique<Bottom>(4.18, 0.5, 50, 0, 0, false, "blue"));
  higgs_decay_particles.push_back(std::make_unique<Bottom>(4.18, 0.5, 50, 0, 0, true, "antiblue"));
  higgs_decay_particles.push_back(std::make_unique<Photon>(0, 0, 0, 0, 0, false));
  higgs_decay_particles.push_back(std::make_unique<Photon>(0, 0, 0, 0, 0, false));
  higgs.add_higgs_decay_particles(higgs_decay_particles); 
  higgs.check_decay_charges(); 
  higgs.print_particle_data();
  // testing invariant mass calculation
  std::cout<<"Checking invariant mass calculation for a test electron."<<std::endl;
  Electron test_electron(0.511, 0, 0, 0, 0, false);
  //invariant mass computed outside of parametrised constructors 
  std::cout<<"Accessing invariant mass function from main()."<<std::endl;
  double test_electron_invariant_mass = test_electron.get_momentum().compute_invariant_mass();
  // instantiating a container of all particles 
  Catalogue catalogue;
  catalogue.add_particle_to_catalogue(std::make_unique<Electron>(0.511, 1.0, 0.1, 0.1, 0.1, false));
  catalogue.add_particle_to_catalogue(std::make_unique<Electron>(0.511, 1.0, 0.1, 0.1, 0.1, true)); 
  catalogue.add_particle_to_catalogue(std::make_unique<Muon>(105.7, 10.0, 0.1, 0.1, 0.1, false, false));
  catalogue.add_particle_to_catalogue(std::make_unique<Muon>(105.7, 10.0, 0.1, 0.1, 0.1, true, false));
  catalogue.add_particle_to_catalogue(std::make_unique<Tau>(1777, 20.0, 0.2, 0.2, 0.2, false));
  catalogue.add_particle_to_catalogue(std::make_unique<Tau>(1777, 20.0, 0.2, 0.2, 0.2, true)); 
  catalogue.add_particle_to_catalogue(std::make_unique<ElectronNeutrino>(0, 0.5, 0.05, 0.05, 0.05, false, false));
  catalogue.add_particle_to_catalogue(std::make_unique<ElectronNeutrino>(0, 0.5, 0.05, 0.05, 0.05, true, false));
  catalogue.add_particle_to_catalogue(std::make_unique<MuonNeutrino>(0, 0.5, 0.05, 0.05, 0.05, false, false));
  catalogue.add_particle_to_catalogue(std::make_unique<MuonNeutrino>(0, 0.5, 0.05, 0.05, 0.05, true, false));
  catalogue.add_particle_to_catalogue(std::make_unique<TauNeutrino>(0, 0.5, 0.05, 0.05, 0.05, false, false));
  catalogue.add_particle_to_catalogue(std::make_unique<TauNeutrino>(0, 0.5, 0.05, 0.05, 0.05, true, false));
  catalogue.add_particle_to_catalogue(std::make_unique<UpQuark>(2.3, 0.3, 0.03, 0.03, 0.03, false, "red"));
  catalogue.add_particle_to_catalogue(std::make_unique<UpQuark>(2.3, 0.3, 0.03, 0.03, 0.03, true, "antired")); 
  catalogue.add_particle_to_catalogue(std::make_unique<DownQuark>(4.8, 0.4, 0.04, 0.04, 0.04, false, "blue"));
  catalogue.add_particle_to_catalogue(std::make_unique<DownQuark>(4.8, 0.4, 0.04, 0.04, 0.04, true, "antiblue")); 
  catalogue.add_particle_to_catalogue(std::make_unique<Charm>(1275, 1.5, 0.15, 0.15, 0.15, false, "green"));
  catalogue.add_particle_to_catalogue(std::make_unique<Charm>(1275, 1.5, 0.15, 0.15, 0.15, true, "antigreen")); 
  catalogue.add_particle_to_catalogue(std::make_unique<Strange>(95, 0.95, 0.095, 0.095, 0.095, false, "red"));
  catalogue.add_particle_to_catalogue(std::make_unique<Strange>(95, 0.95, 0.095, 0.095, 0.095, true, "antired")); 
  catalogue.add_particle_to_catalogue(std::make_unique<Top>(173000, 2.0, 0.2, 0.2, 0.2, false, "red"));
  catalogue.add_particle_to_catalogue(std::make_unique<Top>(173000, 2.0, 0.2, 0.2, 0.2, true, "antired"));
  catalogue.add_particle_to_catalogue(std::make_unique<Bottom>(4180, 1.0, 0.1, 0.1, 0.1, false, "blue"));
  catalogue.add_particle_to_catalogue(std::make_unique<Bottom>(4180, 1.0, 0.1, 0.1, 0.1, true, "antiblue")); 
  catalogue.add_particle_to_catalogue(std::make_unique<Photon>(0, 1.0, 0.1, 0.1, 0.1, false));
  catalogue.add_particle_to_catalogue(std::make_unique<Gluon>(0, 1.0, 0.1, 0.1, 0.1, false, "red", "antired"));
  catalogue.add_particle_to_catalogue(std::make_unique<Gluon>(0, 1.0, 0.1, 0.1, 0.1, true, "green", "antigreen"));
  catalogue.add_particle_to_catalogue(std::make_unique<WPlusBoson>(80.4, 5.0, 0.5, 0.5, 0.5, false));
  catalogue.add_particle_to_catalogue(std::make_unique<WMinusBoson>(80.4, 5.0, 0.5, 0.5, 0.5, false));
  catalogue.add_particle_to_catalogue(std::make_unique<ZBoson>(91.2, 6.0, 0.6, 0.6, 0.6, false));
  catalogue.add_particle_to_catalogue(std::make_unique<HiggsBoson>(126, 7.0, 0.7, 0.7, 0.7));
  //particle report 
  std::cout<<std::endl; 
  std::cout<<std::endl; 
  std::cout<<"Report on particle information for catalogue particles - "<<std::endl; 
  std::cout<<std::endl; 
  catalogue.print_all_catalogue_particles();
  std::cout<<"Total number of particles in the container: "<<catalogue.total_number_catalogue_particles()<<std::endl;
  std::vector<int> individual_particle_count = catalogue.get_each_particle_count();
  const char* particle_names[] =      
  {
    "higgs bosons: ",            //same as catalogue order 
    "electrons: ", "muons: ", "taus: ",
    "electron neutrinos: ", "muon neutrinos: ", "tau neutrinos: ",
    "up quarks: ", "down quarks: ", "charm quarks: ", "strange quarks: ", "top quarks: ", "bottom quarks: ",
    "gluons: ", "photons: ", "w minus bosons: ", "w plus bosons: ", "z bosons: "
  };
  std::cout<<"Numbers of each particle in the container:"<<std::endl;
  for(int i = 0; i < individual_particle_count.size(); ++i)
  {
    std::cout<<particle_names[i]<<individual_particle_count[i]<<std::endl;
  }
  // 4-momentum in container 
  std::vector<double> total_momentum = catalogue.combined_catalogue_momentum();
  std::cout<<"Total 4-Momentum: Energy = "<<total_momentum[0]<<std::endl;
  std::cout<<"Momentum in the x-direction ="<<total_momentum[1]<<std::endl;
  std::cout<<"Momentum in the y-direction = "<<total_momentum[2]<<std::endl;
  std::cout<<"Momentum in the z-direction = "<<total_momentum[3]<<std::endl;
  //testing template 
  Muon muon(105.7, 25, 0, 0, 0, false, false);   
  bool carries_less_energy = particle_energy_comparison(electron, muon);
  if(carries_less_energy == true)
    std::cout<<"The electron has less energy than the muon."<<std::endl;
  else
    std::cout<<"The electron has energy greater than or equal to the muon."<<std::endl;
  //testing lambda function 
  catalogue.particles_sorted_by_descending_energy();
  catalogue.print_catalogue_by_energy();
  //instantiating particles with the wrong characteristics  
  Muon muon_test(105.7, 1, 0, 0, 0, true, true); 
  Gluon gluon(0, 1.0, 0.1, 0.1, 0.1, false, "antired", "antired");
  return 0;
}