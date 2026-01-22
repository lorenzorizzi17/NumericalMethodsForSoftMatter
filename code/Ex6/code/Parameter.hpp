#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

class Parameters {
public:
    // Spatial length are measured in units of sigma (diameter of the molecules!)
    double param_rho;      // Number density
    double param_N;        // Number of particles
    double param_temperature; // Temperature of the system    

    // Verlet algorithm parameters
    double param_r_cut; // Cutoff radius
    double param_skin;  // Skin for Verlet lists
    int param_M;        // Number of cells per dimension for cell list
    double param_max_displacement; // Maximum displacement for MC moves

    // Dynamic parameters. Energies are rescaled with respect to epsilon_SW = 1, k_B = 1
    int param_num_MC_steps; // Number of MC steps to perform
    double param_max_cos_theta; // Maximum cosine of the angle for patch interactions
    double param_energyWell; // Depth of the energy well for patch interactions
    double param_gamma; // Gamma factor

    void readFromFile(std::string const& filename);
};



void Parameters::readFromFile(std::string const& filename){
    std::ifstream infile(filename);
    if(!infile.is_open()){
        throw std::runtime_error("Could not open parameter file: " + filename);
    }
    // Starts parsing the parameter configuration yaml file. Simple key: value pairs
    std::string line;
    while(std::getline(infile, line)){
        std::istringstream iss(line);
        std::string key;
        if(std::getline(iss, key, ':')){
            std::string value_str;
            if(std::getline(iss, value_str)){
                value_str.erase(0, value_str.find_first_not_of(" \t"));
                if(key == "rho"){
                    param_rho = std::stod(value_str);
                } else if(key == "N"){
                    param_N = std::stoi(value_str);
                } else if(key == "r_cut"){
                    param_r_cut = std::stod(value_str);
                } else if(key == "skin"){
                    param_skin = std::stod(value_str);
                } else if(key == "M"){
                    param_M = std::stoi(value_str);
                } else if(key == "max_displacement"){  
                    param_max_displacement = std::stod(value_str);
                } else if(key == "num_MC_steps"){
                    param_num_MC_steps = std::stoi(value_str);
                } else if(key == "temperature"){
                    param_temperature = std::stod(value_str);
                } else if(key=="max_cos_theta"){
                    param_max_cos_theta=std::stod(value_str);
                } else if(key=="energyWell"){
                    param_energyWell=std::stod(value_str);
                } else if(key=="gamma"){
                    param_gamma=std::stod(value_str);
                }
            }
        }
    }
    infile.close();
    // Sanity check on the parameters
    double L = std::pow(static_cast<double>(param_N) / param_rho, 1.0 / 3.0);
    double ax = int(L/(param_r_cut+param_skin));
    if(ax < param_M){
        std::cerr << "WARNING: M is too large" << std::endl;
    }
    if (param_max_displacement > 0.5 * param_skin){
        std::cerr << "WARNING: max_displacement is larger than half the skin" << std::endl;
    }
}