#include "ActiveIsingModel.hpp"
#include<iostream>
#include<fstream>

int main(int argc, char* argv[]){
    // Physical parameters, to be consistent with the saved file
    bool flux = false;
    double beta;
    float rho;
    double epsilon;
    int totalSteps;
    int saveEvery;
    std::string prestring;
    if(argc == 6){
        beta = std::stod(argv[1]);
        rho = std::stod(argv[2]);
        epsilon = std::stod(argv[3]);
        totalSteps = std::stoi(argv[4]);
        saveEvery = std::stoi(argv[5]);
    } else if (argc == 7) {
        beta = std::stod(argv[1]);
        rho = std::stod(argv[2]);
        epsilon = std::stod(argv[3]);
        totalSteps = std::stoi(argv[4]);
        saveEvery = std::stoi(argv[5]);
        flux = true;
        prestring = "f";
    } else {
        std::cout << "Usage: ./player.out <beta> <rho> <epsilon> or ./player.out <beta> <rho> <epsilon> f \n";
        return -1;
    }


    int N_CELL_X = 200;
    int N_CELL_Y = 50;
    // Simulation parameters



    // Open file
    std::string filename = "../AIMData/Nx" + std::to_string(N_CELL_X) + 
                            "_Ny" + std::to_string(N_CELL_Y) + 
                            "_beta" + std::to_string(beta).substr(0,5) + 
                            "_eps" + std::to_string(epsilon).substr(0,5) + 
                            "_rho" + std::to_string(rho).substr(0,5) +
                            "_Ntot" + std::to_string(totalSteps) +
                            "_Nevry" + std::to_string(saveEvery) +
                            prestring + ".dat";
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error{"Couldn't find the file you were looking for"};
    };


    int* mag_buffer = new int[N_CELL_X*N_CELL_Y];
    int* den_buffer = new int[N_CELL_X*N_CELL_Y];
    int* flux_buffer = nullptr;
    if (flux) {
        flux_buffer = new int[N_CELL_X];
    }
    
    
    int burnin = 300; // before 200
    std::vector<double> densityProfile(N_CELL_X, 0.0);
    std::vector<double> magnetizationProfile(N_CELL_X, 0.0);
    double sigma_flux = 0.0;
    double mag_flux = 0.0;
    int current_frame = 0;
    while (current_frame < totalSteps / saveEvery) {

        file.read(reinterpret_cast<char*>(mag_buffer), N_CELL_X * N_CELL_Y * sizeof(int));
        file.read(reinterpret_cast<char*>(den_buffer), N_CELL_X * N_CELL_Y * sizeof(int));
        if (flux) {
            file.read(reinterpret_cast<char*>(flux_buffer), N_CELL_X * sizeof(int));
        }

        double rolling_avg_mag = 0;
        double rolling_avg_den = 0;
        // Start drawing the spins
        for (int col = 0; col < N_CELL_X; ++col) {

            double averageMag = std::accumulate(mag_buffer + N_CELL_Y*col, mag_buffer + N_CELL_Y*(col+1), 0.f) / float(N_CELL_Y);
            double averageDen = std::accumulate(den_buffer+N_CELL_Y*col, den_buffer+N_CELL_Y*(col+1), 0.f)/ float(N_CELL_Y);

            if (current_frame > burnin){
                densityProfile[col] += averageDen;
                magnetizationProfile[col] += averageMag;

                double flux_avg = std::accumulate(flux_buffer, flux_buffer + N_CELL_Y, 0.0) / N_CELL_Y;
                double sq_sum = std::inner_product(flux_buffer, flux_buffer + N_CELL_Y, flux_buffer, 0.0);
                double stdev = std::sqrt(sq_sum / N_CELL_Y - flux_avg * flux_avg);

                sigma_flux += stdev / double(totalSteps / saveEvery - burnin);
            }

        }

        if (current_frame > burnin){
            double current_avg_mag = std::accumulate(magnetizationProfile.begin(), magnetizationProfile.end(), 0.0) / N_CELL_X;
            double current_square_sum_mag = std::inner_product(magnetizationProfile.begin(), magnetizationProfile.end(), magnetizationProfile.begin(), 0.0);
            double current_stddev_mag = std::sqrt(current_square_sum_mag / N_CELL_X - current_avg_mag * current_avg_mag);
            mag_flux += current_stddev_mag / double(totalSteps / saveEvery - burnin);
        }

        current_frame++;
    }
    file.close();
    double final_avg_mag = std::accumulate(magnetizationProfile.begin(), magnetizationProfile.end(), 0.0) / N_CELL_X;
    double final_avg_den = std::accumulate(densityProfile.begin(), densityProfile.end(), 0.0) / N_CELL_X;
    //std::cout << "\n Final average magnetization per column: " << final_avg_mag << "\n";
    //std::cout << " Final average density per column: " << final_avg_den << "\n";
    std::cout << epsilon << " " << rho << " " << sigma_flux << "\n";
    //std::cout << epsilon << " " << rho << " " << std::abs(final_avg_mag) / final_avg_den << "\n";
    //std::cout << " Sigma flux: " << sigma_flux << "\n";
    //std::cout << "End frames\n";
}