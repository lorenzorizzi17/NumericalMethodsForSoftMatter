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
    // Graphical parameters
    int SizeHorizontal = 1200;
    int SizeVertical = 600;
    int sizeWindowPlot = 600;
    int delta = 40;
    double y_max = 15;
    
    double ay = (SizeVertical-2*delta)/float(3*N_CELL_Y);
    double ax = (SizeHorizontal)/float(N_CELL_X);
    double spacing = float(sizeWindowPlot) / float(N_CELL_X);
    int deltay = 10;

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

    // Launch the windows
    sf::RenderWindow window(sf::VideoMode(SizeHorizontal, SizeVertical), "Replay Simulation");
    sf::RenderWindow window_plot(sf::VideoMode(sizeWindowPlot, sizeWindowPlot), "Plotter");
    window.setFramerateLimit(120);
    window_plot.setFramerateLimit(120);

    int* mag_buffer = new int[N_CELL_X*N_CELL_Y];
    int* den_buffer = new int[N_CELL_X*N_CELL_Y];
    int* flux_buffer = nullptr;
    if (flux) {
        flux_buffer = new int[N_CELL_X];
    }
    
    int current_frame = 0;
    sf::Clock clock;
    std::vector<double> densityProfile(N_CELL_X, 0.0);
    std::vector<double> magnetizationProfile(N_CELL_X, 0.0);
    double sigma_flux = 0.0;
    int burnin = 200;
    while (window.isOpen() && current_frame < totalSteps / saveEvery) {
        window.clear(sf::Color::Black);
        window_plot.clear(sf::Color::Black);
        // Poll event
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) window.close();
        }
        if (window_plot.isOpen()){
            while (window_plot.pollEvent(event)){
                if (event.type == sf::Event::Closed)
                window_plot.close();
            }
        }
        
        // Maybe a bottleneck here? Assign to two threads?
        file.read(reinterpret_cast<char*>(mag_buffer), N_CELL_X * N_CELL_Y * sizeof(int));
        file.read(reinterpret_cast<char*>(den_buffer), N_CELL_X * N_CELL_Y * sizeof(int));
        if (flux) {
            file.read(reinterpret_cast<char*>(flux_buffer), N_CELL_X * sizeof(int));
        }


        
        double rolling_avg_mag = 0;
        double rolling_avg_den = 0;
        // Start drawing the spins
        for (int i = 0; i < N_CELL_X * N_CELL_Y; ++i) {
            // column-major order !
            float x = (i / N_CELL_Y) * ax;
            float y = (i % N_CELL_Y) * ay;
            
            // Magnetization
            sf::RectangleShape rectMag(sf::Vector2f(ax, ay));
            rectMag.setPosition(x, y);
            sf::Uint8 mag_int = static_cast<sf::Uint8>(255.0 * (1.0 - mag_buffer[i]) / 2.0);
            rectMag.setFillColor(sf::Color(mag_int, mag_int, mag_int));
            window.draw(rectMag);
            // Density
            sf::RectangleShape rectDen(sf::Vector2f(ax, ay));
            rectDen.setPosition(x, y + N_CELL_Y * ay + delta);
            float scale = std::min(1.0f, (float)den_buffer[i] / float(y_max));
            sf::Uint8 sub = static_cast<sf::Uint8>(255.0f * (1.0f - scale));
            rectDen.setFillColor(sf::Color(255, sub, sub));
            window.draw(rectDen);

            rolling_avg_mag += float(mag_buffer[i]) / N_CELL_Y;
            rolling_avg_den += float(den_buffer[i]) / N_CELL_Y;

            if (i % N_CELL_Y == N_CELL_Y - 1) {//at the end of a column
                if (current_frame > burnin){
                    densityProfile[i / N_CELL_Y] += rolling_avg_den / float(totalSteps / saveEvery - burnin);
                    magnetizationProfile[i / N_CELL_Y] += rolling_avg_mag / float(totalSteps / saveEvery - burnin);

                    // COmpute stddev of the flux
                    double flux_avg = std::accumulate(flux_buffer, flux_buffer + N_CELL_Y, 0.0) / N_CELL_Y;
                    double sq_sum = std::inner_product(flux_buffer, flux_buffer + N_CELL_Y, flux_buffer, 0.0);
                    double stdev = std::sqrt(sq_sum / N_CELL_Y - flux_avg * flux_avg);

                    sigma_flux += stdev / double(totalSteps / saveEvery - burnin);

                    //std::cout

                }
     
                
                float heightDensity = ((-0.5*sizeWindowPlot)/(y_max))*(rolling_avg_den) + sizeWindowPlot / 2;
                float heightMagnetization = ((deltay-0.5*sizeWindowPlot)/(y_max))*(rolling_avg_mag) + sizeWindowPlot / 2;
                float x_coord = spacing * (i / N_CELL_Y );
                
                sf::CircleShape point(2);
                point.setPosition(x_coord, heightDensity); 
                point.setFillColor(sf::Color::Green);
                sf::CircleShape pointMag(2);
                pointMag.setPosition(x_coord, heightMagnetization);
                pointMag.setFillColor(sf::Color::Cyan);
                
                window_plot.draw(pointMag);
                window_plot.draw(point);
                sf::RectangleShape line(sf::Vector2f(SizeHorizontal, 2));
                line.setPosition(0, sizeWindowPlot/2);
                line.setFillColor(sf::Color::White);
                window_plot.draw(line);

                rolling_avg_mag = 0.f;
                rolling_avg_den = 0.f;

                // Now displays the fluxes
                sf::RectangleShape fluxRect(sf::Vector2f(ax,N_CELL_Y*ay));
                fluxRect.setPosition((i / N_CELL_Y)*ax, 2*N_CELL_Y*ay + 2*delta);
                double flux = double(flux_buffer[i / N_CELL_Y]) / 100;
                sf::Uint8 mag_int = static_cast<sf::Uint8>(255.0 * (1.0 - flux) / 2.0);
                fluxRect.setFillColor(sf::Color(mag_int, mag_int, mag_int));
                window.draw(fluxRect);
            }
        }

        window_plot.display();
        window.display();

        current_frame++;
        sf::Time elapsed = clock.getElapsedTime();
        //std::cout << "Running at fps: " << 1.0 / elapsed.asSeconds() << "\r";
        clock.restart();
    }
    file.close();
    double final_avg_mag = std::accumulate(magnetizationProfile.begin(), magnetizationProfile.end(), 0.0) / N_CELL_X;
    double final_avg_den = std::accumulate(densityProfile.begin(), densityProfile.end(), 0.0) / N_CELL_X;
    //std::cout << "\n Final average magnetization per column: " << final_avg_mag << "\n";
    //std::cout << " Final average density per column: " << final_avg_den << "\n";
    std::cout << beta << " " << rho << " " << std::abs(final_avg_mag) / final_avg_den << "\n";
    std::cout << " Sigma flux: " << sigma_flux << "\n";
    //std::cout << "End frames\n";
}