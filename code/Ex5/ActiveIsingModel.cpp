#include "ActiveIsingModel.hpp"
#include <iostream>
#include <fstream>

//Helper function to write to file (no flux)
void writeSingleFrame(std::ofstream& file, const std::vector<Cells>& cells, int Nx, int Ny, double beta, double epsilon, double rho){
    int k = 0;

    std::vector<int> mag_buffer(Nx*Ny);
    std::vector<int> den_buffer(Nx*Ny);
    

    for(const auto& cell : cells) {
        mag_buffer[k] = cell.getMagnetization();
        den_buffer[k] = cell.getDensity();
        k++;
    }
    file.write(reinterpret_cast<const char*>(mag_buffer.data()), mag_buffer.size() * sizeof(int));
    file.write(reinterpret_cast<const char*>(den_buffer.data()), den_buffer.size() * sizeof(int));
}

//Helper function to write to file (flux)
void writeSingleFrame(std::ofstream& file, const std::vector<Cells>& cells, std::vector<int>& fluxes, int Nx, int Ny, double beta, double epsilon, double rho){
    int k = 0;

    std::vector<int> mag_buffer(Nx*Ny);
    std::vector<int> den_buffer(Nx*Ny);

    for(const auto& cell : cells) {
        mag_buffer[k] = cell.getMagnetization();
        den_buffer[k] = cell.getDensity();
        k++;
    }
    file.write(reinterpret_cast<const char*>(mag_buffer.data()), mag_buffer.size() * sizeof(int));
    file.write(reinterpret_cast<const char*>(den_buffer.data()), den_buffer.size() * sizeof(int));
    file.write(reinterpret_cast<const char*>(fluxes.data()), fluxes.size() * sizeof(int));
}

const std::vector<Cells>& ActiveIsingModel::getSpins() const {
    return m_cells;
}

ActiveIsingModel::ActiveIsingModel(int N_CELL_X, int N_CELL_Y, double density, double epsilon, double T, initType IT, bool flux) : N_CELL_X(N_CELL_X), N_CELL_Y(N_CELL_Y), m_epsilon(epsilon), m_T(T), m_rho(density) {
    
    N_SPINS = int(density * N_CELL_X * N_CELL_Y);

    std::random_device rd; 
    m_gen = std::mt19937(rd());

    m_cells.reserve(N_CELL_X * N_CELL_Y);
    if(flux){
        m_saveFlux = true;
        m_fluxes.resize(N_CELL_X, 0); // Initialize fluxes to zero for each column
    }
    std::uniform_int_distribution<int> dist(0, 1);
    switch (IT) {
        case UNIFORM:
            for(int j = 0; j < N_CELL_Y*N_CELL_X; j++){
                m_cells.emplace_back(Cells(density));
                for(int i = 0; i < int(N_SPINS / (N_CELL_X * N_CELL_Y)); i++){
                    int orientation = dist(m_gen) * 2 - 1; 
                    m_cells[j].addSpin(orientation);
                }
            }
            break;
        case PEAKED:
            // All spins in 4 central cells
            for(int j = 0; j < N_CELL_Y*N_CELL_X; j++){
                m_cells.emplace_back(Cells(N_SPINS));   //Introduce empty cells (they will reserve space for 2*N_SPINS each)
                int central_cell_x = N_CELL_X / 2;
                int central_cell_y = N_CELL_Y / 2;
                int central_cell_index = central_cell_y * N_CELL_X + central_cell_x;
                if (j == central_cell_index){
                    for(int _ = 0; _ < N_SPINS; _++){
                        int orientation = +1; 
                        m_cells[j].addSpin(orientation);
                    }
                }
            }
            break;
        default:
            std::cerr << "Invalid initialization type. Using UNIFORM instead." << std::endl;
            break;
    }

    m_deltaT = 2/(4*m_D + std::exp(1/m_T));
}

// In the graphical windows, we will display both the magnetizazion and the local density of each cell
void ActiveIsingModel::runGraphics(){
    //GRAPHICAL PARAMETERS
    int a = 5;
    sf::RenderWindow window(sf::VideoMode(N_CELL_X*a, 2*N_CELL_Y*a+10*a), "Active Ising Model");
    int sizeWindowPlot = 600;
    sf::RenderWindow window_plot(sf::VideoMode(sizeWindowPlot, sizeWindowPlot), "Plotter");

    
    while (window.isOpen()){
        sf::Event event;
        while (window.pollEvent(event)){
            if (event.type == sf::Event::Closed)
            window.close();
        }
        if (window_plot.isOpen()){
            while (window_plot.pollEvent(event)){
                if (event.type == sf::Event::Closed)
                window_plot.close();
            }
        }
        
        
        window.clear(sf::Color::Black);
        window_plot.clear(sf::Color::Black);
        //draw a straight line at avg_density
        sf::RectangleShape avgDensityLine(sf::Vector2f(sizeWindowPlot, 2));
        avgDensityLine.setFillColor(sf::Color::White);
        avgDensityLine.setPosition(0, sizeWindowPlot / 2);
        window_plot.draw(avgDensityLine);
        
        for (auto it = m_cells.begin(); it != m_cells.end(); it++){
            double avg_mag = it->getMagnetization() / it->getDensity();
            sf::RectangleShape rectangle(sf::Vector2f(a, a));
            int idx = std::distance(m_cells.begin(), it);
            int x_coord = (idx / N_CELL_Y)*a;
            int y_coord = (idx % N_CELL_Y)*a;
            rectangle.setPosition(x_coord, y_coord);
            
            sf::Uint8 intensity = static_cast<sf::Uint8>(255.0 * (1.0 - avg_mag) / 2.0);
            
            sf::Color cell_color;
            cell_color.r = intensity;
            cell_color.g = intensity;
            cell_color.b = intensity;
            
            rectangle.setFillColor(cell_color);
            window.draw(rectangle);
        }
        // DENSITY DISPLAY
        const float DENSITY_SCALE_MAX = 2.0f * m_rho;
        for (auto it = m_cells.begin(); it != m_cells.end(); it++){
            int local_den = it->getDensity(); // Accesso diretto tramite iteratore
            sf::RectangleShape rectangle(sf::Vector2f(a, a));
            int idx = std::distance(m_cells.begin(), it);
            int x_coord = (idx / N_CELL_Y)*a;
            int y_coord = (idx % N_CELL_Y)*a;
            rectangle.setPosition(x_coord, y_coord + N_CELL_Y*a + 10*a); // Offset per il display densità
            
            // 1. Calcola il fattore di scala (crampato tra 0 e 1)
            float scale_factor = (DENSITY_SCALE_MAX > 0) ? std::min(1.0f, (float)local_den / DENSITY_SCALE_MAX) : 0.0f;
            
            // 2. Mappatura Invertita (Bianco -> Rosso):
            // Il canale R rimane alto (255) o viene fissato.
            // I canali G e B si riducono man mano che la densità aumenta, partendo da 255.
            
            sf::Uint8 red_channel = 255;
            // La sottrazione deve essere scalata in modo che arrivi a 0 quando scale_factor=1.0
            sf::Uint8 subtractive_channel = static_cast<sf::Uint8>(255.0f * (1.0f - scale_factor));

            sf::Color cell_color;
            cell_color.r = red_channel;        // Rosso fisso a 255 (o vicino)
            cell_color.g = subtractive_channel; // Va da 255 (Bianco) a 0 (Rosso)
            cell_color.b = subtractive_channel; // Va da 255 (Bianco) a 0 (Rosso)
            
            rectangle.setFillColor(cell_color);
            window.draw(rectangle);
        }
        // A full MCMC step:
        for(int _ = 0; _ < N_SPINS; _++){
            this->MCMCStep();
        }
        m_time += m_deltaT;
        std::cout << "Time: " << m_time << "\r";
        window.display();


        // PLOTTING SECTION (to be implemented)
        double spacing = sizeWindowPlot / N_CELL_X;
        double y_max = 3*N_SPINS / (N_CELL_X * N_CELL_Y);
        int k = 0;
        int deltay = 10;
        for(auto it = m_cells.begin(); it!= m_cells.end(); it+=N_CELL_Y){
            /// it punta ad una cella
            double total_mag_col = 0.0;
            double total_rho_col = 0.0;
            for(auto col_it = it; col_it != it + N_CELL_Y; ++col_it) {
                double rho = col_it->getDensity();
                total_rho_col += rho;
                if (rho > 0) {
                    total_mag_col += col_it->getMagnetization(); 
                }
            }
            double avg_rho_y = total_rho_col / N_CELL_Y;
            double avg_mag_y = total_mag_col / N_CELL_Y;
            double x_coord = spacing * k;

            double heightDensity = ((deltay-0.5*sizeWindowPlot)/(y_max))*(avg_rho_y) + sizeWindowPlot / 2;
            double heightMagnetization = ((deltay-0.5*sizeWindowPlot)/(y_max))*(avg_mag_y) + sizeWindowPlot / 2;

            sf::CircleShape point(2); // Raggio 2 è più visibile
            point.setPosition(x_coord, heightDensity); 
            point.setFillColor(sf::Color::Green);
            sf::CircleShape pointMag(2);
            pointMag.setPosition(x_coord, heightMagnetization);
            pointMag.setFillColor(sf::Color::Cyan);

            window_plot.draw(pointMag);
            window_plot.draw(point);
            k++;
        }
        window_plot.display();

    }
}

void ActiveIsingModel::run(int totalSteps, int saveEvery){
    // Prepare the file
    std::string prestring; 
    if(m_saveFlux){
        m_fluxes.resize(N_CELL_X, 0); 
        prestring = "f";
    } else {
        prestring = "";
    }
    
    std::string filename = "../AIMData/Nx" + std::to_string(N_CELL_X) + 
                            "_Ny" + std::to_string(N_CELL_Y) + 
                            "_beta" + std::to_string(1.0/m_T).substr(0,5) + 
                            "_eps" + std::to_string(m_epsilon).substr(0,5) + 
                            "_rho" + std::to_string(m_rho).substr(0,5) +
                            "_Ntot" + std::to_string(totalSteps) +
                            "_Nevry" + std::to_string(saveEvery) +
                            prestring + ".dat";
    std::ofstream fileMag(filename, std::ios::binary);
    for(int step = 0; step < totalSteps; step++){
        if (m_saveFlux){  
            m_fluxes.assign(N_CELL_X, 0); 
        }
        for(int _ = 0; _ < N_SPINS; _++){  
            this->MCMCStep();
        }
        m_time += m_deltaT;
        if (step % saveEvery == 0){
            if(!m_saveFlux){
                writeSingleFrame(fileMag, this->getSpins(), N_CELL_X, N_CELL_Y, 1. / m_T, m_epsilon, m_rho);
            } else {
                writeSingleFrame(fileMag, this->getSpins(), m_fluxes, N_CELL_X, N_CELL_Y, 1. / m_T, m_epsilon, m_rho);
            }
        }

    }
    if(fileMag.is_open()) fileMag.close(); 
}

void ActiveIsingModel::MCMCStep(){
    // Perform a single step of the MCMC algorithm
    std::uniform_int_distribution<int> cell_dist(0, N_CELL_X * N_CELL_Y - 1);
    int cell_index = cell_dist(m_gen);
    int cell_density = m_cells[cell_index].getDensity();
    // If a cell is empty, try again
    while (cell_density <= 0) { // !!!!
        cell_index = cell_dist(m_gen);
        cell_density = m_cells[cell_index].getDensity();
    }
    int cell_mag = m_cells[cell_index].getMagnetization();
    // Select a random spin in that cell
    std::uniform_int_distribution<int> spin_dist(0, cell_density - 1);
    int spin_index = spin_dist(m_gen);
    // Decide whether to flip or hop according to the correct rates
    int current_spin = m_cells[cell_index][spin_index];


    double p_hop = 4 * m_deltaT * m_D;
    double p_spinflip = m_deltaT * std::exp(- float(current_spin) * (1./m_T) * (double(cell_mag)/double(cell_density)));
    
    std::uniform_real_distribution<float> move_dist(0.0, 1.0);
    float rnd = move_dist(m_gen);
    if (rnd < p_spinflip){
        MCMC_SpinFlip(cell_index, spin_index);
    } else if (rnd < p_spinflip + p_hop){ 
        MCMC_Hop(cell_index, spin_index);
    } else {}
}

void ActiveIsingModel::MCMC_SpinFlip(int cell_index, int spin_index){
    m_cells[cell_index][spin_index] = - m_cells[cell_index][spin_index]; // Accept the flip
}

// Here I have to 
void ActiveIsingModel::MCMC_Hop(int cell_index, int spin_index) {
    const int current_spin = m_cells[cell_index][spin_index];
    // Choose a direction to hop to (1->Right, 2->Down, 3->Left, 4->Up)
    std::uniform_real_distribution<float> dir_dist(0,1);
    float rnd = dir_dist(m_gen);

    int target_cell_index;
    int direction;
    if (rnd < 0.25) {
        //up
        target_cell_index = (cell_index % N_CELL_Y == 0) ?
        (cell_index + N_CELL_Y -1) :
        (cell_index - 1); 
        direction = 4;
    } else if (rnd < 0.5) {
        //down
        target_cell_index = (cell_index % N_CELL_Y == N_CELL_Y - 1) ?
        (cell_index - (N_CELL_Y - 1)) :
        (cell_index + 1); 
        direction = 2;
    } else if (rnd < 0.75 - float(current_spin*m_epsilon)/4) {
        //left
        target_cell_index = (cell_index / N_CELL_Y == 0) ?
        (cell_index + N_CELL_Y*(N_CELL_X - 1)) :
        (cell_index - N_CELL_Y); 
        direction = 3;
    } else {
        //right
        target_cell_index = (cell_index / N_CELL_Y == N_CELL_X - 1) ? // PBC
        (cell_index - N_CELL_Y * (N_CELL_X - 1)) :
        (cell_index + N_CELL_Y); 
        direction = 1;
    }
    // Remove the spin from the current cell
    int last_index = m_cells[cell_index].getDensity() - 1;
    // Swap
    m_cells[cell_index][spin_index] = m_cells[cell_index][last_index]; 
    m_cells[cell_index].getData().pop_back(); 
    m_cells[target_cell_index].addSpin(current_spin);

    if (m_saveFlux){
        if (direction == 1) { //particle has moved right
            // Update fluxes. First of all, find the original column of the particle:
            int col_old = cell_index / N_CELL_Y;
            m_fluxes[col_old] += 1;
        } else if (direction == 3) { //particle has moved left
            int col_old = cell_index / N_CELL_Y;
            int channel_to_decrement;
            if (col_old == 0) {
                channel_to_decrement = N_CELL_X - 1; 
            } else {
                channel_to_decrement = col_old - 1;
            }    
            m_fluxes[channel_to_decrement] += -1;
        }
    }
}



