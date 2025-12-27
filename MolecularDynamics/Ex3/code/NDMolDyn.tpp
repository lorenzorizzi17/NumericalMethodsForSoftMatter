#include<cmath>
#include<random>
#include<fstream>
#include<numeric>
#include<assert.h>
#include <string> 
#include <sstream> 

#define DEBUG_PLUS

// Squared distance with MIC
template<int Dim>
double NDMolDyn<Dim>::dist_sq_pbc(const vecd& pos_i, const vecd& pos_j) const {
    double dist_sq = 0.0;
    for(int k = 0; k < Dim; ++k) {
        double dr = pos_i[k] - pos_j[k];
        dr = dr - m_L * std::round(dr / m_L);
        dist_sq += dr * dr;
    }
    return dist_sq;
}

// Helper function, more efficient than std::pow when exponent is an integer
template<typename T, int Dim>
T powerDim(T base){
    T tmp = 1.;
    for(int i = 0; i < Dim; ++i){
        tmp *= base;
    }
    return tmp;
}

// Constructor from explicit parameters
template<int Dim>
NDMolDyn<Dim>::NDMolDyn(int const numParticles, double const rho, double const rcut, double const skin, int const M, double temp): m_rcut(rcut), m_rho(rho), m_skin(skin), m_M(M), m_initialTemperature(temp) {

    m_L = std::pow(static_cast<double>(numParticles) / m_rho, 1.0 / static_cast<double>(Dim));
    m_V = powerDim<double, Dim>(m_L);

    // Next, we have to fill the m_particles vector with numParticles particles. Make sure no overlaps occur (check N^2 distances. Not optimal, but ok for initialization)
    m_particles.reserve(numParticles);
    for(int i = 0; i < numParticles; ++i) {
        vecd temp_pos;
        bool overlap = true;
        while(overlap) {
            for(int k = 0; k < Dim; ++k) {
                temp_pos[k] = (static_cast<double>(rand()) / RAND_MAX) * m_L;  // RANDOM ENGINE!!
            }
            // Check for overlaps
            overlap = false;
            for(const auto& p : m_particles) {
                double d2 = dist_sq_pbc(temp_pos, p.position);
                if(d2 < 1) { // Overlap detected
                    overlap = true;
                    break;
                }
            }
        }
        m_particles.emplace_back(temp_pos); 
    }
    // Particles are now initialed with random positions in the box. Now, the initialization of headOfChain
    int total_cells = powerDim<int, Dim>(m_M);
    m_headOfChain.resize(total_cells);
    std::fill(m_headOfChain.begin(), m_headOfChain.end(), -1); // -1 means empty cell
    // Now init of the verlet lists
    m_verletLists.resize(numParticles);
    m_backupParticles.resize(numParticles);
    buildVerletLists();

    // Initialize forces vector ...
    m_forces.resize(numParticles, {0.0});
    // ... and precompute the forces
    computeForces();

    // Velocity initialization: Maxwell-Boltzmann distribution
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, std::sqrt(m_initialTemperature));
    for(auto& p : m_particles) {
        for(int k = 0; k < Dim; ++k) {
            p.velocity[k] = distribution(generator);
        }
    }
    vecd v_cm = {0.0};
    for(const auto& p : m_particles) {
        for(int k = 0; k < Dim; ++k) {
            v_cm[k] += p.velocity[k];
        }
    }

    for(int k = 0; k < Dim; ++k) {
        v_cm[k] /= static_cast<double>(m_particles.size());
    }

    for(auto& p : m_particles) {
        for(int k = 0; k < Dim; ++k) {
            p.velocity[k] -= v_cm[k];
        }
    }

    #ifdef DEBUG_PLUS
    std::cout << "Initialized system with " << numParticles << " particles in " << Dim << "D." << std::endl;
    #endif
};


template<int Dim>
NDMolDyn<Dim>::NDMolDyn(std::string const& filename, double const L, double const rcut, double const skin, int const M, double const temp) : m_rcut(rcut), m_skin(skin), m_M(M), m_initialTemperature(temp) {
    // Next, we have to fill the m_particles vector with numParticles particles. Make sure no overlaps occur (check N^2 distances. Not optimal, but ok for initialization)
    std::ifstream infile(filename);
    int row = 0;
    double max_coord;
    while(infile){
        std::string line;
        std::getline(infile, line);
        if(line.empty()) break;
        if(row > 1){
            std::istringstream iss(line);
            std::string species;
            vecd position;
            iss >> species >> position[0] >> position[1] >> position[2];
            max_coord = std::max(position[0], max_coord);
            Particle<Dim> p(position);
            m_particles.push_back(p);
        }
        row++;
    }
    infile.close();
    int numParticles = m_particles.size();

    m_L = L;
    m_V = powerDim<double, Dim>(m_L);
    m_rho = double(numParticles) / m_V;
    // Particles are now initialed with random positions in the box. Now, the initialization of headOfChain
    int total_cells = powerDim<int, Dim>(m_M);
    m_headOfChain.resize(total_cells);
    std::fill(m_headOfChain.begin(), m_headOfChain.end(), -1); // -1 means empty cell
    // Now init of the verlet lists
    m_verletLists.resize(numParticles);
    m_backupParticles.resize(numParticles);
    buildVerletLists();

    // Initialize forces vector ...
    m_forces.resize(numParticles, {0.0});
    // ... and precompute the forces
    computeForces();

    // Velocity initialization: Maxwell-Boltzmann distribution
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, std::sqrt(m_initialTemperature));
    for(auto& p : m_particles) {
        for(int k = 0; k < Dim; ++k) {
            p.velocity[k] = distribution(generator);
        }
    }

    #ifdef DEBUG_PLUS
    std::cout << "Initialized system with " << numParticles << " particles in " << Dim << "D." << std::endl;
    #endif
}

template<int Dim>
void NDMolDyn<Dim>::buildVerletLists(){
    for(auto& list : m_verletLists) {
        list.clear(); 
    }
    std::fill(m_headOfChain.begin(), m_headOfChain.end(), -1);

    double cellSize = m_L / static_cast<double>(m_M);
    double threshold_sq = (m_rcut + m_skin) * (m_rcut + m_skin);

    for(auto it = m_particles.begin(); it != m_particles.end(); ++it){
        vecd temp_pos = it->position;
        multidx cell_idx;
        for(int k = 0; k < Dim; ++k){
            int idx_k = static_cast<int>(temp_pos[k] / cellSize);
            idx_k = (idx_k % m_M + m_M) % m_M;
            cell_idx[k] = idx_k; 
        }
        int linear_idx = 0;
        for(int k = 0; k < Dim; ++k){
            linear_idx += cell_idx[k] * std::pow(m_M, k); 
        }
        m_particles[std::distance(m_particles.begin(), it)].next = m_headOfChain[linear_idx];
        m_headOfChain[linear_idx] = std::distance(m_particles.begin(), it);
    }

    int total_cells = powerDim<int, Dim>(m_M);
    
    for(int cell_idx = 0; cell_idx < total_cells; ++cell_idx){
        if(m_headOfChain[cell_idx] == -1) continue;

        multidx current_cell_coords;
        int temp_idx = cell_idx;
        for(int k = 0; k < Dim; ++k) { 
            current_cell_coords[k] = temp_idx % m_M;
            temp_idx /= m_M;
        }

        std::vector<int> unique_neighbors;
        unique_neighbors.reserve(powerDim<int, Dim>(3));

        for(int n = 0; n < powerDim<int, Dim>(3); ++n) {
            int neighbor_cell_linear_idx = 0;
            int stride = 1; 
            int temp_n = n; 
            for(int k = 0; k < Dim; ++k) {
                int offset = (temp_n % 3) - 1; 
                temp_n /= 3;
                int c = current_cell_coords[k] + offset;
                c = (c % m_M + m_M) % m_M;
                neighbor_cell_linear_idx += c * stride;
                stride *= m_M;
            }

            bool already_added = false;
            for(int existing : unique_neighbors){
                if(existing == neighbor_cell_linear_idx){
                    already_added = true;
                    break;
                }
            }
            if(!already_added){
                unique_neighbors.push_back(neighbor_cell_linear_idx);
            }
        }

        int particle_idx = m_headOfChain[cell_idx];
        while(particle_idx != -1){
            auto& verletList = m_verletLists[particle_idx];
            vecd pos_i = m_particles[particle_idx].position;
            
            for(int neighbor_idx : unique_neighbors) {
                int j_idx = m_headOfChain[neighbor_idx];                
                while(j_idx != -1) {
                    if(particle_idx != j_idx) { 
                        vecd pos_j = m_particles[j_idx].position;
                        double d2 = dist_sq_pbc(pos_i, pos_j);
                        if(d2 < threshold_sq) {
                            verletList.push_back(j_idx);
                        }
                    }
                    j_idx = m_particles[j_idx].next;
                }
            }
            particle_idx = m_particles[particle_idx].next;
        }
    }

    for(int i = 0; i < m_particles.size(); ++i){
        m_backupParticles[i] = m_particles[i].position;
    }
};

template<int Dim>
void NDMolDyn<Dim>::computeForces(){
    // Reset forces
    std::fill(m_forces.begin(), m_forces.end(), vecd{0.0, 0.0, 0.0});

    double rcut_sq = m_rcut * m_rcut;

    for (int i = 0; i < m_particles.size(); ++i){
        vecd pos_i = m_particles[i].position;
        
        for(const int j_idx : m_verletLists[i]){
            vecd pos_j = m_particles[j_idx].position;
            double d2 = dist_sq_pbc(pos_i, pos_j);

            // Compute Lennard-Jones force if within cutoff
            if (d2 < rcut_sq) {
                double inv_r2 = 1.0 / d2;
                double inv_r6 = inv_r2 * inv_r2 * inv_r2;
                // La formula 48 * ... è corretta per la forza vettoriale
                double f_scalar = 48.0 * inv_r6 * (inv_r6 - 0.5) * inv_r2; 
                
                for(int k = 0; k < Dim; ++k){
                    double dr = pos_i[k] - pos_j[k];
                    dr = dr - m_L * std::round(dr / m_L);
                    m_forces[i][k] += f_scalar * dr;
                }
            }
        }
    }
}

template<int Dim>
double NDMolDyn<Dim>::computePotEnergy() const {
    double total_potential_energy = 0.0;
    double rcut_sq = m_rcut * m_rcut;
    
    // Rallenta molto, dovrei spostarlo fuori dal ciclo
    double e_shift = 4.0 * (std::pow(1.0 / m_rcut, 12) - std::pow(1.0 / m_rcut, 6));

    for (int i = 0; i < m_particles.size(); ++i){
        vecd pos_i = m_particles[i].position;
    
        for(const int j_idx : m_verletLists[i]) {
            vecd pos_j = m_particles[j_idx].position;
            
            double d2 = dist_sq_pbc(pos_i, pos_j);

            if (d2 < rcut_sq) {
                double inv_r2 = 1.0 / d2;
                double inv_r6 = inv_r2 * inv_r2 * inv_r2;
                double potential_energy_ij = 4.0 * inv_r6 * (inv_r6 - 1.0);
                
                total_potential_energy += (potential_energy_ij - e_shift);
            }
        }
    }
    // Each pair counted twice
    return total_potential_energy * 0.5;
}


template<int Dim>
double NDMolDyn<Dim>::computeKinEnergy() const {
    double total_kinetic_energy = 0.0;

    for (int i = 0; i < m_particles.size(); ++i){
        vecd vel = m_particles[i].velocity;
        double v2 = 0.0;
        for(int k = 0; k < Dim; ++k){
            v2 += vel[k] * vel[k];   // m = 1
        }
        total_kinetic_energy += 0.5 * v2;
    }
    return total_kinetic_energy;   
}
template<int Dim>
void NDMolDyn<Dim>::MDRun(int NumSteps, double stepSize, int saveEvery, std::string const& filename){
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    // Assumo che m_BerendsenTemperature sia la temperatura target anche per Andersen
    std::normal_distribution<double> normal_dist(0.0, std::sqrt(m_BerendsenTemperature));

    std::ofstream outfile("../data/traj.xyz");
    std::ofstream energyfile(filename);

    for(int step = 0; step < NumSteps; ++step){
        
        double max_displacement_sq_in_step = 0.0;

        // --- STEP 1: Velocità (metà) e Posizione (tutta) ---
        for(int i = 0; i < m_particles.size(); ++i){
            for(int k = 0; k < Dim; ++k){
                m_particles[i].velocity[k] += 0.5 * m_forces[i][k] * stepSize;
                m_particles[i].position[k] += m_particles[i].velocity[k] * stepSize;
                
                if (m_particles[i].position[k] >= m_L) m_particles[i].position[k] -= m_L;
                else if (m_particles[i].position[k] < 0.0) m_particles[i].position[k] += m_L;
            }
            double d2_disp = dist_sq_pbc(m_particles[i].position, m_backupParticles[i]);
            if(d2_disp > max_displacement_sq_in_step) max_displacement_sq_in_step = d2_disp;
        }

        // --- CHECK VERLET LISTS ---
        if(max_displacement_sq_in_step > (m_skin * 0.5) * (m_skin * 0.5)){
            buildVerletLists(); 
            m_numVerletCalls++;
        } 
        computeForces();

        // --- STEP 2: Velocità (seconda metà) ---
        for(int i = 0; i < m_particles.size(); ++i){
            for(int k = 0; k < Dim; ++k){
                m_particles[i].velocity[k] += 0.5 * m_forces[i][k] * stepSize;
            }
        }

        // --- TERMOSTATI ---
        if (m_Thermostat == ThermostatType::BERENDSEN){
            double current_kinetic_energy = computeKinEnergy();
            double current_temperature = (2.0 / (Dim * (m_particles.size()-1))) * current_kinetic_energy;
            double lambda = std::sqrt(1.0 + (1.0 / m_tau_berendsen) * (m_BerendsenTemperature / current_temperature - 1.0));
            for(int i = 0; i < m_particles.size(); ++i){
                for(int k = 0; k < Dim; ++k) m_particles[i].velocity[k] *= lambda;
            }
        } else if (m_Thermostat == ThermostatType::ANDERSEN){
            // Andersen thermostat
            double nu = m_omega_Andersen; 
            double prob_collision = nu * stepSize; 
            for(int i = 0; i < m_particles.size(); ++i){
                double rand_num = uniform_dist(generator);
                if(rand_num < prob_collision){
                    for(int k = 0; k < Dim; ++k){
                        m_particles[i].velocity[k] = normal_dist(generator);
                    }
                }
            }
        }

        // --- SAVING ---
        if (step % saveEvery == 0){
            outfile << m_particles.size() << "\nLattice=\"" << m_L << " 0.0 0.0 0.0 " << m_L << " 0.0 0.0 0.0 " << m_L << "\" Properties=species:S:1:pos:R:3:vel:R:3\n";
            for(const auto& p : m_particles){
                outfile << "Ar " << p.position[0] << " " << p.position[1] << " " << p.position[2] 
                        << " " << p.velocity[0] << " " << p.velocity[1] << " " << p.velocity[2] << "\n";
            }
            double KE = computeKinEnergy();
            double PE = computePotEnergy();
            // Kinetic enery of center of mass only.
            // Compute the velocity of the center of mass
            vecd vMass = {0,0,0};
            for (const auto& p : m_particles){
                vMass[0] += p.velocity[0];
                vMass[1] += p.velocity[1];
                vMass[2] += p.velocity[2];
            }
            double KE_COM = vMass[0]*vMass[0] + vMass[1]*vMass[1]+ vMass[2]*vMass[2];
            energyfile << step << " " << KE / double(m_particles.size()) << " " << PE / double(m_particles.size()) << " " << (KE + PE) / double(m_particles.size()) << " " <<  0.5*KE_COM / (double(m_particles.size()) * double(m_particles.size())) << "\n";
        }
        #ifdef DEBUG_PLUS
        if (saveEvery > 0 && step % 1000 == 0) {
            std::cout << "\rStep: " << step << "/" << NumSteps << std::flush;
        }
        #endif
    }
    
    #ifdef DEBUG_PLUS
    std::cout << "\nVerlet rebuilds: " << m_numVerletCalls << std::endl;
    #endif
    outfile.close();
    energyfile.close();
}