#include<cmath>
#include<random>
#include<fstream>
#include<numeric>
#include<assert.h>


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
NDMolDyn<Dim>::NDMolDyn(int const numParticles, double const rho, double const rcut, double const skin, int const M, double const max_disp, double temp, int MCSteps, double cosThetaMax): m_rcut(rcut), m_rho(rho), m_skin(skin), m_M(M), max_displacement(max_disp), m_temperature(temp), m_numMCSteps(MCSteps), m_cosThetaMax(cosThetaMax){

    m_energyWell = -1;
    m_gamma = 0.3;

    m_L = std::pow(static_cast<double>(numParticles) / m_rho, 1.0 / static_cast<double>(Dim));
    m_V = powerDim<double, Dim>(m_L);

    // Next, we have to fill the m_particles vector with numParticles particles. Make sure no overlaps occur (check N^2 distances. Not optimal, but ok for initialization)
    m_particles.reserve(numParticles);
    for(int i = 0; i < numParticles; ++i) {
        vecd temp_pos;
        bool overlap = true;
        while(overlap) {
            for(int k = 0; k < Dim; ++k) {
                temp_pos[k] = m_unif_dist(m_gen) * m_L;
            }
            // Check for overlaps
            overlap = false;
            for(const auto& p : m_particles) {
                double d2 = dist_sq_pbc(temp_pos, p.position);
                if(d2 < 1.0) { // Overlap detected
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
    #ifdef DEBUG_PLUS
    std::cout << "Initialized system with " << numParticles << " particles in " << Dim << "D." << std::endl;
    #endif
};

// Constructor from Parameters struct
template<int Dim>
NDMolDyn<Dim>::NDMolDyn(Parameters const& params, int seed): m_rcut(params.param_r_cut), m_temperature(params.param_temperature), m_rho(params.param_rho), m_skin(params.param_skin), m_M(params.param_M), max_displacement(params.param_max_displacement), m_numMCSteps(params.param_num_MC_steps), m_cosThetaMax(params.param_max_cos_theta), m_energyWell(params.param_energyWell), m_gamma(params.param_gamma), m_unif_dist(0.0, 1.0), m_gen(seed) {
    // Set linear characteristic length and (hyper)volume
    m_L = std::pow(static_cast<double>(params.param_N) / params.param_rho, 1.0 / static_cast<double>(Dim));
    m_V = powerDim<double, Dim>(m_L);
    // Debug stuff
    #ifdef DEBUG_PLUS
    std::cout << "################ STARTING THE RUN ################" << std::endl;
    std::cout << "Initializing system with the following parameters:" << std::endl;
    std::cout << "Dimension: " << Dim << "D" << std::endl;
    std::cout << "Box length L: " << m_L << std::endl;
    std::cout << "Density rho: " << m_rho << std::endl;
    std::cout << "Temperature T: " << m_temperature << std::endl;
    std::cout << "Number of particles N: " << params.param_N << std::endl;
    std::cout << "Cutoff radius rcut: " << m_rcut << std::endl;
    std::cout << "Skin: " << m_skin << std::endl;
    std::cout << "Cells per dimension M: " << m_M << std::endl;
    std::cout << "Max displacement: " << max_displacement << std::endl;
    std::cout << "Number of MC steps: " << m_numMCSteps << std::endl;
    std::cout << "#################################################" << std::endl;
    #endif

    // Next, we have to fill the m_particles vector with N particles. Make sure no overlaps occur (check N^2 distances. Not optimal, but ok for initialization)
    m_particles.reserve(params.param_N);
    for(int i = 0; i < params.param_N; ++i) {
        vecd temp_pos;
        bool overlap = true;
        while(overlap) {
            for(int k = 0; k < Dim; ++k) {
                temp_pos[k] = m_unif_dist(m_gen) * m_L;
            }
            // Check for overlaps
            overlap = false;
            for(const auto& p : m_particles) {
                double d2 = dist_sq_pbc(temp_pos, p.position);
                if(d2 < 1.0) { // Overlap detected
                    overlap = true;
                    break;
                }
            }
        }
        m_particles.emplace_back(temp_pos); 
    }
    // Particles are now initialed with random positions and random patches in the box. Now, the initialization of headOfChain
    int total_cells = powerDim<int, Dim>(m_M);
    m_headOfChain.resize(total_cells);
    std::fill(m_headOfChain.begin(), m_headOfChain.end(), -1); // -1 means empty cell
    // Now init of the verlet lists
    m_verletLists.resize(params.param_N);
    m_backupParticles.resize(params.param_N);
    buildVerletLists();
    #ifdef DEBUG_PLUS
    std::cout << "Initialized system with " << params.param_N << " particles in " << Dim << "D." << std::endl;
    #endif
};

///
// CELL_INDEX: a multi_index (n_1, n_2, ..., n_Dim) where each n_i is at max M
//
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
// HERE WE HAVE DOUBLE COUNTING, BUT FOR MC NO OTHER WAY OUT

template<int Dim>
double NDMolDyn<Dim>::computeAngle(Patch const& patch, vecd  const& pos_i, vecd const& pos_j, double d2) const {
    // First of all, compute the versor
    double x = pos_j[0] - pos_i[0];
    x = x - m_L * std::round(x / m_L);
    double y = pos_j[1] - pos_i[1];
    y = y - m_L * std::round(y / m_L);
    double z = pos_j[2] - pos_i[2];
    z = z - m_L * std::round(z / m_L);

    double xx = x * patch.getX();
    double yy = y * patch.getY();
    double zz = z * patch.getZ();

    return (xx+yy+zz)/std::sqrt(d2);
}


template<int Dim>
double NDMolDyn<Dim>::getEnergy(vecd const& particle_pos, Patch const& particle_patch, int particle_idx) const {
    double energy = 0.0;
    for(auto it = m_verletLists[particle_idx].begin(); it != m_verletLists[particle_idx].end(); ++it){
        vecd pos_j = m_particles[*it].position;
        double d2 = dist_sq_pbc(particle_pos, pos_j);
        // If within interaction range (we've already checked for overlaps)
        if (d2 < (1.0 + 0.5)*(1.0 + 0.5)) {
            double cos_angle_1 = computeAngle(particle_patch, particle_pos, pos_j, d2);
            double cos_angle_2 = computeAngle(m_particles[*it].getPatch(), pos_j, particle_pos, d2);
            if (std::abs(cos_angle_1) > m_cosThetaMax && std::abs(cos_angle_2) >  m_cosThetaMax) { 
                energy += m_energyWell;
            }
        }   
    }
    return energy;
}

template<int Dim>
double NDMolDyn<Dim>::getEnergyNoPatches(vecd const& particle_pos, Patch const& particle_patch, int particle_idx) const {
    double energy = 0.0;
    for(auto it = m_verletLists[particle_idx].begin(); it != m_verletLists[particle_idx].end(); ++it){
        vecd pos_j = m_particles[*it].position;
        double d2 = dist_sq_pbc(particle_pos, pos_j);
        // If within interaction range (we've already checked for overlaps)
        if (d2 < (1.0 + 0.5)*(1.0 + 0.5)) {
                energy += m_energyWell;
        }   
    }
    return energy;
}

Patch perturbPatch(Patch const& old_patch, double moduler, std::mt19937& gen, std::uniform_real_distribution<double>& unif_dist) {
    // Sample a random versor
    double z = unif_dist(gen) * 2.0 - 1.0;
    double theta = unif_dist(gen) * 2.0 * M_PI;
    double r_xy = std::sqrt(1.0 - z*z); 
    double x = r_xy * std::cos(theta);
    double y = r_xy * std::sin(theta);

    double new_x = old_patch.getX() + moduler * x;
    double new_y = old_patch.getY() + moduler * y;
    double new_z = old_patch.getZ() + moduler * z;
    
    std::array<double, 3> new_direction = {new_x, new_y, new_z};
    return Patch(new_direction);
}

template<int Dim>
void NDMolDyn<Dim>::singleMCStep(){ //Run a single step of MC (thus N proposals)
    // First of all, select a random particle
    int particle_idx = m_unif_dist(m_gen) * m_particles.size();
    // Is it going to be a position move or a patch move?
    bool move_position = m_unif_dist(m_gen) < 0.5;
            move_position = true; // For testing purposes; force the position move only if we are not interested in the patch move
    if(move_position){
        // Save old position
        vecd old_pos = m_particles[particle_idx].position;
        // Propose a new position
        vecd new_pos;
        for(int k = 0; k < Dim; ++k){
            double delta = m_unif_dist(m_gen) * 2.0 * max_displacement - max_displacement;
            new_pos[k] = old_pos[k] + delta;
            if(new_pos[k] >= m_L){ //PBC
                new_pos[k] -= m_L;
            } else if(new_pos[k] < 0.0){
                new_pos[k] += m_L;
            }
        }
        // Now that we have the new_pos, verify whether there is an overlap
        // Preliminary step: check if the new position will overlap with any other particle
        for(auto it = m_verletLists[particle_idx].begin(); it != m_verletLists[particle_idx].end(); ++it){
            vecd pos_j = m_particles[*it].position;
            double d2 = dist_sq_pbc(new_pos, pos_j);
            if(d2 < 1.0){ //Overlap detected! Return and don't do anything, no need to check further (regardless of energy)
                return;
            }
        }
        //If no overlap, compute deltaE. First of all, compute the initial energy
        double initialEnergy = getEnergy(old_pos, m_particles[particle_idx].getPatch(), particle_idx);
        // Do not perturb the patch
        // Compute the new energy
        double newEnergy = getEnergy(new_pos, m_particles[particle_idx].getPatch(), particle_idx);
        double deltaEnergy = newEnergy - initialEnergy;
        // If we reach this point, no overlaps were detected, we can accept the move with a probability of exp(-deltaEnergy)
        double acceptanceProbability = std::exp(-deltaEnergy/m_temperature);
        double randomValue = m_unif_dist(m_gen);
        if(randomValue > acceptanceProbability) {
            return;  // reject
        }
        // If we reach this point, we can accept the move
        m_particles[particle_idx].position = new_pos;
        m_acceptCounter++;

        // Finally, check whether we need to rebuild the verlet lists
        double trigger = dist_sq_pbc(m_backupParticles[particle_idx], new_pos);
        if (trigger > m_skin/2 * m_skin/2) {
            // Rebuild the verlet lists
            buildVerletLists();
            #ifdef DEBUG_PLUS
            m_numVerletCalls++;
            #endif
        }
    } else {
        // Now the position is fixed, no need to verify whether an overlap will occur
        //If no overlap, compute deltaE. First of all, compute the initial energy
        double initialEnergy = getEnergy(m_particles[particle_idx].position, m_particles[particle_idx].getPatch(), particle_idx);
        // Now propose a new patch
        //Patch new_patch =  m_particles[particle_idx].getPatch();
        Patch new_patch = perturbPatch(m_particles[particle_idx].getPatch(), m_gamma, m_gen, m_unif_dist);
        // Compute the new energy
        double newEnergy = getEnergy(m_particles[particle_idx].position, new_patch, particle_idx);
        double deltaEnergy = newEnergy - initialEnergy;
        // If we reach this point, no overlaps were detected, we can accept the move with a probability of exp(-deltaEnergy)
        double acceptanceProbability = std::exp(-deltaEnergy/m_temperature);
        double randomValue = m_unif_dist(m_gen);
        if(randomValue > acceptanceProbability) {
            return;  // reject
        }
        // If we reach this point, we can accept the move
        m_particles[particle_idx].patch = new_patch;
        m_acceptCounter++;
    }    
}


template<int Dim>
void NDMolDyn<Dim>::MCRun(int savePositionsEvery, int saveEnergyEvery, std::string const& outfile_string, std::string const& outfile_energy_string){
    std::ofstream outfile(outfile_string);
    std::ofstream outfile_energy(outfile_energy_string);


    for(int step = 0; step < m_numMCSteps; ++step){
        m_acceptCounter = 0;
        // A single MC step consists of N trial moves:
        for(int _ = 0; _ < m_particles.size(); ++_){
            singleMCStep();
        }
        // If we want to save positions, ...
        if(savePositionsEvery){
            if (step % savePositionsEvery == 0){
                outfile << m_particles.size() << "\n";
                // WIP: adjust for number of dimensions arbitrarily
                if (Dim == 2){
                    outfile << "Lattice=\" " << m_L <<" 0.0 0.0 "<< m_L <<"  \" Properties=species:S:1:pos:R:2:dir:R:2 Time=" << step << "\n";
                } else if (Dim ==3){
                    outfile << "Lattice=\" " << m_L <<" 0.0 0.0 0.0 "<< m_L <<" 0.0 0.0 0.0 " << m_L <<" \" Properties=species:S:1:pos:R:3:dir:R:3 Time=" << step << "\n";
                }
                for(const auto& p : m_particles){
                    outfile << "H ";
                    for(int k = 0; k < Dim; ++k){
                        outfile << p.position[k] << " ";
                    }
                    outfile << p.patch.getX() << " " << p.patch.getY() << " " << p.patch.getZ() << "\n";
                }
            }
        }
        // If we want to save energy, ...
        if(saveEnergyEvery){
            if (step % saveEnergyEvery == 0){
                double totalEnergy = 0.0;
                for(int i = 0; i < m_particles.size(); ++i){
                    totalEnergy += getEnergy(m_particles[i].position, m_particles[i].getPatch(), i);
                }
                totalEnergy *= 0.5; // Correct for double counting
                outfile_energy << step << " " << totalEnergy / static_cast<double>(m_particles.size()) << "\n";
                #ifdef DEBUG_PLUS
                std::cout << "\rAcceptance ratio: " << static_cast<double>(m_acceptCounter) / (m_particles.size());
                if (step % (m_numMCSteps/100) == 0){
                    std::cout << "  Step: " << step << "/" << m_numMCSteps;
                }
                #endif
            }
        }

    }

    #ifdef DEBUG_PLUS
    std::cout << "\nNumber of Verlet list rebuilds: " << m_numVerletCalls << " over " << m_numMCSteps* m_particles.size() << " moves." << std::endl;
    #endif
    outfile.close();
    outfile_energy.close();
}

template<int Dim>
void NDMolDyn<Dim>::MCRun(int savePositionsEvery, int saveEnergyEvery, std::ofstream& outfile, std::ofstream& outfile_energy, int offset){
    std::ofstream temp("../data/Path" + std::to_string(m_idx) + ".dat", std::ios::app);
    for(int s = 0; s < m_numMCSteps; ++s){
        // A single MC step consists of N trial moves:
        for(int _ = 0; _ < m_particles.size(); ++_){
            singleMCStep();
        }

        if(savePositionsEvery > 0){
            if ((offset + s) % savePositionsEvery == 0){
                outfile << m_particles.size() << "\n";
                
                if (Dim == 2){
                    outfile << "Lattice=\" " << m_L <<" 0.0 0.0 "<< m_L 
                            <<" \" Properties=species:S:1:pos:R:2:dir:R:2 Time=" << offset +s << "\n";
                } else if (Dim == 3){
                    outfile << "Lattice=\" " << m_L <<" 0.0 0.0 0.0 "<< m_L <<" 0.0 0.0 0.0 " << m_L 
                            <<" \" Properties=species:S:1:pos:R:3:dir:R:3 Time=" << offset + s << "\n";
                }

                // Scrittura particelle
                for(const auto& p : m_particles){
                    outfile << "H ";
                    for(int k = 0; k < Dim; ++k){
                        outfile << p.position[k] << " ";
                    }
                    outfile << p.patch.getX() << " " << p.patch.getY() << " " << p.patch.getZ() << "\n";
                }
            }
        }

        if(saveEnergyEvery > 0){
            if ((offset + s) % saveEnergyEvery == 0){
                double totalEnergy = 0.0;
                for(int i = 0; i < m_particles.size(); ++i){
                    totalEnergy += getEnergy(m_particles[i].position, m_particles[i].getPatch(), i);
                }
                totalEnergy *= 0.5; // Correct for double counting

                outfile_energy << offset + s << " " << totalEnergy / static_cast<double>(m_particles.size()) << "\n";

                //temp << offset+s << " " << m_temperature << "\n";
            }
        }
    }
}

template<int Dim>
double NDMolDyn<Dim>::getTotalEnergy() const {
    double totalEnergy = 0.0;
    for(int i = 0; i < m_particles.size(); ++i){
        totalEnergy += getEnergy(m_particles[i].position, m_particles[i].getPatch(), i);
    }
    totalEnergy *= 0.5; // Correct for double countingù
    return totalEnergy;
}