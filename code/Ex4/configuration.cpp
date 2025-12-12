#include "configuration.hpp"
#include <iostream>
#include <deque>

// Standard constructor
SpinConfiguration::SpinConfiguration(int N, double temperature, double p) : m_N(N), m_T(temperature) {
    // Initialize the random engine
    std::random_device rd;
    m_rng = std::mt19937(rd());
    // Reserve a bit of space in the heap so that we don't have to reallocate memory
    m_spins.reserve(N);
    // Initialize the MCMC engine (default to MetropolisHastings, if the user doesn't set it explicitly)
    m_engine = MCMCengine(this, MCMCType::MetropolisHastings, m_N);
    // Initialize spins randomly
    std::uniform_real_distribution<double> dist(0, 1);
    for (int i = 0; i < N; ++i) {
        m_spins.push_back(dist(m_rng) < p);
    }
}

SpinConfiguration::SpinConfiguration(int N, double temperature, double p, double h) : m_N(N), m_T(temperature), m_h(h) {
    // Initialize the random engine
    std::random_device rd;
    m_rng = std::mt19937(rd());
    // Reserve a bit of space in the heap so that we don't have to reallocate memory
    m_spins.reserve(N);
    // Initialize the MCMC engine (default to MetropolisHastings, if the user doesn't set it explicitly)
    m_engine = MCMCengine(this, MCMCType::MetropolisHastings, m_N);
    // Initialize spins randomly
    std::uniform_real_distribution<double> dist(0, 1);
    for (int i = 0; i < N; ++i) {
        m_spins.push_back(dist(m_rng) < p);
    }
}

// Getters
double SpinConfiguration::getMagnetization() const {
    int sum = 0;
    for (const auto& spin : m_spins) {
        sum += spin ? 1 : -1;
    }
    return static_cast<double>(sum) / double(m_N);
}

std::array<int, 4> SpinConfiguration::getNeighbourhoodSpins(int row, int col) const {
    int L = static_cast<int>(std::sqrt(m_N));
    std::array<int, 4> neighbours;
    // PBC!
    if(m_PBC){
        int right_col = (col + 1) % L;
        int left_col  = (col - 1 + L) % L;
        int up_row    = (row - 1 + L) % L;
        int down_row  = (row + 1) % L;
        // get the values
        neighbours[0] = m_spins[row * L + right_col] ? 1 : -1;  // Right
        neighbours[1] = m_spins[row * L + left_col]  ? 1 : -1;  // Left
        neighbours[2] = m_spins[up_row * L + col]    ? 1 : -1;  // Up
        neighbours[3] = m_spins[down_row * L + col]  ? 1 : -1;  // Down
    } else {
        // Non-PBC case
        neighbours[0] = (col < L - 1) ? (m_spins[row * L + (col + 1)] ? 1 : -1) : 0; // Right
        neighbours[1] = (col > 0) ? (m_spins[row * L + (col - 1)] ? 1 : -1) : 0;     // Left
        neighbours[2] = (row > 0) ? (m_spins[(row - 1) * L + col] ? 1 : -1) : 0;     // Up
        neighbours[3] = (row < L - 1) ? (m_spins[(row + 1) * L + col] ? 1 : -1) : 0; // Down
    }

    return neighbours;
}

void SpinConfiguration::step() {
    m_engine.run();
}

void SpinConfiguration::sweep() {
    m_engine.sweep();
}

void SpinConfiguration::run(int MAX_TIME, int BURNIN, int THINNING) {
    // Run the simulation for a given number of steps
    if (!m_results.getMagnetizations().empty()) {
        m_results.getMagnetizations().clear();
    }
    if (!m_results.getEnergy().empty()) {
        m_results.getEnergy().clear();
    }
    m_results.getMagnetizations().reserve(MAX_TIME);
    m_results.getEnergy().reserve(MAX_TIME);

    //std::cout << "\nProgress over (n. samples): " << MAX_TIME - BURNIN << std::endl;
    //std::cout << "\n-- Burn-in phase --\n";
    bool inBurnin = false;
    for (int time = 0; time < MAX_TIME; ++time) {
        for(int _= 0; _ < THINNING; _++) {
            this->sweep();
        }
        if (time >= BURNIN && m_keepTrackof[0]) {
            m_results.storeMagnetization(this->getMagnetization());
        }
        if (time >= BURNIN && m_keepTrackof[1]) {
            m_results.storeEnergy(this->getEnergy());
        }
        if (time == BURNIN){
            inBurnin = true;
            //std::cout << "\n-- Sampling phase --\n" ;
        }
        if (!inBurnin) {
            //std::cout << " " << time << " ";
        } else {
            //std::cout << " " << time - BURNIN << " ";
        }
    }
}

void SpinConfiguration::runGraphics(int MAX_TIME, int sizeBlur = 0) {
    int time = 0;
    int L = static_cast<int>(std::sqrt(m_N));
    int controllerLength = 400;
    int controllerHeight = 400;
    int upperSliderMargin = 50;
    int downMargin = 30;

    // Disable tracking by default in graphics mode
    this->keepTrack(false, false); 

    sf::RenderWindow window(sf::VideoMode(800, 800), "Ising Model");
    sf::RenderWindow controller(sf::VideoMode(controllerLength, controllerHeight), "Controller");
    window.setFramerateLimit(60);
    controller.setFramerateLimit(60);

    // Load the font
    sf::Font font;
    if (!font.loadFromFile("../res/Ubuntu-Regular.ttf")) {
        throw std::runtime_error("Could not load font");
    }

    //  SLIDERS\\
    // Dragging state for both sliders
    bool draggingT = false;
    bool draggingH = false;
    // Slider parameters
    float sliderWidth = 350.f;
    float sliderHeight = 3.f;
    float knobRadius = 5.f;
    float minH = -1.f, maxH = 1.f;
    float maxTemp = 5.f;
    // Slider positions (one above the other)
    sf::Vector2f sliderPosition((controllerLength - sliderWidth) / 2.f, upperSliderMargin);         // Temperature
    sf::Vector2f sliderHPosition((controllerLength - sliderWidth) / 2.f, 2*upperSliderMargin + sliderHeight); // Magnetic field

    sf::RectangleShape sliderBar(sf::Vector2f(sliderWidth, sliderHeight));
    sliderBar.setPosition(sliderPosition);
    sliderBar.setFillColor(sf::Color(200, 200, 200));
    

    sf::RectangleShape sliderHBar(sf::Vector2f(sliderWidth, sliderHeight));
    sliderHBar.setPosition(sliderHPosition);
    sliderHBar.setFillColor(sf::Color(200, 200, 200));
    controller.draw(sliderHBar);

    sf::RectangleShape zeroLine(sf::Vector2f(sliderWidth, 1.f));
    zeroLine.setPosition((controllerLength - sliderWidth) / 2.f, 3*upperSliderMargin + 0.5*(controllerHeight - 3*upperSliderMargin - downMargin));
    zeroLine.setFillColor(sf::Color(150, 150, 150, 100));

    // Double ended queue to store points for magnetization plot
    std::deque<sf::CircleShape> points;
    std::cout << "Progress over (n. samples): " << MAX_TIME << std::endl;
    while (window.isOpen() && controller.isOpen() && time < MAX_TIME) {
        std::cout << "Step: " << time << "\r";
        // Handle events in main window
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
                controller.close();
            }
        }

        // Handle events in controller window and modify the parameters if needed
        while (controller.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                controller.close();
                window.close();
            } else if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
                // Collect mouse position
                double mx = static_cast<double>(event.mouseButton.x);
                double my = static_cast<double>(event.mouseButton.y);

                // Check if it is inside the temp slider
                if (mx >= sliderPosition.x && mx <= sliderPosition.x + sliderWidth &&
                    my >= sliderPosition.y - 10 && my <= sliderPosition.y + sliderHeight + 10) {
                    draggingT = true;
                }
                // Check if it is inside the magnetic field slider
                else if (mx >= sliderHPosition.x && mx <= sliderHPosition.x + sliderWidth &&
                         my >= sliderHPosition.y - 10 && my <= sliderHPosition.y + sliderHeight + 10) {
                    draggingH = true;
                }
            }
            else if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left) {
                draggingT = false;
                draggingH = false;
            }
            else if (event.type == sf::Event::MouseMoved) {
                float mx = static_cast<float>(event.mouseMove.x);
                if (draggingT) {
                    mx = std::clamp(mx, sliderPosition.x, sliderPosition.x + sliderWidth);
                    float ratio = (mx - sliderPosition.x) / sliderWidth;
                    m_T = 0 + ratio * (maxTemp - 0);
                }
                else if (draggingH) {
                    mx = std::clamp(mx, sliderHPosition.x, sliderHPosition.x + sliderWidth);
                    float ratioH = (mx - sliderHPosition.x) / sliderWidth;
                    m_h = -1 + ratioH * (maxH - (minH));
                }
            }
        }

        // Draw the lattice
        window.clear(sf::Color::Black);
        int reducedL = L / (sizeBlur > 0 ? sizeBlur : 1);
        for (int rrow = 0; rrow < reducedL; ++rrow) {
            for (int rcol = 0; rcol < reducedL; ++rcol) {
                sf::RectangleShape cell(sf::Vector2f(800.f / reducedL, 800.f / reducedL));
                cell.setPosition(rcol * (800.f / reducedL), rrow * (800.f / reducedL));
                double avgSpin = 0.;
                for (int r = rcol * sizeBlur; r < (rcol + 1) * sizeBlur; ++r) {
                    for (int c = rrow * sizeBlur; c < (rrow + 1) * sizeBlur; ++c) {
                        int index = r * L + c;
                        int spin = m_spins[index] ? 1 : -1;
                        avgSpin += double(spin) / (sizeBlur * sizeBlur);
                    }
                }
                int colorValue = static_cast<int>((avgSpin + 1) * 127.5);
                colorValue = std::clamp(colorValue, 0, 255);
                sf::Color color(255 - colorValue, 0, colorValue);
                cell.setFillColor(color);
                window.draw(cell);
            }
        }
        
        // Perform the sweep (each frame is THINNING sweeps)
        for (int _ = 0; _ < 1; ++_) {
            this->sweep();
        }
        time++;
        // Main window is ready to be displayed
        window.display();

        // Update the controller window
        controller.clear(sf::Color(40, 40, 40));
        controller.draw(sliderBar);
        controller.draw(sliderHBar);
        
        // Draw the knob for the temperature slider
        double normalizedT = (m_T - 0.) / (maxTemp - 0.);
        sf::CircleShape sliderKnob(knobRadius);
        sliderKnob.setFillColor(sf::Color::White);
        sliderKnob.setOrigin(knobRadius, knobRadius);
        sliderKnob.setPosition(sliderPosition.x + normalizedT * sliderWidth, sliderPosition.y + sliderHeight / 2.f);
        controller.draw(sliderKnob);
        sf::Text textT;
        textT.setFont(font);
        textT.setCharacterSize(16);
        textT.setFillColor(sf::Color::White);
        textT.setString("T = " + std::to_string(m_T).substr(0, 4));
        textT.setPosition(sliderPosition.x, sliderPosition.y - 30.f);
        controller.draw(textT);

        // Draw the knob for the magnetic field slider
        double normalizedH = (m_h - (minH)) / (maxH - (minH));
        sf::CircleShape sliderHKnob(knobRadius);
        sliderHKnob.setFillColor(sf::Color::White);
        sliderHKnob.setOrigin(knobRadius, knobRadius);
        sliderHKnob.setPosition(sliderHPosition.x + normalizedH * sliderWidth, sliderHPosition.y + sliderHeight / 2.f);
        controller.draw(sliderHKnob);
        sf::Text textH;
        textH.setFont(font);
        textH.setCharacterSize(16);
        textH.setFillColor(sf::Color::White);
        textH.setString("H = " + std::to_string(m_h).substr(0, 4));
        textH.setPosition(sliderHPosition.x, sliderHPosition.y - 30.f);
        controller.draw(textH);


        // Magnetization plot
        sf::RectangleShape graphPlate(sf::Vector2f(sliderWidth, controllerHeight - 3*upperSliderMargin - downMargin));
        graphPlate.setPosition((controllerLength - sliderWidth) / 2.f, 2*upperSliderMargin + 2*sliderHeight + 40.f);
        graphPlate.setFillColor(sf::Color(70, 70, 70));
        controller.draw(graphPlate);    
        int pointsToBeDisplayed = 30;
        if (points.size() >= pointsToBeDisplayed){
            points.pop_front();
        }

        sf::CircleShape newPoint;
        newPoint.setRadius(5);
        newPoint.setOutlineThickness(1);
        newPoint.setOutlineColor(sf::Color(0,0,0));
        newPoint.setOrigin(0,newPoint.getRadius());
        sf::Color c = sf::Color::White;
        newPoint.setFillColor(c);
        newPoint.setPosition(sf::Vector2f((controllerLength + sliderWidth)/2, 3*upperSliderMargin + 0.5*graphPlate.getSize().y*(1 - this->getMagnetization())));
        points.push_back(newPoint);
            
        for(auto& p : points){
            double new_x = p.getPosition().x - (graphPlate.getSize().x)/pointsToBeDisplayed;
            p.setPosition(new_x, p.getPosition().y); 
            controller.draw(p);
        }    

        
        controller.draw(zeroLine);

        controller.display();
    }
}



double SpinConfiguration::getEnergy() const {
    double energy = 0.0;
    int L = static_cast<int>(std::sqrt(m_N));
    for (int r = 0; r < L; ++r) {
        for (int c = 0; c < L; ++c) {
            int i = r*L + c;
            int spin = m_spins[i] ? 1 : -1;
            
            // Considera solo due direzioni: destra e sotto
            int right = m_spins[r*L + ((c+1)%L)] ? 1 : -1;
            int down  = m_spins[((r+1)%L)*L + c] ? 1 : -1;
            
            energy -= spin * right;
            energy -= spin * down;
        }
    }
    return energy / static_cast<double>(m_N);
}