#include <SFML/Graphics.hpp>
#include "NDMolDyn.hpp"
#include <iostream>

// Parametri Grafici
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;

int main() {

    int N = 200;
    double rho = 0.05;
    double L = std::sqrt(static_cast<double>(N) / rho);
    double r_cut = 1;
    double skin = 0.5;
    int M = 10;
    double max_disp = 0.2;
    NDMolDyn<2> sim(N, L, r_cut, skin, M, max_disp);

    
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Hard Spheres MC Simulation");
    window.setFramerateLimit(60); 

    // Fattore di scala: Mondo Fisico -> Pixel
    // Es: se L=10 e Window=800, allora 1 unità fisica = 80 pixel
    float scale = static_cast<float>(WINDOW_WIDTH) / L;
    
    // Raggio grafico (sigma è il diametro!)
    float radiusPixel = (1.0 / 2.0) * scale;

    // Creiamo una forma base da riutilizzare
    sf::CircleShape circle(radiusPixel);
    circle.setFillColor(sf::Color(0, 255, 255, 200)); // Ciano semi-trasparente
    circle.setOutlineThickness(1.0f);
    circle.setOutlineColor(sf::Color::White);

    // 3. Loop Principale
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        for(int i=0; i < N; ++i) { // A single MC step every frame!
            sim.singleMCStep();
        }
        std::cout << "\rAcceptance ratio: " << static_cast<double>(sim.m_acceptCounter)/double(N);
        sim.m_acceptCounter = 0;


        // --- RENDER ---
        window.clear(sf::Color(20, 20, 30)); // Sfondo scuro

        // Ottieni posizioni attuali
        std::vector<std::array<double, 2>> positions = sim.getPositions();

        for (const auto& pos : positions) {
            // Conversione coordinate: Fisiche -> Grafiche
            float x_screen = pos[0] * scale;
            float y_screen = pos[1] * scale;

            // SFML disegna le sfere partendo dall'angolo in alto a sinistra del cerchio.
            // Dobbiamo sottrarre il raggio per centrarle sulla coordinata fisica.
            circle.setPosition(x_screen - radiusPixel, y_screen - radiusPixel);
            
            window.draw(circle);

            // --- GESTIONE GRAFICA PBC (Effetto Pac-Man) ---
            // Se una sfera è sul bordo, disegniamo le sue "copie fantasma" 
            // per vedere l'attraversamento fluido.
            
            bool drawGhost = false;
            float ghostX = x_screen;
            float ghostY = y_screen;

            // Check X axis
            if (pos[0] < 1.0) { ghostX += WINDOW_WIDTH; drawGhost = true; }
            else if (pos[0] > L - 1.0) { ghostX -= WINDOW_WIDTH; drawGhost = true; }

            // Check Y axis
            if (pos[1] < 1.0) { ghostY += WINDOW_HEIGHT; drawGhost = true; }
            else if (pos[1] > L - 1.0) { ghostY -= WINDOW_HEIGHT; drawGhost = true; }

            if (drawGhost) {
                // Disegna solo se necessario per risparmiare risorse
                circle.setPosition(ghostX - radiusPixel, ghostY - radiusPixel);
                window.draw(circle);
            }
        }

        window.display();
    }

    return 0;
}