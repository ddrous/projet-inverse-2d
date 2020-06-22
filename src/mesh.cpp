#include <string>

#include "include/mesh.hpp"

using namespace std;


Mesh::Mesh(const Config &cfg){
    x_min = atof(cfg.values.at("x_min").c_str());
    x_max = atof(cfg.values.at("x_max").c_str());
    N = atoi(cfg.values.at("N").c_str());

    if (x_min >= x_max)
        throw std::string("ERREUR: Vérifiez que x_min < x_max");
    if (N <= 0)
        throw std::string("ERREUR: Vérifiez que N > 0");

    dx = (x_max - x_min) / N;

    cells = new double*[N+2];
    for (int j = 0; j < N+2; j++)
        cells[j] = new double[3];
}


void Mesh::create_cells(){
    // Les mailles fantomes sont indicees 0 et N+1
    for (int j = 0; j < N+2; j++){
        cells[j][0] = (j-1)*dx;             // extremite gauche de la maille j
        cells[j][1] = (j-1)*dx + dx/2.;     // centre de la maille j
        cells[j][2] = (j)*dx;               // extremite droite de la maille j
    }
}


Mesh::~Mesh(){
    for (int j = 0; j < N+2; j++)
        delete[] cells[j];
    delete[] cells;
}