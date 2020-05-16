#include <string>

#include "mesh.hpp"

using namespace std;


Mesh::Mesh(double x_left, double x_right, int size){
    x_min = x_left; x_max = x_right; N = size;
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