#include <string>

#include "maille.hpp"

using namespace std;


Mesh::Mesh(double a_hat, double b_hat, int N_hat){
    if (a > b || N < 0)
        throw std::string("ERREUR: Vérifier que a < b et N > 0");
    
    a = a_hat; b = b_hat; N = N_hat;

    // dx correspond a \Delta x
    dx = (b-a)/N;       

    cells = new double*[N+2];
    for (int j = 0; j < N+2; j++)
        cells[j] = new double[3];
}


void Mesh::create_cells(){
    // Les mailles fantomes sont indicees 0 et N+1
    for (int j = 0; j < N+2; j++){
        cells[j][0] = (j-1)*dx;
        cells[j][1] = (j-1)*dx + dx/2.;
        cells[j][2] = (j)*dx;
    }
}


Mesh::~Mesh(){
    for (int j = 0; j < N+2; j++)
        delete[] cells[j];
    delete[] cells;
}