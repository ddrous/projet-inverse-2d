#include <iostream>
#include <iomanip>
#include <string>

#include "maille.hpp"

// using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;


/************************************************
 * ClasseMesh
 */ 
Mesh::Mesh(double a_hat, double b_hat, int N_hat){
    assert (a_hat < b_hat && N_hat > 0);
        // assert ( a_hat >= b_hat || N_hat <= 0 && "ERREUR: Vérifier que a<b et N>0");)
        // CHK(throw string("ERREUR: Vérifier que a<b et N>0"););
    
    a = a_hat; b = b_hat; N = N_hat;                    // n'oublie pas les asserts !!!
    cells = MatrixXd(N+2, 3);

}

void Mesh::create_cells(){
        
        // dx correspond a \Delta x
        dx = (b-a)/N;       

        for (int j = 0; j < N+2; j++){     // Les mailles fantomes sont indicees 0 et N
            cells(j, 0) = (j-1)*dx;
            cells(j, 1) = (j-1)*dx + dx/2.;
            cells(j, 2) = (j)*dx;         // Optimiser ceci, creer une 
        }

}

Mesh::~Mesh(){

}