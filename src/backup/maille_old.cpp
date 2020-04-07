#include <iostream>
#include <iomanip>
#include <string>

#include "maille.hpp"

// using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;

// /************************************************
//  * Macro pour vérifier les appels aux fonctions
//  * 
//  */
// #define CHK(v)         \
//     try {v}           \
//     catch(const std::exception& e){                  \
//         std::cerr << #v << std::endl << e.what() << std::endl;       \
//         throw e;        \
//     } while (0)


// /************************************************
//  * Vérifie les appels aux fonctions qui renvoient -1 en cas d'erreur
//  * 
//  */
// #define CHK_PS(v)         \
//     do {                  \
//         if ((v) == -1)    \
//             msg_erreur(1, #v); \
//     } while (0)

// /************************************************
//  * Vérifie les allocations de memoire
//  * 
//  */
// #define CHK_MA(v)         \
//     do {                  \
//         if ((v) == NULL)    \
//             msg_erreur(1, #v); \
//     } while (0)

// /************************************************
//  * Fonction pour indiquer une erreur et exit()
//  */
// void msg_erreur(int syserr, std::string msg) {
//     std::cerr << msg << std::endl;
//     if (syserr)
//         perror("");
//     exit(1);            // Apres il ne faudra pas exit ....
// }



/************************************************
 * ClasseMesh
 */ 
Mesh::Mesh(double a_hat, double b_hat, int N_hat){
    assert (a_hat < b_hat && N_hat > 0);
        // assert ( a_hat >= b_hat || N_hat <= 0 && "ERREUR: Vérifier que a<b et N>0");)
        // CHK(throw string("ERREUR: Vérifier que a<b et N>0"););
    
    a = a_hat; b = b_hat; N = N_hat;                    // n'oublie pas les asserts !!!
}



/**************************************************
 * Classe UniformMesh
 */ 
UniformMesh::UniformMesh(double a, double b, int N): Mesh(a, b, N){
    /**
     * Pour les N volumes, Chaque volumes possede 3 characteristiques
     * l'extremite gauche -> cells[j][0]
     * Le centre -> cells[j][1]
     * l'extremite droite -> cells[j][2]
     */ 
    cells = MatrixXd(N, 3);
}

void UniformMesh::create_cells(){
        
        // dx correspond a \Delta x
        dx = (b-a)/N;       

        for (size_t j = 0; j < N; j++){
            cells(j, 0) = j*dx;
            cells(j, 1) = j*dx + dx/2.;
            cells(j, 2) = (j+1)*dx;         // Optimiser ceci, creer une 
        }

}

UniformMesh::~UniformMesh(){
    // if (cells != NULL){
    // //     for (size_t j = 0; j < N; j++)
    // //         delete[] cells[j];
    // //     delete[] cells;   
    // }   
}