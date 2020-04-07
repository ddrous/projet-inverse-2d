/***********************************************
 * Attention si le fichier maille.hpp est deja inclu
 */ 
#ifndef DEJA_INCLU_MAILLE     
#define DEJA_INCLU_MAILLE

#include <Eigen/Dense>
// using Eigen::MatrixXd;
// using namespace Eigen;


/************************************************
 * Classe pour créer le maillage, apres on aura solver
 * 
 */
class Mesh{
    // protected:                                       // Accesseurs mutateurs
    public:
        // Parametres de depart
        double a;
        double b;
        int N;
        double dx;       
        Eigen::MatrixXd cells;

    // public:

        // Constructeur
        Mesh(double a_hat, double b_hat, int N_hat);

        void create_cells();
        
        // Destructeur
        virtual ~ Mesh();

};

#endif
