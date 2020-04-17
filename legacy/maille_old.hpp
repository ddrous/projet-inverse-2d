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

    // public:
        // Constructeur
        Mesh(double a_hat, double b_hat, int N_hat);
        
        // Destructeur
        virtual ~ Mesh(){};

};


/**************************************************
 * Classe derivees pour les maillages uniformes
 */ 
class UniformMesh: public Mesh{
    // private:                                         // protected.. ceci apres!
    public:
        // Le volume de la maille/cellule
        double dx;       
        // Les differents mailles du maillage                    
        // double **cells;
        // try{
        // Eigen::MatrixXd cells = Eigen::MatrixXd(N, 3);
        Eigen::MatrixXd cells;
        // }catch(std::exception& e ){}
    // public:
        UniformMesh(double a, double b, int N);

        void create_cells();

        virtual ~UniformMesh();
};


#endif
