#include "maille.hpp"

/***********************************************
 * Attention si le fichier maille.hpp est deja inclu
 */ 
#ifndef DEJA_INCLU_SOLVEUR     
#define DEJA_INCLU_SOLVEUR

#include <Eigen/Dense>


/************************************************
 * Classe pour créer le solveur
 * 
 */
class Solver{
    public:
        // Maillage du probleme
        Mesh *mesh;
        // Destination des donnes
        std::string export_1;            // Signaux en tout temps aux bords
        std::string export_2;            // Signaux au temps final en tout x
        // Parametres physiques du probleme
        double c;
        double a;
        double C_v;
        std::string rho_expr;         // Suppposons constant pour la premiere partie
        std::string sigma_a_expr;     // Faire une fonction(rho, T) pour chaque sigma
        std::string sigma_c_expr;
        // Autres Parametres
        double CFL;                 // Valeur de la CFL; dt=CFL*dx/c
        double epsilon;              // Precision sur E-Theta
        double Temps_max;               // Temps final
        // Variables du solveur
        std::string E_x_0_expr;                // Energie des photons
        std::string E_0_t_expr;                // Sur le bord droit
        std::string E_N_t_expr;                // Sur le bord gauche
        Eigen::VectorXd E;                
        std::string F_x_0_expr;                // Flux
        std::string F_0_t_expr;                // Flux
        std::string F_N_t_expr;                // Flux
        Eigen::VectorXd F;                
        std::string T_x_0_expr;
        std::string T_0_t_expr;
        std::string T_N_t_expr;
        Eigen::VectorXd T;                // Temperature des photons


        /***************
         * Les fonctions-parametres utilsees
         */
        double rho(double x);
        double sigma_a(double rho, double T);
        double sigma_c(double rho, double T);
        double E_x_0(double x);
        double E_0_t(double t);
        double E_N_t(double t);
        double F_x_0(double x);
        double F_0_t(double t);
        double F_N_t(double t);
        double T_x_0(double x);
        double T_0_t(double t);
        double T_N_t(double t);

        /***************
         * Constructeur
         */
        Solver(Mesh *new_mesh, double *double_values, std::string *str_values);


        /***************
         * Resout le probleme et exporte les resultats au fur et a mesure
         */
        void solve();


        /***************
         * Affiche les resultats sur la console
         */
        void dislay();


        /***************
         * Exporte les resultats au temps final
         */
        void export_final();


        /***************
         * Destructeur
         */
        virtual ~ Solver(){};

};

#endif
