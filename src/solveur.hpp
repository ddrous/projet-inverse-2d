#include "maille.hpp"

/***********************************************
 * Attention si le fichier maille.hpp est deja inclu
 */ 
#ifndef DEJA_INCLU_SOLVEUR     
#define DEJA_INCLU_SOLVEUR

#include <Eigen/Dense>
// #include <Eigen/Sparse>


/************************************************
 * Classe pour créer le solveur
 * 
 */
class Solver{
    public:
        // Maillage du probleme
        Mesh *mesh;
        // Parametres physiques du probleme
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
        // // Coefficients pour la resolutions iterative de l'etape 1
        // double alpha, beta, gamma, delta;
        // // Matrices pour la resolutions iterative de l'etape 2
        // Eigen::SparseMatrix<double> A, B, C, D;
        // Eigen::VectorXd G_1, G_2, G_3, G_4;


        /***************
         * Constructeur
         */
        Solver(Mesh *new_mesh,
                double new_a,
                double new_C_v,
                std::string new_rho,
                std::string new_sigma_a,
                std::string new_sigma_c,
                double new_CFL,
                double new_epsilon,
                double new_Temps_max,
                std::string new_E_x_0,
                std::string new_E_0_t,
                std::string new_E_N_t,
                std::string new_F_x_0,
                std::string new_F_0_t,
                std::string new_F_N_t,
                std::string new_T_x_0,        // tests sur la validite des donnees
                std::string new_T_0_t,        // tests sur la validite des donnees
                std::string new_T_N_t);        // tests sur la validite des donnees

        // /***************
        //  * Initialise les talbeaux du solveur
        //  */
        // void init_solve();

        // /***************
        //  * retourne sigma_a a tout instant
        //  */
        // double sigma_a(double rho, double T_n);

        // /***************
        //  * retourne sigma_c a tout instant
        //  */
        // double sigma_c(double rho, double T_n);

        /***************
         * Resout de facon iterative
         */
        void solve(std::string nom_fichier);

        /***************
         * Affiche les resultats sur la console
         */
        void dislay();


        /***************
         * Exporte les resultats au format CSV
         */
        void export_csv(std::string nom_fichier);


        /***************
         * Destructeur
         */
        virtual ~ Solver(){};

};

#endif