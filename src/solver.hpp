
#ifndef INCLUDED_SOLVER
#define INCLUDED_SOLVER

#include <vector>
#include <fstream>

#include "mesh.hpp"

/************************************************
 * Classe pour resoudre le probleme
 */
class Solver{
    public:
        // Maillage des volumes finis pour le probleme
        Mesh *mesh;

        // Parametres physiques
        double c;                   // Vitesse de la lumiere
        double a;                   // Constante de Stefen-Boltzmann

        // Propriete du domaine
        double C_v;                 // Capacite calorifique a volume constant

        // Parametres du probleme
        double CFL;                 // Condition CFL: dt=CFL*dx/c
        double precision;           // Precision sur les calculs
        double t_init;              // Temps de depart de la simulation
        double t_simu;              // Temps de la simulation

        // Autres proprietes du domaine
        std::string rho_expr;           // Expression de la densite du domaine
        std::string sigma_a_expr;       // Expression de l'opacite d'absorption
        std::string sigma_c_expr;       // Expression de l'opacite de scatering

        // Energies des photons
        std::vector<double> E;              // Solution numerique du problme: E(t, x)
        std::string E_exact_expr;           // Solution exacte: E_exact(t, x)
        std::string E_init_expr;            // Au temps initial: E(t_init, x)
        std::vector<double> E_left;         // Sur le bord gauche: E(t, x_0)
        std::vector<double> E_right;        // Sur le bord droit: E(t, n_N)

        // Flux de photons
        std::vector<double> F;              // Solution numerique du problme: F(t, x)
        std::string F_exact_expr;           // Solution exacte: F_exact(t, x)
        std::string F_init_expr;            // Au temps initial: F(t_init, x)
        std::vector<double> F_left;         // Sur le bord gauche: F(t, x_0)
        std::vector<double> F_right;        // Sur le bord droit: F(t, x_N)

        // Temperatures du milieux
        std::vector<double> T;              // Solution numerique du problme: T(t, x)
        std::string T_exact_expr;           // Solution exacte: T_exact(t, x)
        std::string T_init_expr;            // Au temps initial: T(t_init, x)
        std::vector<double> T_left;         // Sur le bord gauche: T(t, x_0)
        std::vector<double> T_right;        // Sur le bord droit: T(t, x_N)

        // Autres parametres du probleme
        double dt;                      // Pas de temps
        int step_count;                 // Nombre d'iterations en temps
        std::vector<double> time_steps; // Les temps a chaque pas

        /***************
         * Constructeur
         * @param params: les parametres du solveur autres que le maillage
         */
        Solver(Mesh *new_mesh, std::string *params);

        /***************
         * Fonction pour calculer rho a partir de son expression rho_expr
         */
        double rho(double x);

        /***************
         * Fonction pour calculer sigma_a
         */
        double sigma_a(double rho, double T);

        /***************
         * Fonction pour calculer sigma_c
         */        
        double sigma_c(double rho, double T);

        /***************
         * Calcule E(t, x), solution exacte
         */ 
        double E_exact(double t, double x);

        /***************
         * Calcule E(t_init, x), energie a la position x au temps initial
         */ 
        double E_init(double x);

        /***************
         * Calcule F(t, x), solution exacte
         */ 
        double F_exact(double t, double x);

        /***************
         * Calcule F(t_init, x)
         */ 
        double F_init(double x);

        /***************
         * Calcule T(t, x), solution exacte
         */ 
        double T_exact(double t, double x);

        /***************
         * Calcule T(t_init, x)
         */ 
        double T_init(double x);

        /***************
         * Resout le probleme et exporte les resultats au fur et a mesure
         */
        void solve();

        /***************
         * Affiche les resultats sur la console
         */
        void display();

        /***************
         * Exporte les resultats en tout temps aux bords du domaine
         */
        void export_temporal(std::string file_name);

        /***************
         * Exporte les resultats sur tout le domaine au temps final
         */
        void export_spatial(std::string file_name);

        /***************
         * Destructeur vide
         */
        virtual ~ Solver(){};
};

#endif
