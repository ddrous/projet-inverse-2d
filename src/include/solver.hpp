
#ifndef INCLUDED_SOLVER
#define INCLUDED_SOLVER

#include <vector>
#include <map>

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
        double t_0;                 // Temps de depart de la simulation
        double t_f;                 // Temps de la simulation

        // Autres proprietes du domaine
        std::string rho_expr;           // Expression de la densite du domaine
        std::string sigma_a_expr;       // Expression de l'opacite d'absorption
        std::string sigma_c_expr;       // Expression de l'opacite de scatering

        // Solution numerique
        std::vector<double> E;              // Solution numerique du problme: E(t, x)
        std::vector<double> F;              // Solution numerique du problme: F(t, x)
        std::vector<double> T;              // Solution numerique du problme: T(t, x)

        // Au temps initial
        std::string E_0_expr;               // Au temps initial: E(t_0, x)
        std::string F_0_expr;               // Au temps initial: F(t_0, x)
        std::string T_0_expr;               // Au temps initial: T(t_0, x)

        // Conditions imposees sur le bord gauche
        std::string E_l_expr;               // E(t, x_min)
        std::string F_l_expr;               // F(t, x_min)
        std::string T_l_expr;               // T(t, x_min)

        // Conditions imposees sur le bord droit
        std::string E_r_expr;               // E(t, x_max)
        std::string F_r_expr;               // F(t, x_max)
        std::string T_r_expr;               // T(t, x_max)

        // Solution exacte
        std::string E_exact_expr;           // Solution exacte: E_exact(t, x)
        std::string F_exact_expr;           // Solution exacte: F_exact(t, x)
        std::string T_exact_expr;           // Solution exacte: T_exact(t, x)

        // Resultats a exporter a gauche
        std::vector<double> E_left;         // Sur le bord gauche: E(t, x_0)
        std::vector<double> F_left;         // Sur le bord gauche: F(t, x_0)
        std::vector<double> T_left;         // Sur le bord gauche: T(t, x_0)

        // Resultats a exporter a gauche
        std::vector<double> E_right;        // Sur le bord droit: E(t, n_N)
        std::vector<double> F_right;        // Sur le bord droit: F(t, x_N)
        std::vector<double> T_right;        // Sur le bord droit: T(t, x_N)

        // Evolutions a observer
        double **E_evol;            // E sur tout les bords et au centre en tout temps
        double **T_evol;            // T sur tout le domaine en 5 temps precis

        // Autres parametres du probleme
        double dt;                      // Pas de temps
        int step_count;                 // Nombre d'iterations en temps
        std::vector<double> time_steps; // Les temps a chaque pas

        /***************
         * Constructeur
         * @param params: Tous les parametres du probleme
         */
        Solver(Mesh *new_mesh, std::map<std::string, std::string> &params);

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
         * Calcule E(t_0, x), energie a la position x au temps initial
         */ 
        double E_0(double x);

        /***************
         * Calcule F(t_0, x)
         */ 
        double F_0(double x);

        /***************
         * Calcule T(t_0, x)        // Eviter T_0 = 0 !
         */ 
        double T_0(double x);

        /***************
         * Calcule E(t, x_min), energie a la position x_min en tout temps
         */ 
        double E_l(double t);

        /***************
         * Calcule F(t, x_min), energie a la position x_min en tout temps
         */ 
        double F_l(double t);

        /***************
         * Calcule T(t, x_min), energie a la position x_min en tout temps
         */ 
        double T_l(double t);

        /***************
         * Calcule E(t, x_max), energie a la position x_min en tout temps
         */ 
        double E_r(double t);

        /***************
         * Calcule F(t, x_max), energie a la position x_min en tout temps
         */ 
        double F_r(double t);

        /***************
         * Calcule T(t, x_max), energie a la position x_min en tout temps
         */ 
        double T_r(double t);

        /***************
         * Calcule E(t, x), solution exacte
         */ 
        double E_exact(double t, double x);

        /***************
         * Calcule F(t, x), solution exacte
         */ 
        double F_exact(double t, double x);

        /***************
         * Calcule T(t, x), solution exacte
         */ 
        double T_exact(double t, double x);

        /***************
         * Resout le probleme et exporte les resultats au fur et a mesure
         */
        void solve();

        /***************
         * Affiche les resultats sur la console
         */
        void display();

        /***************
         * Destructeur vide
         */
        virtual ~ Solver();
};

#endif
