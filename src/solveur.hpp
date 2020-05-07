
#ifndef DEJA_INCLU_SOLVEUR     
#define DEJA_INCLU_SOLVEUR

#include <vector>
#include <fstream>

#include "maille.hpp"
#include "exprtk.hpp"       // Pour transformer les expressions en fonctions

/************************************************
 * Classe pour resoudre le probleme
 * 
 */
class Solver{
    public:
        // Maillage des volumes finis pour le probleme
        Mesh *mesh;

        // Parametres physiques
        double c;                   // Vitesse de la lumiere
        double a;                   // Constante de Boltzmann

        // Propriete du domaine
        double C_v;                 // Capacite calorifique a volume constant

        // Parametres du probleme
        double CFL;                 // Condition CFL: dt=CFL*dx/c
        double epsilon;             // Precision sur les calculs
        double Temps_max;           // Temps de la simulation

        // Autres proprietes du domaine
        std::string rho_expr;           // Densite du domaine
        std::string sigma_a_expr;       // Opacite d'absorption
        std::string sigma_c_expr;       // Opacite de scatering

        // Energies des photons
        std::vector<double> E;          // Solution du probleme: E(x, t)
        std::string E_x_0_expr;         // Au temps initial E(x, 0)
        std::string E_0_t_expr;         // Sur le bord droit: E(x_0, t)
        std::string E_N_t_expr;         // Sur le bord gauche: E(x_N, t)

        // Flux de photons
        std::vector<double> F;          // Solution du problme: F(x, t)
        std::string F_x_0_expr;         // Au temps initial F(x, 0)
        std::string F_0_t_expr;         // Sur le bord droit: F(x_0, t)
        std::string F_N_t_expr;         // Sur le bord gauche: F(x_N, t)

        // Temperatures du milieux
        std::vector<double> T;          // Solution du problme: T(x, t)
        std::string T_x_0_expr;         // Au temps initial T(x, 0)
        std::string T_0_t_expr;         // Sur le bord droit: T(x_0, t)
        std::string T_N_t_expr;         // Sur le bord gauche: T(x_N, t)

        // Fichiers ou sont exportes les donnees
        std::string export_1;            // Signaux en tout temps aux bords
        std::string export_2;            // Signaux au temps final en tout x

        /***************
         * Constructeur
         * @param doubles: les parametres qui sont doubles dans l'ordre ci-haut
         * @param strings: les autres parametres (expressions pour les fonctions
         * et noms des fichiers pour les resultats  a exporter)
         */
        Solver(Mesh *new_mesh, double *doubles, std::string *strings);

        /***************
         * Fonction pour calculer rho a partir de rho_expr
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
         * Calcule E(x, 0), energie a la position x au temps initial
         */ 
        double E_x_0(double x);

        /***************
         * Calcule E(x_0, t), energie sur le bord gauche, au temps t
         * Correspond a l'energie sur la maille fantome de gauche au temps t
         */ 
        double E_0_t(double t);

        /***************
         * Calcule E(x_N, t), energie sur le bord droit, au temps t
         * Correspond a l'energie sur la maille fantome de droite au temps t
         */ 
        double E_N_t(double t);

        /***************
         * Calcule F(x, 0)
         */ 
        double F_x_0(double x);

        /***************
         * Calcule F(x_0, t)
         */ 
        double F_0_t(double t);

        /***************
         * Calcule F(x_N, t)
         */ 
        double F_N_t(double t);

        /***************
         * Calcule T(x, 0)      // Eviter T(x) == 0, mu_q devient inf 
         */ 
        double T_x_0(double x);

        /***************
         * Calcule T(x_0, t)
         */ 
        double T_0_t(double t);

        /***************
         * Calcule T(x_N, t)
         */ 
        double T_N_t(double t);

        /***************
         * Resout le probleme et exporte les resultats au fur et a mesure
         */
        void solve();

        /***************
         * Affiche les resultats sur la console
         */
        void dislay();

        /***************
         * Exporte les resultats au temps courant aux bords du domaine
         * @param file: le fichier (deja ouvert) dans lequel on ecrit
         */
        void export_current(double t, std::ofstream &file);

        /***************
         * Exporte les resultats au temps final sur tout le domaine
         */
        void export_final();

        /***************
         * Destructeur vide
         */
        virtual ~ Solver(){};
};

#endif
