#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>

#include "solveur.hpp"

// using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;

/************************************************
 * Constructeur
 */ 
Solver::Solver(Mesh *new_mesh,
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
                std::string new_T_N_t){
    mesh = new_mesh;
    a = new_a;
    C_v = new_C_v;
    // rho = VectorXd(mesh.N);
    rho_expr = new_rho;
    // sigma_a = VectorXd(mesh.N);
    sigma_a_expr = new_sigma_a;
    // sigma_c = VectorXd(mesh.N);
    sigma_c_expr = new_sigma_c;
    CFL = new_CFL;
    epsilon = new_epsilon;
    Temps_max = new_Temps_max;
    E_x_0_expr = new_E_x_0;                // Energie des photons
    E_0_t_expr = new_E_0_t;                // Sur le bord droit
    E_N_t_expr = new_E_N_t;                // Sur le bord gauche
    E = VectorXd(mesh->N+2);
    F_x_0_expr = new_F_x_0;                // Flux
    F_0_t_expr = new_F_0_t;                // Flux
    F_N_t_expr = new_F_N_t;                // Flux
    F = VectorXd(mesh->N+2);
    T_x_0_expr = new_T_x_0;
    T_0_t_expr = new_T_0_t;
    T_N_t_expr = new_T_N_t;
    T = VectorXd(mesh->N+2);
    // a l'nistant t
} 

/***************
 * Fonction pour calculer rho
 */
double rho(double x){
    // parse rho_str
    // return 1062;       // Densite moyenne du corps humain
    return 3821.4;       // Cas test de Olson-Auer-hall
}

/***************
 * Fonction pour calculer jour sigma_a
 */
double sigma_a(double rho, double T){
    // parse sigma_a_str
    return 299792458;       // Approximation de diffusion
}

/***************
 * Fonction pour calculer jour sigma_c
 */
double sigma_c(double rho, double T){
    // parse sigma_c_str
    return 299792458;       // Approximation de diffusion
}

/**
 * Calcule sigma_c_j_plus_un_demi
 */ 
double flux_sigma_c(double sigma_c_left, double sigma_c_right){
    return 0.5 * (sigma_c_left + sigma_c_right);
}

/**
 * Calcule M_j_plus_un_demi
 */ 
double flux_M(double Delta_x, double flux_sigma_c){
    return 2 / (2 + Delta_x * flux_sigma_c);
}

/**
 * Calcule E_j_plus_un_demi
 */ 
double flux_E(double flux_M, double E_left, double E_right, double F_left, double F_right){
    return flux_M * ((E_right + E_left)/2 - (F_right - F_left)/2);
}

/**
 * Calcule F_j_plus_un_demi
 */ 
double flux_F(double flux_M, double F_left, double F_right, double E_left, double E_right){
    return flux_M * ((F_right + F_left)/2 - (E_right - E_left)/2);
}

/**
 * Calcule E(x, 0), energie a la position x au temps initial
 */ 
double E_x_0(double x){
    return 0;
}

/**
 * Calcule E(x_0, t), energie sur le bord gauche, au temps t
 * Correspond a l'energie sur la maille fantome de gauche au temps t
 */ 
double E_0_t(double t){
    return 0;
}

/**
 * Calcule E(x_N, t), energie sur le bord droit, au temps t
 * Correspond a l'energie sur la maille fantome de droite au temps t
 */ 
double E_N_t(double t){
    return 0;
}

/**
 * Calcule F(x, 0)
 */ 
double F_x_0(double x){
    return 0;
}

/**
 * Calcule F(x_0, t)
 */ 
double F_0_t(double t){
    return 0;
}

/**
 * Calcule F(x_N, t)
 */ 
double F_N_t(double t){
    return 0;
}

/**
 * Calcule T(x, 0)      // Eviter T(x) == 0, mu_q devient inf 
 */ 
double T_x_0(double x){
    // return 0.56234 * 11.6*1e6;       // Cas test de Olson-Auer-hall
    // return 300 ;       // Cas test de Marshak lineaire 
    return (0.4<=x && x<=0.6)? 310:297 ;       // Interieur/exterieur du corps humain
}

/**
 * Calcule T(x_0, t)
 */ 
double T_0_t(double t){
    return 0;
}

/**
 * Calcule T(x_N, t)
 */ 
double T_N_t(double t){
    return 0;
}

/***************************************************
 * Utilise les etape 1 et 2 de facon iterative
 */
void Solver::solve(std::string nom_fichier){
    // Les constantes
    int N = mesh->N;
    double dx = mesh->dx;
    double c = 299792458;

    // Les varaibles pour l'etape 1
    double E_n, E_next, T_n, F_n, F_next, Theta, Theta_n, Theta_next;
    // Les variables pour l'etape 2
    VectorXd E_etoile(N+2), F_etoile(N+2), T_etoile(N+2);
    VectorXd E_suiv(N+2), F_suiv(N+2), T_suiv(N+2);

    /**
     * Initialisation de la boucle
     */
    // Maille fantome a gauche ==> meme valeur que la maille d'a cote
    E(0) = E_x_0(mesh->cells(1, 1));
    F(0) = F_x_0(mesh->cells(1, 1));
    T(0) = T_x_0(mesh->cells(1, 1));
    // Mailles du milieux
    for (int j = 1; j < N+1; j++){
        E(j) = E_x_0(mesh->cells(j, 1));
        F(j) = F_x_0(mesh->cells(j, 1));
        T(j) = T_x_0(mesh->cells(j, 1));
    }
    // Maille fantome a droite ==> meme valeur que la maille d'a cote
    E(N+1) = E_x_0(mesh->cells(N, 1));
    F(N+1) = F_x_0(mesh->cells(N, 1));
    T(N+1) = T_x_0(mesh->cells(N, 1));
    
    // Temps courant de delta t
    double t = 0;
    double dt = CFL*dx/c;

    // Pour l'criture dans le fichier
    ofstream log_file;
    // log_file.open(nom_fichier, ios_base::app);      // Ajout
    log_file.open(nom_fichier);      // Ecraser
    if(!log_file)
        throw string ("Erreur d'ouverture du fichier dat");
    log_file << "t, " << "E_0, " << "E_N, " << "F_0, " << "F_N, "  << "T_0, " << "T_N\n"; 
    log_file << t << ", " << E(1) << ", " << E(N) << ", " << F(1) << ", " << F(N) << ", " << T(1) << ", " << T(N) << "\n";

    /**
     * Boucle de resolution
     */
    while (t < Temps_max)
    {
        t += dt;

        /***********
         * Etape 1
         */
        // E = E_0;
        for (size_t j = 1; j < N+1; j++){
            // Initialisation etape 1
            E_n = E(j);
            F_n = F(j);
            T_n = T(j);
            Theta_n = a * pow(T_n, 4);
            Theta = Theta_n;
            // while (abs(E(j) - Theta) > epsilon){
            //     double rho_tmp = rho(mesh->cells(j, 1));
            //     double sigma_a_tmp = sigma_a(rho_tmp, T(j));
            //     double tmp_1 = (1/dt) + c*sigma_a_tmp;
            //     double alpha = 1/dt/tmp_1;
            //     double beta = c*sigma_a_tmp/tmp_1;
            //     double mu_q = 1/ (pow(T_n, 3) + T_n*pow(T(j), 2) + T(j)*pow(T_n, 2) + pow(T(j), 3));
            //     double tmp_2 = (rho_tmp*C_v*mu_q/dt) + c*sigma_a_tmp;
            //     double gamma = rho_tmp*C_v*mu_q/dt/tmp_2;
            //     double delta = c*sigma_a_tmp/tmp_2;

            //     E(j) = alpha*E_n + beta*Theta;
            //     F(j) = F_n;
            //     Theta = gamma*Theta_n + delta*E(j);

            //     T(j) = pow(Theta/a, 0.25);      // Necessaire pour les calculs
            // }

            E_next = E(j);
            F_next = F(j);
            Theta_next = Theta;
            do{
                E(j) = E_next;
                F(j) = F_next;
                Theta = Theta_next;
                T(j) = pow(Theta/a, 0.25);      // Necessaire pour les calculs

                double rho_tmp = rho(mesh->cells(j, 1));
                double sigma_a_tmp = sigma_a(rho_tmp, T(j));
                double tmp_1 = (1/dt) + c*sigma_a_tmp;
                double alpha = 1/dt/tmp_1;
                double beta = c*sigma_a_tmp/tmp_1;
                double mu_q = 1/ (pow(T_n, 3) + T_n*pow(T(j), 2) + T(j)*pow(T_n, 2) + pow(T(j), 3));
                double tmp_2 = (rho_tmp*C_v*mu_q/dt) + c*sigma_a_tmp;
                double gamma = rho_tmp*C_v*mu_q/dt/tmp_2;
                double delta = c*sigma_a_tmp/tmp_2;

                // E_next = alpha*E_n + beta*Theta;
                // F_next = F_n;
                // Theta_next = gamma*Theta_n + delta*E_next;

                E_next = (alpha*E_n + gamma*beta*Theta_n) / (1 - beta*delta);
                F_next = F_n;
                Theta_next = (gamma*Theta_n + alpha*delta*E_next) / (1-beta*delta);

            } while (abs(E_next-E(j)) > epsilon && abs(Theta_next-Theta) > epsilon);
            
        }
        
        /**
         * Remplissage des mailles fantomes
         */
        // for (int i = 1; i < N+1; i++)
        //     E(0) = E_0_t(t);
        E(0) = E(1);
        F(0) = F(1);
        T(0) = T(1);
        E(N+1) = E(N);
        F(N+1) = F(N);
        T(N+1) = T(N);

        // Initialisation etape 2
        E_etoile = E;
        F_etoile = F;
        T_etoile = T;

        /***********
         * Etape 2
         */
        for (size_t j = 1; j < N+1; j++){
            double x_left = mesh->cells(j-1, 1);
            double x_center = mesh->cells(j, 1);
            double x_right = mesh->cells(j+1, 1);
            double sigma_c_left = sigma_c(rho(x_left), T(j-1));
            double sigma_c_center = sigma_c(rho(x_center), T(j));
            double sigma_c_right = sigma_c(rho(x_right), T(j+1));

            double flux_sigma_c_left = flux_sigma_c(sigma_c_left, sigma_c_center);
            double flux_sigma_c_right = flux_sigma_c(sigma_c_center, sigma_c_right);
            double flux_M_left = flux_M(dx, flux_sigma_c_left);
            double flux_M_right = flux_M(dx, flux_sigma_c_right);

            double tmp = (1/dx) + (c/2)*(flux_M_right*flux_sigma_c_right + flux_M_left*flux_sigma_c_left);
            double Alpha = 1/dt/tmp;
            double Beta = c/dx/tmp;

            E_suiv(j) = E_etoile(j) - (c*dt/dx)*(flux_F(flux_M_right, F(j), F(j+1), E(j), E(j+1)) - flux_F(flux_M_left, F(j-1), F(j), E(j-1), E(j)));

            F_suiv(j) = Alpha*F_etoile(j) - Beta*(flux_E(flux_M_right, E(j), E(j+1), F(j), F(j+1)) - flux_E(flux_M_left, E(j-1), E(j), F(j-1), F(j)));

            T_suiv(j) = T_etoile(j);
        }

        E = E_suiv;
        F = F_suiv;
        T = T_suiv;

        // Ajout des solution aux bords a cet instant dans le fichier dat
        log_file << t << ", " << E(1) << ", " << E(N) << ", " << F(1) << ", " << F(N) << ", " << T(1) << ", " << T(N) << "\n";
    }

    // fermeture du fichier  
    log_file.close();

};

/***************************************************
 * Affiche sur la console
 */
void Solver::dislay(){
    cout << "\nE:\t" ;
    for (int j = 1; j < mesh->N+1; j++)
        cout << E(j) << "  ";

    cout << "\nF:\t" ;
    for (int j = 1; j < mesh->N+1; j++)
        cout << F(j) << "  ";

    cout << "\nT:\t" ;
    for (int j = 1; j < mesh->N+1; j++)
        cout << T(j) << "  ";

    cout << "\n";

    // cout << endl << E << endl << F << endl << T << endl;
};

/***************************************************
 * Export
 */
void Solver::export_csv(std::string nom_fichier){

    ofstream log_file(nom_fichier);
    if(log_file){
        log_file << "x" << ", " << "E" << ", " << "F"  << ", " << "T\n"; 
        for (int j = 1; j < mesh->N+1; j++)
            log_file << mesh->cells(j, 1) << ", " << E(j) << ", " << F(j) << ", " << T(j) << "\n"; 
        log_file.close();
    }else {throw string ("Erreur d'ouverture de fichier");}
};
