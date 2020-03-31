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
Solver::Solver(UniformMesh *new_mesh,
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
    E = VectorXd(mesh->N);
    F_x_0_expr = new_F_x_0;                // Flux
    F_0_t_expr = new_F_0_t;                // Flux
    F_N_t_expr = new_F_N_t;                // Flux
    F = VectorXd(mesh->N);
    T_x_0_expr = new_T_x_0;
    T_0_t_expr = new_T_0_t;
    T_N_t_expr = new_T_N_t;
    T = VectorXd(mesh->N);
    // a l'nistant t
} 


/***************
 * Fonction pour mettrecalculer rho
 */
double rho(double x){
    // parse rho_str
    return 10;       // calcule siga_a pour la maille j a l'instant T_n+1
}


/***************
 * Fonction pour mettre a jour sigma_a pour le temps T_n+1
 */
double sigma_a(double rho, double T){
    // parse sigma_a_str
    return 299792458;       // calcule siga_a pour la maille j a l'instant T_n+1
}

/***************
 * Fonction pour mettre a jour sigma_c pour le temps T_n+1
 */
double sigma_c(double rho, double T){
    // parse sigma_c_str
    return 299792458;       // Proche de c
}


/**
 * Calcule alpha pour la maille j
 */ 
double update_alpha(double dt, double c, double sigma_a){
    return (1/dt)*1/((1/dt)+c*sigma_a);
}

/**
 * Calcule alpha pour la maille j
 */ 
double update_beta(double dt, double c, double sigma_a){
    return 1;
}

/**
 * Calcule alpha pour la maille j
 */ 
double update_gamma(double dt, double c, double sigma_a, double rho, double mu_q){
    return 1;
}

/**
 * Calcule alpha pour la maille j
 */ 
double update_delta(double dt, double c, double sigma_a, double rho, double mu_q){
    return 1;
}

/**
 * Calcule mu_q
 */ 
double update_mu_q(double T_n){
    return 1;
}

// /**
//  * Mets a jour A
//  */ 
// void update_A_C(SparseMatrix<double> *A, int N, double alpha){
//     for (size_t i = 0; i < N-1; i++)
//         A->insert(i, i+1) = alpha;
//     for (size_t i = 1; i < N; i++)
//         A->insert(i, i-1) = alpha;
//     for (size_t i = 0; i < N; i++)
//         A->insert(i, i) = alpha*(-2+1/alpha);
// }

// /**
//  * Mets a jour A
//  */ 
// void update_B_D(SparseMatrix<double> *B, int N, double alpha){
//     for (size_t i = 0; i < N-1; i++)
//         B->insert(i, i+1) = -alpha;
//     for (size_t i = 1; i < N; i++)
//         B->insert(i, i-1) = alpha;
//     for (size_t i = 0; i < N; i++)
//         B->insert(i, i) = 0;
// }

// void update_Sigma_a(SparseMatrix<double> &Sigma_a, VectorXd &sigma_a, int N){
//     for (size_t j = 0; j < N; i++)
//         Sigma_a.coeffRef(j, j) = sigma_a(j);    
// }

/**
 * Calcule T_0_t
 */ 
double update_M(double dx, double sigma_c){
    return 1;
}

/**
 * Calcule T_0_t
 */ 
double update_Alpha(double M, double c, double dt, double dx){
    return 1;
}

/**
 * Update M
 */ 
double update_Beta(double M, double c, double dt, double dx, double sigma_c){
    return 1;
}


/**
 * Calcule E_x_0 au temps initial
 */ 
double E_x_0(double x){
    return 1;
}


/**
 * Calcule E_0_t, energie a la position initiale, au temps t
 */ 
double E_0_t(double t){
    return 1;
}

/**
 * Calcule E_N_t,
 */ 
double E_N_t(double t){
    return 1;
}

/**
 * Calcule 
 */ 
double F_x_0(double x){
    return 1;
}

/**
 * Calcule F_0_t
 */ 
double F_0_t(double t){
    return 1;
}

/**
 * Calcule F_N_t
 */ 
double F_N_t(double t){
    return 1;
}


/**
 * Calcule E_0_t, energie a la position initiale, au temps t
 */ 
double T_x_0(double x){
    return 1;
}

/**
 * Calcule T_0_t
 */ 
double T_0_t(double t){
    return 1;
}

/**
 * Calcule T_0_t
 */ 
double T_N_t(double t){
    return 1;
}



// /**
//  * Mets a jour G
//  */ 
// void update_G(VectorXd *G, double E_0_t){
//     G->coeffRef(0) = E_0_t;
// }


// /***************
//  * retourne sigma_a a tout instant
//  */
// double Solver::sigma_a(double rho, double T_n){
//     return 299792458;       // Proche de c
// }

// /***************
//  * retourne sigma_c a tout instant
//  */
// double Solver::sigma_c(double rho, double T_n){
//     return 299792458;       // Proche de c
// }

/***************************************************
 * Utilise les etape 1 et 2 de facon iterative
 */
void Solver::solve(){
    // Nombre de mailles
    int N = mesh->N;
    double dx = mesh->dx;
    // Coefficients pour la resolutions iterative de l'etape 1
    double Theta, rho_tmp, sigma_a_tmp, sigma_c_tmp, mu_q, alpha, beta, gamma, delta;
    // Matrices pour la resolutions iterative de l'etape 2
    double M, Alpha, Beta;
    // SparseMatrix<double> M(N, N), Sigma_a(N, N), Sigma_c(N, N), Alpha(N, N), Beta(N, N), A(N, N), B(N, N), C(N, N), D(N, N);
    VectorXd E_suiv(N), F_suiv(N), T_suiv(N);
    
    // Assemblage initial des matrices
    // M.reserve(VectorXi::Constant(N,1));
    // Sigma_a.reserve(VectorXi::Constant(N,1));
    // Sigma_c.reserve(VectorXi::Constant(N,1));
    // Alpha.reserve(VectorXi::Constant(N,1));
    // Beta.reserve(VectorXi::Constant(N,1));
    // A.reserve(VectorXi::Constant(N,3));
    // B.reserve(VectorXi::Constant(N,3));
    // C.reserve(VectorXi::Constant(N,3));
    // D.reserve(VectorXi::Constant(N,3));

    // Autres parametres
    double c = 299792458;       //m.s-1

    // Initialisation
    for (int j = 0; j < N; j++){
        E(j) = E_x_0(mesh->cells(j, 1));
        F(j) = F_x_0(mesh->cells(j, 1));
        T(j) = T_x_0(mesh->cells(j, 1));
    }
    double dt = CFL*dx/c;
    double t = 0;
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
        for (size_t j = 0; j < N; j++){
            rho_tmp = rho(mesh->cells(j, 1));
            Theta = a * std::pow(T(j), 4);
            while (std::abs(E(j) - Theta) > epsilon){
                sigma_a_tmp = sigma_a(rho_tmp, T(j));
                sigma_c_tmp = sigma_c(rho_tmp, T(j));       // Inutile
                mu_q = update_mu_q(t);
                alpha = update_alpha(dt, c, sigma_a_tmp);
                beta = update_beta(dt, c, sigma_a_tmp);
                gamma = update_gamma(dt, c, sigma_a_tmp, rho_tmp, mu_q);
                delta = update_delta(dt, c, sigma_a_tmp, rho_tmp, mu_q);

                E(j) = alpha*E(j) + beta*Theta;
                Theta = gamma*Theta + delta*E(j);

                T(j) = pow(Theta/4, 0.25);
            }

        }
        
        /***********
         * Etape 2
         */
        for (size_t j = 0; j < N; j++){
            sigma_a_tmp = sigma_a(rho(mesh->cells(j, 1)), T(j));       //Inutile
            sigma_c_tmp = sigma_c(rho(mesh->cells(j, 1)), T(j));
            M = update_M(sigma_c_tmp, dx);
            Alpha = update_Alpha(M, c, dt, dx);
            Beta = update_Beta(M, c, dt, dx, sigma_c_tmp);

            if (j == 0){
                E_suiv(j) = Alpha*(E_0_t(t) + (-2+1/Alpha)*E(j) + E(j+1)) + Alpha*(F_0_t(t) - F(j+1));
                F_suiv(j) = Beta*(F_0_t(t) + (-2+1/Beta)*F(j) + F(j+1)) + Beta*(E_0_t(t) - E(j+1));
            }else if(j == N-1){
                E_suiv(j) = Alpha*(E(j-1) + (-2+1/Alpha)*E(j) + E_N_t(t)) + Alpha*(F(j-1) - F_N_t(t));
                F_suiv(j) = Beta*(E(j-1) + (-2+1/Beta)*E(j) + F_N_t(t)) + Beta*(E(j-1) - E_N_t(t));
            } else{       
                E_suiv(j) = Alpha*(E(j-1) + (-2+1/Alpha)*E(j) + E(j+1)) + Alpha*(F(j-1) - F(j+1));  
                E_suiv(j) = Beta*(F(j-1) + (-2+1/Beta)*F(j) + F(j+1)) + Beta*(E(j-1) - E(j+1));  
            }

            T_suiv(j) = T(j);
        }
        // update_Sigma_a(Sigma_a, sigma_a, N);         // Place les sigma sur la diag
        // update_Sigma_c(Sigma_c, sigma_c, N);
        // update_M(&M, &Sigma_c, dx, N);
        // update_Alpha(&Alpha, &M, &Sigma_c, c, dt, dx, N);
        // update_Beta(&Beta, &M, &Sigma_c, N);

        E = E_suiv;
        F = F_suiv;
        T = T_suiv;
    }


};

/***************************************************
 * Affiche sur la console
 */
void Solver::dislay(){
    cout << endl << E << endl << F << endl << T << endl;
};

/***************************************************
 * Export
 */
void Solver::export_csv(std::string nom_fichier){

    ofstream log_file(nom_fichier);
    if(log_file){
        log_file << "x\t" <<  "E\t"  << "F\t"  << "T\n"; 
        for (int j = 0; j < mesh->N; j++)
            log_file << mesh->cells(j, 1) << "\t" << E(j) << "\t" << F(j) << "\t" << T(j) << "\n"; 
        log_file.close();
    }else {throw string ("Erreur d'ouverture de fichier");}
};