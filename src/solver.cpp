#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>

#include "solver.hpp"

#include "muParser.h"       // Pour transformer des expressions (chaines de caracteres) en fonctions

using namespace mu;
using namespace std;

Solver::Solver(Mesh *new_mesh, double *doubles, std::string *strings){
    mesh = new_mesh;

    c = doubles[0];
    a = doubles[1];

    C_v = doubles[2];

    CFL = doubles[3];
    epsilon = doubles[4];
    t_final = doubles[5];

    dt = CFL * mesh->dx/c;
    time_steps = (int)(t_final / dt) + 1;

    rho_expr = strings[0];
    sigma_a_expr = strings[1];
    sigma_c_expr = strings[2];

    E = vector<double>(mesh->N+2);
    E_x_0_expr = strings[3];
    E_0 = vector<double>(time_steps);
    E_N = vector<double>(time_steps);

    F = vector<double>(mesh->N+2);
    F_x_0_expr = strings[4];
    F_0 = vector<double>(time_steps);
    F_N = vector<double>(time_steps);

    T = vector<double>(mesh->N+2);
    T_x_0_expr = strings[5];
    T_0 = vector<double>(time_steps);
    T_N = vector<double>(time_steps);

    export_1 = strings[6];
    export_2 = strings[7];

    // Validite des proprietes du probleme
    if (CFL <= 0)
        throw std::string("ERREUR: Vérifiez que CFL > 0");
    if (epsilon <= 0)
        throw std::string("ERREUR: Vérifiez que epsilon > 0");
    if (t_final <= 0)
        throw std::string("ERREUR: Vérifiez que t_final > 0");
} 


double Solver::rho(double x){
    static int first_call = 1;
    static Parser p;

    if (first_call == 1){           // Seulement au 1er appel
        p.DefineVar("x", &x); 
        p.SetExpr(rho_expr);
    }

    return p.Eval();
}


double Solver::sigma_a(double rho, double T){
    static int first_call = 1;
    static Parser p;

    if (first_call == 1){
        p.DefineVar("rho", &rho); 
        p.DefineVar("T", &T); 
        p.SetExpr(sigma_a_expr);
    }

    return p.Eval();
}


double Solver::sigma_c(double rho, double T){
    static int first_call = 1;
    static Parser p;

    if (first_call == 1){
        p.DefineVar("rho", &rho); 
        p.DefineVar("T", &T); 
        p.SetExpr(sigma_c_expr);
    }

    return p.Eval();
}


double Solver::E_x_0(double x){
    static int first_call = 1;
    static Parser p;

    if (first_call == 1){
        p.DefineVar("x", &x); 
        p.SetExpr(E_x_0_expr);
    }

    return p.Eval();
}


double Solver::F_x_0(double x){
    static int first_call = 1;
    static Parser p;

    if (first_call == 1){
        p.DefineVar("x", &x); 
        p.SetExpr(F_x_0_expr);
    }

    return p.Eval();
}


double Solver::T_x_0(double x){
    static int first_call = 1;
    static Parser p;

    if (first_call == 1){
        p.DefineVar("x", &x); 
        p.SetExpr(T_x_0_expr);
    }

    return p.Eval();
}


/**
 * Calcule les flux sigma_c_{j+1/2} et sigma_c_{j-1/2}
 */ 
double flux_sigma_c(double sigma_c_left, double sigma_c_right){
    return 0.5 * (sigma_c_left + sigma_c_right);
}


/**
 * Calcule les flux M_{j+1/2} et M_{j-1/2}
 */ 
double flux_M(double Delta_x, double flux_sigma_c){
    return 2 / (2 + Delta_x * flux_sigma_c);
}


/**
 * Calcule les flux E_{j+1/2} et E_{j-1/2}
 */ 
double flux_E(double flux_M, double E_left, double E_right, double F_left, double F_right){
    return flux_M * ((E_right + E_left)/2 - (F_right - F_left)/2);
}


/**
 * Calcule les flux F_{j+1/2} et F_{j-1/2}
 */ 
double flux_F(double flux_M, double F_left, double F_right, double E_left, double E_right){
    return flux_M * ((F_right + F_left)/2 - (E_right - E_left)/2);
}


void Solver::solve(){
    // Raccourcissons les noms des constantes
    int N = mesh->N;            // Nombre effectif de mailles
    double dx = mesh->dx;       // Delta x

    // Les variables pour l'etape 1
    double E_n, E_next, T_n, F_n, F_next, Theta, Theta_n, Theta_next;
    
    // Les variables pour l'etape 2
    vector<double> E_etoile(N+2), F_etoile(N+2), T_etoile(N+2);
    vector<double> E_suiv(N+2), F_suiv(N+2), T_suiv(N+2);

    // Initialisation de la doucle de resolution
    for (int j = 1; j < N+1; j++){
        double x = mesh->cells[j][1];           // Centre de la maille
        E[j] = E_x_0(x);
        F[j] = F_x_0(x);
        T[j] = T_x_0(x);
    }

    // Temps courant et indice de l'iteration
    double t = 0;
    int step = -1;      // -1 pour une initialisation a 0

    /**
     * Boucle de resolution
     */
    while (t < t_final){
        t += dt;
        step += 1;

        // Signaux aux bords pour ce pas d'iteration en vue de l'export
        E_0[step] = E[1];
        F_0[step] = F[1];
        T_0[step] = T[1];
        E_N[step] = E[N];
        F_N[step] = F[N];
        T_N[step] = T[N];

        for (int j = 1; j < N+1; j++){
            // Initialisation etape 1
            E_n = E[j];
            F_n = F[j];
            T_n = T[j];
            Theta_n = a * pow(T_n, 4);
            Theta = Theta_n;

            E_next = E[j];
            F_next = F[j];
            Theta_next = Theta;

            /*********************************
             * Etape 1
             */
            do{
                E[j] = E_next;
                F[j] = F_next;
                Theta = Theta_next;
                T[j] = pow(Theta/a, 0.25);      // Necessaire pour les calculs

                double rho_tmp = rho(mesh->cells[j][1]);
                double sigma_a_tmp = sigma_a(rho_tmp, T[j]);
                double tmp_1 = (1/dt) + c*sigma_a_tmp;
                double alpha = 1/dt/tmp_1;
                double beta = c*sigma_a_tmp/tmp_1;
                double mu_q = 1/ (pow(T_n, 3) + T_n*pow(T[j], 2) + T[j]*pow(T_n, 2) + pow(T[j], 3));
                double tmp_2 = (rho_tmp*C_v*mu_q/dt) + c*sigma_a_tmp;
                double gamma = rho_tmp*C_v*mu_q/dt/tmp_2;
                double delta = c*sigma_a_tmp/tmp_2;

                E_next = (alpha*E_n + gamma*beta*Theta_n) / (1 - beta*delta);
                F_next = F_n;
                Theta_next = (gamma*Theta_n + alpha*delta*E_n) / (1-beta*delta);

            } while (abs(E_next-E[j]) > epsilon && abs(Theta_next-Theta) > epsilon);
        }
        
        // Initialisation etape 2
        E_etoile = E;
        F_etoile = F;
        T_etoile = T;

        // Remplissage des mailles fantomes -> neumann naturel
        E[0] = E[1];
        F[0] = F[1];
        T[0] = T[1];
        E[N+1] = E[N];
        F[N+1] = F[N];
        T[N+1] = T[N];

        /***********************************
         * Etape 2
         */
        for (int j = 1; j < N+1; j++){
            double x_left = mesh->cells[j-1][1];
            double x_center = mesh->cells[j][1];
            double x_right = mesh->cells[j+1][1];

            double sigma_c_left = sigma_c(rho(x_left), T[j-1]);
            double sigma_c_center = sigma_c(rho(x_center), T[j]);
            double sigma_c_right = sigma_c(rho(x_right), T[j+1]);

            double flux_sigma_c_left = flux_sigma_c(sigma_c_left, sigma_c_center);
            double flux_sigma_c_right = flux_sigma_c(sigma_c_center, sigma_c_right);
            double flux_M_left = flux_M(dx, flux_sigma_c_left);
            double flux_M_right = flux_M(dx, flux_sigma_c_right);

            double tmp = (1/dt) + (c/2)*(flux_M_right*flux_sigma_c_right + flux_M_left*flux_sigma_c_left);
            double Alpha = 1/dt/tmp;
            double Beta = c/dx/tmp;

            E_suiv[j] = E_etoile[j] - (c*dt/dx)*(flux_F(flux_M_right, F[j], F[j+1], E[j], E[j+1]) - flux_F(flux_M_left, F[j-1], F[j], E[j-1], E[j]));

            F_suiv[j] = Alpha*F_etoile[j] - Beta*(flux_E(flux_M_right, E[j], E[j+1], F[j], F[j+1]) - flux_E(flux_M_left, E[j-1], E[j], F[j-1], F[j]));

            T_suiv[j] = T_etoile[j];
        }

        E = E_suiv;
        F = F_suiv;
        T = T_suiv;
    }
};


void Solver::display(){
    cout << "E:\t" ;
    for (int j = 1; j < mesh->N+1; j++)
        cout << E[j] << "  ";

    cout << "\nF:\t" ;
    for (int j = 1; j < mesh->N+1; j++)
        cout << F[j] << "  ";

    cout << "\nT:\t" ;
    for (int j = 1; j < mesh->N+1; j++)
        cout << T[j] << "  ";

    cout << "\n";
};


void Solver::export_temporal(){
    ofstream file;
    file.open(export_1, ios_base::app);             // Ajout dans le fichier
    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + export_1 + "'");

    file << mesh->x_min << "," << mesh->x_max << "," << mesh->N << "," << c << "," << a << "," << C_v << ","<< CFL << "," << epsilon << "," << t_final << ",\"" << rho_expr << "\",\"" << sigma_a_expr << "\",\"" << sigma_c_expr << "\",\"" << E_x_0_expr << "\",\"" << F_x_0_expr << "\",\"" << T_x_0_expr << "\"," << dt << "," << time_steps << ",";

    file << "\"[";
    for (int n = 0; n < time_steps; n++){
        file << E_0[n];
        if (n != time_steps-1)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < time_steps; n++){
        file << E_N[n];
        if (n != time_steps-1)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < time_steps; n++){
        file << F_0[n];
        if (n != time_steps-1)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < time_steps; n++){
        file << F_N[n];
        if (n != time_steps-1)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < time_steps; n++){
        file << T_0[n];
        if (n != time_steps-1)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < time_steps; n++){
        file << T_N[n];
        if (n != time_steps-1)
            file << ", ";
    }
    file << "]\"\n";

    file.close();
};


void Solver::export_spatial(){
    ofstream file;
    file.open(export_2, ios_base::app);             // Ajout dans le fichier
    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + export_2 + "'");

    file << mesh->x_min << "," << mesh->x_max << "," << mesh->N << "," << c << "," << a << "," << C_v << ","<< CFL << "," << epsilon << "," << t_final << ",\"" << rho_expr << "\",\"" << sigma_a_expr << "\",\"" << sigma_c_expr << "\",\"" << E_x_0_expr << "\",\"" << F_x_0_expr << "\",\"" << T_x_0_expr << "\"," << dt << "," << time_steps << ",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << mesh->cells[j][1];
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << rho(mesh->cells[j][1]);
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << sigma_a(rho(mesh->cells[j][1]), T[j]);
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << sigma_c(rho(mesh->cells[j][1]), T[j]);
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";
     
    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << E_x_0(mesh->cells[j][1]);
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << E[j];
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << F_x_0(mesh->cells[j][1]);
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << F[j];
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << T_x_0(mesh->cells[j][1]);
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << T[j];
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\"\n";

    file.close();
};
