#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>

#include "solver.hpp"

#include "muParser.h"       // Pour transformer des expressions (chaines de caracteres) en fonctions

using namespace mu;
using namespace std;

Solver::Solver(Mesh *new_mesh, std::string *params){
    mesh = new_mesh;

    c = atof(params[0].c_str());
    a = atof(params[1].c_str());

    C_v = atof(params[2].c_str());

    CFL = atof(params[3].c_str());
    precision = atof(params[4].c_str());
    t_init = atof(params[5].c_str());
    t_simu = atof(params[6].c_str());

    // Verifications preliminaires
    if (c <= 0)
        throw string("ERREUR: Vérifiez que c > 0");
    if (a <= 0)
        throw string("ERREUR: Vérifiez que a > 0");
    if (C_v <= 0)
        throw string("ERREUR: Vérifiez que C_v > 0");
    if (CFL <= 0)
        throw string("ERREUR: Vérifiez que CFL > 0");
    if (precision <= 0)
        throw string("ERREUR: Vérifiez que la precision est > 0");
    if (t_init < 0)
        throw string("ERREUR: Vérifiez que t_init >= 0");
    if (t_simu <= 0)
        throw string("ERREUR: Vérifiez que t_simu > 0");

    dt = CFL * mesh->dx/c;
    step_count = (int)(t_simu / dt) + 1;
    time_steps = vector<double>(step_count);

    rho_expr = params[7];
    sigma_a_expr = params[8];
    sigma_c_expr = params[9];

    E = vector<double>(mesh->N+2);
    E_exact_expr = params[10];
    E_init_expr = params[11];
    E_left = vector<double>(step_count);
    E_right = vector<double>(step_count);

    F = vector<double>(mesh->N+2);
    F_exact_expr = params[12];
    F_init_expr = params[13];
    F_left = vector<double>(step_count);
    F_right = vector<double>(step_count);

    T = vector<double>(mesh->N+2);
    T_exact_expr = params[14];
    T_init_expr = params[15];
    T_left = vector<double>(step_count);
    T_right = vector<double>(step_count);
} 


double Solver::rho(double x){
    static int first_call = 1;
    static Parser p;

    p.DefineVar("x", &x); 
    if (first_call == 1){           // Un essai d'optimisation
        p.SetExpr(rho_expr);
        first_call = 0;
    }

    return p.Eval();
}


double Solver::sigma_a(double rho, double T){
    static int first_call = 1;
    static Parser p;

    p.DefineVar("rho", &rho); 
    p.DefineVar("T", &T); 
    if (first_call == 1){ 
        p.SetExpr(sigma_a_expr); 
        first_call = 0;
    }

    return p.Eval();
}


double Solver::sigma_c(double rho, double T){
    static int first_call = 1;
    static Parser p;

    p.DefineVar("rho", &rho); 
    p.DefineVar("T", &T); 
    if (first_call == 1){
        p.SetExpr(sigma_c_expr);
        first_call = 0;
    }

    return p.Eval();
}


double Solver::E_exact(double t, double x){
    static int first_call = 1;
    static Parser p;

    p.DefineVar("t", &t);
    p.DefineVar("x", &x);
    p.DefineVar("t_init", &t_init);
    if (first_call == 1){
        p.SetExpr(E_exact_expr);
        first_call = 0;
    }

    return p.Eval();
}


double Solver::E_init(double x){
    static int first_call = 1;
    static Parser p;

    p.DefineVar("x", &x);
    p.DefineVar("t_init", &t_init);
    if (first_call == 1){
        p.SetExpr(E_init_expr);
        first_call = 0;
    }

    return p.Eval();
}


double Solver::F_exact(double t, double x){
    static int first_call = 1;
    static Parser p;

    p.DefineVar("t", &t);
    p.DefineVar("x", &x);
    p.DefineVar("t_init", &t_init);
    if (first_call == 1){
        p.SetExpr(F_exact_expr);
        first_call = 0;
    }

    return p.Eval();
}


double Solver::F_init(double x){
    static int first_call = 1;
    static Parser p;

    p.DefineVar("x", &x); 
    p.DefineVar("t_init", &t_init);
    if (first_call == 1){
        p.SetExpr(F_init_expr);
        first_call = 0;
    }

    return p.Eval();
}


double Solver::T_exact(double t, double x){
    static int first_call = 1;
    static Parser p;

    p.DefineVar("t", &t);
    p.DefineVar("x", &x);
    p.DefineVar("t_init", &t_init);
    if (first_call == 1){
        p.SetExpr(T_exact_expr);
        first_call = 0;
    }

    return p.Eval();
}


double Solver::T_init(double x){
    static int first_call = 1;
    static Parser p;

    p.DefineVar("x", &x); 
    p.DefineVar("t_init", &t_init);
    if (first_call == 1){
        p.SetExpr(T_init_expr);
        first_call = 0;
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
    int N = mesh->N;            // Nombre de mailles interieures
    double dx = mesh->dx;       // Delta x

    // Les variables pour l'etape 1
    double E_n, E_next, T_n, F_n, F_next, Theta, Theta_n, Theta_next;
    
    // Les variables pour l'etape 2
    vector<double> E_etoile(N+2), F_etoile(N+2), T_etoile(N+2);
    vector<double> E_suiv(N+2), F_suiv(N+2), T_suiv(N+2);

    // Initialisation de la doucle de resolution
    for (int j = 1; j < N+1; j++){
        double x = mesh->cells[j][1];           // Centre de la maille
        E[j] = E_init(x);
        F[j] = F_init(x);
        T[j] = T_init(x);
    }

    // Temps courant (translate de t_init) et indice de l'iteration
    double t = 0;
    int n = 0;

    /**
     * Boucle de resolution
     */
    while (t <= t_simu){
        // Signaux aux bords du domaine pour ce pas d'iteration en vue de l'export
        E_left[n] = E[1];
        F_left[n] = F[1];
        T_left[n] = T[1];
        E_right[n] = E[N/2];                // ERREUR VOLONTAIRE!
        F_right[n] = F[N];
        T_right[n] = T[N];

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
                Theta_next = (gamma*Theta_n + alpha*delta*E_n) / (1 - beta*delta);

            } while (abs(E_next-E[j]) > precision && abs(Theta_next-Theta) > precision);
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
            double x_left = mesh->cells[j-1][1];        // Centre de la maille de gauche
            double x_center = mesh->cells[j][1];        // Centre de cette maille
            double x_right = mesh->cells[j+1][1];       // Centre de la maille de droite

            double sigma_c_left = sigma_c(rho(x_left), T[j-1]);
            double sigma_c_center = sigma_c(rho(x_center), T[j]);
            double sigma_c_right = sigma_c(rho(x_right), T[j+1]);

            double flux_sigma_c_left = flux_sigma_c(sigma_c_left, sigma_c_center);
            double flux_sigma_c_right = flux_sigma_c(sigma_c_center, sigma_c_right);
            double flux_M_left = flux_M(dx, flux_sigma_c_left);
            double flux_M_right = flux_M(dx, flux_sigma_c_right);

            double tmp = (1/dt) + (c/2)*(flux_M_right*flux_sigma_c_right + flux_M_left*flux_sigma_c_left);
            double Alpha = c*dt/dx;
            double Beta = 1/dt/tmp;
            double Gamma = c/dx/tmp;

            E_suiv[j] = E_etoile[j] - Alpha*(flux_F(flux_M_right, F[j], F[j+1], E[j], E[j+1]) - flux_F(flux_M_left, F[j-1], F[j], E[j-1], E[j]));

            F_suiv[j] = Beta*F_etoile[j] - Gamma*(flux_E(flux_M_right, E[j], E[j+1], F[j], F[j+1]) - flux_E(flux_M_left, E[j-1], E[j], F[j-1], F[j]));

            T_suiv[j] = T_etoile[j];
        }

        E = E_suiv;
        F = F_suiv;
        T = T_suiv;

        time_steps[n] = t;
        t += dt;
        n += 1;
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


void Solver::export_temporal(std::string file_name){
    ofstream file;
    file.open(file_name, ios_base::app);             // Ajout dans le fichier
    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_name + "'");

    file << mesh->x_min << "," << mesh->x_max << "," << mesh->N << "," << c << "," << a << "," << C_v << ","<< CFL << "," << precision << "," << t_init << "," << t_simu << ",\"" << rho_expr << "\",\"" << sigma_a_expr << "\",\"" << sigma_c_expr << "\",\"" << E_init_expr << "\",\"" << F_init_expr << "\",\"" << T_init_expr << "\"," << dt << "," << step_count << ",";

    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << time_steps[n];
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    //************************** Sol numerique
    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << E_left[n];
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << E_right[n];
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << F_left[n];
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << F_right[n];
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << T_left[n];
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << T_right[n];
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    //************************** Sol exacte
    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << E_exact(time_steps[n], mesh->cells[1][1]);
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << E_exact(time_steps[n], mesh->cells[mesh->N][1]);
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << F_exact(time_steps[n], mesh->cells[1][1]);
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << F_exact(time_steps[n], mesh->cells[mesh->N][1]);
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << T_exact(time_steps[n], mesh->cells[1][1]);
        if (n != step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < step_count; n++){
        file << T_exact(time_steps[n], mesh->cells[mesh->N][1]);
        if (n != step_count-1) file << ", ";
    }
    file << "]\"\n";

    file.close();
};


void Solver::export_spatial(std::string file_name){
    ofstream file;
    file.open(file_name, ios_base::app);             // Ajout dans le fichier
    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_name + "'");

    file << mesh->x_min << "," << mesh->x_max << "," << mesh->N << "," << c << "," << a << "," << C_v << ","<< CFL << "," << precision << "," << t_init << "," << t_simu << ",\"" << rho_expr << "\",\"" << sigma_a_expr << "\",\"" << sigma_c_expr << "\",\"" << E_init_expr << "\",\"" << F_init_expr << "\",\"" << T_init_expr << "\"," << dt << "," << step_count << ",";

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
     
    //************************** Sol initiale
    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << E_init(mesh->cells[j][1]);
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << F_init(mesh->cells[j][1]);
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << T_init(mesh->cells[j][1]);
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    //************************** Sol numerique finale
    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << E[j];
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
        file << T[j];
        if (j != mesh->N)
            file << ", ";
    }
    file << "]\",";

    //************************** Sol exacte finale
    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << E_exact(time_steps[step_count-1], mesh->cells[j][1]);
        if (j != mesh->N) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << F_exact(time_steps[step_count-1], mesh->cells[j][1]);
        if (j != mesh->N) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < mesh->N+1; j++){
        file << T_exact(time_steps[step_count-1], mesh->cells[j][1]);
        if (j != mesh->N) file << ", ";
    }
    file << "]\"\n";

    file.close();
};
