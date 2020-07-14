#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>

#include "include/solver.hpp"

#include "muParser.h"       // Pour transformer des expressions (chaines de caracteres) en fonctions

using namespace mu;
using namespace std;

/* Type utilisé pour le stockage des signaux */
typedef std::vector<double> vector_t;


Solver::Solver(const Mesh *new_mesh, const Config &cfg){
    mesh = new_mesh;

    c = atof(cfg.values.at("c").c_str());
    a = atof(cfg.values.at("a").c_str());

    C_v = atof(cfg.values.at("C_v").c_str());

    CFL = atof(cfg.values.at("CFL").c_str());
    precision = atof(cfg.values.at("precision").c_str());
    t_0 = atof(cfg.values.at("t_0").c_str());
    t_f = atof(cfg.values.at("t_f").c_str());

    rho_expr = cfg.values.at("rho");
    sigma_a_expr = cfg.values.at("sigma_a");
    sigma_c_expr = cfg.values.at("sigma_c");

    E = vector_t(mesh->n_cells);
    F = vector_t(mesh->n_cells);
    T = vector_t(mesh->n_cells);

    E_0_expr = cfg.values.at("E_0");
    F_0_expr = cfg.values.at("F_0");
    T_0_expr = cfg.values.at("T_0");

    E_u_expr = cfg.values.at("E_u");
    F_u_expr = cfg.values.at("F_u");
    T_u_expr = cfg.values.at("T_u");

    E_d_expr = cfg.values.at("E_d");
    F_d_expr = cfg.values.at("F_d");
    T_d_expr = cfg.values.at("T_d");

    E_l_expr = cfg.values.at("E_l");
    F_l_expr = cfg.values.at("F_l");
    T_l_expr = cfg.values.at("T_l");

    E_r_expr = cfg.values.at("E_r");
    F_r_expr = cfg.values.at("F_r");
    T_r_expr = cfg.values.at("T_r");

    E_exact_expr = cfg.values.at("E_exact");
    F_exact_expr = cfg.values.at("F_exact");
    T_exact_expr = cfg.values.at("T_exact");

    /* Verifications preliminaires */
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
    if (t_0 < 0)
        throw string("ERREUR: Vérifiez que t_0 >= 0");
    if (t_f <= 0)
        throw string("ERREUR: Vérifiez que t_f > 0");

    dt = CFL * mesh->dx/c;
    double tmp = t_f / dt;
    if (tmp == floor(tmp))        // si entier
        step_count = floor(tmp);
    else
        step_count = floor(tmp) + 1;
    time_steps = vector_t(step_count);

    /* A exporter */
    allocate(E_up, step_count, mesh->N);
    allocate(F_up, step_count, mesh->N);
    allocate(T_up, step_count, mesh->N);

    allocate(E_down, step_count, mesh->N);
    allocate(F_down, step_count, mesh->N);
    allocate(T_down, step_count, mesh->N);

    allocate(E_left, step_count, mesh->M);
    allocate(F_left, step_count, mesh->M);
    allocate(T_left, step_count, mesh->M);

    allocate(E_right, step_count, mesh->M);
    allocate(F_right, step_count, mesh->M);
    allocate(T_right, step_count, mesh->M);
} 

/**
 * Allocate space for 2D matrix
 */ 
void allocate(double **matrix, int nrows, int ncols){
    matrix = new double*[nrows];
    for (int i = 0; i < nrows; i++)
        matrix[i] = new double[ncols]; 
}

/**
 * Free space from 2D matrix
 */ 
void free(double **matrix, int nrows){
    for (int j = 0; j < nrows; j++)
        delete[] matrix[j];
    delete[] matrix;
}


/**
 * Laplacian smoothing (k-means) d'un vecteur, avec k=3
 */ 
vector_t smooth(vector_t& input){
    int size = input.size();
    vector_t output(size);
    
    output[0] = (input[0]+input[1]+input[2]) / 3;
    for (int j = 1; j < size-1; j++){
        output[j] = (input[j-1]+input[j]+input[j+1]) / 3;
    }
    output[size-1] = (input[size-3]+input[size-2]+input[size-1]) / 3;

    return output;
}


vector_t Solver::niche(int n_niche, double z_min, double z_max, int n_smooth){
    /* Vecteur qui va contenir le signal en crenaux */
    vector_t signal((mesh->N+2)*(mesh->M+2));

    /* Les attributs du signal */
    double ** attr = new double*[n_niche];
    for (int i = 0; i < n_niche; i++)
        attr[i] = new double[4];

    srand(time(NULL));
    for (int l = 0; l < n_niche; l++){
        /* Attributs des crenaux pris au hazard */
        // attr[k][0] = rand() % (N-1) + 1;                     // position
        // attr[k][1] = rand() % (N/3) + 5;                  // largeur
        // attr[k][2] = ((double) rand() / (RAND_MAX)) * (y_max-1.) + 1;    // hauteur
        
        /* Les memes attributs a chaque fois */
        attr[l][0] = (0.7*mesh->N);                     // abcisse
        attr[l][0] = (0.4*mesh->M);                     // ordonnee
        attr[l][1] = (0.1*mesh->n_cells);                    // diametre
        attr[l][2] = z_max;                                 // hauteur

    }

    /* Placement des crenaux */
    for (int k = 0; k < mesh->N+2; k++){
        signal[k] = z_min;
        int i = mesh->coord[k][0];
        int j = mesh->coord[k][1];
        for (int l = 0; l < n_niche; l++){
            if (sqrt(pow(i - attr[l][0], 2) + pow(j - attr[l][1], 2)) <= attr[l][2]/2){
                signal[j] = attr[l][3];
                break;
            }
        }
    }

    /* Suppression des attributs devenus inutiles */
    for (int j = 0; j < n_niche; j++)
        delete[] attr[j];
    delete[] attr;

    /* Lissage du signal */
    for (int i = 0; i < n_smooth; i++)
        signal = smooth(signal);

    return signal;
}


double Solver::rho(double x, double y){
    static int first_call = 1;
    static int rho_niche = rho_expr.compare("custom");

    if (rho_niche == 0){
        static vector_t signal(mesh->N+2);
        if (first_call == 1){signal = niche(1, 0.01, 1.0, (int)(0.05*mesh->n_cells)); first_call = 0;}
        int i = int((x - mesh->x_min + mesh->dx/2.) / mesh->dx);
        int j = int((y - mesh->y_min + mesh->dy/2.) / mesh->dy);
        int k = i + j*(mesh->M+2);
        return signal[j];
    }
    else{
        static Parser p;
        p.DefineVar("x", &x);
        p.DefineVar("y", &y);
        if (first_call == 1){p.SetExpr(rho_expr); first_call = 0;}
        return p.Eval();
    }
}


double Solver::sigma_a(double rho, double T){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("rho", &rho); 
    p.DefineVar("T", &T); 
    if (first_call == 1){ p.SetExpr(sigma_a_expr); first_call = 0; }

    return p.Eval();
}


double Solver::sigma_c(double rho, double T){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("rho", &rho);
    p.DefineVar("T", &T); 
    if (first_call == 1){ p.SetExpr(sigma_c_expr); first_call = 0; }

    return p.Eval();
}


double Solver::E_0(double x, double y){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("x", &x);
    p.DefineVar("y", &y);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(E_0_expr); first_call = 0; }

    return p.Eval();
}


double Solver::F_0(double x, double y){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("x", &x); 
    p.DefineVar("y", &y);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(F_0_expr); first_call = 0; }

    return p.Eval();
}


double Solver::T_0(double x, double y){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("x", &x); 
    p.DefineVar("y", &y);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(T_0_expr); first_call = 0; }

    return p.Eval();
}


double Solver::E_u(double t, double x){
    static int first_call = 1;
    static int neumann = E_u_expr.compare("neumann");

    if (neumann == 0){
        int i = int((x - mesh->x_min + mesh->dx/2.) / mesh->dx);
        int k = i + (mesh->M)*(mesh->M+2);
        return E[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("x", &x);
        if (first_call == 1){ p.SetExpr(E_u_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::F_u(double t, double x){
    static int first_call = 1;
    static int neumann = F_u_expr.compare("neumann");

    if (neumann == 0){
        int i = int((x - mesh->x_min + mesh->dx/2.) / mesh->dx);
        int k = i + (mesh->M)*(mesh->M+2);
        return F[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("x", &x);
        if (first_call == 1){ p.SetExpr(F_u_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::T_u(double t, double x){
    static int first_call = 1;
    static int neumann = T_u_expr.compare("neumann");

    if (neumann == 0){
        int i = int((x - mesh->x_min + mesh->dx/2.) / mesh->dx);
        int k = i + (mesh->M)*(mesh->M+2);
        return T[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("x", &x);
        if (first_call == 1){ p.SetExpr(T_u_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::E_d(double t, double x){
    static int first_call = 1;
    static int neumann = E_d_expr.compare("neumann");

    if (neumann == 0){
        int i = int((x - mesh->x_min + mesh->dx/2.) / mesh->dx);
        int k = i + 0*(mesh->M+2);
        return E[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("x", &x);
        if (first_call == 1){ p.SetExpr(E_d_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::F_d(double t, double x){
    static int first_call = 1;
    static int neumann = F_d_expr.compare("neumann");

    if (neumann == 0){
        int i = int((x - mesh->x_min + mesh->dx/2.) / mesh->dx);
        int k = i + 0*(mesh->M+2);
        return F[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("x", &x);
        if (first_call == 1){ p.SetExpr(F_d_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::T_d(double t, double x){
    static int first_call = 1;
    static int neumann = T_d_expr.compare("neumann");

    if (neumann == 0){
        int i = int((x - mesh->x_min + mesh->dx/2.) / mesh->dx);
        int k = i + 0*(mesh->M+2);
        return T[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("x", &x);
        if (first_call == 1){ p.SetExpr(T_d_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::E_l(double t, double y){
    static int first_call = 1;
    static int neumann = E_l_expr.compare("neumann");

    if (neumann == 0){
        int j = int((y - mesh->y_min + mesh->dy/2.) / mesh->dy);
        int k = 0 + j*(mesh->M+2);
        return E[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("y", &y);
        if (first_call == 1){ p.SetExpr(E_l_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::F_l(double t, double y){
    static int first_call = 1;
    static int neumann = F_l_expr.compare("neumann");

    if (neumann == 0){
        int j = int((y - mesh->y_min + mesh->dy/2.) / mesh->dy);
        int k = 0 + j*(mesh->M+2);
        return F[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("y", &y);
        if (first_call == 1){ p.SetExpr(F_l_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::T_l(double t, double y){
    static int first_call = 1;
    static int neumann = T_l_expr.compare("neumann");

    if (neumann == 0){
        int j = int((y - mesh->y_min + mesh->dy/2.) / mesh->dy);
        int k = 0 + j*(mesh->M+2);
        return T[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("y", &y);
        if (first_call == 1){ p.SetExpr(T_l_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::E_r(double t, double y){
    static int first_call = 1;
    static int neumann = E_r_expr.compare("neumann");

    if (neumann == 0){
        int j = int((y - mesh->y_min + mesh->dy/2.) / mesh->dy);
        int k = mesh->N + j*(mesh->M+2);
        return E[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("y", &y);
        if (first_call == 1){ p.SetExpr(E_r_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::F_r(double t, double y){
    static int first_call = 1;
    static int neumann = F_r_expr.compare("neumann");

    if (neumann == 0){
        int j = int((y - mesh->y_min + mesh->dy/2.) / mesh->dy);
        int k = mesh->N + j*(mesh->M+2);
        return F[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("y", &y);
        if (first_call == 1){ p.SetExpr(F_r_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::T_r(double t, double y){
    static int first_call = 1;
    static int neumann = T_r_expr.compare("neumann");

    if (neumann == 0){
        int j = int((y - mesh->y_min + mesh->dy/2.) / mesh->dy);
        int k = mesh->N + j*(mesh->M+2);
        return T[k];
    }
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("y", &y);
        if (first_call == 1){ p.SetExpr(T_r_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::E_exact(double t, double x, double y){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("t", &t);
    p.DefineVar("x", &x);
    p.DefineVar("y", &y);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(E_exact_expr); first_call = 0; }

    return p.Eval();
}


double Solver::F_exact(double t, double x, double y){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("t", &t);
    p.DefineVar("x", &x);
    p.DefineVar("y", &y);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(F_exact_expr); first_call = 0; }

    return p.Eval();
}


double Solver::T_exact(double t, double x, double y){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("t", &t);
    p.DefineVar("x", &x);
    p.DefineVar("y", &y);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(T_exact_expr); first_call = 0; }

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
double flux_M(double dx, double flux_sigma_c){
    return 2 / (2 + dx * flux_sigma_c);
}


/**
 * Calcule les flux E_{j+1/2} et E_{j-1/2}
 */ 
double flux_E(double flux_M, double E_left, double E_right, double F_left, double F_right){
    return flux_M * ((E_right + E_left) / 2 - (F_right - F_left) / 2);
}


/**
 * Calcule les flux F_{j+1/2} et F_{j-1/2}
 */ 
double flux_F(double flux_M, double F_left, double F_right, double E_left, double E_right){
    return flux_M * ((F_right + F_left) / 2 - (E_right - E_left) / 2);
}


void Solver::save_animation(int time_step){
    string file_name = "data/anim/animation." + to_string(time_step) + ".csv";
    ofstream file;
    file.open(file_name, ios_base::trunc);

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_name + "'");

    file << "x,rho,E,F,T,Tr\n";

    for (int j = 1; j < mesh->N+1; j++){
        file << mesh->cells[j][1] << "," << rho(mesh->cells[j][1]) << "," << E[j] << "," << F[j] << "," << T[j] << "," << pow(E[j]/a, 0.25) << "\n";
    }

    file.close();
}


void Solver::phase_1(){
    /* Des variables necessaires pour cette etape */
    double Theta;       // Theta = a*T^4 
    double E_n, F_n, T_n, Theta_n; 
    double E_next, F_next, Theta_next;
    
    for (int j = 1; j < mesh->N+1; j++){
        // Initialisation etape 1 
        E_n = E[j];
        F_n = F[j];
        T_n = T[j];
        Theta_n = a * pow(T_n, 4);
        Theta = Theta_n;

        E_next = E[j];
        F_next = F[j];
        Theta_next = Theta;

        do{
            E[j] = E_next;
            F[j] = F_next;
            Theta = Theta_next;
            
            T[j] = pow(Theta/a, 0.25);
            double mu_q = 1/ (pow(T_n, 3) + T_n*pow(T[j], 2) + T[j]*pow(T_n, 2) + pow(T[j], 3));
            if (isnan(mu_q))
                cerr << "ATTENTION! mu = nan" << endl;

            double rho_tmp = rho(mesh->cells[j][1]);
            double sigma_a_tmp = sigma_a(rho_tmp, T[j]);

            double tmp_1 = (1/dt) + c*sigma_a_tmp;
            double alpha = 1/dt/tmp_1;
            double beta = c*sigma_a_tmp/tmp_1;

            double tmp_2 = (rho_tmp*C_v*mu_q/dt) + c*sigma_a_tmp;
            double gamma = rho_tmp*C_v*mu_q/dt/tmp_2;
            double delta = c*sigma_a_tmp/tmp_2;

            E_next = (alpha*E_n + gamma*beta*Theta_n) / (1 - beta*delta);
            F_next = F_n;
            Theta_next = (gamma*Theta_n + alpha*delta*E_n) / (1 - beta*delta);

        } while (abs(E_next-E[j]) > precision && abs(Theta_next-Theta) > precision);
    }
};


void Solver::phase_1_equilibre(){
    /* Des variables necessaires pour cette etape */
    double Theta;       // Theta = a*T^4
    double E_n, F_n, T_n, Theta_n; 
    double E_next, F_next, Theta_next;
    
    for (int j = 1; j < mesh->N+1; j++){
        // Initialisation
        E_n = E[j];
        F_n = F[j];
        T_n = T[j];
        Theta_n = a * pow(T_n, 4);
        Theta = Theta_n;

        E_next = E[j];
        F_next = F[j];
        Theta_next = Theta;
            // cout << "E = "<< E_next << " Theta = " << Theta_next<< endl;

        do{
            E[j] = E_next;
            F[j] = F_next;
            Theta = Theta_next;

            T[j] = pow(Theta/a, 0.25);
            double mu_q = 1/ (pow(T_n, 3) + T_n*pow(T[j], 2) + T[j]*pow(T_n, 2) + pow(T[j], 3));
            if (isnan(mu_q))
                cerr << "ATTENTION! mu = nan" << endl;

            double rho_tmp = rho(mesh->cells[j][1]);
            double alpha = c * sigma_a(rho_tmp, T[j]) * dt;
            double beta = rho_tmp * C_v * mu_q;

            double X_n = E_n - Theta_n;
            double Y_n = E_n + beta*Theta_n;

            double X = X_n / (1 + alpha*(1 + (1./beta)));
            double Y = Y_n;

            E_next = (beta*X + Y) / (1 + beta);
            Theta_next = (-X + Y) / (1 + beta);

        } while (abs(E_next-E[j]) > precision && abs(Theta_next-Theta) > precision);
    }
};


void Solver::phase_2(){
    /* Vecteurs necessaires pour cette etape */
    vector_t E_etoile(mesh->N+2), F_etoile(mesh->N+2), T_etoile(mesh->N+2);
    vector_t E_suiv(mesh->N+2), F_suiv(mesh->N+2), T_suiv(mesh->N+2);

    /* Initialisation de l'etape */
    E_etoile = E;
    F_etoile = F;
    T_etoile = T;

    for (int j = 1; j < mesh->N+1; j++){
        double x_left = mesh->cells[j-1][1];        // Centre de la maille de gauche
        double x_center = mesh->cells[j][1];        // Centre de cette maille
        double x_right = mesh->cells[j+1][1];       // Centre de la maille de droite

        double sigma_c_left = sigma_c(rho(x_left), T[j-1]);
        double sigma_c_center = sigma_c(rho(x_center), T[j]);
        double sigma_c_right = sigma_c(rho(x_right), T[j+1]);

        double flux_sigma_c_left = flux_sigma_c(sigma_c_left, sigma_c_center);
        double flux_sigma_c_right = flux_sigma_c(sigma_c_center, sigma_c_right);
        double flux_M_left = flux_M(mesh->dx, flux_sigma_c_left);
        double flux_M_right = flux_M(mesh->dx, flux_sigma_c_right);

        double tmp = (1/dt) + (c/2)*(flux_M_right*flux_sigma_c_right + flux_M_left*flux_sigma_c_left);
        double alpha = -c*dt/mesh->dx;
        double beta = 1/dt / tmp;
        double gamma = c*(flux_M_right-flux_M_left)/mesh->dx / tmp;
        double delta = -c/mesh->dx / tmp;

        E_suiv[j] = E_etoile[j] + alpha*(flux_F(flux_M_right, F[j], F[j+1], E[j], E[j+1]) - flux_F(flux_M_left, F[j-1], F[j], E[j-1], E[j]));

        F_suiv[j] = beta*F_etoile[j] + gamma*E[j] + delta*(flux_E(flux_M_right, E[j], E[j+1], F[j], F[j+1]) - flux_E(flux_M_left, E[j-1], E[j], F[j-1], F[j]));

        T_suiv[j] = T_etoile[j];
    }

    E = E_suiv;
    F = F_suiv;
    T = T_suiv;
};


void Solver::solve(){
    /* Initialisation de la doucle de resolution */
    for (int j = 1; j < mesh->N+1; j++){
        double x = mesh->cells[j][1];           // Centre de la maille j
        E[j] = E_0(x);
        F[j] = F_0(x);
        T[j] = T_0(x);
    }

    /* Temps courant (translaté de t_0) et indice de l'iteration */
    double t = 0;
    int n = 0;

    /**
     * Boucle de resolution
     */
    while (t <= t_f){
        /* Enregistrement des signaux pour ce temps */
        // save_animation(n);

        /* Signaux exportés */
        E_left[n] = E[1];
        F_left[n] = F[1];
        T_left[n] = T[1];
        E_right[n] = E[mesh->N];
        F_right[n] = F[mesh->N];
        T_right[n] = T[mesh->N];

        /* etape 1 */ 
        // phase_1();
        phase_1_equilibre();

        /* Remplissage des mailles fantomes */
        E[0] = E_l(t);
        F[0] = F_l(t);
        T[0] = T_l(t);
        E[mesh->N+1] = E_r(t);
        F[mesh->N+1] = F_r(t);
        T[mesh->N+1] = T_r(t);

        /* etape 2 */
        phase_2();

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

Solver::~Solver(){
    free(E_up, step_count);
    free(F_up, step_count);
    free(T_up, step_count);

    free(E_down, step_count);
    free(F_down, step_count);
    free(T_down, step_count);

    free(E_left, step_count);
    free(F_left, step_count);
    free(T_left, step_count);

    free(E_right, step_count);
    free(F_right, step_count);
    free(T_right, step_count);
}