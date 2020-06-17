#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>

#include "include/solver.hpp"

#include "muParser.h"       // Pour transformer des expressions (chaines de caracteres) en fonctions

using namespace mu;
using namespace std;

Solver::Solver(Mesh *new_mesh, std::map<std::string, std::string> &params){
    mesh = new_mesh;

    c = atof(params["c"].c_str());
    a = atof(params["a"].c_str());

    C_v = atof(params["C_v"].c_str());

    CFL = atof(params["CFL"].c_str());
    precision = atof(params["precision"].c_str());
    t_0 = atof(params["t_0"].c_str());
    t_f = atof(params["t_f"].c_str());

    rho_expr = params["rho"];
    sigma_a_expr = params["sigma_a"];
    sigma_c_expr = params["sigma_c"];

    E = vector<double>(mesh->N+2);
    F = vector<double>(mesh->N+2);
    T = vector<double>(mesh->N+2);

    E_0_expr = params["E_0"];
    F_0_expr = params["F_0"];
    T_0_expr = params["T_0"];

    E_l_expr = params["E_l"];
    F_l_expr = params["F_l"];
    T_l_expr = params["T_l"];

    E_r_expr = params["E_r"];
    F_r_expr = params["F_r"];
    T_r_expr = params["T_r"];

    E_exact_expr = params["E_exact"];
    F_exact_expr = params["F_exact"];
    T_exact_expr = params["T_exact"];

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
    time_steps = vector<double>(step_count);

    // A exporter
    E_left = vector<double>(step_count);
    T_left = vector<double>(step_count);
    F_left = vector<double>(step_count);

    F_right = vector<double>(step_count);
    E_right = vector<double>(step_count);
    T_right = vector<double>(step_count);

    E_evol = new double*[3];
    for (int i = 0; i < 3; i++)
        E_evol[i] = new double[step_count];
    
    T_evol = new double*[5];
    for (int i = 0; i < 5; i++)
        T_evol[i] = new double[mesh->N];
} 

/**
 * Laplacian smoothing (k-means) d'un vecteur
 */ 
vector<double> smooth(vector<double>& input){
    int size = input.size();
    vector<double> output(size);
    
    output[0] = (input[0]+input[1]+input[2]) / 3;
    for (int j = 1; j < size-1; j++){
        output[j] = (input[j-1]+input[j]+input[j+1]) / 3;
    }
    output[size-1] = (input[size-3]+input[size-2]+input[size-1]) / 3;

    return output;
}

/**
 * Calcule rho sous forme de fonction crenaux
 */ 
double niche(double x, int nb, double y_min, double y_max, Mesh* mesh, int smooth_nb){
    static int first_call = 1;

    // Vecteur contenant le signal
    static vector<double> values(mesh->N+2);

    if (first_call==1){
        // Les attributs des crenaux pris au hazard
        double ** attr = new double*[nb];
        for (int i = 0; i < nb; i++)
            attr[i] = new double[3];

        srand(time(NULL));
        for (int k = 0; k < nb; k++){
            // attr[k][0] = rand() % mesh->N + 1;                     // position
            // attr[k][1] = rand() % mesh->N/20 + 1;                  // largeur
            // attr[k][2] = ((double) rand() / (RAND_MAX)) * (y_max-1) + 1;    // hauteur
            
            /* pour 1 creneau et N = 500 */
            attr[k][0] = (int)(0.7*mesh->N);                     // position
            attr[k][1] = (int)(0.05*mesh->N);                  // largeur
            attr[k][2] = y_max;    // hauteur
        }
        
        // Recherche des pics dans le vecteur d'entre
        for (int j = 0; j < mesh->N+2; j++){
            values[j] = y_min;
            for (int k = 0; k < nb; k++){
                if (abs(j - attr[k][0]) <= attr[k][1]){
                    values[j] = attr[k][2];
                    break;
                }
            }
        }

        // Suppression des attributs maintenatn inutiles
        for (int j = 0; j < nb; j++)
            delete[] attr[j];
        delete[] attr;

        // Lissage du signal
        for (int i = 0; i < smooth_nb; i++)
            values = smooth(values);

    first_call = 0;
    }

    int index = int((x - mesh->x_min) * mesh->N / (mesh->x_max - mesh->x_min));
    return values[index + 1];
}


double Solver::rho(double x){
    static int first_call = 1;
    static int rho_niche = rho_expr.compare("crenaux");

    if (rho_niche == 0)
        return niche(x, 1, 1, 10.0, mesh, (int)(0.1*mesh->N));
    else{
        static Parser p;
        p.DefineVar("x", &x);
        if (first_call == 1){ p.SetExpr(rho_expr); first_call = 0; }
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


double Solver::E_0(double x){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("x", &x);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(E_0_expr); first_call = 0; }

    return p.Eval();
}


double Solver::F_0(double x){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("x", &x); 
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(F_0_expr); first_call = 0; }

    return p.Eval();
}


double Solver::T_0(double x){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("x", &x); 
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(T_0_expr); first_call = 0; }

    return p.Eval();
}


double Solver::E_l(double t){
    static int first_call = 1;
    static int neumann = E_l_expr.compare("neumann");

    if (neumann == 0)
        return E[1];
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("x_min", &mesh->cells[1][0]);
        if (first_call == 1){ p.SetExpr(E_l_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::F_l(double t){
    static int first_call = 1;
    static int neumann = F_l_expr.compare("neumann");

    if (neumann == 0)
        return F[1];
    else{ 
        static Parser p;
        p.DefineVar("t", &t); 
        p.DefineVar("x_min", &mesh->cells[1][0]);
        if (first_call == 1){ p.SetExpr(F_l_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::T_l(double t){
    static int first_call = 1;
    static int neumann = T_l_expr.compare("neumann");

    if (neumann == 0)
        return T[1];
    else{ 
        static Parser p;
        p.DefineVar("t", &t); 
        p.DefineVar("x_min", &mesh->cells[1][0]);
        if (first_call == 1){ p.SetExpr(T_l_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::E_r(double t){
    static int first_call = 1;
    static int neumann = E_r_expr.compare("neumann");

    if (neumann == 0)
        return E[mesh->N];
    else{ 
        static Parser p;
        p.DefineVar("t", &t);
        p.DefineVar("x_max", &mesh->cells[mesh->N][2]);
        if (first_call == 1){ p.SetExpr(E_r_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::F_r(double t){
    static int first_call = 1;
    static int neumann = F_r_expr.compare("neumann");
    if (neumann == 0)
        return F[mesh->N];
    else{ 
        static Parser p;
        p.DefineVar("t", &t); 
        p.DefineVar("x_max", &mesh->cells[mesh->N][2]);
        if (first_call == 1){ p.SetExpr(F_r_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::T_r(double t){
    static int first_call = 1;
    static int neumann = T_r_expr.compare("neumann");

    if (neumann == 0)
        return T[mesh->N];
    else{ 
        static Parser p;
        p.DefineVar("t", &t); 
        p.DefineVar("x_max", &mesh->cells[mesh->N][2]);
        if (first_call == 1){ p.SetExpr(T_r_expr); first_call = 0; }
        return p.Eval();
    }
}


double Solver::E_exact(double t, double x){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("t", &t);
    p.DefineVar("x", &x);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(E_exact_expr); first_call = 0; }

    return p.Eval();
}


double Solver::F_exact(double t, double x){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("t", &t);
    p.DefineVar("x", &x);
    p.DefineVar("t_0", &t_0);
    if (first_call == 1){ p.SetExpr(F_exact_expr); first_call = 0; }

    return p.Eval();
}


double Solver::T_exact(double t, double x){
    static int first_call = 1;

    static Parser p;
    p.DefineVar("t", &t);
    p.DefineVar("x", &x);
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
        E[j] = E_0(x);
        F[j] = F_0(x);
        T[j] = T_0(x);
    }

    // Temps courant (translate de t_0) et indice de l'iteration
    double t = 0;
    int n = 0;

    /**
     * Boucle de resolution
     */
    while (t <= t_f){
        // Sauvegarde de l'animation
        this->save_animation(n);

        // Signaux aux bords du domaine pour ce pas d'iteration en vue de l'export
        E_left[n] = E[1];
        F_left[n] = F[1];
        T_left[n] = T[1];
        E_right[n] = E[N];
        F_right[n] = F[N];
        T_right[n] = T[N];

        // Evolution de l'energie en 3 points
        E_evol[0][n] = E[1];
        E_evol[1][n] = E[(N+1)/2];
        E_evol[2][n] = E[N];

        // Evolution de la temperature sur tout le domaine
        if (n == 1) 
            for (int j = 1; j < N+1; j++) T_evol[0][j] = T[j];
        else if (n == 1*step_count/4)
            for (int j = 1; j < N+1; j++) T_evol[1][j] = T[j];
        else if (n == 2*step_count/4)
            for (int j = 1; j < N+1; j++) T_evol[2][j] = T[j];
        else if (n == 3*step_count/4)
            for (int j = 1; j < N+1; j++) T_evol[3][j] = T[j];
        else if (n == step_count-1)
            for (int j = 1; j < N+1; j++) T_evol[4][j] = T[j];

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
                if (isnan(mu_q))
                    cerr << "ATTENTION! mu = " << mu_q << endl;
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

        // Remplissage des mailles fantomes
        E[0] = E_l(t);
        F[0] = F_l(t);
        T[0] = T_l(t);
        E[N+1] = E_r(t);
        F[N+1] = F_r(t);
        T[N+1] = T_r(t);

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


Solver::~ Solver(){
    for (int i = 0; i < 3; i++)
        delete[] E_evol[i];
    delete[] E_evol;

    for (int i = 0; i < 5; i++)
        delete[] T_evol[i];
    delete[] T_evol;
};