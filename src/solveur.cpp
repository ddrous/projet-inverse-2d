#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>

#include "solveur.hpp"

using namespace exprtk;
using namespace std;

Solver::Solver(Mesh *new_mesh, double *doubles, std::string *strings){
    mesh = new_mesh;

    c = doubles[0];
    a = doubles[1];

    C_v = doubles[2];

    CFL = doubles[3];
    epsilon = doubles[4];
    Temps_max = doubles[5];

    rho_expr = strings[0];
    sigma_a_expr = strings[1];
    sigma_c_expr = strings[2];

    E = vector<double>(mesh->N+2);
    E_x_0_expr = strings[3];
    E_0_t_expr = strings[4];
    E_N_t_expr = strings[5];

    F = vector<double>(mesh->N+2);
    F_x_0_expr = strings[6];
    F_0_t_expr = strings[7];
    F_N_t_expr = strings[8];

    T = vector<double>(mesh->N+2);
    T_x_0_expr = strings[9];
    T_0_t_expr = strings[10];
    T_N_t_expr = strings[11];

    export_1 = strings[12];
    export_2 = strings[13];
} 


double Solver::rho(double x){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){                             // seulement au 1er appel
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("x", x);
        expression.register_symbol_table(symbol_table);
        parser.compile(rho_expr,expression);
        first_call = 0;
    }

    return expression.value();
}


double Solver::sigma_a(double rho, double T){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("rho", rho);
        symbol_table.add_variable("T", T);
        expression.register_symbol_table(symbol_table);
        parser.compile(sigma_a_expr,expression);
        first_call = 0;
    }

    return expression.value();
}


double Solver::sigma_c(double rho, double T){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("rho", rho);
        symbol_table.add_variable("T", T);
        expression.register_symbol_table(symbol_table);
        parser.compile(sigma_c_expr,expression);
        first_call = 0;
    }

    return expression.value();
}


double Solver::E_x_0(double x){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("x", x);
        expression.register_symbol_table(symbol_table);
        parser.compile(E_x_0_expr,expression);
        first_call = 0;
    }

    return expression.value();
}


double Solver::E_0_t(double t){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("t", t);
        expression.register_symbol_table(symbol_table);
        parser.compile(E_0_t_expr,expression);
        first_call = 0;
    }

    // return expression.value();
    return E[1];
}


double Solver::E_N_t(double t){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("t", t);
        expression.register_symbol_table(symbol_table);
        parser.compile(E_N_t_expr,expression);
        first_call = 0;
    }

    // return expression.value();
    return E[mesh->N];
}


double Solver::F_x_0(double x){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("x", x);
        expression.register_symbol_table(symbol_table);
        parser.compile(F_x_0_expr,expression);
        first_call = 0;
    }

    return expression.value();
}


double Solver::F_0_t(double t){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("t", t);
        expression.register_symbol_table(symbol_table);
        parser.compile(F_0_t_expr,expression);
        first_call = 0;
    }

    // return expression.value();
    return F[1];
}


double Solver::F_N_t(double t){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("t", t);
        expression.register_symbol_table(symbol_table);
        parser.compile(F_N_t_expr,expression);
        first_call = 0;
    }

    // return expression.value();
    return F[mesh->N];
}


double Solver::T_x_0(double x){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("x", x);
        expression.register_symbol_table(symbol_table);
        parser.compile(T_x_0_expr,expression);
        first_call = 0;
    }

    return expression.value();
}


double Solver::T_0_t(double t){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("t", t);
        expression.register_symbol_table(symbol_table);
        parser.compile(T_0_t_expr,expression);
        first_call = 0;
    }

    // return expression.value();
    return T[1];
}


double Solver::T_N_t(double t){
    static int first_call = 1;
    static expression<double> expression;

    if (first_call == 1){
        symbol_table<double> symbol_table;
        parser<double> parser;
        symbol_table.add_variable("t", t);
        expression.register_symbol_table(symbol_table);
        parser.compile(T_N_t_expr,expression);
        first_call = 0;
    }

    // return expression.value();
    return T[mesh->N];
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
    // Les constantes
    int N = mesh->N;            // Nombre effectif de mailles
    double dx = mesh->dx;       // Delta x
    double dt = CFL*dx/c;       // Delta t

    // Les varaibles pour l'etape 1
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

    // Temps courant a l'interieur de la boucle
    double t = 0;

    // Pour l'export de la solution au temps t aux bords du domaine
    ofstream file_1;
    file_1.open(export_1, ios_base::app);          // Ajouter
    // file_1.open(export_1);                            // Ecraser
    if(!file_1)
        throw string ("Erreur d'ouverture du fichier " + export_1);
    // file_1 << "t, " << "E_0, " << "E_N, " << "F_0, " << "F_N, "  << "T_0, " << "T_N\n"; 
    export_current(t, file_1);

    /**
     * Boucle de resolution
     */
    while (t < Temps_max){
        t += dt;

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

        // Remplissage des mailles fantomes
        E[0] = E_0_t(t);
        F[0] = F_0_t(t);
        T[0] = T_0_t(t);
        E[N+1] = E_N_t(t);
        F[N+1] = F_N_t(t);
        T[N+1] = T_N_t(t);

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

        // Ajout des solution aux bords a cet instant dans le fichier export_1
        export_current(t, file_1);
    }

    // fermeture du fichier
    file_1 << "\n"; 
    file_1.close();
};


void Solver::dislay(){
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


void Solver::export_current(double t, std::ofstream &file){
    file << t << ", " << E[1] << ", " << E[mesh->N] << ", " << F[1] << ", " << F[mesh->N] << ", " << T[1] << ", " << T[mesh->N] << "\n";
};


void Solver::export_final(){
    ofstream file_2;
    file_2.open(export_2, ios_base::app);        // Pour l'ajout
    // ofstream file_2(export_2);                      // Pour ecraser

    if(file_2){
        // file_2 << "x" << ", " << "E" << ", " << "F"  << ", " << "T\n"; 
        for (int j = 1; j < mesh->N+1; j++)
            file_2 << mesh->cells[j][1] << ", " << E[j] << ", " << F[j] << ", " << T[j] << "\n"; 
        file_2 << "\n"; 
        file_2.close();
    }else {throw string ("Erreur d'ouverture de fichier " + export_2);}
};
