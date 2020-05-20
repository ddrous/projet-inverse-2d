#include <iostream>
#include <fstream>

#include "include/exporter.hpp"

using namespace std;

Exporter::Exporter(Solver *new_solver){
    if (new_solver == NULL)
        throw string("ERREUR: Pas de solveur a exporter");

    solver = new_solver;
}

void Exporter::case_1(std::string file_spatial, std::string file_temporal){
    // Export Spatial
    ofstream file;
    file.open(file_spatial, ios_base::trunc);

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_spatial + "'");

    Solver &s = *solver;        // Pou raccourcir

    file << "x, E_0, F_0, E+F, E-F\n";

    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.mesh->cells[j][1] << "," << s.E_0(s.mesh->cells[j][1]) << "," << s.F_0(s.mesh->cells[j][1]) << "," << s.E[j]+s.F[j] << "," << s.E[j]-s.F[j] << "\n" ;
    }

    file.close();

    // Export Temporelle
    file.open(file_temporal, ios_base::trunc);

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_spatial + "'");

    file << "t, E_left+F_left, E_left-F_left, E_right+F_right, E_right-F_right\n";

    for (int n = 0; n < s.step_count; n++){
        file << s.time_steps[n] << "," << s.E_left[n]+s.F_left[n] << "," << s.E_left[n]-s.F_left[n] << "," << s.E_right[n]+s.F_right[n] << "," << s.E_right[n]-s.F_right[n] << "\n" ;
    }

    file.close();
}

void Exporter::case_2(std::string file_spatial, std::string file_temporal){
    // Export Spatial
    ofstream file;
    file.open(file_spatial, ios_base::trunc);

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_spatial + "'");

    Solver &s = *solver;        // Pou raccourcir

    file << "x, E_0, E, E_exact\n";

    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.mesh->cells[j][1] << "," << s.E_0(s.mesh->cells[j][1]) << "," << s.E[j] << "," << s.E_exact(s.time_steps[s.step_count-1], s.mesh->cells[j][1]) << "\n" ;
    }

    file.close();

    // Export Temporelle
    file.open(file_temporal, ios_base::trunc);

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_spatial + "'");

    file << "t, E_left, E_center, E_right, E_exact_left, E_exact_center, E_exact_right\n";

    for (int n = 0; n < s.step_count; n++){
        file << s.time_steps[n] << "," << s.E_evol[0][n] << "," << s.E_evol[1][n] << "," << s.E_evol[2][n] << "," << s.E_exact(s.time_steps[n], s.mesh->cells[1][1]) << "," << s.E_exact(s.time_steps[n], s.mesh->cells[(s.mesh->N+1)/2][1]) << "," << s.E_exact(s.time_steps[n], s.mesh->cells[s.mesh->N][1]) <<  "\n" ;
    }

    file.close();

}

void Exporter::case_3(std::string file_spatial, std::string file_temporal){
    // Export Spatial
    ofstream file;
    file.open(file_spatial, ios_base::trunc);

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_spatial + "'");

    Solver &s = *solver;        // Pou raccourcir
    int n = s.step_count;

    file << "x, t=" << 0 << ", t=t_f/4, t=t_f/2, t=3*t_f/4, t=t_f" << "\n";

    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.mesh->cells[j][1] << "," << s.T_evol[0][j] << "," << s.T_evol[1][j] << "," << s.T_evol[2][j] << "," << s.T_evol[3][j] << "," << s.T_evol[4][j] << "\n" ;
    }

    file.close();
}


void Exporter::spatial(std::string file_name, std::string mode){
    ofstream file;
    if ((mode == "append") == 0)
        file.open(file_name, ios_base::app);
    else
        file.open(file_name, ios_base::trunc);

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_name + "'");

    Solver &s = *solver;        // Pou raccourcir

    //**************************** Parametres du probleme
    file << s.mesh->x_min << "," << s.mesh->x_max << "," << s.mesh->N << "," << s.c << "," << s.a << "," << s.C_v << ","<< s.CFL << "," << s.precision << "," << s.t_0 << "," << s.t_f << ",\"" << s.rho_expr << "\",\"" << s.sigma_a_expr << "\",\"" << s.sigma_c_expr << "\",\"" << s.E_0_expr << "\",\"" << s.F_0_expr << "\",\"" << s.T_0_expr << "\"," << s.dt << "," << s.step_count << ",";

    //************************** Abcisse 
    file << "\"[";
    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.mesh->cells[j][1];
        if (j != s.mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.rho(s.mesh->cells[j][1]);
        if (j != s.mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.sigma_a(s.rho(s.mesh->cells[j][1]), s.T[j]);
        if (j != s.mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.sigma_c(s.rho(s.mesh->cells[j][1]), s.T[j]);
        if (j != s.mesh->N)
            file << ", ";
    }
    file << "]\",";
     
    //************************** Sol initiale
    file << "\"[";
    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.E_0(s.mesh->cells[j][1]);
        if (j != s.mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.F_0(s.mesh->cells[j][1]);
        if (j != s.mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.T_0(s.mesh->cells[j][1]);
        if (j != s.mesh->N)
            file << ", ";
    }
    file << "]\",";

    //************************** Sol numerique finale
    file << "\"[";
    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.E[j];
        if (j != s.mesh->N)
            file << ", ";
    }
    file << "]\",";


    file << "\"[";
    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.F[j];
        if (j != s.mesh->N)
            file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int j = 1; j < s.mesh->N+1; j++){
        file << s.T[j];
        if (j != s.mesh->N)
            file << ", ";
    }
    file << "]\"\n";

    file.close();
}


void Exporter::temporal(std::string file_name, std::string mode){
    ofstream file;
    if ((mode == "append") == 0)
        file.open(file_name, ios_base::app);
    else
        file.open(file_name, ios_base::trunc);

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_name + "'");

    Solver &s = *solver;        // Pou raccourcir

    //**************************** Parametres du probleme
    file << s.mesh->x_min << "," << s.mesh->x_max << "," << s.mesh->N << "," << s.c << "," << s.a << "," << s.C_v << ","<< s.CFL << "," << s.precision << "," << s.t_0 << "," << s.t_f << ",\"" << s.rho_expr << "\",\"" << s.sigma_a_expr << "\",\"" << s.sigma_c_expr << "\",\"" << s.E_0_expr << "\",\"" << s.F_0_expr << "\",\"" << s.T_0_expr << "\"," << s.dt << "," << s.step_count << ",";

    //**************************Abcisse 
    file << "\"[";
    for (int n = 0; n < s.step_count; n++){
        file << s.time_steps[n];
        if (n != s.step_count-1) file << ", ";
    }
    file << "]\",";

    //************************** Sol numerique
    file << "\"[";
    for (int n = 0; n < s.step_count; n++){
        file << s.E_left[n];
        if (n != s.step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < s.step_count; n++){
        file << s.E_right[n];
        if (n != s.step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < s.step_count; n++){
        file << s.F_left[n];
        if (n != s.step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < s.step_count; n++){
        file << s.F_right[n];
        if (n != s.step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < s.step_count; n++){
        file << s.T_left[n];
        if (n != s.step_count-1) file << ", ";
    }
    file << "]\",";

    file << "\"[";
    for (int n = 0; n < s.step_count; n++){
        file << s.T_right[n];
        if (n != s.step_count-1) file << ", ";
    }
    file << "]\",";

    file.close();
}
