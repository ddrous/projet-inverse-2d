#include <iostream>
#include <fstream>
#include <map>
#include <cmath>

#include "include/exporter.hpp"

using namespace std;


Exporter::Exporter(const Solver *new_solver){
    if (new_solver == NULL)
        throw string("ERREUR: Pas de solveur a exporter");

    solver = new_solver;
}


// void Exporter::case_1(std::string file_spatial, std::string file_temporal){
//     /* Export Spatial */
//     ofstream file;
//     file.open(file_spatial, ios_base::trunc);

//     if(!file)
//         throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_spatial + "'");

//     Solver s = *solver;        // Pou raccourcir les notations

//     file << "x,E_0,F_0,E+F,E-F\n";

//     for (int j = 1; j < s.mesh->N+1; j++){
//         file << s.mesh->cells[j][1] << "," << s.E_0(s.mesh->cells[j][1]) << "," << s.F_0(s.mesh->cells[j][1]) << "," << s.E[j]+s.F[j] << "," << s.E[j]-s.F[j] << "\n" ;
//     }

//     file.close();

//     /* Export Temporelle */
//     file.open(file_temporal, ios_base::trunc);

//     if(!file)
//         throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_spatial + "'");

//     file << "t,E_left+F_left,E_left-F_left,E_right+F_right,E_right-F_right\n";

//     for (int n = 0; n < s.step_count; n++){
//         file << s.time_steps[n] << "," << s.E_left[n]+s.F_left[n] << "," << s.E_left[n]-s.F_left[n] << "," << s.E_right[n]+s.F_right[n] << "," << s.E_right[n]-s.F_right[n] << "\n" ;
//     }

//     file.close();
// }


// void Exporter::case_2(std::string file_spatial, std::string file_temporal){
//     /* Export Spatial */
//     ofstream file;
//     file.open(file_spatial, ios_base::trunc);

//     if(!file)
//         throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_spatial + "'");

//     Solver s = *solver;

//     if (s.E_exact_expr.empty() == false){
//         file << "x,E_0, E,E_exact\n";

//         for (int j = 1; j < s.mesh->N+1; j++){
//             file << s.mesh->cells[j][1] << "," << s.E_0(s.mesh->cells[j][1]) << "," << s.E[j] << "," << s.E_exact(s.time_steps[s.step_count-1], s.mesh->cells[j][1]) << "\n" ;
//         }

//         file.close();

//         /* Export Temporelle */
//         file.open(file_temporal, ios_base::trunc);

//         if(!file)
//             throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_spatial + "'");

//         file << "t,E_left,E_right,E_exact_left,E_exact_right\n";

//         for (int n = 0; n < s.step_count; n++){
//             file << s.time_steps[n] << "," << s.E_left[n] << "," << s.E_right[n] << "," << s.E_exact(s.time_steps[n], s.mesh->cells[1][1]) << "," << s.E_exact(s.time_steps[n], s.mesh->cells[s.mesh->N][1]) <<  "\n";
//         }

//         file.close();
//     } 
//     else
//         throw string ("ERREUR: E_exact non fourni");

// }


// void Exporter::case_3(std::string file_spatial, std::string file_temporal){
//     /* Export Spatial */
//     ofstream file;
//     file.open(file_spatial, ios_base::trunc);

//     if(!file)
//         throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_spatial + "'");

//     Solver s = *solver;
//     int n = s.step_count;

//     file << "x,T_radiation,T_matter" << "\n";

//     for (int j = 1; j < s.mesh->N+1; j++){
//         file << s.mesh->cells[j][1] << "," << s.T[j] << "," << pow(s.E[j]/s.a, 0.25) << "\n" ;
//     }

//     file.close();
// }



/**
 * Ecris un signal spatial a n_rows lignes et n_cols colones dans le fichier file
 */
void write_spatial(ofstream &file, vector<double> const &signal, int n_rows, int n_cols){
    file << "\"[";
    for (int j = 0; j < n_cols; j++){
        file << "[";
        for (int i = 0; i < n_rows; i++){
            int k = cell_id(i+1, j+1, n_rows+2, n_cols+2);
            file << signal[k];
            if (i != n_rows-1) file << ", ";
        }
        file << "]";
        if (j != n_cols-1) file << ", ";
    }
    file << "]\"";
}


/**
 * Ecris un signal temporel a n_rows lignes et n_cols colones dans le fichier file
 */
void write_temporal(ofstream &file, double **signal, int n_rows, int n_cols){
    file << "\"[";
    for (int n = 0; n < n_rows; n++){
        file << "[";
        for (int i = 0; i < n_cols; i++){
            file << signal[n][i];
            if (i != n_cols-1) file << ", ";
        }
        file << "]";
        if (n != n_rows-1) file << ", ";
    }
    file << "]\"";
}


void Exporter::write_dataframe(std::string file_name, std::string mode){
    ofstream file;
    if (mode.compare("append") == 0)
        file.open(file_name, ios_base::app);
    else
        file.open(file_name, ios_base::trunc);

    if(!file)
        throw string ("ERREUR: Erreur d'ouverture du fichier '" + file_name + "'");


    //**************************** En mode trunc, il faut rajouter l'en tete
    if (mode.compare("append") != 0)
        file << "x_min,x_max,y_min,y_max,N,M,c,a,C_v,CFL,precision,t_0,t_f,step_count,rho_expr,sigma_a_expr,sigma_c_expr,E_0_expr,F_0_x_expr,F_0_y_expr,T_0_expr,E_u_expr,F_u_x_expr,F_u_y_expr,T_u_expr,E_d_expr,F_d_x_expr,F_d_y_expr,T_d_expr,E_l_expr,F_l_x_expr,F_l_y_expr,T_l_expr,E_r_expr,F_r_x_expr,F_r_y_expr,T_r_expr,rho_attr,rho,E_f,F_f_x,F_f_y,T_f,E_u,F_u,T_u,E_d,F_d,T_d,E_l,F_l,T_l,E_r,F_r,T_r" << "\n";

// ---------------------------------------------------------------------- Avec une reference
    // Solver s = *solver;

    // //**************************** Parametres du probleme
    // file << s.mesh->x_min << "," << s.mesh->x_max << "," << s.mesh->y_min << "," << s.mesh->y_max << "," << s.mesh->N << "," << s.mesh->M << "," << s.c << "," << s.a << "," << s.C_v << ","<< s.CFL << "," << s.precision << "," << s.t_0 << "," << s.t_f << "," << s.step_count << ",\"" << s.rho_expr << "\",\"" << s.sigma_a_expr << "\",\"" << s.sigma_c_expr << "\",\"" << s.E_0_expr << "\",\"" << s.F_0_x_expr << "\",\"" << s.F_0_y_expr << "\",\"" << s.T_0_expr << "\",\"" << s.E_u_expr << "\",\"" << s.F_u_x_expr << "\",\"" << s.F_u_y_expr << "\",\"" << s.T_u_expr << "\",\"" << s.E_d_expr << "\",\"" << s.F_d_x_expr << "\",\"" << s.F_d_y_expr << "\",\"" << s.T_d_expr << "\",\"" << s.E_l_expr << "\",\"" << s.F_l_x_expr << "\",\"" << s.F_l_y_expr << "\",\"" << s.T_l_expr << "\",\"" << s.E_r_expr << "\",\"" << s.F_r_x_expr << "\",\"" << s.F_r_y_expr << "\",\"" << s.T_r_expr << "\",";

    // //************************** Attributs du crenau de densite
    // file << "\"[";
    // for (int l = 0; l < s.n_niche; l++){
    //     if (l != 0) file << ", ";
    //     file << "(" << s.attr[l][0] << ", " << s.attr[l][1] << ", " << s.attr[l][2] << ", " << s.attr[l][3] << ")";
    // }
    // file << "]\",";

    // //************************** Densite
    // file << "\"[";
    // for (int j = 1; j <= s.mesh->M; j++){
    //     file << "[";
    //     for (int i = 1; i <= s.mesh->N; i++){
    //         file << s.rho(s.mesh->x[i], s.mesh->y[j]);
    //         if (i != s.mesh->N) file << ", ";
    //     }
    //     file << "]";
    //     if (j != s.mesh->M) file << ", ";
    // }
    // file << "]\",";

    // write_spatial(file, s.E, s.mesh->N, s.mesh->M);    file << ",";

    // vector<double> F_x(s.mesh->n_cells);
    // vector<double> F_y(s.mesh->n_cells);
    // for (int k = 0; k < s.mesh->n_cells; k++){
    //     F_x[k] = s.F[k][0];
    //     F_y[k] = s.F[k][1];
    // }    
    // write_spatial(file, F_x, s.mesh->N, s.mesh->M);    file << ",";
    // write_spatial(file, F_y, s.mesh->N, s.mesh->M);    file << ",";

    // write_spatial(file, s.T, s.mesh->N, s.mesh->M);    file << ",";

    // //************************** Signaux sur les bords
    // // file << "\"[";
    // // for (int n = 0; n < s.step_count; n++){
    // //     file << "[";
    // //     for (int i = 1; i < s.mesh->N+1; i++){
    // //         // int k = cell_id(i, s.mesh->M, s.mesh->N+2, s.mesh->M+2);
    // //         file << s.E_up[n][i-1];
    // //         if (i != s.mesh->N) file << ", ";
    // //     }
    // //     file << "]";
    // //     if (n != s.step_count-1) file << ", ";
    // // }
    // // file << "]\",";

    // write_temporal(file, s.E_up, s.step_count, s.mesh->N);    file << ",";
    // write_temporal(file, s.F_up, s.step_count, s.mesh->N);    file << ",";
    // write_temporal(file, s.T_up, s.step_count, s.mesh->N);    file << ",";

    // write_temporal(file, s.E_down, s.step_count, s.mesh->N);    file << ",";
    // write_temporal(file, s.F_down, s.step_count, s.mesh->N);    file << ",";
    // write_temporal(file, s.T_down, s.step_count, s.mesh->N);    file << ",";

    // write_temporal(file, s.E_left, s.step_count, s.mesh->M);    file << ",";
    // write_temporal(file, s.F_left, s.step_count, s.mesh->M);    file << ",";
    // write_temporal(file, s.T_left, s.step_count, s.mesh->M);    file << ",";

    // write_temporal(file, s.E_right, s.step_count, s.mesh->M);    file << ",";
    // write_temporal(file, s.F_right, s.step_count, s.mesh->M);    file << ",";
    // write_temporal(file, s.T_right, s.step_count, s.mesh->M);    file << "\n";

// ---------------------------------------------------------------------- Avec un pointeur

    const Solver *s = solver;           // Pour raccourcir

    //**************************** Parametres du probleme
    file << s->mesh->x_min << "," << s->mesh->x_max << "," << s->mesh->y_min << "," << s->mesh->y_max << "," << s->mesh->N << "," << s->mesh->M << "," << s->c << "," << s->a << "," << s->C_v << ","<< s->CFL << "," << s->precision << "," << s->t_0 << "," << s->t_f << "," << s->step_count << ",\"" << s->rho_expr << "\",\"" << s->sigma_a_expr << "\",\"" << s->sigma_c_expr << "\",\"" << s->E_0_expr << "\",\"" << s->F_0_x_expr << "\",\"" << s->F_0_y_expr << "\",\"" << s->T_0_expr << "\",\"" << s->E_u_expr << "\",\"" << s->F_u_x_expr << "\",\"" << s->F_u_y_expr << "\",\"" << s->T_u_expr << "\",\"" << s->E_d_expr << "\",\"" << s->F_d_x_expr << "\",\"" << s->F_d_y_expr << "\",\"" << s->T_d_expr << "\",\"" << s->E_l_expr << "\",\"" << s->F_l_x_expr << "\",\"" << s->F_l_y_expr << "\",\"" << s->T_l_expr << "\",\"" << s->E_r_expr << "\",\"" << s->F_r_x_expr << "\",\"" << s->F_r_y_expr << "\",\"" << s->T_r_expr << "\",";

    //************************** Attributs du crenau de densite
    file << "\"[";
    for (int l = 0; l < s->n_niche; l++){
        if (l != 0) file << ", ";
        file << "(" << s->attr[l][0] << ", " << s->attr[l][1] << ", " << s->attr[l][2] << ", " << s->attr[l][3] << ")";
    }
    file << "]\",";

    /* Donnes spatialles */
    write_spatial(file, s->rho_vec, s->mesh->N, s->mesh->M);    file << ",";

    write_spatial(file, s->E, s->mesh->N, s->mesh->M);    file << ",";

    vector<double> F_x(s->mesh->n_cells);
    vector<double> F_y(s->mesh->n_cells);
    for (int k = 0; k < s->mesh->n_cells; k++){
        F_x[k] = s->F[k][0];
        F_y[k] = s->F[k][1];
    }    
    write_spatial(file, F_x, s->mesh->N, s->mesh->M);    file << ",";
    write_spatial(file, F_y, s->mesh->N, s->mesh->M);    file << ",";

    write_spatial(file, s->T, s->mesh->N, s->mesh->M);    file << ",";

    /* Donnees temporelles */
    write_temporal(file, s->E_up, s->step_count, s->mesh->N);    file << ",";
    write_temporal(file, s->F_up, s->step_count, s->mesh->N);    file << ",";
    write_temporal(file, s->T_up, s->step_count, s->mesh->N);    file << ",";

    write_temporal(file, s->E_down, s->step_count, s->mesh->N);    file << ",";
    write_temporal(file, s->F_down, s->step_count, s->mesh->N);    file << ",";
    write_temporal(file, s->T_down, s->step_count, s->mesh->N);    file << ",";

    write_temporal(file, s->E_left, s->step_count, s->mesh->M);    file << ",";
    write_temporal(file, s->F_left, s->step_count, s->mesh->M);    file << ",";
    write_temporal(file, s->T_left, s->step_count, s->mesh->M);    file << ",";

    write_temporal(file, s->E_right, s->step_count, s->mesh->M);    file << ",";
    write_temporal(file, s->F_right, s->step_count, s->mesh->M);    file << ",";
    write_temporal(file, s->T_right, s->step_count, s->mesh->M);    file << "\n";

    file.close();
}
