#include <iostream>
#include <fstream>

#include "config.hpp"

using namespace std;


Config::Config(string file_path){
    file_name = file_path;

    values["x_min"] = "";
    values["x_max"] = "";
    values["y_min"] = "";
    values["y_max"] = "";
    values["N"] = "";

    values["c"] = "";
    values["a"] = "";
    values["C_v"] = "";

    values["CFL"] = "";
    values["precision"] = "";
    values["t_0"] = "";
    values["t_f"] = "";

    values["rho"] = "";
    values["sigma_a"] = "";
    values["sigma_c"] = "";

    values["E_exact"] = "";
    values["F_exact_x"] = "";
    values["F_exact_y"] = "";
    values["T_exact"] = "";

    values["E_0"] = "";
    values["F_0_x"] = "";
    values["F_0_y"] = "";
    values["T_0"] = "";

    values["E_u"] = "";
    values["F_u_x"] = "";
    values["F_u_y"] = "";
    values["T_u"] = "";

    values["E_d"] = "";
    values["F_d_x"] = "";
    values["F_d_y"] = "";
    values["T_d"] = "";

    values["E_l"] = "";
    values["F_l_x"] = "";
    values["F_l_y"] = "";
    values["T_l"] = "";

    values["E_r"] = "";
    values["F_r_x"] = "";
    values["F_r_y"] = "";
    values["T_r"] = "";

    values["export_file"] = "";
    values["write_mode"] = "";

    size = values.size();
}


void Config::read(){
    map<string, string> :: iterator it;     // parcours des valeurs

    map<string, int> count;                 // compte les occurence de chaque parametre dans le fichier
    for (it = values.begin(); it != values.end(); ++it)
        count[it->first] = 0;

    string name;                        // nom du parametre en cours de traitement

    bool read_unkown = false;           // lecture d'un parametre inconnue
    string unkown_name;                 // nom de l'inconnu

    ifstream file(file_name);
    if(file){
        while(!file.eof()){
            file >> name;
            it = values.find(name);
            if (it != values.end()){
                count[it->first] ++ ;
                file >> values[it->first];
                cout << "   -- " << name << " : " << values[it->first] << endl;
            } 
            else{
                read_unkown = true;
                unkown_name = name;
                break;
            }
        }
        
        file.close();
    
    } else
        throw string ("ERREUR: Erreur d'ouverture du fichier de configuration '" + file_name + "'");

    if (read_unkown == true)
        throw string ("ERREUR: Parametre inconnu '" + unkown_name + "' dans le fichier de configuration");

    for (it = values.begin(); it != values.end(); ++it){
        if (count[it->first] < 1 && !(it->first == "E_exact") && !(it->first == "F_exact_x") && !(it->first == "F_exact_y") && !(it->first == "T_exact")){
            throw string ("ERREUR: Paramètre '"+it->first+"' manquant dans le fichier de configuration");}
        if (count[it->first] > 1)
            throw string ("ERREUR: Paramètre '"+it->first+"' dupliqué dans le fichier de configuration");
    }
}


Config::~Config(){}