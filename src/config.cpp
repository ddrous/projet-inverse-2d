#include <iostream>
#include <fstream>

#include "config.hpp"

using namespace std;


Config::Config(string file_path){
    file_name = file_path;

    values["x_min"] = "";
    values["x_max"] = "";
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
    values["F_exact"] = "";
    values["T_exact"] = "";

    values["E_0"] = "";
    values["F_0"] = "";
    values["T_0"] = "";

    values["E_l"] = "";
    values["F_l"] = "";
    values["T_l"] = "";

    values["E_r"] = "";
    values["F_r"] = "";
    values["T_r"] = "";
    
    values["export_spatial"] = "";
    values["export_temporal"] = "";

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
        if (count[it->first] < 1){
            throw string ("ERREUR: Parametre '"+it->first+"' manquant dans le fichier de configuration");}
        if (count[it->first] > 1)
            throw string ("ERREUR: Parametre '"+it->first+"' duplique dans le fichier de configuration");
    }
}


Config::~Config(){}