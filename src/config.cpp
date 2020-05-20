#include <iostream>
#include <fstream>

#include "config.hpp"

using namespace std;


Config::Config(string file_path){
    file_name = file_path;

    names["x_min"] = 0;
    names["x_max"] = 1;
    names["N"] = 2;
    names["c"] = 3;
    names["a"] = 4;
    names["C_v"] = 5;
    names["CFL"] = 6;
    names["precision"] = 7;
    names["t_0"] = 8;
    names["t_f"] = 9;
    names["rho"] = 10;
    names["sigma_a"] = 11;
    names["sigma_c"] = 12;
    names["E_exact"] = 13;
    names["E_0"] = 14;
    names["F_exact"] = 15;
    names["F_0"] = 16;
    names["T_exact"] = 17;
    names["T_0"] = 18;
    names["export_temporal"] = 19;
    names["export_spatial"] = 20;

    size = names.size();
    values = new string[size];
}


void Config::read(){
    map<string, int> :: iterator it;    // par cours des maps

    string name;                        // nom du parametre en cours de traitement

    bool read_unkown = false;           // lecture d'un parametre inconnue
    string unkown_name;                 // nom de l'inconnu

    map<string, int> count = names;    // compte les occurence de chaque parametre dans le fichier
    for (it = count.begin(); it != count.end(); it++) 
        count[it->first] = 0;

    ifstream file(file_name);
    if(file){
        while(!file.eof()){
            file >> name;
            it = names.find(name);
            if (it != names.end()){
                count[it->first] ++ ;
                file >> values[it->second];
                cout << "   -- " << name << " : " << values[it->second] << endl;
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

    for (it = count.begin(); it != count.end(); it++){
        if (count[it->first] < 1)
            throw string ("ERREUR: Parametre '"+it->first+"' manquant dans le fichier de configuration");
        if (count[it->first] > 1)
            throw string ("ERREUR: Parametre '"+it->first+"' duplique dans le fichier de configuration");
    }
}


Config::~Config(){
    delete[] values;
}