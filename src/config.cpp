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
    names["epsilon"] = 7;
    names["t_final"] = 8;
    names["rho"] = 9;
    names["sigma_a"] = 10;
    names["sigma_c"] = 11;
    names["E_init"] = 12;
    names["F_init"] = 13;
    names["T_init"] = 14;
    names["export_temporal"] = 15;
    names["export_spatial"] = 16;

    size = names.size();
    doubles = new double[9];
    strings = new string[8];
}


void Config::read(){
    map<string, int> :: iterator it;

    string name;                        // nom du parametre en cours de traitement

    bool read_unkown = false;           // lecture d'un parametre inconnue
    string unkown_name;                 // nom de l'inconnu

    int *read_count = new int[size];    // compte les occurence de chaque parametre dans le fichier
    for (int i = 0; i < size; i++) 
        read_count[i] = 0;
    bool not_found = false;             // indique s'il manque un ou plusieurs parametre
    bool duplicate = false;             // indique si un parametre est definie plus d'une fois

    ifstream file(file_name);
    if(file){
        while(!file.eof()){
            file >> name;
            it = names.find(name);
            if (it != names.end()){

                read_count[it->second] ++ ;

                if (it->second < 9){
                    file >> doubles[it->second];
                    cout << "   -- " << name << " : " << doubles[it->second] << endl;
                }
                else{
                    file >> strings[it->second - 9];
                    cout << "   -- " << name << " : " << strings[it->second - 9] << endl;
                }
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

    for (int i = 0; i < size; i++){
        if (read_count[i] < 1)
            not_found = true;
        if (read_count[i] > 1)
            duplicate = true;
    }
    if (not_found == true)
        throw string ("ERREUR: Parametre(s) manquant(s) dans le fichier de configuration");
    if (duplicate == true)
        throw string ("ERREUR: Parametre(s) duplique(s) dans le fichier de configuration");

    delete[] read_count;
}


Config::~Config(){
    delete[] doubles;
    delete[] strings;
}