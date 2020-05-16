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
    names["E_x_0"] = 12;
    names["E_0_t"] = 13;
    names["E_N_t"] = 14;
    names["F_x_0"] = 15;
    names["F_0_t"] = 16;
    names["F_N_t"] = 17;
    names["T_x_0"] = 18;
    names["T_0_t"] = 19;
    names["T_N_t"] = 20;
    names["export_1"] = 21;
    names["export_2"] = 22;

    size = names.size();
    doubles = new double[9];
    strings = new string[14];
}


void Config::read(){
    map<string, int> :: iterator it;

    string name;                        // nom de la variable en cours de traitement

    bool read_unkown = false;           // lecture d'une varaible inconnue
    string unkown_name;                 // nom de l'inconnu

    int *read_count = new int[size];    // compte les occurence de chaque variable dans le fichier
    for (int i = 0; i < size; i++) 
        read_count[i] = 0;
    bool not_found = false;             // indique s'il manque une ou plusieurs variables
    bool duplicate = false;             // indique si une variables est definie plus d'une fois

    ifstream file(file_name);
    if(file){
        while(!file.eof()){
            file >> name;
            it = names.find(name);
            if (it != names.end()){

                read_count[it->second] ++ ;

                if (it->second < 9){
                    file >> doubles[it->second];
                    cout << "    " << name << " = " << doubles[it->second] << endl;
                }
                else{
                    file >> strings[it->second - 9];
                    cout << "    " << name << " = " << strings[it->second - 9] << endl;
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
        throw string ("ERREUR: Erreur d'ouverture du fichier de configuration");

    if (read_unkown == true)
        throw string ("ERREUR: Variable inconnue '" + unkown_name + "' dans le fichier de configuration");

    for (int i = 0; i < size; i++){
        if (read_count[i] < 1)
            not_found = true;
        if (read_count[i] > 1)
            duplicate = true;
    }
    if (not_found == true)
        throw string ("ERREUR: Variable(s) manquante(s) dans le fichier de configuration");
    if (duplicate == true)
        throw string ("ERREUR: Variable(s) dupliquee(s) dans le fichier de configuration");

    delete[] read_count;
}


Config::~Config(){
    delete[] doubles;
    delete[] strings;
}