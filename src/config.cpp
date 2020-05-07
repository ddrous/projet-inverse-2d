#include <iostream>
#include <fstream>

#include "config.hpp"

using namespace std;


Config::Config(string nom_du_fichier){
    nom_fichier = nom_du_fichier;
    valeurs = new double[9];
    fonctions = new string[14];
}


void Config::read(){
    // iterateur pour la "map" des variables
    map<string, int> :: iterator it;

    noms["xmin"] = 0;
    noms["xmax"] = 1;
    noms["N"] = 2;

    noms["c"] = 3;
    noms["a"] = 4;
    noms["C_v"] = 5;
    noms["CFL"] = 6;
    noms["epsilon"] = 7;
    noms["Temps_max"] = 8;

    noms["rho"] = 9;
    noms["sigma_a"] = 10;
    noms["sigma_c"] = 11;
    noms["E_x_0"] = 12;
    noms["E_0_t"] = 13;
    noms["E_N_t"] = 14;
    noms["F_x_0"] = 15;
    noms["F_0_t"] = 16;
    noms["F_N_t"] = 17;
    noms["T_x_0"] = 18;
    noms["T_0_t"] = 19;
    noms["T_N_t"] = 20;

    noms["export_1"] = 21;
    noms["export_2"] = 22;

    bool read_success = true;   // indique le succes de la lecture du fichier
    string nom;    // nom de la variable en cours de traitement

    ifstream file(nom_fichier);

    if(file){
        cout << nom << "Lecture des parametres en cours ... ";
        
        int i = 0;
        while(!file.eof()){
            file >> nom;
            it = noms.find(nom);
            if (it != noms.end()){
                if (it->second < 9){
                    file >> valeurs[it->second];
                    cout << "\n" << nom << " = " << valeurs[it->second] << "";
                }
                else{
                    file >> fonctions[it->second - 9];
                    cout << "\n" << nom << " = " << fonctions[it->second - 9] << "";   
                }
            }
            else
                read_success = false;
            i++;
        }
        
        file.close();
    
    }else
        throw string ("Erreur d'ouverture de fichier du fichier de configurations");

    if (read_success == true){
        cout << "\nLecture des 23 parametres OK!" << endl;
    } else
        throw string ("Fichier config mal ecrit");

}


Config::~Config(){
    delete[] valeurs;
    delete[] fonctions;
}