#include <iostream>
#include <fstream>

#include "mesh.hpp"
#include "solver.hpp"
#include "config.hpp"

#include "muParser.h"

using namespace std;
using namespace mu;

int main(int argc, char * argv[]){
    try {
        cout << "\n======== Equation du transfer radiatif ========" << endl;

        if(argc < 2 || argc > 2)
            throw string("ERREUR: Fournissez un (et un seul) fichier de configuration");

        // Lecture du fichier config
        Config cfg = Config(argv[1]);
        cout << "\nLecture des " << cfg.size << " parametres en cours ... " << endl;
        cfg.read();
        cout << "Lecture OK !" << endl;

        // Creation d'un maillage uniforme
        double x_min = cfg.doubles[0];
        double x_max = cfg.doubles[1];
        int N = (int)cfg.doubles[2];
        Mesh m = Mesh(x_min, x_max, N);
        cout << "\nCreation des " << m.N+2 << " mailles en cours ... " << endl;
        m.create_cells();
        cout << "Creation OK !" << endl;
    
        // Resolution du probleme
        Solver s = Solver(&m, &cfg.doubles[3], cfg.strings);
        cout << "\nResolution (" << s.step_count << " iterations) en cours ..." << endl;
        s.solve();
        cout << "Resolution OK !" << endl;

        cout << "\nExportation des resultats en cours ..." << endl;
        s.export_temporal();
        s.export_spatial();
        cout << "Export OK !" << endl;

        // cout << "\nSignaux au temps final:" << endl;
        // s.display();

        cout << "\nSignaux en tout temps aux bords du domaine exportes dans '" << s.export_1 << "'"  << endl;
        cout << "Signaux au temps final sur tout le domaine exportes dans '" << s.export_2 << "'"  << endl << endl;

        cout << "================================================"  << endl << endl;
    }
    catch(const string &e){
        cerr << endl << e << endl << endl;
        exit(1);
    }
    catch (Parser::exception_type &e){
        cerr << endl << "ERREUR: Dans l'analyseur de fonctions: " << e.GetMsg() << endl << endl;
        exit(1);
    }
    catch(const exception &e){
        cerr << endl << "ERREUR: \n" << e.what() << endl << endl;
        exit(1);
    }

    return 0;
}