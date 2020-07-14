#include <iostream>
#include <fstream>

#include "include/config.hpp"
#include "include/mesh.hpp"
#include "include/solver.hpp"
#include "include/exporter.hpp"

#include "muParser.h"

using namespace std;
using namespace mu;


int main(int argc, char * argv[]){
    try {
        cout << "\n======== Equation du transfer radiatif en 2D ========" << endl;

        if(argc < 2 || argc > 2)
            throw string("ERREUR: Fournissez un (et un seul) fichier de configuration");

        /* Lecture du fichier config */
        Config cfg = Config(argv[1]);
        cout << "\nLecture des " << cfg.size << " parametres en cours ... " << endl;
        cfg.read();
        cout << "Lecture OK !" << endl;

        /* Creation d'un maillage uniforme */
        Mesh m = Mesh(cfg);
        cout << "\nCreation des " << m.n_cells << " mailles (M = " << m.M << ") en cours ... " << endl;
        m.create_cells();
        // for (int k = 0; k < m.n_cells; k++)
        // {
        //     cout << "k = " << k << " i = " << m.coord[k][0] << " j = " << m.coord[k][1] << endl;
        //     cout << "up = " << m.neighb[k][0] << endl;
        //     cout << "down = " << m.neighb[k][1] << endl;
        //     cout << "left = " << m.neighb[k][2] << endl;
        //     cout << "right = " << m.neighb[k][3] << endl;
        //     cout << " " << endl;
        // }

        cout << "Creation OK !" << endl;
    
        /* Resolution du probleme */
        // Solver s = Solver(&m, cfg);
        // cout << "\nResolution (" << s.step_count << " iterations) en cours ..." << endl;
        // s.solve();
        // cout << "Resolution OK !" << endl;
        // // cout << "\nSignaux au temps final:" << endl;
        // // s.display();

        /* Exportation */
        // Exporter ex = Exporter(&s);
        // cout << "\nExportation des resultats en cours ..." << endl;
        // ex.spatial(cfg.values["export_spatial"], "append");
        // ex.temporal(cfg.values["export_temporal"], "append");

        // ex.case_1("data/case_1_spatial.csv", "data/case_1_temporal.csv");
        // // ex.case_2("data/case_2_spatial.csv", "data/case_2_temporal.csv");
        // ex.case_3("data/case_3_spatial.csv", "data/case_3_temporal.csv");
        // cout << "Export OK !" << endl;

        // cout << "\nResultats du cas test 1 dans:  -- 'data/case_1_spatial.csv'  -- 'data/case_1_temporal.csv'" << endl;
        // if(s.E_exact_expr.empty() == false)
        //     cout << "Resultats du cas test 2 dans:  -- 'data/case_2_spatial.csv'  -- 'data/case_2_temporal.csv'" << endl;
        // cout << "Resultats du cas test 3 dans:  -- 'data/case_3_spatial.csv'" << endl;

        // cout << "\nSignaux en tout temps aux bords du domaine exportes dans '" << cfg.values["export_spatial"] << "'"  << endl;
        // cout << "Signaux au temps final sur tout le domaine exportes dans '" << cfg.values["export_temporal"] << "'"  << endl;

        cout << "\n================================================"  << endl << endl;
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