#include <iostream>
#include <fstream>

#include "maille.hpp"
#include "solveur.hpp"
#include "config.hpp"

using namespace std;

int main(int argc, char * argv[]){
    cout << endl;

    try {
        if(argc < 2)
        throw string("Fournissez un fichier de configuration");

        // Lecture du fichier config
        Config cfg1 = Config(argv[1]);
        cfg1.read();

        // Creation d'un maillage uniforme
        double x_min = cfg1.valeurs[0];
        double x_max = cfg1.valeurs[1];
        int N = (int)cfg1.valeurs[2];
        Mesh m1 = Mesh(x_min, x_max, N);
        m1.create_cells();
        cout << "\nLes configurations du maillage:" << "\ta = "<< m1.a << "\tb = "<< m1.b << "\tN = "<< m1.N << "\tdx = " << m1.dx << endl;
    
        // Resolution du probleme
        Solver s1 = Solver(&m1, &cfg1.valeurs[3], cfg1.fonctions);
        cout << "\nResolution en cours ..." << endl;
        s1.solve();
        cout << "Resolution OK!" << endl;
        // cout << "\nSignaux au temps final:" << endl;
        // s1.dislay();
        s1.export_final();

        cout << "\nSignaux en tout temps aux bords du domaine exportes dans '" << s1.export_1 << "'"  << endl;
        cout << "Signaux au temps final sur tout le domaine exportes dans '" << s1.export_2 << "'"  << endl;
    }
    catch(const string& e){
        cerr << e << endl;
        exit(1);
    }
    catch(const exception& e){
        cerr << e.what() << endl;
        exit(1);
    }

    cout << endl;
    return 0;
}