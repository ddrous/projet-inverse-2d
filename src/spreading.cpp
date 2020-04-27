#include <iostream>
#include <fstream>

#include "maille.hpp"
#include "solveur.hpp"
#include "config.hpp"

/************************************************
 * Macro pour vérifier les exceptions
 * 
 */
#define CHK_EX(v)         \
    do{                 \
        try {v}           \
        catch(const std::string& e){        \
            std::cerr << e << std::endl;        \
            exit(1);                            \
        }                           \
        catch(const std::exception& e){                  \
            std::cerr << #v << std::endl << e.what() << std::endl;       \
            exit(1);        \
        }                       \
    } while (0)


using namespace Eigen;
using namespace std;

int main(int argc, char * argv[]){


    if(argc < 2)
        CHK_EX(throw string("Fournissez un fichier de configuration"););

    // Lecture du fichier config
    Config cfg1 = Config(argv[1]);
    CHK_EX(cfg1.read(););

    // Config du maillage
    double x_min = cfg1.valeurs[0];
    double x_max = cfg1.valeurs[1];
    int N = (int)cfg1.valeurs[2];

    // Affichage du maillage
    Mesh m1 = Mesh(x_min, x_max, N);
    CHK_EX(
        m1.create_cells();
        std::cout << "\n\nLes config du maillage:" << "\ta = "<< m1.a << "\tb = "<< m1.b << "\tN = "<< m1.N << "\tdx = " << m1.dx << std::endl;
        std::cout << "Les volumes: left - center - right \n" << m1.cells << std::endl;
    );
   
    // Resolution du probleme
    Solver s1 = Solver(&m1, cfg1.valeurs, cfg1.fonctions);

    CHK_EX(
        s1.solve();
        cout << "\nSignaux au temps final:" << endl;
        s1.dislay();
        s1.export_final();
    
    );
    
    return 0;
}