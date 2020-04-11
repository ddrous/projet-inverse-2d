#include <iostream>
#include <fstream>

#include "maille.hpp"
#include "solveur.hpp"

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




// using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;

int main(int argc, char * argv[]){

    double x_min = 0;
    double x_max = 1;
    int N = 10;

    
    Mesh m1 = Mesh(x_min, x_max, N);
 
    CHK_EX(
        
        m1.create_cells();

        std::cout << "my mesh" << "\ta = "<< m1.a << "\tb = "<< m1.b << "\tN = "<< m1.N << "\tdx = " << m1.dx << std::endl;
        std::cout << "Mailles: \n" << m1.cells << std::endl;

    );

    double a = 5.670374419*1e-8;            // Constante de Boltzmann
    double C_v = 3500;          // Capacite thermique moyenne du coprs humain
    string rho = "rho(x)";      // Fonction a definir
    string sigma_a = "sigma_a(rho, T)";
    string sigma_c = "sigma_c(rho, T)";
    double CFL = 0.99;
    double epsilon = 1e-6;
    double Temps_max = 1e-6;
    string E_x_0 = "E(0, x)";
    string E_0_t = "E(t, 0)";
    string E_N_t = "E(t, N)";   // A mesurer a droite du domaine
    string F_x_0 = "F(0, x)";
    string F_0_t = "F(t, 0)";
    string F_N_t = "F(t, N)";   // A mesurer
    string T_x_0 = "T(0, x)";
    string T_0_t = "T(t, 0)";
    string T_N_t = "T(t, N)";   // A mesurer

    
    Solver s1 = Solver(&m1, a, C_v, rho, sigma_a, sigma_c, CFL, epsilon, Temps_max, E_x_0, E_0_t, E_N_t, F_x_0, F_0_t, F_N_t, T_x_0, T_0_t, T_N_t);

    CHK_EX(

        s1.solve("../data/log.dat");
        s1.dislay();
        s1.export_csv("../data/log.csv");
    
    );
    

    return 0;
}