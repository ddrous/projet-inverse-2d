#ifndef INCLUDED_MESH
#define INCLUDED_MESH

#include "config.hpp"

/*****************************************
 * Classe pour créer le maillage
 */
class Mesh{
    public:
        // Parametres du maillage
        double x_min;            // borne gauche du domaine
        double x_max;            // borne droite du domaine
        int N;                  // nombre de mailles/volumes intermediares
        double dx;              // Delta x
        double **cells;         // gauche-centre-droite pour chaque maille (y compris les 2 mailles fantomes)

        // Constructeur
        Mesh(const Config &cfg);

        // Creation des differents volumes uniformes
        void create_cells();
        
        // Destructeur
        virtual ~ Mesh();
};

#endif
