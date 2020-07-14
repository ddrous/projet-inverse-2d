#ifndef INCLUDED_MESH
#define INCLUDED_MESH

#include <vector>

#include "config.hpp"

/*****************************************
 * Classe pour créer le maillage
 */
class Mesh{
    public:
        // Parametres du maillage
        double x_min;            // borne gauche du domaine
        double x_max;            // borne droite du domaine
        double y_min;            // borne du haut du domaine
        double y_max;            // borne du bas du domaine
        int N;                  // nombre de mailles/volumes en verticales
        int M;                  // nombre de mailles/volumes en horizontale
        int n_cells;            // nombre de mailles/volumes totales

        double dx;              // Delta x
        double dy;              // Delta y

        std::vector<double> x;         // abscisses des centres des mailles
        std::vector<double> y;         // ordonnees des centres des mailles

        int **coord;       // indices i,j identifiants chaque maille
        int **neighb;        // numero des 4 mailles voisines

        // Constructeur
        Mesh(const Config &cfg);

        // Creation des differents volumes uniformes
        void create_cells();
        
        // Destructeur
        virtual ~ Mesh();
};

#endif
