#ifndef DEJA_INCLU_MAILLE     
#define DEJA_INCLU_MAILLE

/*****************************************
 * Classe pour créer le maillage
 */
class Mesh{
    public:
        // Parametres de depart
        double a;
        double b;
        int N;
        double dx;       
        double **cells;

        // Constructeur
        Mesh(double a_hat, double b_hat, int N_hat);

        // Creation des differents volumes uniformes
        void create_cells();
        
        // Destructeur
        virtual ~ Mesh();
};

#endif
