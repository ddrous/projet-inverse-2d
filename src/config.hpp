#ifndef DEJA_INCLU_CONFIG     
#define DEJA_INCLU_CONFIG

#include <map>

/************************************************
 * Classe pour lire le fichier de configuration
 */
class Config{
    public:
        // nom du fichier
        std::string nom_fichier;
        // noms des variables
        std::map<std::string, int> noms;
        // valeurs des parametres qui sont des flotants
        double *valeurs;
        // valeurs des parametres qui sont des expressions
        std::string *fonctions;


        /***************
         * Constructeur
         */
        Config(std::string nom_du_fichier);


        /**************
         * Fonction pour lire le fichier.cfg
         */
        void read();

        /***************
         * Destructeur
         */
        virtual ~Config();
};

#endif
