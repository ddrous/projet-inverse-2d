#ifndef INCLUDED_CONFIG
#define INCLUDED_CONFIG

#include <map>

/************************************************
 * Classe pour lire le fichier de configuration
 */
class Config{
    public:
        // nom du fichier config
        std::string file_name;
        // noms des variables et leurs identifiants
        std::map<std::string, int> names;
        // nombre de variables a lire
        int size;
        // valeurs des parametres qui sont des nombres
        double *doubles;
        // valeurs des parametres qui sont des chaines de caracteres (expressons de fonctions)
        std::string *strings;

        /***************
         * Constructeur
         */
        Config(std::string file_path);

        /**************
         * Fonction pour lire le fichier .cfg
         */
        void read();

        /***************
         * Destructeur
         */
        virtual ~Config();
};

#endif
