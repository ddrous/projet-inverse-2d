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
        // noms des parametres et leurs identifiants
        std::map<std::string, int> names;
        // nombre de parametres
        int size;
        // valeurs des parametres (tous comme des chaines de caracteres)
        std::string *values;

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
