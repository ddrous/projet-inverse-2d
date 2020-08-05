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
        // noms des parametres et leurs valeurs
        std::map<std::string, std::string> values;
        // nombre maximal de parametres
        int size;
        // nombre effectif de parametres
        int n_param;

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
