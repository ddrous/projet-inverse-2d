#ifndef INCLUDED_EXPORT
#define INCLUDED_EXPORT

#include <map>
#include <fstream>

#include "solver.hpp"

/************************************************
 * Classe pour exporter les donnees
 */
class Exporter{
    public:
        // solveur dont on eporte les resultats
        Solver *solver;
        /***************
         * Constructeur
         */
        Exporter(Solver *new_solver);

        /**************
         * Exporte une dataframe de donnees spatiales
         */
        void spatial(std::string file_name, std::string mode);

        /**************
         * Exporte une dataframe de donnees temportelles
         */
        void temporal(std::string file_name, std::string mode);

        /**************
         * Exporte le cas test 1
         */
        void case_1(std::string file_spatial, std::string file_temporal);

        /**************
         * Exporte le cas test 2
         */
        void case_2(std::string file_spatial, std::string file_temporal);

        /**************
         * Exporte le cas test 3
         */
        void case_3(std::string file_spatial, std::string file_temporal);

        /***************
         * Destructeur
         */
        virtual ~Exporter(){};
};

#endif
