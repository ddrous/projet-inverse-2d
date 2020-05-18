# Transfer Radiatif Inverse (MOCO)

### Objectif
Résolution du problème inverse de propagation de la lumière par l'utilisation de la méthode des volumes finis et des réseaux de neurones.

### Documentation
- Indications: `doc/projetM1.pdf`
- Rapport version 0: `doc/rapport_v0.pdf`
- Rapport version 1: `doc/rapport_v1.pdf`

## __1ère partie: Résolution de l'EDP__    
_Les commandes sont à exécuter à partir du répertoire racine du projet._

### Compilation
- Avec Cmake: `cmake --build build`
- Avec g++: `g++ -I thirdparty/muparser/include src/main.cpp src/solver.cpp src/mesh.cpp src/config.cpp thirdparty/muparser/src/muParser.cpp thirdparty/muparser/src/muParserBase.cpp thirdparty/muparser/src/muParserBytecode.cpp thirdparty/muparser/src/muParserCallback.cpp thirdparty/muparser/src/muParserError.cpp thirdparty/muparser/src/muParserTokenReader.cpp -o build/transfer`   

### Exécution
- Pour une seule exécution: __`build/transfer src/config/dump.cfg`__ [Voir fichiers de configuration](https://github.com/feelpp/csmi-m1-2020-moco-inverse/tree/master/src/config)
- Pour générer les données a étudier: __`bash src/simu/data_dump.sh`__

Les résultats sont exportés dans le répertoire `data`:
- __`df_temporal.csv`__ pour les signaux aux bords du domaine en tous temps. 
Ses colonnes sont `x_min, x_max, N, c, a, C_v, CFL, epsilon, t_final, rho_expr, sigma_a_expr, sigma_c_expr, E_init_expr, F_init_expr, T_init_expr, dt, step_count, t, E_left, E_right, F_left, F_right, T_left, T_right`
- __`df_spatial.csv`__ pour les signaux sur tout le domaine au temps final. 
Ses colonnes sont `x_min, x_max, N, c, a, C_v, CFL, epsilon, t_final, rho_expr, sigma_a_expr, sigma_c_expr, E_init_expr, F_init_expr, T_init_expr, dt, step_count, x, rho, sigma_a, sigma_c, E_init, E_final, F_init, F_final, T_init, T_final`

## __2ème partie: Analyse des données__   
- Analyse des données: `src/notebook/analyse_des_donnees.ipynb` [Version Colab](https://colab.research.google.com/drive/17eqqFvVzvzFqB8URGFR9-YQmqDNxU5Ax?usp=sharing).  
- Réseaux de neurones: `src/notebook/reseaux_de_neurones.ipynb` [Version Colab](https://colab.research.google.com/drive/1DXee80oz_6OqLDHdnO00VjK62TdKSE5O?usp=sharing).

## Ressources utilisées:
- __muParser__ pour transformer des expressions en fonctions: [Example](https://beltoforion.de/article.php?a=muparser&s=idExample#idExample) - [Instructions](https://beltoforion.de/article.php?a=muparser&p=building)
