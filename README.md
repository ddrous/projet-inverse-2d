# csmi-m1-2020-moco-inverse

Pour la résolution du problème inverse de propagation de la lumière par l'utilisation de la méthode des volumes finis et des réseaux de neurones.  

Les différent rapports se trouvent dans le répertoire `doc`.


## __1ere partie: Résolution de l'EDP__    
Les commandes indiquées sont a exécuter à partir du répertoire racine du projet.   

Les fichiers sources sont dans le répertoire `src`. Pour compiler le programme, on pourra: 
- Utiliser g++: `g++ -I thirdparty/exprtk src/main.cpp src/solver.cpp src/mesh.cpp src/config.cpp -o build/transfer`   
- Utiliser Cmake: `cmake --build build`

L'exécutable `transfer` est généré dans le répertoire `build`. Les fichiers de configuration `*.cfg` se trouvent dans le répertoire `src/config`. Pour lancer le programme, on pourra executer: __`build/transfer src/config/rho.cfg`__     

Les résultats sont exportés dans le répertoire `data`. Par exemple:
- __`df_temporal.csv`__ pour les signaux aux bords du domaine en tous temps. 
Ses colonnes sont `x_min, x_max, N, c, a, C_v, CFL, epsilon, t_final, rho_exp, sigma_a_expr, sigma_c_expr, E_x_0_expr, F_x_0_expr, T_x_0_expr, dt, time_steps, E_0, E_N, F_0, F_N, T_0, T_N`

- __`df_spatial.csv`__ pour les signaux sur tout le domaine au temps final. 
Ses colonnes sont `x_min, x_max, N, c, a, C_v, CFL, epsilon, t_final, x, rho_exp, sigma_a_expr, sigma_c_expr, E_x_0_expr, F_x_0_expr, T_x_0_expr, dt, time_steps, rho, sigma_a, sigma_c, E_init, E_final, F_init, F_final, T_init, T_final`

Pour generer plusieurs lignes dans les fichier d'export d'un seul coup, on peut executer le script `bash src/simu/rho.sh`. Ce script lance une simulation qui fait varier `rho` entre `rho_min` et `rho_max` (voir `src/simu/rho_1.txt`) tout en gardant le reste des parametres constants (voir `src/simu/rho_2.txt`). 


## __2ème partie: Analyse des données__   
- L'analyse des donnes se fait dans le notebook `src/notebook/analyse_des_donnees.ipynb` [Visualiser le notebook sur Colab](https://colab.research.google.com/drive/17eqqFvVzvzFqB8URGFR9-YQmqDNxU5Ax?usp=sharing).  
- La prediction grace aux reseaux de neuronnes se fera dans `src/notebook/reseaux_de_neurones.ipynb` [Visualiser le notebook sur Colab](https://colab.research.google.com/drive/1DXee80oz_6OqLDHdnO00VjK62TdKSE5O?usp=sharing).


## Lien vers les ressources utilisees:
- Pour transformer des expressions en fonctions: [ExprTk](http://www.partow.net/programming/exprtk/)
