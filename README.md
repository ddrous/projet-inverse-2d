# csmi-m1-2020-moco-inverse

Pour la résolution du problème inverse de propagation de la lumière par l'utilisation de la méthode des volumes finis et des réseaux de neurones.  

Les différent rapports se trouvent dans le répertoire `doc`.

## __1ere partie: Résolution de l'EDP__    
Les commandes indiquées sont a exécuter à partir du répertoire racine du projet.   

Pour lancer le programme: `build/transfer src/config/rho.cfg`     

Les fichiers de configuration `*.cfg` se trouvent dans le répertoire `src/config`.   

Les fichiers sources sont dans le répertoire `src`.   

Pour compiler le programme: 
`g++ -I thirdparty/exprtk src/main.cpp src/solver.cpp src/mesh.cpp src/config.cpp -o build/transfer`   

L'exécutable __`transfer`__ est généré dans le répertoire `build`.    

Les résultats sont exportés dans le répertoire `data`:
> `df_temporal.csv` pour les signaux aux bords du domaine en tous temps. 
Les colonnes sont `x_min, x_max, N, c, a, C_v, CFL, epsilon, t_final, rho_exp, sigma_a_expr, sigma_c_expr, E_x_0_expr, F_x_0_expr, T_x_0_expr, E_0, E_N, F_0, F_N, T_0, T_N`
> `df_spatial.csv` pour les signaux sur tout le domaine au temps final. Les colonnes sont `x_min, x_max, N, c, a, C_v, CFL, epsilon, t_final, x, rho_exp, sigma_a_expr, sigma_c_expr, E_x_0_expr, F_x_0_expr, T_x_0_expr, rho, sigma_a, sigma_c, E_x_0, E, F_x_0, F, T_x_0, T`


## __2ème partie: Analyse des données__   
A partir du 13 mai.


### Ce que le programme peut faire
- Résoudre numériquement l'EDP du transfert radiatif
- Exporter les données   
- Admettre un fichier de configuration
- Prendre en entrée des fonction a t=0 pour E, F, et T 
- Compiler avec CMAKE

### Ce que le programme ne peut pas faire  
- Statistique descriptive sur les données   

### Links to some ressources:
- For function parsing: [ExprTk](http://www.partow.net/programming/exprtk/)

