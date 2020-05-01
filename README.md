# csmi-m1-2020-moco-inverse

Pour la résolution du problème inverse de propagation de la lumière par l'utilisation de la méthode des volumes finis et des réseaux de neurones.  

Les différent rapports se trouvent dans le repertoire `doc`.

## __1ere partie: Résolution de l'EDP__    
Les commandes indiquées sont a exécuter à partir du repertoire racine du projet.   

Pour lancer le programme: `build/spreading tools/config_1.cfg`     

Les fichiers de configuration `config_*.cfg` se trouvent dans le repertoire `tools`.   

Les fichiers sources sont dans le repertoire `src`.   

Pour compiler le programme: 
`g++ -I thirdparty/eigen-3.3.7 -I thirdparty/exprtk src/spreading.cpp src/solveur.cpp src/maille.cpp src/config.cpp -o build/spreading`   

Les exécutables sont générés dans le repertoire `build`.    

Les résultats sont exportés dans les fichiers `data/log1.csv` pour les signaux aux bords du domaine en tous temps; et `data/log2.csv` pour les signaux sur tout le domaine au temps final.


## __2ème partie: Analyse des données__   
A partir du 5 mai.


### Ce que le programme peut faire
- Résoudre numériquement l'EDP du transfert radiatif
- Exporter les données   
- Admettre un fichier de configuration
- Prendre en entrée des fonction a t=0 pour E, F, et T 

### Ce que le programme ne peut pas faire  
- Compiler avec CMAKE
- Statistique descriptive sur les données   

### Links to some ressources:
- For vector operations: [Eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html) 
- For function parsing: [ExprTk](http://www.partow.net/programming/exprtk/)

