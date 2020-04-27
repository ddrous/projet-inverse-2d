# csmi-m1-2020-moco-inverse

Pour la résolution du problème inverse de propagation de la lumière par l'utilisation de la méthode des volumes finis et des réseaux de neurones.  

Les différent rapports se trouvent dans le repertoire `doc`.

## __1ere partie: Résolution de l'EDP__    
Pour lancer le programme, on pourra exécuter la commande: `g++ -I thirdparty/eigen-3.3.7 src/spreading.cpp src/solveur.cpp src/maille.cpp src/config.cpp -o build/spreading && build/spreading src/config.cfg` a partir du repertoire racine du projet.    

Les fichiers source et les fichiers de configurations sont dans le repertoire `src`. 

Les exécutables sont générés dans le repertoire `build`.  

Les résultats sont exportés dans les fichiers `data/log1.csv` pour les signaux aux bords du domaine a tous temps; et `data/log2.csv` pour les signaux en tout x au temps final.


## __2ème partie: Analyse des données__   
A partir du 29 avril


### Ce que le programme peut faire en ce moment
- Résoudre numériquement l'EDP du transfert radiatif
- Exporter les données   
- Admettre un fichier de configuration


### Ce que le programme ne peut pas encore faire  
- Prendre en entrée des fonction a t=0 pour E, F, et T 
- Compiler avec CMAKE
- Prédire a laide des réseaux de neurones   

### Links to some ressources:
- For vector operations: [Eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html) 
- For function parsing: [cparse](https://github.com/cparse/cparse/wiki/Getting-Started)


