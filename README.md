# csmi-m1-2020-moco-inverse

Pour la résolution du problème inverse de propagation de la lumière par l'utilisation de la méthode des volumes finis et des réseaux de neurones.

## __1ere partie: Résolution de l'EDP__    
Pour lancer le programme, on pourra exécuter la commande: `g++ -I ./eigen-3.3.7 spreading.cpp solveur.cpp maille.cpp -o spreading && ./spreading` a partir du répertoire `src`.  

Les résultats sont exportés dans les fichiers `data/log.dat` pour les signaux aux bords du domaine a tous temps; et `data/log.csv` pour les signaux en tout x au temps final.


## __2ème partie: Réseaux de neurones__
Pas encore commencée.   
A partir du 29 avril


### Ce que le programme peut faire en ce moment
- Résoudre numériquement l'EDP du transfert radiatif
- Exporter les données   


### Ce que le programme ne peut pas encore faire  
- Prendre en entrée des fonction a t=0 pour E, F, et T 
- Admettre un fichier de configuration
- Compiler avec CMAKE
- Prédire a laide des réseaux de neurones   

### Links to some ressources:
- For vector operations: [Eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html) 
- For function parsing: [cparse](https://github.com/cparse/cparse/wiki/Getting-Started)


