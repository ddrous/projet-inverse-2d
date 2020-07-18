# Transfer Radiatif Inverse 2D (MOCO)

### Objectif
Résolution du problème inverse de propagation de la lumière par la méthode des volumes finis en 2D et les réseaux de neurones.

### Documents
- Indications et travail à faire: `doc/projetM1.pdf`  
- Rapport de stage: `doc/rapport.pdf`  

## __1ère partie: Résolution de l'EDP__    
_Les commandes indiquées sont à exécuter à partir du répertoire racine du projet._

### Compilation
- Avec Cmake: `cmake --build build`

### Exécution
- Pour une simple exécution: __`build/transfer src/config/simu.cfg`__ (voir ci dessous comment écrite un fichier de configuration)

__Format des fichiers de configuration:__  
![Instructions for configuration](data/img/config.png)

- Pour générer un tas de données à étudier: __`bash src/simu/data_dump.sh`__ 

### Résultats
- __`data/df_simu.csv`__ pour les signaux a exporter (densité et ses attributs eventuels, énergie, flux et température).

## __2ème partie: Reseaux de neurones__   
- Analyse des données: `src/notebook/analyse_des_donnees.ipynb` [Version Colab](https://colab.research.google.com/drive/17eqqFvVzvzFqB8URGFR9-YQmqDNxU5Ax?usp=sharing).  
- Réseaux de neurones: `src/notebook/reseaux_de_neurones.ipynb` [Version Colab](https://colab.research.google.com/drive/1DXee80oz_6OqLDHdnO00VjK62TdKSE5O?usp=sharing).
- Visualisation et animation: `src/notebook/visualisation_2d.ipynb`.

## Ressources utilisées:
- __muParser__ pour transformer des expressions en fonctions: [Example](https://beltoforion.de/article.php?a=muparser&s=idExample#idExample) - [Instructions](https://beltoforion.de/article.php?a=muparser&p=building)
