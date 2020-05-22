# Transfer Radiatif Inverse (MOCO)

### Objectif
Résolution du problème inverse de propagation de la lumière par la méthode des volumes finis et les réseaux de neurones.

### Documents
- Indications et travail à faire: `doc/projetM1.pdf`
- Rapport version 0: `doc/rapport_v0.pdf`
- Rapport version 1: `doc/rapport_v1.pdf`

## __1ère partie: Résolution de l'EDP__    
_Les commandes indiquées sont à exécuter à partir du répertoire racine du projet._

### Compilation
- Avec Cmake: `cmake --build build`

### Exécution
- Pour une simple exécution: __`build/transfer src/config/case_1.cfg`__ [Voir case_2 et case_3 pour les autres cas test](https://github.com/feelpp/csmi-m1-2020-moco-inverse/tree/master/src/config)
Explication du format des fichiers de configuration:  
![Instructions for configuration](data/img/config.png)

- Pour générer un tas de données à étudier: __`bash src/simu/gauss_dump.sh`__ 

### Résultats
- __`data/df_temporal.csv`__ pour les signaux aux bords du domaine en tous temps.
- __`data/df_spatial.csv`__ pour les signaux sur tout le domaine au temps final.

## __2ème partie: Analyse des données__   
- Analyse des données: `src/notebook/analyse_des_donnees.ipynb` [Version Colab](https://colab.research.google.com/drive/17eqqFvVzvzFqB8URGFR9-YQmqDNxU5Ax?usp=sharing).  
- Réseaux de neurones: `src/notebook/reseaux_de_neurones.ipynb` [Version Colab](https://colab.research.google.com/drive/1DXee80oz_6OqLDHdnO00VjK62TdKSE5O?usp=sharing).

## Ressources utilisées:
- __muParser__ pour transformer des expressions en fonctions: [Example](https://beltoforion.de/article.php?a=muparser&s=idExample#idExample) - [Instructions](https://beltoforion.de/article.php?a=muparser&p=building)
