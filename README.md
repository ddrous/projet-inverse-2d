# Transfer Radiatif Inverse 2D (MOCO)

### Objectif
Résolution du problème de propagation de la lumière par la méthode des volumes finis en 2D et apprentissage par un réseau de neurones.

### Documents
- Indications et travail à faire: `doc/projetM1.pdf`  
- Rapport de stage: `doc/rapport.pdf`  

## __1ère partie: Résolution de l'EDP__    
_Les commandes indiquées sont à exécuter à partir du répertoire racine du projet._

### Compilation
Avec Cmake (utiliser __GCC 7.4.0__ ou __Clang 6.0.0__ en cas d'erreurs):   
- `cmake -S . -B build`  
- `cmake --build build`  

### Exécution
- Pour une simple exécution: __`build/transfer src/config/simu.cfg`__ 
- Pour générer un tas de données à étudier: __`bash src/simu/data_dump.sh`__ 

__Format des fichiers de configuration (`data/img/config.pdf`):__   

![Instructions for configuration](data/img/config.png)


### Résultats
![Quelques resultats](data/img/energie_flux.png)

- __`data/df_simu.csv`__: fichier CSV pour les signaux à exporter (densité, énergie, flux et température).
- __`data/part_1.sds`__: fichier binaire SDS (source-densité-signal) pour les les sauvegardes.

## __2ème partie: Apprentissage__   
- Regression: `src/notebook/regression_1d.ipynb` [(version Colab)](https://colab.research.google.com/drive/1kyPV7in4heCWPRQZh4iYFJSlx2lMZqmp?usp=sharing)  
- Classification: `src/notebook/classification.ipynb` [(version Colab)](https://colab.research.google.com/drive/1acJNt3krkJ0rK-RZHzDmxEcn6zAlAv92?usp=sharing)  
- Animation des résultats: `src/notebook/visualisation_2d.ipynb`
- Lecture du format CSV: `src/notebook/sauvegarde_2d.ipynb`
- Lecture du format SDS: `src/notebook/format_binaire.ipynb` [(version Colab)](https://colab.research.google.com/drive/1pbuw_aORnOMEFs824BE2jOz4rCTevC3K?usp=sharing).  

## Ressources utilisées:
- __muParser__: pour transformer des expressions en fonctions: [Example](https://beltoforion.de/article.php?a=muparser&s=idExample#idExample) - [Instructions](https://beltoforion.de/article.php?a=muparser&p=building)
