# Transfer Radiatif Inverse 2D (MOCO)

### Objectif
Résolution du problème de propagation de la lumière par la méthode des volumes finis en 2D et apprentissage par un réseau de neurones.

### Documents
- Indications et travail à faire: `doc/guidelines/projetM1.pdf`  
- Rapport de stage en français: `doc/fr/rapport.pdf`  
- Rapport en anglais: `doc/eng/report.pdf`  

## __1ère partie: Résolution de l'EDP__    
_Les commandes indiquées sont à exécuter à partir du répertoire racine du projet._

### Compilation
- `cmake -S . -B build`  
- `cmake --build build`  

### Exécution
- Pour une simple exécution: __`build/transfer src/config/simu.cfg`__ 
- Pour faire une série de simulations: __`bash src/simu/data_dump.sh`__ 

#### Format des fichiers de configuration (`data/img/config.pdf`)   

![Instructions for configuration](data/img/config.png)

### Formats possibles de sauvegarde
- __`data/df_simu.csv`__: exemple de fichier CSV pour les signaux exportés.
- __`data/part_1.sds`__: exemple de fichier binaire SDS (source-densité-signal) pour la sauvegarde d'une série de simulation.

### Visualisation des resultats
- __`src/notebook/visualisation_2d.ipynb`__[(version Jupyter)](https://github.com/desmond-rn/projet-inverse-2d/blob/master/src/notebook/visualisation_2d.ipynb)
  
![Quelques resultats](data/img/energie_flux.png)


## __2ème partie: Apprentissage__   
- Regression: `src/notebook/regression_1d.ipynb` [(version Colab)](https://colab.research.google.com/drive/1kyPV7in4heCWPRQZh4iYFJSlx2lMZqmp?usp=sharing)  
- Classification: `src/notebook/classification.ipynb` [(version Colab)](https://colab.research.google.com/drive/1acJNt3krkJ0rK-RZHzDmxEcn6zAlAv92?usp=sharing)  
- Lecture du format CSV: `src/notebook/sauvegarde_2d.ipynb` [(version Jupyter)](https://github.com/desmond-rn/projet-inverse-2d/blob/master/src/notebook/sauvegarde_2d.ipynb)
- Lecture du format SDS: `src/notebook/format_binaire.ipynb` [(version Colab)](https://colab.research.google.com/drive/1pbuw_aORnOMEFs824BE2jOz4rCTevC3K?usp=sharing).  

## Ressources utilisées:
- __muParser__: pour transformer des expressions en fonctions: [Example](https://beltoforion.de/article.php?a=muparser&s=idExample#idExample) - [Instructions](https://beltoforion.de/article.php?a=muparser&p=building)
