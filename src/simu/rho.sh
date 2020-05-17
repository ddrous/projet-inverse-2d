#!/bin/bash

# Script pour generer les donnees en varaint la densite 'rho'
#   la variable 'rho' ainsi que les nons des fichiers ou exporter sont ecrits dans 'rho_1.txt'
#   les constantes pour cette simulation sont ecrites dans 'rho_1.txt'
#   les fichiers 'rho_1.txt' et 'rho_2.txt' forment le fichie de configuaration dans 'src/config/rho.cfg'

# fichiers d'exportation
file_1=data/df_temporal.csv
file_2=data/df_spatial.csv

# Pour vider les fichiers d'exportation. Et remplacer le contenu par les en-tetes du csv
overwrite=true
if [ "$overwrite" = true ] ; then
    cat src/simu/headers_temporal.txt > "$file_1"
    cat src/simu/headers_spatial.txt > "$file_2"
fi

# Pour generer 'rho' entre rho_min et rho_max
rho_min=1
rho_max=10
nbiter=1
i=0
a=$(awk -v i="${i}" -v nbiter="${nbiter}" -v rho_min="${rho_min}" -v rho_max="${rho_max}" 'BEGIN{print ((rho_max-rho_min)/(nbiter-i))}')
b=$(awk -v rho_min="${rho_min}" 'BEGIN{print (rho_min)}')
rho=0
while [ $i -lt $nbiter ]; do
    rho=$(awk -v i="${i}" -v a="${a}" -v b="${b}" 'BEGIN{print (a*i+b)}')
    echo -e "rho $rho\n" > src/simu/rho_1.txt
    # echo -e "rho exp(x)\n" > src/simu/rho_1.txt
    cat src/simu/rho_1.txt src/simu/rho_2.txt > src/config/rho.cfg
    build/transfer src/config/rho.cfg > /dev/null
	i=$((i+1))
done