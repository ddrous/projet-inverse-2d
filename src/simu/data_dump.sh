#!/bin/bash

# Script pour generer les donnees en variant de facon aleatoire nos parametres

nbiter=50
i=0
for (( i = 0 ; i < $nbiter; i++ )); do
    build/transfer src/config/tmp.cfg > /dev/null
done
