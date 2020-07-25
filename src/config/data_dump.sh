#!/bin/bash

# Script pour generer les donnees en variant 'rho'
#   Ces variables ainsi sont ecrites dans 'src/config/part1.txt'
#   Les constantes pour ces simulation sont ecrites dans 'src/config/part2.txt'
#   Les fichiers 'part1' et 'part2' forment le fichier de configuaration 'src/config/tmp.cfg'

# Pour generer 'rho' de position entre pos_min et pos_max et de hauteur entre h_min et h_max 
pos_min=0.2
pos_max=0.8
h_min=1
h_max=10
nbiter_rho_pos=4
nbiter_rho_h=2
h=1
a1=$(awk -v h="${h}" -v nbiter="${nbiter_rho_h}" -v min="${h_min}" -v max="${h_max}" 'BEGIN{print ((max-min)/(nbiter-h))}')
b1=$(awk -v min="${h_min}" 'BEGIN{print (min)}')

i=1
j=1
a2=$(awk -v i="${i}" -v nbiter="${nbiter_rho_pos}" -v min="${pos_min}" -v max="${pos_max}" 'BEGIN{print ((max-min)/(nbiter-i))}')
b2=$(awk -v min="${pos_min}" 'BEGIN{print (min)}')

# Pour generer une source E entre la position min et max
pos_min=0.1
pos_max=0.7
nbiter_source=4
k=1
a3=$(awk -v k="${k}" -v nbiter="${nbiter_source}" -v min="${pos_min}" -v max="${pos_max}" 'BEGIN{print ((max-min)/(nbiter-k))}')
b3=$(awk -v min="${pos_min}" 'BEGIN{print (min)}')

# Indique quelle source activer (haut=0, gauche=1, ou bas=2)
l=0

# Pour omptabiliser les iterations totales
m=1
nbiter_rho_h=1
nbiter_total=$(awk -v nbiter1="${nbiter_rho_h}" -v nbiter2="${nbiter_rho_pos}" -v nbiter3="${nbiter_source}" 'BEGIN{print (nbiter1*nbiter2*nbiter2*nbiter3*3)}')

simulate () {
    echo "simulation $m sur $nbiter_total en cours ..."
    m=$((m+1))

    echo -e "$rho$source$write_mode" > src/config/part1.txt
    cat src/config/part1.txt src/config/part2.txt > src/config/tmp.cfg
    # build/transfer src/config/tmp.cfg > /dev/null
}

# Boucle des simulations
rho_x=0.15
rho_y=0.15
rho_h=0.1
source_pos_min=0.0
source_pos_max=0.1
for (( h = 1 ; h <= $nbiter_rho_h; h++ )); do
    rho_h=$(awk -v h="${h}" -v a="${a1}" -v b="${b1}" 'BEGIN{print (a*(h-1)+b)}')
    for (( i = 1 ; i <= $nbiter_rho_pos; i++ )); do
        rho_x=$(awk -v i="${i}" -v a="${a2}" -v b="${b2}" 'BEGIN{print (a*(i-1)+b)}')
        for (( j = 1 ; j <= $nbiter_rho_pos; j++ )); do
            rho_y=$(awk -v j="${j}" -v a="${a2}" -v b="${b2}" 'BEGIN{print (a*(j-1)+b)}')
            rho=$(echo "rho crenau($rho_x,$rho_y,0.1,$rho_h)\n\n")
            for (( l = 0 ; l <= 2; l++ )); do
                for (( k = 1 ; k <= $nbiter_source; k++ )); do
                    source_pos_min=$(awk -v k="${k}" -v a="${a3}" -v b="${b3}" 'BEGIN{print (a*(k-1)+b)}')
                    source_pos_max=$(awk -v pos_min="${source_pos_min}" 'BEGIN{print (pos_min + 0.2)}')

                    if [[ $m -eq 1 ]]
                    then
                    write_mode=$(echo "write_mode trunc\n")
                    else
                    write_mode=$(echo "write_mode append\n")
                    fi

                    if [[ $l -eq 0 ]]
                    then
                        source=$(echo "E_u ponctuel($source_pos_min,$source_pos_max)\nF_u_x 0\nF_u_y 0\nT_u 5\n\nE_d neumann\nF_d_x neumann\nF_d_y neumann\nT_d neumann\n\nE_l neumann\nF_l_x neumann\nF_l_y neumann\nT_l neumann\n\n")
                        simulate
                    elif [[ $l -eq 1 ]]
                    then
                        source=$(echo "E_l ponctuel($source_pos_min,$source_pos_max)\nF_l_x 0\nF_l_y 0\nT_l 5\n\nE_d neumann\nF_d_x neumann\nF_d_y neumann\nT_d neumann\n\nE_u neumann\nF_u_x neumann\nF_u_y neumann\nT_u neumann\n\n")
                        simulate
                    else
                        source=$(echo "E_d ponctuel($source_pos_min,$source_pos_max)\nF_d_x 0\nF_d_y 0\nT_d 5\n\nE_l neumann\nF_l_x neumann\nF_l_y neumann\nT_l neumann\n\nE_u neumann\nF_u_x neumann\nF_u_y neumann\nT_u neumann\n\n")
                        simulate
                    fi
                done
            done
        done
    done
done
