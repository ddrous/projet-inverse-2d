#!/bin/bash

# Script pour generer les donnees en variant 'rho'
#   Ces variables ainsi sont ecrites dans 'src/config/part1.txt'
#   Les constantes pour ces simulation sont ecrites dans 'src/config/part2.txt'
#   Les fichiers 'part1' et 'part2' forment le fichier de configuaration 'src/config/tmp.cfg'

# Pour generer un crenau sur 'rho' situe entre pos_min et pos_max et de hauteur entre h_min et h_max 
pos_x_min=0.692
pos_x_max=0.738
pos_y_min=0.3
pos_y_max=0.7
h_min=1
h_max=2
nbiter_rho_pos_x=4
nbiter_rho_pos_y=40
nbiter_rho_h=2
h=1
a1=$(awk -v h="${h}" -v nbiter="${nbiter_rho_h}" -v min="${h_min}" -v max="${h_max}" 'BEGIN{print ((max-min)/(nbiter-h))}')
b1=$(awk -v min="${h_min}" 'BEGIN{print (min)}')

i=1
a20=$(awk -v i="${i}" -v nbiter="${nbiter_rho_pos_x}" -v min="${pos_x_min}" -v max="${pos_x_max}" 'BEGIN{print ((max-min)/(nbiter-i))}')
b20=$(awk -v min="${pos_x_min}" 'BEGIN{print (min)}')

j=1
a21=$(awk -v j="${j}" -v nbiter="${nbiter_rho_pos_y}" -v min="${pos_y_min}" -v max="${pos_y_max}" 'BEGIN{print ((max-min)/(nbiter-j))}')
b21=$(awk -v min="${pos_y_min}" 'BEGIN{print (min)}')
# Pour generer une source E entre la position min et max
pos_min=0.1
pos_max=0.7
nbiter_source=4
k=1
a3=$(awk -v k="${k}" -v nbiter="${nbiter_source}" -v min="${pos_min}" -v max="${pos_max}" 'BEGIN{print ((max-min)/(nbiter-k))}')
b3=$(awk -v min="${pos_min}" 'BEGIN{print (min)}')

# Indique quelle source activer (haut=0, gauche=1, ou bas=2)
l=1

# Pour comptabiliser les iterations totales
m=1
nbiter_total=$(awk -v nbiter1="${nbiter_rho_pos_x}" -v nbiter2="${nbiter_rho_pos_y}" -v nbiter3="${nbiter_source}" 'BEGIN{print (nbiter1*nbiter2*nbiter3)}')
simu_count=$(echo "simu_count $nbiter_total\n")

simulate () {
    echo "simulation $m sur $nbiter_total en cours ..."
    m=$((m+1))

    echo -e "$rho$source$write_mode$simu_count" > src/config/part14.txt
    cat src/config/part14.txt src/config/part24.txt > src/config/tmp4.cfg
    build/transfer src/config/tmp4.cfg > /dev/null
}

# Boucle des simulations
rho_x=0.15
rho_y=0.15
rho_h=0.1
source_pos_min=0.0
source_pos_max=0.1
# for (( h = 1 ; h <= $nbiter_rho_h; h++ )); do
for (( h = 1 ; h <= 1; h++ )); do
    rho_h=$(awk -v h="${h}" -v a="${a1}" -v b="${b1}" 'BEGIN{print (a*(h-1)+b)}')
    for (( i = 1 ; i <= $nbiter_rho_pos_x; i++ )); do
        rho_x=$(awk -v i="${i}" -v a="${a20}" -v b="${b20}" 'BEGIN{print (a*(i-1)+b)}')
        for (( j = 1 ; j <= $nbiter_rho_pos_y; j++ )); do
            rho_y=$(awk -v j="${j}" -v a="${a21}" -v b="${b21}" 'BEGIN{print (a*(j-1)+b)}')
            rho=$(echo "rho crenau($rho_x,$rho_y,0.1,$rho_h)\n\n")
            for (( l = 1 ; l <= 1; l++ )); do
                for (( k = 1 ; k <= $nbiter_source; k++ )); do
                    source_pos_min=$(awk -v k="${k}" -v a="${a3}" -v b="${b3}" 'BEGIN{print (a*(k-1)+b)}')
                    source_pos_max=$(awk -v pos_min="${source_pos_min}" 'BEGIN{print (pos_min + 0.2)}')

                    if [[ $m -eq 1 ]]
                    then
                    write_mode=$(echo "write_mode truncate\n")
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
