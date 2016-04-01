#!/bin/bash
#networks reconstruction
echo -n "You are running Gene4x"
echo -n "Multiplex reconstruction... "
mkdir temp
echo -n "Enter the path to the mRNA dataset file > "
read text1
./ARACNE/aracne2 -i $text1 -o exp_net.adj -H ./ARACNE #exp_net.adj es un output de aracne vainilla... aquí podríamos meter parallel ARACNE ... also, un echo de que esto ya acabo...
awk 'NR>=18' ./exp_net.adj > ./exp_net_1.txt #exp_net_1.txt es el adj sin las lineas de comentarios
perl ./code/conversion_network1.pl -g1 exp_net_1.txt -g2 exp_net.txt #exp_net.txt es un "sif-like" del tipo v1 v2 mi
rm exp_net.adj
rm exp_net_1.txt
cp exp_net.txt ./temp/
#nodes selection
cat  $text1 | awk '{ print $1 }' > ./temp/exp1.txt #lista de genes
awk 'NR>=2' ./temp/exp1.txt > ./temp/exp.txt #resto de la matriz de expresion
rm ./temp/exp1.txt
split -l $(( $( wc -l < ./temp/exp_net.txt ) / 2 + 1 )) ./temp/exp_net.txt ./temp/exp_net.txt
mv ./temp/exp_net.txtaa ./temp/exp_net_primaparte.txt #output de split
mv ./temp/exp_net.txtab ./temp/exp_net_secondaparte.txt #output de split
Rscript ./code/unique_nodes.R #este codigo necesita que todo lo que esta en el script este en data, NO en la carpeta interna. Also, podemos optimizar con un fread de data.table, aunque perdamos universalidad. Todos deberian tener data.table 
cat ./temp/net_exp_unique_nodes_primaparte.txt ./temp/net_exp_unique_nodes_secondaparte.txt > ./temp/net_exp_unique_nodes.txt
#integer conversion .... aqui saca unos sif-likes que parece pega con el sed en el multiplex...
Rscript ./code/conversion_integer_exp.R
Rscript ./code/conversion_integer_ppi.R
Rscript ./code/conversion_integer_mirna.R
Rscript ./code/conversion_integer_tf.R
sed -i 's/"//g' ./temp/net_*_unique_nodes_int.txt #quitar las comillas del output de los rscript... se pueden pulir para que el output no los tenga
echo -n "Multiplex reconstruction finished "
#filtering
echo -n "Layers filtering... "
python ./code/SBV_filter_exp.py
python ./code/SBV_filter_mirna.py
python ./code/SBV_filter_tf.py
Rscript ./code/SBV_threshold_exp.R
Rscript ./code/SBV_threshold_mirna.R
Rscript ./code/SBV_threshold_tf.R
echo -n "Layers filtering... finished"
echo -n "Enter the integer corrisponding to the optimal alpha value for the expression network (0=0.005,1=0.01,2=0.02,3=0.03,4=0.05,5=0.1,6=0.2,7=0.3,8=0.4,9=0.5)> "
read alpha_exp
echo -n "Enter the integer corrisponding to the optimal alpha value for the transcription factor co-taregting network (0=0.005,1=0.01,2=0.02,3=0.03,4=0.05,5=0.1,6=0.2,7=0.3,8=0.4,9=0.5)> "
read alpha_tf
echo -n "Enter the integer corrisponding to the optimal alpha value for the microRNA co-taregting network (0=0.005,1=0.01,2=0.02,3=0.03,4=0.05,5=0.1,6=0.2,7=0.3,8=0.4,9=0.5)> "
read alpha_mirna
sed -i "s/{'weight': //g" "./temp/test.edgelist_exp"$alpha_exp".txt" 
sed -i "s/{'weight': //g" "./temp/test.edgelist_tf"$alpha_tf".txt" 
sed -i "s/{'weight': //g" "./temp/test.edgelist_mirna"$alpha_mirna".txt" 
sed  "s/}//g" "./temp/test.edgelist_exp"$alpha_exp".txt" > "./temp/net_exp_final.txt"
sed  "s/}//g" "./temp/test.edgelist_tf"$alpha_tf".txt" > "./temp/net_tf_final.txt"
sed  "s/}//g" "./temp/test.edgelist_mirna"$alpha_mirna".txt"  > "./temp/net_mirna_final.txt"
rm ./temp/test.edgelist*.txt
printf '%s\n' "../temp/net_exp_final.txt" "../temp/net_tf_final.txt" "../temp/net_mirna_final.txt" "../temp/net_ppi_unique_nodes_int.txt" >> ./temp/slice.txt
#consensus clustering
echo -n "Community detection... "
echo -n "Enter the enter the algorithm you want to use for network clustering (0: oslom undirected, 2: infomap undirected, 4: louvain, 5: label propagation method, 8: modularity optimization)> "
read alg
cd clustering_programs_5_2
python select.py -slices ../temp/slice.txt -p $alg -f ../temp/slice_out -c 100
cp ../temp/slice_out/results_consensus/tp ../temp/
cd ..
awk 'NR%2==0' temp/tp > temp/comm.txt
awk 'NF>=4' temp/comm.txt > temp/comm_filt.txt
sed -i 's/ /	/g' temp/comm_filt.txt
Rscript ./code/conversion_communities.R
rm -r temp
echo -n "Community detection finsihed, look at the output: comm_final_genesymbol.txt"
