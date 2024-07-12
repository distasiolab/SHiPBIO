declare -a arr=("","","")

for ((i=0; i<${#arr[@]}; i++)); do python ReadAnnotationsToAnnData.py -f "${arr[$i]}" -c "${arr[$i]}".FullSizeImage_LeidenClusters.tif_Annotations.json ; done
