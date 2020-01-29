cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/"

sp_array=("BELLII" "BILINEATA" "BRUNNEICAPILLUS" "FLAVICEPS" "FUSCA" "CRISSALE" "CURVIROSTRE" "MELANURA" "NITENS" "SINUATUS")
#ch_array=("1A" "1B" "4A" "LG2" "LG5" "LGE22" "mtDNA" "Z")
#ch_array=("mtDNA" "Z")

for chr in {26..28}; do ## done through 7, up to 28 
#for chr in ${ch_array[@]}; do
#echo $chr >> gdm_results_jan2020_full_$spp.txt 2>&1
for spp in ${sp_array[@]}; do
echo $chr >> gdm_results_jan2020_full_${spp}_${chr}.txt 2>&1 
echo $spp >> gdm_results_jan2020_full_${spp}_${chr}.txt 2>&1 
nohup Rscript scripts/run_GDMS.R $spp $chr >> gdm_results_jan2020_full_${spp}_${chr}.txt 2>&1 2>&1 &
done
done

mv *png ./GDM_results/univariate/
cd GDM_results/univariate/; mv *chr* ./CHRS/ONLY_VARIABLE/; mv *pc1* ./PC1/ONLY_VARIABLE/; mv *pc2* ./PC2/ONLY_VARIABLE/; mv *pc3* ./PC3/ONLY_VARIABLE/; mv *morph* ./ALL_MORPH/ONLY_VARIABLE/; mv *gene* ./GENES/ONLY_VARIABLE/; cd ../../;
