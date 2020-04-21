cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/"

sp_array=("BELLII" "BILINEATA" "BRUNNEICAPILLUS" "FLAVICEPS" "FUSCA" "CRISSALE" "CURVIROSTRE" "MELANURA" "NITENS" "SINUATUS")
#sp_array=("BILINEATA")
#ch_array=("1A" "1B" "4A" "LG2" "LG5" "LGE22" "mtDNA" "Z")

## all species independently
for spp in ${sp_array[@]}; do
nohup Rscript "/Users/kprovost/Documents/Github/GDM_pipeline/scripts/GDM Analyses/run_GDMS.R" $spp >> gdm_results_feb2020_full_${spp}_morphpcs.txt 2>&1 &
done

# for chr in {3..4}; do ## done through 4, up to 28 
# #for chr in ${ch_array[@]}; do
# #echo $chr >> gdm_results_jan2020_full_$spp.txt 2>&1
# for spp in ${sp_array[@]}; do
# echo $chr >> gdm_results_feb2020_full_${spp}_${chr}.txt 2>&1 
# echo $spp >> gdm_results_feb2020_full_${spp}_${chr}.txt 2>&1 
# nohup Rscript "/Users/kprovost/Documents/Github/GDM_pipeline/scripts/GDM Analyses/run_GDMS.R" $spp $chr >> gdm_results_jan2020_full_${spp}_${chr}.txt 2>&1 &
# done
# done

cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/"
mv *png ./GDM_results/
mv *spline* ./GDM_results/
mv *_variable_* ./GDM_results/
mv *modelparam* ./GDM_results/
jobs


#cd GDM_results/univariate/; mv *chr* ./CHRS/ONLY_VARIABLE/; mv *pc1* ./PC1/ONLY_VARIABLE/; mv *pc2* ./PC2/ONLY_VARIABLE/; mv *pc3* ./PC3/ONLY_VARIABLE/; mv *morph* ./ALL_MORPH/ONLY_VARIABLE/; mv *gene* ./GENES/ONLY_VARIABLE/; cd ../../;


## run everything simultaneously
# nohup Rscript "/Users/kprovost/Documents/Github/GDM_pipeline/scripts/GDM Analyses/run_GDMS.R" >> gdm_results_jan2020_full_allspp_allchr.txt 2>&1 &
