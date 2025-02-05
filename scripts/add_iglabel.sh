## add new iglabels no novel alleles

## step 1 copy the novel allele file new_alleles_without_iglabel_*Date*

## change the file name in process_R24_novels.py

## cd into gldb-macaque/macaca_mulatta

cd gldb-macaque/macaca_mulatta

## step 2 run process_R24_novels.py
python ../python/process_R24_novels.py 'bosinger_watson_R24' 'MUSA_set_alleles_without_iglabel_bosinger_watson_2025-02-04.csv'

## cd into IGH/bosinger_watson_R24 to process the new sequences

cd IGH/bosinger_watson_R24

## step 3 run IgLabel query

python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igh_with_trio.csv novels.fasta results.csv actions.csv

## step 4 change the action to new_label to anything that is not 'exact' in reason. delete any that are 'exact'
## save the file as actions_as_processed_1.csv

echo "$(grep -E '("?exact"?)' actions.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions.csv | wc -l)"

grep -v -E ',("?exact"?)$' actions.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_1.csv

## step 5 add the new sequences to the database

python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igh_with_trio.csv actions_as_processed_1.csv "Bosinger_Watson R24 after bulk processing"

## step 6 repeat the query on the postponed sequences

python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igh_with_trio.csv novels_postponed.fasta results.csv actions_as_processed_2.csv

echo "$(grep -E '("?exact"?)' actions_as_processed_2.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions_as_processed_2.csv | wc -l)"

grep -v -E '("?exact"?)' actions_as_processed_2.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_3.csv

## step 7 repeat steps 4 and 5 and add the sequences

python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igh_with_trio.csv actions_as_processed_3.csv "Bosinger_Watson R24 after bulk processing"

## cd to gldb-macaque/macaca_mulatta

cd ../../

## step 8 check that all novel alleles were added

python ../python/process_R24_novels.py 'bosinger_watson_R24' 'MUSA_set_alleles_without_iglabel_bosinger_watson_2025-02-04.csv'

## repeat the process for all loci

# IGK

cd IGK/bosinger_watson_R24
python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igk_with_trio.csv novels.fasta results.csv actions.csv
echo "$(grep -E '("?exact"?)' actions.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions.csv | wc -l)"
grep -v -E '("?exact"?)' actions.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_1.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igk_with_trio.csv actions_as_processed_1.csv "Bosinger_Watson R24 after bulk processing"
python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igk_with_trio.csv novels_postponed.fasta results.csv actions_as_processed_2.csv
echo "$(grep -E '("?exact"?)' actions_as_processed_2.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions_as_processed_2.csv | wc -l)"
grep -v -E '("?exact"?)' actions_as_processed_2.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_3.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igk_with_trio.csv actions_as_processed_3.csv "Bosinger_Watson R24 after bulk processing"
cd ../../
python ../python/process_R24_novels.py 'bosinger_watson_R24' 'MUSA_set_alleles_without_iglabel_bosinger_watson_2025-02-04.csv'

# IGL

cd IGL/bosinger_watson_R24
python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igl_with_trio.csv novels.fasta results.csv actions.csv
echo "$(grep -E '("?exact"?)' actions.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions.csv | wc -l)"
grep -v -E '("?exact"?)' actions.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_1.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igl_with_trio.csv actions_as_processed_1.csv "Bosinger_Watson R24 after bulk processing"
python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igl_with_trio.csv novels_postponed.fasta results.csv actions_as_processed_2.csv
echo "$(grep -E '("?exact"?)' actions_as_processed_2.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions_as_processed_2.csv | wc -l)"
grep -v -E '("?exact"?)' actions_as_processed_2.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_3.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igl_with_trio.csv actions_as_processed_3.csv "Bosinger_Watson R24 after bulk processing"
cd ../../
python ../python/process_R24_novels.py 'bosinger_watson_R24' 'MUSA_set_alleles_without_iglabel_bosinger_watson_2025-02-04.csv'



### do the same for imgt, vrc, and guo

## imgt

#cd gldb-macaque/macaca_mulatta
python ../python/process_R24_novels.py 'imgt_dec_28' 'MUSA_set_alleles_without_iglabel_imgt_2025-02-04.csv'
cd IGH/imgt_dec_28
python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igh_with_trio.csv novels.fasta results.csv actions.csv
echo "$(grep -E '("?exact"?)' actions.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions.csv | wc -l)"
grep -v -E '("?exact"?)' actions.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_1.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igh_with_trio.csv actions_as_processed_1.csv "imgt"
cd ../../
python ../python/process_R24_novels.py 'imgt_dec_28' 'MUSA_set_alleles_without_iglabel_imgt_2025-02-04.csv'
cd IGK/imgt_dec_28
python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igk_with_trio.csv novels.fasta results.csv actions.csv
echo "$(grep -E '("?exact"?)' actions.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions.csv | wc -l)"
grep -v -E '("?exact"?)' actions.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_1.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igk_with_trio.csv actions_as_processed_1.csv "imgt"
cd ../../
python ../python/process_R24_novels.py 'imgt_dec_28' 'MUSA_set_alleles_without_iglabel_imgt_2025-02-04.csv'
cd IGL/imgt_dec_28
python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igl_with_trio.csv novels.fasta results.csv actions.csv
echo "$(grep -E '("?exact"?)' actions.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions.csv | wc -l)"
grep -v -E '("?exact"?)' actions.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_1.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igl_with_trio.csv actions_as_processed_1.csv "imgt"
cd ../../
python ../python/process_R24_novels.py 'imgt_dec_28' 'MUSA_set_alleles_without_iglabel_imgt_2025-02-04.csv'

## vrc
python ../python/process_R24_novels.py 'vrc' 'MUSA_set_alleles_without_iglabel_vrc_2025-02-04.csv'
cd IGH/vrc
python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igh_with_trio.csv novels.fasta results.csv actions.csv
echo "$(grep -E '("?exact"?)' actions.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions.csv | wc -l)"
grep -v -E '("?exact"?)' actions.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_1.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igh_with_trio.csv actions_as_processed_1.csv "VRC"
cd ../../
python ../python/process_R24_novels.py 'vrc' 'MUSA_set_alleles_without_iglabel_vrc_2025-02-04.csv'
cd IGK/vrc
python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igk_with_trio.csv novels.fasta results.csv actions.csv
echo "$(grep -E '("?exact"?)' actions.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions.csv | wc -l)"
grep -v -E '("?exact"?)' actions.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_1.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igk_with_trio.csv actions_as_processed_1.csv "VRC"

python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igk_with_trio.csv novels_postponed.fasta results.csv actions_as_processed_2.csv
echo "$(grep -E '("?exact"?)' actions_as_processed_2.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions_as_processed_2.csv | wc -l)"
grep -v -E '("?exact"?)' actions_as_processed_2.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_3.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igk_with_trio.csv actions_as_processed_3.csv "VRC"

cd ../../
python ../python/process_R24_novels.py 'vrc' 'MUSA_set_alleles_without_iglabel_vrc_2025-02-04.csv'
cd IGL/vrc
python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igl_with_trio.csv novels.fasta results.csv actions.csv
echo "$(grep -E '("?exact"?)' actions.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions.csv | wc -l)"
grep -v -E '("?exact"?)' actions.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_1.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igl_with_trio.csv actions_as_processed_1.csv "VRC"
cd ../../
python ../python/process_R24_novels.py 'vrc' 'MUSA_set_alleles_without_iglabel_vrc_2025-02-04.csv'

## guo
python ../python/process_R24_novels.py 'guo' 'MUSA_set_alleles_without_iglabel_guo_2025-02-04.csv'
cd IGH/guo
python ../../../../IgLabel/iglabel.py query ../db/macaca_mulatta_igh_with_trio.csv novels.fasta results.csv actions.csv
echo "$(grep -E '("?exact"?)' actions.csv | wc -l) exact matches, expected remaining lines $(grep -v -E '("?exact"?)' actions.csv | wc -l)"
grep -v -E '("?exact"?)' actions.csv | awk -F',' 'NR==1 {print; next} { $3 = "new_label"; print }' OFS=',' > actions_as_processed_1.csv
python ../../../../IgLabel/iglabel.py add ../db/macaca_mulatta_igh_with_trio.csv actions_as_processed_1.csv "guo"
cd ../../
python ../python/process_R24_novels.py 'guo' 'MUSA_set_alleles_without_iglabel_guo_2025-02-04.csv'
