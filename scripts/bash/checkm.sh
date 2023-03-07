

module load CheckM/1.1.2-foss-2018a-Python-3.6.4

checkm lineage_wf \
  -t 50 \
  -x fa \
  /user_data/sh/differential_digestion/mmlong2/results/bins \
  /user_data/sh/differential_digestion/mmlong2/results/checkm \
  --tmpdir /user_data/sh/tmp

  module purge