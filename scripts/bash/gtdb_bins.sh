

module load GTDBTk/1.5.0-foss-2020b-Python-3.8.6

gtdbtk classify_wf \
  --genome_dir /user_data/sh/differential_digestion/mmlong2/results/bins \
  --out_dir /user_data/sh/differential_digestion/mmlong2/results/gtdbtk \
  --extension fa \
  --cpus 50

module purge