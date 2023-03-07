
wd=/user_data/sh/differential_digestion
mmlong2-lite -p 70 \
  -o  $wd/results/mmlong2 \
  -np $wd/data/processed/NP/NP_R1041_400bps_sub.fastq \
  -cov $wd/data/processed/coverage_samples.csv \
  -tmp /user_data/sh/tmp

