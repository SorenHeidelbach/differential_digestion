WD="/user_data/sh/differential_digestion"
cd $WD

# Subset 10.4.1 data to MboI data bases

rasusa -b 7700MB \
  -i data/processed/NP/NP_R1041_400bps_sub.fastq \
  -o data/processed/NP/NP_R1041_400bps_sub_7.7GB.fastq


