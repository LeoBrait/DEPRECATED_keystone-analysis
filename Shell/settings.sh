#!/bin/bash

# Analysis title
analysis_frame="phyla_august21"

# Package manager path
package_manager="miniconda3"
source ~/$package_manager/etc/profile.d/conda.sh

# Inputs
annotated_table_absolute="kraken_biomedb_absolute_phyla.csv"
annotated_table_relative="kraken_biomedb_relative_phyla.csv"
metadata_table="biome_classification.csv"

# Preprocessing inputs
multiplicative_const=1
minimum_samples=5

# Computational resources
parallel=25


############################# Performance analysis #############################
iterations=(300 400 500 600 700 800 900 1000
            1000 2000 3000 4000 5000 6000 7000
            8000 9000 10000)

seeds=(1 2 3 4 5 6 7 8 9 10
11 12 13 14 15 16 17 18 19 20
21 22 23 24 25 26 27 28 29 30
31 32 33 34 35 36 37 38 39 40
41 42 43 44 45 46 47 48 49 50)

# Test subsets for performance analysis
communities_path="data/${analysis_frame}/community_subsets"
iterations_test_subsets=( 
    "${communities_path}/animal_host-associated.chyme.tsv"          #N=11
    "${communities_path}/animal_host-associated.animal_feces.tsv"   #N=342
    "${communities_path}/human_host-associated.human-gut.tsv"       #N=16
    "${communities_path}/saline_water.coastal_seawater.tsv"         #N=151
    "${communities_path}/saline_water.hypersaline_water.tsv"        #N=8
    "${communities_path}/soil.savanna_soil.tsv"                     #N=18
    "${communities_path}/freshwater.freshwater_biofilm.tsv"         #N=39
    "${communities_path}/groundwater.porous_contaminated.tsv"       #N=50
    "${communities_path}/sediment.oil_affected_sediment.tsv")       #N=16

############################ Final Fastspar ####################################
definitive_iter=10000
remove=10

############################ Bootstrapping settings ############################
synthetic_communities=1000


############################# Taxa to exclude ##################################
#This is still nos implemented, but planned to be a dictionary given to the
#python preprocessing script.



# Taxa to exclude
# Settings for taxa to drop based on habitats
# declare -A habitats_to_drop_taxa
# habitats_to_drop_taxa=(["estuarine_seawater"]="Myxococcota")
# habitats_to_drop_taxa["salt_lake"]="Candidatus Azambacteria"
# habitats_to_drop_taxa["salt_lake"]="Thermomicrobiota"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Liptonbacteria"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Tagabacteria"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Terrybacteria"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Spechtbacteria"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Microgenomates"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Wallbacteria"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Coatesbacteria"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Blackallbacteria"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Calescamantes"
# habitats_to_drop_taxa["salt_lake"]="Chrysiogenota"
# habitats_to_drop_taxa["salt_lake"]="Coprothermobacterota"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Diapherotrites"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Aenigmarchaeota"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Thermoplasmatota"
# habitats_to_drop_taxa["salt_lake"]="Candidatus Thorarchaeota"
