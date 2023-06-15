#!/bin/bash
#PBS -N SparCC_000
#PBS -V
#PBS -l nodes=1:ppn=48
#PBS -q exec_8_600

cd /mnt/d/documentos/leonardo/01_universidade/01_biome_lab/1_paper_production/keystone-analysis/Python/keystone-fastspar-correlation-master


./run_fastspar.sh ls &> outs/phyla_000.out
./run_fastspar.sh -Sr &> outs/phyla_001.out
./run_fastspar.sh community_matrix/agricultural_soil.csv &> outs/phyla_002.out
./run_fastspar.sh community_matrix/alkaline_environment.csv &> outs/phyla_003.out
./run_fastspar.sh community_matrix/anaerobic_sediment.csv &> outs/phyla_004.out
./run_fastspar.sh community_matrix/animal_feces.csv &> outs/phyla_005.out
./run_fastspar.sh community_matrix/animal_gut.csv &> outs/phyla_006.out
./run_fastspar.sh community_matrix/animal_manure.csv &> outs/phyla_007.out
./run_fastspar.sh community_matrix/beef.csv &> outs/phyla_008.out
./run_fastspar.sh community_matrix/bodily_fluid.csv &> outs/phyla_009.out
./run_fastspar.sh community_matrix/bone.csv &> outs/phyla_010.out
./run_fastspar.sh community_matrix/chyme.csv &> outs/phyla_011.out
./run_fastspar.sh community_matrix/coastal_marine_sediment.csv &> outs/phyla_012.out
./run_fastspar.sh community_matrix/coastal_seawater.csv &> outs/phyla_013.out
./run_fastspar.sh community_matrix/compost.csv &> outs/phyla_014.out
./run_fastspar.sh community_matrix/contaminated_soil.csv &> outs/phyla_015.out
./run_fastspar.sh community_matrix/coral.csv &> outs/phyla_016.out
./run_fastspar.sh community_matrix/coral_reef_biofilm.csv &> outs/phyla_017.out
./run_fastspar.sh community_matrix/coral_reef_seawater.csv &> outs/phyla_018.out
./run_fastspar.sh community_matrix/desert_soil.csv &> outs/phyla_019.out
./run_fastspar.sh community_matrix/estuarine_seawater.csv &> outs/phyla_020.out
./run_fastspar.sh community_matrix/farm_soil.csv &> outs/phyla_021.out
./run_fastspar.sh community_matrix/forest_soil.csv &> outs/phyla_022.out
./run_fastspar.sh community_matrix/freshwater_biofilm.csv &> outs/phyla_023.out
./run_fastspar.sh community_matrix/garden_soil.csv &> outs/phyla_024.out
./run_fastspar.sh community_matrix/geyser.csv &> outs/phyla_025.out
./run_fastspar.sh community_matrix/grassland.csv &> outs/phyla_026.out
./run_fastspar.sh community_matrix/human-gut.csv &> outs/phyla_027.out
./run_fastspar.sh community_matrix/human-skin.csv &> outs/phyla_028.out
./run_fastspar.sh community_matrix/human_associated_biofilm.csv &> outs/phyla_029.out
./run_fastspar.sh community_matrix/human_feces.csv &> outs/phyla_030.out
./run_fastspar.sh community_matrix/hydrothermal_vent.csv &> outs/phyla_031.out
./run_fastspar.sh community_matrix/hypersaline_water.csv &> outs/phyla_032.out
./run_fastspar.sh community_matrix/karst-porous.csv &> outs/phyla_033.out
./run_fastspar.sh community_matrix/lake.csv &> outs/phyla_034.out
./run_fastspar.sh community_matrix/lung.csv &> outs/phyla_035.out
./run_fastspar.sh community_matrix/mangrove_sediment.csv &> outs/phyla_036.out
./run_fastspar.sh community_matrix/meadow_soil.csv &> outs/phyla_037.out
./run_fastspar.sh community_matrix/oceanic_seawater.csv &> outs/phyla_038.out
./run_fastspar.sh community_matrix/oil_affected_sediment.csv &> outs/phyla_039.out
./run_fastspar.sh community_matrix/park_soil.csv &> outs/phyla_040.out
./run_fastspar.sh community_matrix/pasture_soil.csv &> outs/phyla_041.out
./run_fastspar.sh community_matrix/peat_soil.csv &> outs/phyla_042.out
./run_fastspar.sh community_matrix/permafrost.csv &> outs/phyla_043.out
./run_fastspar.sh community_matrix/plant-associated.csv &> outs/phyla_044.out
./run_fastspar.sh community_matrix/plant-associated_hypersaline_water.csv &> outs/phyla_045.out
./run_fastspar.sh community_matrix/porous_contaminated.csv &> outs/phyla_046.out
./run_fastspar.sh community_matrix/prairie.csv &> outs/phyla_047.out
./run_fastspar.sh community_matrix/rhizosphere.csv &> outs/phyla_048.out
./run_fastspar.sh community_matrix/river.csv &> outs/phyla_049.out
./run_fastspar.sh community_matrix/river_sediment.csv &> outs/phyla_050.out
./run_fastspar.sh community_matrix/saline_evaporation_pond.csv &> outs/phyla_051.out
./run_fastspar.sh community_matrix/saline_marsh.csv &> outs/phyla_052.out
./run_fastspar.sh community_matrix/saline_marsh_sediment.csv &> outs/phyla_053.out
./run_fastspar.sh community_matrix/saliva.csv &> outs/phyla_054.out
./run_fastspar.sh community_matrix/savanna_soil.csv &> outs/phyla_055.out
./run_fastspar.sh community_matrix/sewage.csv &> outs/phyla_056.out
./run_fastspar.sh community_matrix/shrubland_soil.csv &> outs/phyla_057.out
./run_fastspar.sh community_matrix/sludge.csv &> outs/phyla_058.out
./run_fastspar.sh community_matrix/sponge.csv &> outs/phyla_059.out
./run_fastspar.sh community_matrix/stream_biofilm.csv &> outs/phyla_060.out
./run_fastspar.sh community_matrix/wastewater_treatment_plant.csv &> outs/phyla_061.out
./run_fastspar.sh community_matrix/wood_fall.csv &> outs/phyla_062.out
wait
