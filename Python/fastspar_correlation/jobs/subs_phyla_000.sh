#!/bin/bash
#PBS -N SparCC_000
#PBS -V
#PBS -l nodes=1:ppn=48
#PBS -q exec_8_600



bash run_fastspar.sh community_matrix/agricultural_soil.phyla.csv &> outs/phyla_002.out
bash run_fastspar.sh community_matrix/alkaline_environment.phyla.csv &> outs/phyla_003.out
bash run_fastspar.sh community_matrix/anaerobic_sediment.phyla.csv &> outs/phyla_004.out
bash run_fastspar.sh community_matrix/animal_feces.phyla.csv &> outs/phyla_005.out
bash run_fastspar.sh community_matrix/animal_gut.phyla.csv &> outs/phyla_006.out
bash run_fastspar.sh community_matrix/animal_manure.phyla.csv &> outs/phyla_007.out
bash run_fastspar.sh community_matrix/beef.phyla.csv &> outs/phyla_008.out
bash run_fastspar.sh community_matrix/bodily_fluid.phyla.csv &> outs/phyla_009.out
bash run_fastspar.sh community_matrix/bone.phyla.csv &> outs/phyla_010.out
bash run_fastspar.sh community_matrix/chyme.phyla.csv &> outs/phyla_011.out
bash run_fastspar.sh community_matrix/coastal_marine_sediment.phyla.csv &> outs/phyla_012.out
bash run_fastspar.sh community_matrix/coastal_seawater.phyla.csv &> outs/phyla_013.out
bash run_fastspar.sh community_matrix/compost.phyla.csv &> outs/phyla_014.out
bash run_fastspar.sh community_matrix/contaminated_soil.phyla.csv &> outs/phyla_015.out
bash run_fastspar.sh community_matrix/coral.phyla.csv &> outs/phyla_016.out
bash run_fastspar.sh community_matrix/coral_reef_biofilm.phyla.csv &> outs/phyla_017.out
bash run_fastspar.sh community_matrix/coral_reef_seawater.phyla.csv &> outs/phyla_018.out
bash run_fastspar.sh community_matrix/desert_soil.phyla.csv &> outs/phyla_019.out
bash run_fastspar.sh community_matrix/estuarine_seawater.phyla.csv &> outs/phyla_020.out
bash run_fastspar.sh community_matrix/farm_soil.phyla.csv &> outs/phyla_021.out
bash run_fastspar.sh community_matrix/forest_soil.phyla.csv &> outs/phyla_022.out
bash run_fastspar.sh community_matrix/freshwater_biofilm.phyla.csv &> outs/phyla_023.out
bash run_fastspar.sh community_matrix/garden_soil.phyla.csv &> outs/phyla_024.out
bash run_fastspar.sh community_matrix/geyser.phyla.csv &> outs/phyla_025.out
bash run_fastspar.sh community_matrix/grassland.phyla.csv &> outs/phyla_026.out
bash run_fastspar.sh community_matrix/human-gut.phyla.csv &> outs/phyla_027.out
bash run_fastspar.sh community_matrix/human-skin.phyla.csv &> outs/phyla_028.out
bash run_fastspar.sh community_matrix/human_associated_biofilm.phyla.csv &> outs/phyla_029.out
bash run_fastspar.sh community_matrix/human_feces.phyla.csv &> outs/phyla_030.out
bash run_fastspar.sh community_matrix/hydrothermal_vent.phyla.csv &> outs/phyla_031.out
bash run_fastspar.sh community_matrix/hypersaline_water.phyla.csv &> outs/phyla_032.out
bash run_fastspar.sh community_matrix/karst-porous.phyla.csv &> outs/phyla_033.out
bash run_fastspar.sh community_matrix/lake.phyla.csv &> outs/phyla_034.out
bash run_fastspar.sh community_matrix/lung.phyla.csv &> outs/phyla_035.out
bash run_fastspar.sh community_matrix/mangrove_sediment.phyla.csv &> outs/phyla_036.out
bash run_fastspar.sh community_matrix/meadow_soil.phyla.csv &> outs/phyla_037.out
bash run_fastspar.sh community_matrix/oceanic_seawater.phyla.csv &> outs/phyla_038.out
bash run_fastspar.sh community_matrix/oil_affected_sediment.phyla.csv &> outs/phyla_039.out
bash run_fastspar.sh community_matrix/park_soil.phyla.csv &> outs/phyla_040.out
bash run_fastspar.sh community_matrix/pasture_soil.phyla.csv &> outs/phyla_041.out
bash run_fastspar.sh community_matrix/peat_soil.phyla.csv &> outs/phyla_042.out
bash run_fastspar.sh community_matrix/permafrost.phyla.csv &> outs/phyla_043.out
bash run_fastspar.sh community_matrix/plant-associated.phyla.csv &> outs/phyla_044.out
bash run_fastspar.sh community_matrix/plant-associated_hypersaline_water.phyla.csv &> outs/phyla_045.out
bash run_fastspar.sh community_matrix/porous_contaminated.phyla.csv &> outs/phyla_046.out
bash run_fastspar.sh community_matrix/prairie.phyla.csv &> outs/phyla_047.out
bash run_fastspar.sh community_matrix/rhizosphere.phyla.csv &> outs/phyla_048.out
bash run_fastspar.sh community_matrix/river.phyla.csv &> outs/phyla_049.out
bash run_fastspar.sh community_matrix/river_sediment.phyla.csv &> outs/phyla_050.out
bash run_fastspar.sh community_matrix/saline_evaporation_pond.phyla.csv &> outs/phyla_051.out
bash run_fastspar.sh community_matrix/saline_marsh.phyla.csv &> outs/phyla_052.out
bash run_fastspar.sh community_matrix/saline_marsh_sediment.phyla.csv &> outs/phyla_053.out
bash run_fastspar.sh community_matrix/saliva.phyla.csv &> outs/phyla_054.out
bash run_fastspar.sh community_matrix/savanna_soil.phyla.csv &> outs/phyla_055.out
bash run_fastspar.sh community_matrix/sewage.phyla.csv &> outs/phyla_056.out
bash run_fastspar.sh community_matrix/shrubland_soil.phyla.csv &> outs/phyla_057.out
bash run_fastspar.sh community_matrix/sludge.phyla.csv &> outs/phyla_058.out
bash run_fastspar.sh community_matrix/sponge.phyla.csv &> outs/phyla_059.out
bash run_fastspar.sh community_matrix/stream_biofilm.phyla.csv &> outs/phyla_060.out
bash run_fastspar.sh community_matrix/wastewater_treatment_plant.phyla.csv &> outs/phyla_061.out
bash run_fastspar.sh community_matrix/wood_fall.phyla.csv &> outs/phyla_062.out
wait
