# Archetype_SOM

A script for illustrating the analysis in doi TBC. Creates archetypes of landscapes, farmed landscapes and farm management systems in Great Britain.

Script name: Self-organising Maps for the production of landscape, farmed landscape or farm management system archetypes 

Purpose of script: To decide on the formulation of the SOM grid and to run the 1000 iterations of the SOM on Tier 1 inputs. Then Sort the results of SOM iterations to derive outputs for archetypes

If you use the information/code in this script, please reference the ERL paper.

Author: Cecily Goodwin

Date Created: 2022-06-22

Copyright (c) Cecily Goodwin, 2021

---------------------------

Notes: An example of the analysis for the run of SOMs and the production of SOM results and maps

Input data cannot be made available due to Third Party licensing constraints, so scripts are provided for an illustration of the analytical process. 
Description of the input data required to replicate these analyses can be found at Goodwin et al. 2022, doi TBC. and EIDC data collection (https://doi.org/10.5285/3b44375a-cbe6-468c-9395-41471054d0f3)
   
- Produces and stores a csv file storing all results of the clustering analysis: ClusterNumber_Tier1s_all_indicies.csv
- Produces and stores 10 R data files of 100 SOM iterations each: SOMs 100 - 1000.Rdata

Produces the following Openly accessible outputs, which can be found at EIDC (https://doi.org/10.5285/3b44375a-cbe6-468c-9395-41471054d0f3):
      - 1 raster object: Tier1.tif: the assignment of each 1km square in GB to an archetype (Band 1) and the distance of each 1km square from its assigned archetype (Band 2)
      - 1 csv files: Tier1_code_table.csv: the average codebook vectors of each archetype based on the iterations which produce the most consistent set of archetypes

Runs a naming procedure at the end to ascribe descriptors to each archetype based on the variables they are associated with.
