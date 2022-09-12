## ---------------------------
##
## Script name: Self-organising Maps for the production of landscape, farmed landscape or farm management system archetypes 
##
## Purpose of script: To decide on the formulation of the SOM grid and to run the 1000 iterations of the SOM on Tier 1 inputs. Then Sort the results of SOM iterations to derive outputs for archetypes
##
## Author: Cecily Goodwin
##
## Date Created: 2022-06-22
##
## Copyright (c) Cecily Goodwin, 2021
##
## ---------------------------
##
## Notes: An example of the analysis for the run of SOMs and the production of SOM results and maps
##
## Input data cannot be made available due to Third Party licensing constraints, so scripts are provided for an illustration of the analytical process. 
## Description of the input data required to replicate these analyses can be found at Goodwin et al. 2022, doi TBC. and EIDC data collection TBC
##   
##  Produces and stores a csv file storing all results of the clustering analysis: ClusterNumber_Tier1s_all_indicies.csv
##  Produces and stores 10 R data files of 100 SOM iterations each: SOMs 100 - 1000.Rdata
##
##  Produces the following Openly accessible outputs, which can be found at: EIDC TBC 
##      1 raster object: Tier1.tif: the assignment of each 1km square in GB to an archetype (Band 1) and the distance of each 1km square from its assigned archetype (Band 2)
##      1 csv files: Tier1_code_table.csv: the average codebook vectors of each archetype based on the iterations which produce the most consistent set of archetypes
##
##  Runs a naming procedure at the end to ascribe descriptors to each archetype based on the variables they are associated with.
##
## ---------------------------

options(scipen = 999, digits = 4) # 
memory.limit(10000000)     # For PC - memory limit increase

## ---------------------------

## load up data via scripts 

# The data needed to run this script is detailed in Goodwin et al. 2022, doi TBC. 
file_path <-  "P:/NEC07065_AGLAND/WP1/CG/" 

source(paste0(file_path, "Scripts/Tier1/Tier1_VariableExtraction (1).R"))

file_path <- "P:/NEC07065_AGLAND/WP1/CG/"

## ---------------------------

## load up packages:  

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)

library(raster)
library(flextable)
library(cluster)
library(clusterSim)
library(clValid)
library(kohonen)

## load up functions and proj4string:

sortpaste <- function(x){
  sorted <- sort(x)
  return(paste0(sorted, collapse="\   "))
}

calc_arch_per_1km <- function(cell_assignment) {
  a <- sort(table(cell_assignment), decreasing = TRUE)
  if(length(a)==1) {
    tc <- names(a)[1]
    ct <- as.numeric(a)[1]/length(cell_assignment) 
    return(data.frame(topC=tc, cert=ct, nextC=NA, cert2=NA))
  } else {
    tc <- names(a)[1]
    ct <- as.numeric(a)[1]/length(cell_assignment)  
    nc <- names(a)[2]
    ct2 <- as.numeric(a)[2]/length(cell_assignment) 
    return(data.frame(topC=tc, cert=ct,  nextC=nc, cert2=ct2))
  }
}

OSGB <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs" 

## ---------------------------

#write.csv(ExVars_L, file=paste0(file_path, "/Archetype_output_data/Tier1_data/Tier1_input.csv"), row.names=FALSE)

############ Clusters ####

# Takes quite a while to run over all grid configurations

# set up dataframe - for cluster combinations, number, distances and indices. 
clus_results <- data.frame(clust_combi = as.character(), clustnum = as.numeric(), clustDist = as.numeric(), 
                           iDB = as.numeric(), iG1 = as.numeric())

# set up range of number of grid dimensions want to look at e.g. 2 to 12 in either dimension
dims <- 2:12
combinations <- expand.grid(dims,dims)
combinations$c_no <- combinations$Var1*combinations$Var2
combo_runs <- combinations[combinations$c_no>8 & combinations$c_no<20,] # Run all cluster numbers for a subset if desired e.g. a total of 8 to 20 clusters

# Just data to go in without coordinates
Input_drivers <- ExVars_L[,-which(names(ExVars_L) %in% c('POINT_X', 'POINT_Y', 'SqID'))]

# Run over each cluster grid configuration
for (i in 1:length(combo_runs$Var1)){
  
  distancesmean <- c() 
  
  for (j in 1:30){
    current <-supersom(list(scale(Input_drivers)),
                       rlen = 500, mode="pbatch",
                       cores = 3, normalizeDataLayers = TRUE, grid = somgrid(combo_runs$Var1[i], combo_runs$Var2[i], "hexagonal"),maxNA.fraction = 0.99)
    distancesmean[j] <- mean(current$distances)
    print(j)
  }
  
  curNum <- data.frame(clust_combi = paste0(combo_runs$Var1[i], '_', combo_runs$Var2[i]), 
                       clustnum = rep((combo_runs$Var1[i] * combo_runs$Var2[i]),length(distancesmean)), 
                       clustDist = distancesmean)
  
  clus_results <- rbind(clus_results, curNum)
  print(paste0(i, "/", length(combo_runs$Var1)))
  print(as.numeric(combo_runs[i,]))
}

plot(current, type="dist.neighbours")
plot(current, type="quality")

#write.csv(clus_results, file="C:/Users/cecgoo/OneDrive - UKCEH/Archetype_output_data/Tier1_data/Cluster_results/ClusterNumber_Tier1s_all_indicies.csv")

# Load in data and plot

#clus_results <- read.csv(paste0(file_path, "Archetype_output_data/Tier1_data/Cluster_results/ClusterNumber_Tier1s_all_indicies.csv"))

ST <- 'Tier1 clusters'

clus_results$clust_combi <- as.factor(clus_results$clust_combi)

Clus_num <- data.frame(clustnum=seq(min(clus_results$clustnum),max(clus_results$clustnum),1)) # dummy vector of cluster number so values can plot in sequence
CN <- full_join(clus_results, Clus_num, by='clustnum') # join to all numbers 

# plot distances of each cluster configuration
png(paste0(file_path, "Archetype_output_data/Tier1_data/Cluster_results/Cluster_config_all.png"), width=20, height =10, units='in', res=300)
clus_results %>%
  mutate(clus_CN = fct_reorder(clust_combi, clustnum)) %>%
  ggplot(aes(x=as.factor(clus_CN), y=clustDist)) +
  ggtitle(ST) +
  xlab('\nCluster configuration (ordered by cluster number)') +
  ylab('Mean within cluster distance\n') +
  geom_boxplot() +
  stat_summary(fun=median, geom="line", aes(group=1)) +
  stat_summary(fun=median, geom="point") +
  theme(axis.text.x = element_text(angle=90)) 
dev.off()

# plot distances of each cluster number (grouped configuration)
CN %>%
  ggplot( aes(x=as.factor(clustnum), y=clustDist)) +
  ggtitle(ST) +
  xlab('\nCluster number') +
  ylab('Mean within cluster distance\n') +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90)) 

# filter out the best grid configurations for each cluster number
CN %>%
  group_by(clust_combi, clustnum) %>%
  summarise(mcd = round(median(clustDist),2), mcd_cons = max(table(round(clustDist,2)))) %>% # work out the distance and distance consistency of each grid formation 
  group_by(clustnum) %>%
  mutate(min_mcd=min(mcd), cons_mcd=max(mcd_cons)) %>% # work out the minimum distance of each cluster number
  ungroup() %>%
  filter(mcd==min_mcd) %>% # select out the grid formation for each cluster number with the minimum distance
  full_join(., Clus_num) %>% # join back in all possible numbers so has gaps in sequence on plot
  arrange(clustnum) %>% # order by cluster number
  data.frame() -> min_mcd_A 

min_mcd_A %>%
  group_by(clustnum) %>%
  filter(mcd_cons==max(mcd_cons)) %>% # if there are multiple groups with the same median distance. Select out the one with most consistent runs. 
  ungroup() %>%
  full_join(., Clus_num) %>%
  arrange(clustnum) %>%
  data.frame() -> min_mcd

# select out best grid combos for each cluster number and plot distances
png(paste0(file_path, "Archetype_output_data/Tier1_data/Cluster_results/Cluster_config_best.png"), width=20, height =10, units='in', res=300)
CN %>%
  filter(clust_combi %in% min_mcd$clust_combi) %>%
  ggplot(aes(x=as.factor(clustnum), y=clustDist)) +
  ggtitle(paste(ST, 'min_MCD')) +
  xlab('\nCluster number') +
  ylab('distance\n') +
  geom_boxplot() +
  stat_summary(fun=median, geom="line", aes(group=1)) +
  stat_summary(fun=median, geom="point") +
  theme(axis.text.x = element_text(angle =90)) 
dev.off()

############ Run 1000 SOM iterations ####

for (i in 1:10){
  results <- list()
  
  for (j in 1:100){
    current <-supersom(list(scale(Input_drivers)),
                       rlen = 500, mode="pbatch",
                       cores = 3, normalizeDataLayers = TRUE, grid = somgrid(8,2, "hexagonal"),maxNA.fraction = 0.99) #, keep.data=FALSE
    results[[j]] <- current
    print(j)
  }
  
  #save(results,file=paste0(file_path, "R_files/Tier1s/SOMs_",(i*100),".rdata"))
  rm(results)
  gc()
}


# Explanation of SOM elements ##

str(results[[1]]) # each run has 13 objects. 
results[[1]]$data
# data = values for each grid cell. One object in list - [[1]]. Values for each var (cols) for each grid cell (rows)

results[[1]]$codes
# 16 x 23 matrix of variable values (cols) for each cluster (rows)

results[[1]]$distances
# distances of each cell from assigned cluster

results[[1]]$unit.classif
# assigned cluster for each cell

############ Load up Archetype colour pallettes ####

grey_pal <- colorRampPalette(c('lightgrey', 'darkgrey'))
Tier1_cols <- data.frame(colrs=c('yellow', 'orange', 'red',  
                               'darkgreen', 
                               'darkblue', 
                               'lightsalmon3','brown', 'beige', 
                               'darkgrey', 'grey',
                               'lightblue',
                               'lightgreen', 'yellowgreen', 'olivedrab4', 
                               'seagreen', 'darkorange4'), Arch_no=c(16,15,6,13,11,3,10,9,8,7,12,5,14,4,2,1))

Tier1_cols <- Tier1_cols[order(Tier1_cols$Arch_no),]

Tier1_pal <- colorRampPalette(Tier1_cols$colrs)

############ Set up data and load in SOM Rdata files ####

no_nodes <- 16 # Number of nodes in SOM
no_iters <- 100 # Number of iterations within each repblock

repblocks <- 10 # Number of repblocks

Input_drivers <- ExVars_L[,-which(names(ExVars_L) %in% c('POINT_X', 'POINT_Y', 'SqID'))] # SOM input data
Input_drivers_names <- c(names(Input_drivers))

Assigned_cluster <- matrix(NA, nrow=length(Input_drivers[,1]), ncol=(repblocks*no_iters)) # empty dataframe to fill with result of each iterations for each cell
Cluster_dist <- matrix(NA, nrow=length(Input_drivers[,1]), ncol=(repblocks*no_iters)) # same for distances
Cluster_var_load <- Input_drivers[0,] # empty codebook data frame 

for (rb in seq(no_iters, repblocks*no_iters, no_iters)) {
  
  load(paste0(file_path, "R_files/Tier1s/SOMs_",rb,".rdata")) # load each block in turn - created above under 'Run 1000 SOM iterations'
  
  for (r in 1:no_iters){
    
    curcodes <- data.frame(results[[r]]$codes) # get the variable values for each cluster centre. One value per input variable
    iter <- (rb - 100)+r # record the iteration 1:1000.
    
    curcodes$iter <- iter
    curcodes$dist <- tapply(results[[r]]$distances, results[[r]]$unit.classif, mean) # also record the mean distance to assigned centre for each archetype
    curcodes$meandist <- mean(results[[r]]$distances) # mean distances of all archetypes in this run
    
    curassigned <- results[[r]]$unit.classif # which node each grid cell is assigned to
    curdist <- results[[r]]$distances # distance of each cell to its assigned node
    
    for (i in 1:nrow(curcodes)){ # label each node with it's associated variables
      curcodes_vect <- curcodes[i,which(names(curcodes) %in% Input_drivers_names)]
      over0.6 <- paste(sort(names(curcodes_vect[which(round(curcodes_vect,1)>0.6)])), collapse = "_")
      under0.6 <- paste(sort(names(curcodes_vect[which(round(curcodes_vect,1)<(-0.6))])), collapse = "_")
      curcodes$thres_0.6_name[i] <- paste(over0.6, under0.6, collapse=" ")
    }
    
    # A number vector in row order to match up to labelled names
    curcodes$cnum <- 1:length(unique(curassigned))
    
    # Relabel the node numbers with their corresponding variable names
    clusterassigned <- recode(curassigned, !!! setNames(curcodes$thres_0.6_name, curcodes$cnum)) 
    
    Cluster_var_load <- rbind(Cluster_var_load, curcodes) # bind this iterations codebook to the end of the codebook data
    
    Assigned_cluster[,iter]<- clusterassigned # node assigned to each km relabelled with it's associated variables.
    Cluster_dist[,iter]<- curdist 
  }
  print(rb)
  rm(results)
}

############ Hierachical clustering of all nodes across all iterations to remove non-consistent iterations ####

## Hierachical clustering of all values of all codesbooks across iterations
Cluster_var_load %>%
  dplyr::select(names(Input_drivers)) -> load_data

ds <- dist(load_data, method='euclidean')

HC <- hclust(ds, method='ward.D2') 

groups <- cutree(HC, k=no_nodes) # cut the tree into the number of node groups

# Plot to see distance - how defined are groups?
plot(HC)
rect.hclust(HC, k=no_nodes)

# match the codebook row to the iteration and its cluster grouping
group_iters <- left_join(data.frame(row=names(groups), group=groups),
                         data.frame(row=row.names(Cluster_var_load), iter=Cluster_var_load$iter), by='row')

# calculate the frequency with which each iteration appears in each group
group_iters_freq <- data.frame(table(group_iters$iter, group_iters$group))

# remove the iterations which don't appear once in each group
other_iters <- unique(group_iters_freq$Var1[group_iters_freq$Freq!=1])

main_iters_index <- unique(group_iters$iter[!group_iters$iter %in% other_iters])

## Hierachical clustering of all codesbooks of uniquely named clusters - to see where differences are
Cluster_var_load %>%
  group_by(thres_0.6_name) %>%
  dplyr::select(names(Input_drivers)) %>%
  summarise(count=n(), across(.cols=everything(), .fns=mean)) -> LDU

LDU %>%
  dplyr::select(!c('count', 'thres_0.6_name')) -> load_data_uni

# Hierachical clustering
dsu <- dist(load_data_uni, method='euclidean')

HCU <- hclust(dsu, method='ward.D2') 

# Plot with 16 groups
plot(HCU, labels=paste(unique(LDU$thres_0.6_name), LDU$count), hang=-1, cex=0.7)
rect.hclust(HCU, k=no_nodes)

############ Selecting out top cluster set and recode number ####

Cluster_var_load_main <- Cluster_var_load[which(Cluster_var_load$iter %in% main_iters_index),] 

# label each archetype with its cluster group number 
Cluster_var_load_main <- cbind(Cluster_var_load_main, Arch_no=groups[match(rownames(Cluster_var_load_main), names(groups))])

# Extract out the iterations for assigning grid cells 
Assigned_cluster_main <- Assigned_cluster[,main_iters_index]
Cluster_dist_main <- Cluster_dist[,main_iters_index]

# recode the assigned cluster with the group
name_an <- data.frame(table(Cluster_var_load_main$thres_0.6_name, Cluster_var_load_main$Arch_no))
names(name_an) <- c('thres_0.6_name', 'Arch_no', 'Freq')
name_an$thres_0.6_name <- as.character(name_an$thres_0.6_name)
name_an$Arch_no <- as.character(name_an$Arch_no)
name_arch_no_combos <- name_an[name_an$Freq!=0,]

Assigned_cluster_main <- apply(Assigned_cluster_main, 2, function(x) recode(x, !!! setNames(name_arch_no_combos$Arch_no, name_arch_no_combos$thres_0.6_name)))

## Hierachical clustering of all codesbooks of consistent iterations
Cluster_var_load_main %>%
  group_by(thres_0.6_name) %>%
  dplyr::select(names(Input_drivers)) %>%
  summarise(count=n(), across(.cols=everything(), .fns=mean)) -> LDU_m

LDU_m %>%
  dplyr::select(!c('thres_0.6_name', 'count')) -> load_data_uni

# Hierachical clustering
dsum <- dist(load_data_uni, method='euclidean')

HCUm <- hclust(dsum, method='ward.D2') 

# Plot 
plot(HCUm, labels=paste(unique(LDU_m$thres_0.6_name), LDU_m$count), hang=-1, cex=0.7)
rect.hclust(HCUm, k=no_nodes)

############ Calculate the certainty of classification into each cluster and distance for top cluster set ####

# a) ALL ITERATIONS calculate the mean distance of each cell from its top and second assigned cluster
tabs_all <- apply(Assigned_cluster, 1, calc_arch_per_1km) # topC = the archetype label applied most to each cell. nextC - the second commonly assigned archeype 
Cluster_df_all <- do.call(rbind.data.frame, tabs_all)

# b) MAIN CLUSTER SET ITERATIONS calculate the mean distance of each cell from its top and second assigned cluster
tabs <- apply(Assigned_cluster_main, 1, calc_arch_per_1km) # topC = the archetype label applied most to each cell. nextC - the second commonly assigned archeype 
Cluster_df <- do.call(rbind.data.frame, tabs)

for (i in 1:(nrow(Assigned_cluster_main))) { # takes a little while to go through all cells
  if (i==1) {
    meandists <-  data.frame(mean_1_dist=NA, mean_2_dist=NA)[0,]
  } 
  GS_dists <- Cluster_dist_main[i,]
  GSclust <- Assigned_cluster_main[i,]
  a <- sort(table(GSclust), decreasing = TRUE)
  if(length(a)==1) { # if only one cluster assigned just put in the mean distance of this
    mx <- mean(GS_dists)
    meandists <- rbind(meandists, data.frame(mean_1_dist=mx, mean_2_dist=NA))
  } else { # if two calculate mean distance of first and second assigned cluster. 
    mx <- mean(GS_dists[which(GSclust %in% names(a)[1])])
    mx2 <- mean(GS_dists[which(GSclust %in% names(a)[2])])
    meandists <- rbind(meandists, data.frame(mean_1_dist=mx, mean_2_dist=mx2))
  }
}

All_1km_main <- cbind(Input_drivers, Cluster_df, meandists, 
                      data.frame('cert_total'=Cluster_df_all$cert, 'cert2_total'=Cluster_df_all$cert2), # this line only works if each cell is majority assigned the same archeype across all iterations and iterations which produce main set. Can check all cert>cert_total 
                      X=ExVars_L$POINT_X, Y=ExVars_L$POINT_Y)
# All input data is in order of output data of SOMs

############ Write map and distance data ####

Tier1_rast <- rasterFromXYZ(All_1km_main[,c('X', 'Y', 'topC')],crs=OSGB)
Tier1_dist <- rasterFromXYZ(All_1km_main[,c('X', 'Y', 'mean_1_dist')],crs=OSGB)

Tier1 <- stack(Tier1_rast, Tier1_dist)

names(Tier1) <- c('Tier1_GB', 'Tier1_dists_GB')

# check plot
plot(Tier1_rast, col=Tier1_pal(no_nodes), legend=FALSE, box=FALSE, axes=FALSE)
legend(legend=sort(unique(Tier1_rast$topC)), fill=Tier1_pal(no_nodes), 
       cex=0.8, x=69e+04, y=1000000, bty='n', y.intersp=2)

# check distances
plot(Tier1_dist, box=FALSE, axes=FALSE)

writeRaster(Tier1, "C:/Users/cecgoo/OneDrive - UKCEH/EIDC/Data/Tier1ArchetypesGB.tif", format='GTiff')

############ Write codebook data ####

Cluster_var_load_main %>%
  group_by(Arch_no) %>%
  dplyr::select(c(names(Input_drivers), 'dist')) %>%
  summarise(across(.cols=everything(), .fns=mean)) -> Cluster_codes 

Tier1_ct <- Cluster_codes[,c('Arch_no', 'dist', names(Input_drivers))]

write.csv(Tier1_ct, file=paste0(file_path, "Archetype_output_data/Tier1_data/Tier1Codes.csv"), row.names=FALSE)

############ Run names script ####

# Run naming script based on codebooks - produces a seperate csv with names. Names can then be translated to a shortened or more readable form. But a useful starting point 

# Keep text as strings
options(stringsAsFactors = FALSE)

file_path <- "P:/NEC07065_AGLAND/WP1/CG/" #"C:/Users/cecgoo/OneDrive - UKCEH/"  #

naming_conv <- paste0(file_path, "Archetype_output_data/LSA_NamingConventions.csv") # A naming convention table which links high, low and mid levels of variables to text e.g. low 'woody linear features' weighting is 'Open'
codes <- Tier1_ct #paste0(file_path, "Archetype_output_data/LSA_data/LSA_code_table.csv")

# Read csv of variable names and definitions (set up for FLAs at present, need to add LSA variables too)
vardefs <- read.csv(naming_conv)

# Set up soil texture triangle
soiltex <- expand.grid(c("SoilSand", "SoilSilt", "SoilClay"),c(-1,1),c("SoilSand", "SoilSilt", "SoilClay"),c(-1,1))
soiltex <- soiltex[soiltex[,1] != soiltex[,3],]
soiltex <- soiltex[order(soiltex$Var2 + soiltex$Var4),]
soiltex$Tex <- c("Clay", "Silty", "Clay", "Sandy", "Silty", "Sandy",
                 "Silty loam", "Clay loam", "Sandy loam", "Clay loam", "Sandy loam","Silty loam",
                 "Sandy loam","Sandy loam","Silty loam","Silty loam","Clay loam","Clay loam",	
                 "Silty sand", "Sandy clay", "Silty sand", "Silty clay", "Sandy clay", "Silty clay")

# Get list of table names of mean (across iterations) variable loadings for LSAs
tab <- read.csv(codes)

# Set the number of variables on which we want to define names # testing with 5
high_thres <- 3
mid_thres <- 1.5
low_thres <- 0.75

# Get variable names
vnams <- names(tab)[-c(1,2)]

# Using an apply here to do the same thing for each row (i.e. each archetype)
archenams <- apply(tab, 1, function(vloadings) {
  
  # Extract loadings
  loads <- as.numeric(vloadings[-c(1,2)])	
  
  # identify variables above a certain threshold of loading
  vars1 <- vnams[which(abs(round(loads,2)) >= high_thres)]
  vars2 <- vnams[which(abs(round(loads,2)) >= mid_thres & abs(round(loads,2)) < high_thres)]
  vars3 <- vnams[which(abs(round(loads,2)) >= low_thres & abs(round(loads,2)) < mid_thres)]
  
  allvars <- vnams[which(abs(round(loads,2)) >= low_thres)]
  
  # Identify in which direction the loading applies
  varsigns <- sign(loads[which(abs(round(loads,2)) >= low_thres)])
  
  # Match to the definitions table
  Vardefs <- vardefs[match(allvars, vardefs$Variable.Name),]
  
  # create data frame of 'highly/very' for variables with a loading above/below high threshold
  HL <- ifelse(Vardefs$Variable.Name %in% vars1 & Vardefs$Category %in% c('Habitats'), 'Highly',
               ifelse(Vardefs$Variable.Name %in% vars1 & Vardefs$Category %in% c('Soil', 'SoilTexture', 'Climate', 'Topography', 'Infrastructure', 'Features', 'FeaturesCount'), 'Very', NA))
  
  #LL <- ifelse(Vardefs$Variable.Name %in% vars3 & Vardefs$Category %in% c('Habitats', 'Soil', 'SoilTexture'), 'Predominantly',
  #             ifelse(Vardefs$Variable.Name %in% vars3 & Vardefs$Category %in% c('Climate', 'Topography', 'Features'), 'Moderately', 
  #                    ifelse(Vardefs$Variable.Name %in% vars3 & Vardefs$Category %in% c('FeaturesCount'), '', NA)))
  
  # Get defintions from low/high columns accoring to signs of loadings
  Vardefs_name <- ifelse(varsigns > 0, Vardefs$High, Vardefs$Low)
  
  # Add on 'highly/very' to high loaded variables 
  Vardefs_names1 <- ifelse(!is.na(HL), paste(HL, Vardefs_name), Vardefs_name)
  #Vardefs_names1 <- ifelse(!is.na(LL), paste(LL, Vardefs_namesA), Vardefs_namesA)
  
  # Hierarchical naming by category
  
  # Get habitat names
  if("Habitats" %in% Vardefs$Category) {
    Vardefs_names1[Vardefs$Category == "Habitats" #& !is.na(LL)
                   ][-1] <- gsub('Predominantly ', '', Vardefs_names1[Vardefs$Category == "Habitats" #& !is.na(LL)
                                                                      ][-1])
    #Vardefs_names1[Vardefs$Category == "Habitats" & !is.na(HL)][-1] <- gsub('Highly ', '', Vardefs_names1[Vardefs$Category == "Habitats" & !is.na(LL)][-1]) # highly should be there twice
    Hab <- Vardefs_names1[Vardefs$Category == "Habitats"]		
    names(Hab) <- allvars[Vardefs$Category == "Habitats"]
  } else {
    Hab <- vardefs$Neutral[vardefs$Variable.Name=='None' & vardefs$Category =='Habitats'] # Adds in 'neutral landscape' definition if no habitats highly loaded
  }
  
  # Get First features names
  if("Features" %in% Vardefs$Category) {
    Fea <- Vardefs_names1[Vardefs$Category == "Features"]	
    names(Fea) <- allvars[Vardefs$Category == "Features"]
  } else {
    Fea <- NULL
  }
  
  # Get Features names
  if("FeaturesCount" %in% Vardefs$Category) {
    FeaC <- Vardefs_names1[Vardefs$Category == "FeaturesCount"]	
  } else {
    FeaC <- NULL
  }
  
  # Get Climate names
  if("Climate" %in% Vardefs$Category) {
    Vardefs_names1[Vardefs$Category == "Climate" #& !is.na(LL)
                   ][-1] <- gsub('Moderately ', '', Vardefs_names1[Vardefs$Category == "Climate" #& !is.na(LL)
                                                                   ][-1])
    Cli <- Vardefs_names1[Vardefs$Category == "Climate"]	
  } else {
    Cli <- NULL
  }
  
  # Get Topography names
  if("Topography" %in% Vardefs$Category) {
    Top <- Vardefs_names1[Vardefs$Category == "Topography"]	
  } else {
    Top <- NULL
  }
  
  # Soil properties
  if("Soil" %in% Vardefs$Category){
    Vardefs_names1[Vardefs$Category == "Soil" #& !is.na(LL)
                   ][-1] <- gsub('Predominantly ', '', Vardefs_names1[Vardefs$Category == "Soil" #& !is.na(LL)
                                                                      ][-1])
    Soi <- Vardefs_names1[Vardefs$Category == "Soil"]	
  } else {
    Soi <- NULL
  }
  
  # Soil texture - This gets handled a bit differently, using the soil texture triangle		
  if(sum(Vardefs$Category == "SoilTexture") == 1){
    Sot <- Vardefs_names1[Vardefs$Category == "SoilTexture"]
    
  } else if(sum(Vardefs$Category == "SoilTexture") >= 2){
    Sot <- soiltex[soiltex$Var1 == allvars[Vardefs$Category == "SoilTexture"][1] &
                     soiltex$Var2 == varsigns[Vardefs$Category == "SoilTexture"][1] &
                     soiltex$Var3 == allvars[Vardefs$Category == "SoilTexture"][2] &
                     soiltex$Var4 == varsigns[Vardefs$Category == "SoilTexture"][2], 5]
  } else {
    Sot <-  NULL
  }
  
  # Now assemble all these components coherently!
  
  archenam <- paste(
    
    ifelse(is.null(Fea), "", paste(paste(Fea[1], collapse=', '), ', ')),
    
    paste(Hab, collapse = ", "), "landscapes",
    
    ifelse(is.null(FeaC), "", paste("with", paste(FeaC, collapse=", "))),
    
    ifelse(!is.null(Top), paste("on", paste(Top, collapse=", "), "land"), ""),
    
    ifelse(is.null(Soi) & is.null(Sot), "",
           ifelse(!is.null(Soi) & !is.null(Sot), paste("with", paste(Soi, collapse = ", "), Sot, "soils"),
                  ifelse(!is.null(Soi), paste("with", paste(Soi,collapse = ", "), "soils"), paste("with", Sot, "soils")))),
    
    ifelse(!is.null(Cli) & is.null(FeaC) & is.null(Soi) & is.null(Sot), paste("with", paste(Cli,collapse = ", "), "climate"),
           ifelse(is.null(Cli), "", paste("and", paste(Cli,collapse = ", "), "climate"))),
    sep = " ")	
  
  # Tidy names
  archenam <- gsub("\\s+", " ",archenam)  # Replace consecutuve spaces with single space
  archenam <- gsub("\\s,", ",",archenam)  # Remove post-comma space
  archenam <- trimws(archenam)		    # Remove leading/trailing whitespace 
  archenam <- tolower(archenam)
  
  return(archenam)
}) # end apply

# Append new names to table
tab$Archetype.name <- archenams

