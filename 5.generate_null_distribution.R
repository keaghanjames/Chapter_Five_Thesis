setwd("~/Dropbox/GALAH/") #you'll need to alter the filepath depending on where you store the Galah data
source("functions.R")
setwd("GALAH_the_canonII/")


#read in empirical dataset
dataset <- 'GALAH'

load('full_dataset_GALAH_CANON.Rdata')
df
# df <- read.csv(paste0(dataset, ".csv"))
# df <- as.tibble(df)
# df$Star_identification <- as.character(df$Star_identification)
# 
# #split the ratios and error into two seperate dataframes
# df_error <- as.tibble(df[,grep("e_", colnames(df))])
# df_ratios <- as.tibble(df[,-grep("e_", colnames(df))])
# 
# dd_ratios <- df_ratios[,grep("_H", colnames(df_ratios))]
# num_na <- vector()
# for(i in 1:ncol(dd_ratios)){
#   num_na <- c(num_na, length(which(is.na(dd_ratios[,i]))))
# }
# num_na <- data_frame(element = colnames(dd_ratios), nNAs = num_na)
# num_na$nNAs_prop = num_na$nNAs/nrow(dd_ratios)
# head(num_na)
# 
# 
# #drop elements with missing data greater than ten per cent
# df_ratios[,num_na$element[which(num_na$nNAs_prop > 0.1)]] <- NULL
# df_error[,paste0(num_na$element[which(num_na$nNAs_prop > 0.1)], "e_")] <- NULL
# 
# #lets look at the error
# dd_error <- df_error
# dd_error <- dd_error[,-(1:8)]
# dd_error <- dd_error %>% gather()
# dd_error$value <- as.numeric(dd_error$value)
# dd_error
# df_ratios
# 
# df_ratios <- df_ratios[,-(2:6)]
# 
# colnames(df_ratios) <- str_remove(colnames(df_ratios), '_H_dex')
# colnames(df_ratios) <- str_remove(colnames(df_ratios), '_H_dex')
# colnames(df_ratios) <- str_to_lower(colnames(df_ratios))
# colnames(df_error) <- str_remove(colnames(df_error), '_H_dex')
# colnames(df_error) <- str_remove(colnames(df_error), '_H_dex')
# colnames(df_error) <- str_remove(colnames(df_error), 'e_')
# colnames(df_error) <- str_to_lower(colnames(df_error))
# colnames(df_error) <- paste0(colnames(df_error), "_e")
# 
# df <- list()
# df[[1]] <- df_ratios
# df[[2]] <- df_error

df <- remove_stars_missing(df)
df[[1]]

volume_beddell <- read.csv("~/Dropbox/GALAH/dynamics.csv")
volume_beddell <- tibble(volume_beddell)
head(volume_beddell)
beddell_R <- range(volume_beddell$R)
beddell_R
beddell_z <- range(volume_beddell$z)
beddell_z
df <- trim_by_volume(df[[1]], R = beddell_R, z = beddell_z, error = df[[2]])#next we want to sample by volume
df


traits <- colnames(df[[1]])[2:17]
#we need to ensure we only include the traits in the Castalli et al. paper


traits <- combn(traits, 2)
traits <- matrix(traits, byrow = T, ncol = 2)
load('MAPtraits.Rdata')
head(MAP_traits)

ratio_maker <- function(element1, element2, data, error_measure = T, scaler = 3, id_col) {
  if(error_measure == T){
    error <- data[[2]]
    ratios <- data[[1]]
    names <- c(ratios[,1])
    traits_e <- vector()
    for(i in 1:dim(error)[1]){
      traits_e <- c(traits_e,
                    #(runif(1, 0, sqrt((error[i, paste0(element1, "_e")])^2 + (error[i, paste0(element2, "_e")])^2))*sample(c(-1,1), 1))
                    #(runif(1, 0, sqrt((ratios[i, paste0(element1)])^2 + (ratios[i, paste0(element2)])^2))*sample(c(-1,1), 1))
                    rnorm(1, 0, ((sqrt((error[i, paste0(element1, "_e")])^2 + (error[i, paste0(element2, "_e")])^2))/scaler)[[1]])
      )
    }
    
  } else {
    names <- data[,1]
    ratios <- data
  }
  ratio <- (ratios[element1] - ratios[element2])
  if(error_measure == T){
    ratio <- ratio + traits_e
  }
  ratio <- data.frame(ratio)
  colnames(ratio) <- paste0(element1, "/", element2)
  rownames(ratio) <- unlist(names)
  return(ratio)
}

nsim <- 100

#we need to create a synthetic database

synthesise_chems <- function(chem) {
  chem <- runif(length(chem), range(chem)[1], range(chem)[2])
}

synthetic_data <- apply(X = as.data.frame(df[[1]])[,-1], MARGIN = 2, FUN = synthesise_chems)

df_synthetic <- df
df_synthetic[[1]][,2:ncol(df_synthetic[[1]])] <- synthetic_data

make_consensus_trees <- function(traits){
  ratio_maker_simple <- function(tr){
    ratio_maker(tr[1], tr[2], data = df_synthetic,
                error_measure = T, scaler = 3, id_col = 'star_identification')
  }
  
  trees <- list()
  class(trees) <- 'multiPhylo'
  pb = txtProgressBar(min = 0, max = nsim, initial = 0) 
  print('Making ratios and building trees')
    for(i in 1:nsim){
    by_cor <- apply(traits, 1, ratio_maker_simple)
    by_cor <- do.call(cbind, by_cor)
    t <- nj(dist(by_cor)) 
    trees[[i]] <- t
    setTxtProgressBar(pb,i)
  }
  
  print('Creating consensus tree')
  maj <- consensus.edges(trees, consensus.tree = consensus(trees, p = 0.5))
  maj$clade_support <- prop.clades(maj,  trees)
  print('Creating mcc tree')
  mcc <- maxCladeCred(trees, rooted = F)
  mcc$clade_support <- prop.clades(mcc,  trees)
  return(list(maj, mcc))
}


make_mcc_trees <- function(traits){
  ratio_maker_simple <- function(tr){
    ratio_maker(tr[1], tr[2], data = df_synthetic,
                error_measure = T, scaler = 3, id_col = 'star_identification')
  }

  trees <- list()
  class(trees) <- 'multiPhylo'
  pb = txtProgressBar(min = 0, max = nsim, initial = 0)
  print('Making ratios and building trees')
  for(i in 1:nsim){
    by_cor <- apply(traits, 1, ratio_maker_simple)
    by_cor <- do.call(cbind, by_cor)
    t <- nj(dist(by_cor))
    trees[[i]] <- t
    setTxtProgressBar(pb,i)
  }

  print('Creating consensus tree')
  mcc <- maxCladeCred(trees, rooted = F)
  mcc$clade_support <- prop.clades(mcc,  trees)
  return(mcc)
}

system.time(maj <- make_consensus_trees(traits = MAP_traits))

make_synthetic_consensus <- function(run){
  synthetic_data <- apply(X = as.data.frame(df[[1]])[,-1], MARGIN = 2, FUN = synthesise_chems)
  df_synthetic <- df
  df_synthetic[[1]][,2:ncol(df_synthetic[[1]])] <- synthetic_data
  maj <- make_consensus_trees(traits = MAP_traits)
  write.tree(maj[[1]], paste0('NULL/' ,run, 'consensus.tree'))
  write.tree(maj[[2]], paste0('NULL/' ,run, 'mcc.tree'))
  return(maj)
}

#run = 2

nsim <- 100
system.time(maj <- make_synthetic_consensus(run = 1))
require(pbmcapply)
detectCores()
pbmclapply(X = c(2:100), make_synthetic_consensus, mc.cores = 4)


make_synthetic_populations <- function(pops, df = df){

  #get synthetic data for each population
  synthetic_data <- apply(X = as.data.frame(df[[1]])[,-1], MARGIN = 2, FUN = synthesise_chems)[sample(seq(1:nrow(df[[1]])), pops), ]
  synthetic_data <- as.data.frame(synthetic_data)

  #draw sampling probablities from uniform distribution 
  probs <- runif(pops, 0, 1)
  probs <- probs/sum(probs)
  probs
  
  #create synthetic data
  population <- sample(seq(1, pops), size = nrow(df[[1]]), replace = T, prob = probs)
  
  df_synthetic <- tibble()
  for(i in 1:length(population)){
  
    df_synthetic <- bind_rows(df_synthetic, synthetic_data[population[i],])
    #print(i)
  
  }
  df_synthetic <- bind_cols(df[[1]]$star_id, df_synthetic)
  colnames(df_synthetic)[1] <- 'star_id'
  df_synthetic <- list(df_synthetic, df[[2]], pops)
  # synthesised_data[[1]] <- df_synthetic
  # synthesised_data[[2]] <- df[[2]]
  return(df_synthetic)
}


df_synthetic <- make_synthetic_populations(5, df = df)

#maj <- make_consensus_trees(traits = MAP_traits)

populations <- sample(seq(2,10,1), 100, replace = T)

make_population_consensus <- function(run){
  df_synthetic <- make_synthetic_populations(populations[run], df = df)
  maj <- make_consensus_trees(traits = MAP_traits)
  # maj[[1]]$clade_support <- prop.clades(maj,  trees)
  # mcc <- maxCladeCred(trees, rooted = F)
  # mcc$clade_support <- prop.clades(mcc,  trees)
  write.tree(maj[[1]], paste0('POP_NULL/' ,run, 'consensus.tree'))
  write.tree(maj[[2]], paste0('POP_NULL/' ,run, 'mcc.tree'))
  return(maj)
  }

make_population_consensus(run)

#test_tree <- make_consensus_trees(traits = MAP_traits)
detectCores()
pbmclapply(X = c(2:100), make_population_consensus, mc.cores = 4)


#analyse output

require(RPANDA)
?RPANDA
?RPANDA::JSDtree()

mrc_list <- list()
mcc_list <- list()
p_mrc_list <- list()
p_mcc_list <- list()

for(i in 1:100){
  mrc_list[[i]] <- read.tree(paste0('NULL/', i, 'consensus.tree'))
  mcc_list[[i]] <- read.tree(paste0('NULL/', i, 'mcc.tree'))
  p_mrc_list[[i]] <- read.tree(paste0('POP_NULL/', i, 'consensus.tree'))
  p_mcc_list[[i]] <- read.tree(paste0('POP_NULL/', i, 'mcc.tree'))
  print(i)
}

class(mrc_list) <- 'multiPhylo'
class(mcc_list) <- 'multiPhylo'
class(p_mrc_list) <- 'multiPhylo'
class(p_mcc_list) <- 'multiPhylo'

spectSmcc <- list()
spectSmrc <- list()
spectSp_mcc <- list()
spectSp_mrc <- list()

for(i in 1:length(mcc_list)){
  spectSmcc[[i]] <- spectR(mcc_list[[i]])
  spectSmrc[[i]] <- spectR(mrc_list[[i]])
  spectSp_mcc[[i]] <- spectR(p_mcc_list[[i]])
  spectSp_mrc[[i]] <- spectR(p_mrc_list[[i]])
  print(i)
}

null_df <- tibble()

for(i in 1:length(spectSmcc)) {
  null_df <- bind_rows(null_df, c(unlist(spectSmcc[[i]][2:5]), 'MCC', 'Null'))
  null_df <- bind_rows(null_df, c(unlist(spectSmrc[[i]][2:5]), 'MRC', 'Null'))
  null_df <- bind_rows(null_df, c(unlist(spectSp_mcc[[i]][2:5]), 'MCC', 'Population Null'))
  null_df <- bind_rows(null_df, c(unlist(spectSp_mrc[[i]][2:5]), 'MRC', 'Population Null'))
  print(i)
}

empirical_mcc <- read.tree('~/Dropbox/GALAH/HARPS/MCC.tree')
empirical_maj <- read.tree('~/Dropbox/GALAH/HARPS/MAJ.tree')
mcc_spect <- spectR(empirical_mcc)
maj_spect <- spectR(empirical_maj)

null_df <- bind_rows(null_df,
                     bind_rows(c(unlist(mcc_spect[2:5]), 'MCC', 'Empirical'),
                     c(unlist(maj_spect[2:5]), 'MRC', 'Empirical'))
  )
colnames(null_df)[5:6] <- c('Summary_tree', 'Model')

null_df[,1:4] <- null_df[,1:4] %>% mutate_if(is.character,as.numeric)
null_df <- null_df %>% mutate_if(is.character,as.factor)

tail(null_df)

write.csv(null_df, 'null_spectra.csv')

require(ggplot2)
require(ggdist) 
require(ggridges)


p1 <- ggplot(data = null_df, aes(y = Summary_tree, x = asymmetry, fill = Model)) +
  stat_dist_halfeye()
p1
