# Written by Keaghan Yaxley

setwd("~/Dropbox/GALAH/")

library(phytools)
library(ape)
library(maps)
library(phangorn)
library(tidyverse)
library(ggtree)
library(adephylo)
library(paleotree)
library(CollessLike)
require(parallel)
require(ggridges)
require(biotools)


options(scipen = 999) #we need to turn off scientific notation because the names of these stars are long


trim_by_volume <- function(data, R, z, error = NULL){
  RZ <- which(data$R >= R[1] & data$R <= R[2] & data$z >= z[1] & data$z <= z[2])
  data <- data[RZ,]
  if(is.null(error) == F){
  error <- error[RZ, ]
  list <- list()
  list[[1]] <- data
  list[[2]] <- error
  data <- list
  }
  return(data)
}



remove_stars_missing <- function(data){
  d1 <- data[[1]]
  d2 <- data[[2]]
  d2 <- d2[-which(rowSums(is.na(d1)) > 0),]
  d1 <- d1[-which(rowSums(is.na(d1)) > 0),]
  list <- list()
  list[[1]] <- d1
  list[[2]] <- d2
  return(list)
}



###################################################
### Calculating Ratios and Accounting for Error ###
###################################################

#we need a function that takes two elements and a solar data as an input, calculates the ratios and can add random error
# colnames(new_df[[1]])
# element1 <- "o_fe"
# element2 <- "na_fe"
# data <- new_df

get_ratio <- function(element1, element2, data, error_measure = T, scaler = 3, id_col) {
  if(error_measure == T){
    error <- data[[2]]
    ratios <- data[[1]]
    names <- c(ratios[id_col])
    # if(is.na(error[i, paste0(element1, "_e")]) == T){
    #   ratios[i, element1] <- ratios[i, element1]
    # } else {
    traits_e <- vector()
    for(i in 1:dim(error)[1]){
      traits_e <- c(traits_e, 
                   #(runif(1, 0, sqrt((error[i, paste0(element1, "_e")])^2 + (error[i, paste0(element2, "_e")])^2))*sample(c(-1,1), 1))
                   #(runif(1, 0, sqrt((ratios[i, paste0(element1)])^2 + (ratios[i, paste0(element2)])^2))*sample(c(-1,1), 1))
                   rnorm(1, 0, ((sqrt((error[i, paste0(element1, "_e")])^2 + (error[i, paste0(element2, "_e")])^2))/scaler)[[1]])
      )
    }
      # if(is.na(error[i, paste0(element1, "_e")]) == T){
      #   ratios[i, element1] <- ratios[i, element1]
      # } else {
      #   ratios[i, element1] <- ratios[i, element1] + runif(1, 0, error[i, paste0(element1, "_e")])*(sample(c(-1,1),1))
      # }
      # if(is.na(error[i, paste0(element2, "_e")]) == T){
      #   ratios[i, element2] <- ratios[i, element2]
      # } else {
      #   ratios[i, element2] <- ratios[i, element2] + runif(1, 0, error[i, paste0(element2, "_e")])*(sample(c(-1,1),1))
      # }
    } else {
    names <- data[id_col]
    ratios <- data
  }
  ratio <- (ratios[element1] - ratios[element2])
  if(error_measure == T){
    ratio <- ratio + traits_e
  }
  #colnames(ratio) <- paste0(element1, "/", element2)  
  ratio <- tibble(id_col = (names)[[1]], ratio)
  #ratio$id_col <- c(names)
  colnames(ratio)[2] <- paste0(element1, "/", element2)
  return(ratio)
}


#df_traits <- apply(traits, 1, simple_get_ratio)

#function for cleaning up output of simple_get_ratio when run with apply
clean_up_ratios <- function(ratio_output) { #this takes as an input the output of simple_get_ratio and the original data
  new_ratios <- ratio_output[[1]][2]
  for(i in 1:length(ratio_output)){
    new_ratios <- cbind(new_ratios, ratio_output[[i]][1])
  }
  rownames(new_ratios) <- unlist(as.list(new_ratios[,1]))
  new_ratios <- new_ratios[,-1]
  return(new_ratios)
}

#df_traits <- clean_up_ratios(df_traits)

#colnames(df_traits)[1] <- "galah_id"
#now we need to get our samples
# nsim <- 100
# stars <- df_ratios$galah_id
# length_sample <- 100
# replace <- T

get_samples <- function(nsim, stars, length_sample, replace = F, outgroup = NULL) {
  sampled_stars <- list()
  for(i in 1:nsim){
    sam <- sample(stars, length_sample, replace = replace)
    if(replace == F){
      stars <- stars[-which(stars %in% sam == T)]
    }
    if(is.null(outgroup) == F){
      sam <- c(sam, outgroup)
    }
    sampled_stars <- append(sampled_stars, list(sam))
  }
  return(sampled_stars)
}

#the above function creates mulitple samples of our stars with or without replacement
#and we can specify the names of our outgroup stars so they can be included in each sample

#sampled_stars <- get_samples(nsim = nsim, stars = stars, length_sample = length_sample, outgroup = df_outgroup$galah_id)
#sampled_stars #it gives back a list of lists


#####################################
### Distance NJ Trees for Samples ###
#####################################

#we need to specify the elements we want to include in the analysis
# elements <- c("c_fe", "o_fe", "na_fe","mg_fe","al_fe","si_fe","ca_fe","sc_fe","ti_fe","v_fe", "cr_fe","mn_fe","co_fe","ni_fe", 
#               "cu_fe","zn_fe","sr_fe","y_fe", "ba_fe","la_fe","ce_fe","nd_fe","eu_fe")

# #stars <- sampled_stars[[1]]
# data <- cbind(df_ratios, df_traits[,2:ncol(df_traits)])
# 
# get_some_stellar_trees <- function(data, elements, stars, id_col){
#   
#   el_d <- data[which(data[id_col] %in% stars),c(elements, id_col)]
#   rownames(el_d) <- el_d[id_col]
#   el_d$id_col <- NULL
#   
#   distance <- dist(el_d)
#   tr <- nj(distance)
# #  tr$outgroup <- data[which(data$galah_id %in% stars), "outgroup"]
# #  names(tr$outgroup) <- data[which(data$galah_id %in% stars), "galah_id"]
#   return(tr)
# }
#(
similarity_tree <- function(traits){
  tr <- nj(dist(traits))
  return(tr)
}


neighbour_net_tree <- function(traits){
  tr <- neighborNet(dist(traits))
  return(tr)
}


Mahalanobis_dist_tree <- function(traits){
  cov <- cov(traits)
  dist <- D2.dist(traits, cov)
  tr <- nj(dist)
}

go_get_it <- function(stars){
  get_some_stellar_trees(data = data, elements = elements, stars = unlist(stars))
}

ratio_maker <- function(element1, element2, data, error_measure = T, scaler = 3, id_col) {
  if(error_measure == T){
    error <- data[[2]]
    ratios <- data[[1]]
    names <- c(ratios[id_col])
    traits_e <- vector()
    for(i in 1:dim(error)[1]){
      traits_e <- c(traits_e, 
                    #(runif(1, 0, sqrt((error[i, paste0(element1, "_e")])^2 + (error[i, paste0(element2, "_e")])^2))*sample(c(-1,1), 1))
                    #(runif(1, 0, sqrt((ratios[i, paste0(element1)])^2 + (ratios[i, paste0(element2)])^2))*sample(c(-1,1), 1))
                    rnorm(1, 0, ((sqrt((error[i, paste0(element1, "_e")])^2 + (error[i, paste0(element2, "_e")])^2))/scaler)[[1]])
      )
    }
    
  } else {
    names <- data[id_col]
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


tree_build <- function(traits, data, id_col, truncate_thresh = NULL, is.Mahalanobis.dist = F, scaler = 3){
  ratio_maker_simple <- function(element){
    ratio_maker(element1 = element[1], element2 = element[2], data = data, error_measure = T, scaler = scaler, id_col = id_col)
  }
  
  # simple_get_ratio <- function(element){
  #   get_ratio(element1 = element[1], element2 = element[2], data = data, error_measure = T, scaler = scaler, id_col = id_col)
  # }
  df <- apply(traits, 1, ratio_maker_simple)
  df <- do.call(cbind, df)
    # df <- clean_up_ratios(apply(traits, 1, simple_get_ratio))
  if(is.Mahalanobis.dist == F){
    t <- similarity_tree(df)
  } else {
    t <- Mahalanobis_dist_tree(df)
  }
  if(!is.null(truncate_thresh)){
    #hresh <- calc_median_uncertain(data = data, traits = traits)  #sigma
    print(paste0('trunctaion threshold ', truncate_thresh))
    print(truncate_thresh)
    t <- di2multi(t, tol = truncate_thresh, random = F)
  }
  return(t)
}





#run the nj process in parallel
#stellar_trees <- mclapply(sampled_stars, go_get_it, mc.cores = 4)

###############################
### Measuring Tree Distance ###
###############################

#CLI <- unlist(lapply(stellar_trees, colless.like.index))

gamma_ultra <- function(tree){
  t <- gammaStat(force.ultrametric(tree, method = 'extend'))
  return(t)
}



calc_uncertain <- function(data, traits, median = F){
  error <- data[[2]]
  traits_e <- vector()
  for(j in 1:nrow(traits)){
    element1 <- traits[j, 1]
    element2 <- traits[j, 2]
    for(i in 1:dim(error)[1]){
      traits_e <- c(traits_e, sqrt((error[i, paste0(element1, "_e")])^2 + (error[i, paste0(element2, "_e")])^2)
      )
    }
  }
  
  traits_e <- matrix(traits_e, ncol = nrow(traits))
  trait_names <- vector()
  for(i in 1:nrow(traits)){
    trait_names <- c(trait_names, paste0(traits[i,1], "/", traits[i,2]))
  }
  colnames(traits_e) <- trait_names
  
  if(median == T){
    traits_e <- median(traits_e)
  }
  return(traits_e)
}


#Gamma_ultra <- unlist(lapply(stellar_trees, gamma_ultra))


# plot(CLI, Gamma_ultra)
# index <- data_frame(CLI = scale(CLI), Gamma = scale(Gamma_ultra))
# hist(unlist(index))


#############################
### Creating a null model ###
#############################


# range <- range(na.omit(unlist(df_ratios[elements])))
# 
# data_null <- data
# 
# data_null[elements] <- matrix(runif(length(data_null[elements]), min = range[1], max = range[2]),
#                               nrow = nrow(data_null), ncol = ncol(data_null[elements]))
# 
# 
# 
# 
# go_get_it_null <- function(stars){
#   get_some_stellar_trees(data = data_null, elements = elements, stars = unlist(stars))
# }
# 
# 
# stellar_trees_null <- mclapply(sampled_stars, go_get_it_null, mc.cores = 10)
# 
# CLI_null <- unlist(lapply(stellar_trees_null, colless.like.index))
# Gamma_ultra_null <- unlist(lapply(stellar_trees_null, gamma_ultra))
# 
# #need to merge the dataframes and plot
# 
# indices <- as_data_frame(cbind(CLI, Gamma_ultra))
# indices$model <- "H1"
# 
# indices_null <- as_data_frame(cbind(CLI_null, Gamma_ultra_null))
# indices_null$model <- "H0"
# colnames(indices_null) <- colnames(indices)
# 
# indices <- rbind(indices, indices_null)
# head(indices)
# 
# ggplot(indices, aes(x = scale(CLI), y = scale(Gamma_ultra), color = model)) +
#   geom_point()
# 
# #maybe we can try permuting it
# 
# dist(x = index)
# 
# hist(CLI)
# hist(Gamma_ultra)
# 
# ltt.plot(force.ultrametric(t, method = 'extend'))
# 
# stellar_trees[[1]]$outgroup
# 
# getMRCA(stellar_trees[[1]], tip = names(which(stellar_trees[[1]]$outgroup == T)))
# 
# 
# 
# stellar_trees[[1]]<-paintBranches(stellar_trees[[1]],edge=sapply(b,match,stellar_trees[[1]]$tip.label),
#                                   state="b",anc.state="a")
# 
# cols<-setNames(c("black","blue"),c(T,F))
# plot(stellar_trees[[1]],colors=cols,lwd=4,split.vertical=TRUE, show.tip.label = F)
# 
# 
# dotTree(ladderize(stellar_trees[[1]]), stellar_trees[[1]]$outgroup)
# 
# t <- stellar_trees[[1]]
# t$outgroup[which(t$outgroup == T)] <- 'indianred'
# t$outgroup[which(t$outgroup == F)] <- 'lightblue'
# 
# require(ggtree)
# ggtree(t) + 
#   geom_tiplab(color = t$outgroup) 
# 
# geom_hilight(node=getMRCA(stellar_trees[[1]], tip = names(which(stellar_trees[[1]]$outgroup == T)))
#              , fill="gold") 
# 
# 
# grp <- stellar_trees[[1]]$outgroup
# 
# stellar_trees[[1]]$outgroup
# 
# st
# 
# go_get_it(sampled_stars[5])
# t <- stellar_trees[[1]]
# tr$outgroup <- data[which(data$galah_id %in% stars), "outgroup"]
# names(tr$outgroup) <- data[which(data$galah_id %in% stars), "galah_id"]
