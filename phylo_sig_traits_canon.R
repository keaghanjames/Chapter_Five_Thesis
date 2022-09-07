setwd("~/Dropbox/GALAH/") #you'll need to alter the filepath depending on where you store the Galah data
source("functions.R")
setwd("GALAH_the_canonII")


dataset <- 'GALAH_CANON'
df <- read.csv('TheCannon_solar_twin_full_set_predict_14xfe_39548spectra_wVAC.csv')
dim(df)
df$X  <- NULL
#wow so big

df <- as.tibble(df)
df$star_id <- as.character(df$star_id)

R <- df$dynamics_vac_R
z <- df$dynamics_vac_z
#these are the volumes of our stars so we need these for later

df
#split the ratios and error into two seperate dataframes
df_error <- as.tibble(df[,grep("e_", colnames(df))])
df_ratios <- as.tibble(df[,-grep("e_", colnames(df))])
df_ratios
df_error

#lets look at the error 
dd_error <- df_error
dd_error <- dd_error[,-(1:4)]

dd_error <- dd_error %>% gather()
dd_error$value <- as.numeric(dd_error$value)
ggplot(dd_error, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7, scale = 3,
  ) +
  ylab("") +
  xlab("Error") +
  scale_y_discrete(labels = unique(dd_error$key))+
  theme_classic() +
  xlim(c(0,1)) +
  theme(axis.text.y = element_text(face = 'bold', size = 14))

df_ratios
#df_ratios <- df_ratios[,-(2:6)]

colnames(df_ratios) <- str_remove(colnames(df_ratios), '_fe')
colnames(df_error) <- str_remove(colnames(df_error), '_fe')

colnames(df_error) <- str_remove(colnames(df_error), 'e_')
colnames(df_error) <- paste0(colnames(df_error), "_e")
df_error


df <- list()
df[[1]] <- df_ratios[,-c(2:16, 31:75)]
df[[2]] <- df_error[,-(1:4)]
#df <- remove_stars_missing(df)

df[[1]]

rm(df_error, df_ratios, dd_error)

#now we will add volume
df[[1]]$R <- R
df[[1]]$z <- z

save(df, file = paste0('full_dataset_', dataset, '.Rdata'))


elements <- colnames(df[[1]])[2:15]
traits <- combn(elements, 2)
traits <- matrix(traits, byrow = T, ncol = 2)
traits


#now we are going to sample just stars that fall within the same volume as Beddell
volume_beddell <- read.csv("~/Dropbox/GALAH/dynamics.csv")
volume_beddell <- tibble(volume_beddell)
head(volume_beddell)
beddell_R <- range(volume_beddell$R)
beddell_R
beddell_z <- range(volume_beddell$z)
beddell_z
df <- trim_by_volume(df[[1]], R = beddell_R, z = beddell_z, error = df[[2]])#next we want to sample by volume
df


# star_sample <- sample(seq(1, nrow(df[[1]]), 1), 200)
# df[[1]] <- df[[1]][star_sample,]
# df[[2]] <- df[[2]][star_sample,]


#plot the error in our sample
df_error_plot <- df[[2]][,which(colnames(df[[2]]) %in% paste0(elements, '_e'))]
colnames(df_error_plot) <- str_remove(colnames(df_error_plot), '_e')
df_error_plot <- pivot_longer(df_error_plot, cols = Na:Ba)
#df_error_plot$name <- str_to_title(df_error_plot$name)


colours <- c('#1B1016', '#5669A6', '#89AFDC', '#98B1D3', '#6C2B16')

p1 <- ggplot(df_error_plot, aes(x = value, y = name, fill = stat(x))) +
  geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.2, scale = 3,
  ) +
  ylab("") +
  xlab("Error") +
  scale_y_discrete(labels = unique(df_error_plot$name))+
  theme_classic() +
  scale_fill_viridis_c(name = "Error", option = "C") +
  theme(axis.text.y = element_text(face = 'bold', size = 14))
p1

require(ggdist)
p2 <- ggplot(df_error_plot, aes(x = value, y = name, color = name)) +
  stat_dots(fill = NA, scale = 3, quantiles = 50) +
  geom_density_ridges_gradient(quantile_lines = F, quantiles = c(0.025, 0.975), alpha = 0, scale = 2,
                               fill = NA, color = rgb(0,0,0,.5)) +
  ylab("") +
  xlab("Error") +
  theme_classic() +
  theme(legend.position = "None") 
p2

ggsave('Figures/error_distributions.pdf', scale = .5)


### MCMC STARTS HERE ###

setwd('chain_2/')
#if you need to star mid chain
# require(prodlim)
# trace <- read.csv('trace2.csv')
# head(trace)
# trace$X <- NULL
# load('traits_considered.Rdata')
# load('HARPStrees2.Rdata')
# clade_support <- trace$clade_support
# Nnode <- trace$Nnode
# Ntraits <- trace$Ntraits
# pois <- trace$pois
# poisson_shape <- pois[length(pois)]
# trace <- trace$clade_support
# cores <- 4
# sample_traits <- traits_considered[[length(traits_considered)]]
# traits <- traits[-(which(do.call(paste0, as.data.frame(traits)) %in% do.call(paste0, as.data.frame(sample_traits)) == T)),]
# nsim <- 6000
#now jump to the mcmc function

# chain <- 1
# 
# setwd(paste0('chain_', chain))
# 
poisson_shape <- 5
Ntraits_sampled <- rpois(1, poisson_shape) + 1
Ntraits_sampled
sample <- sample(seq(1, nrow(traits), 1), Ntraits_sampled)
sample_traits <-traits[sample,]
traits <- traits[-sample,]
traits
sample_traits

tree_build_mcmc <- function(tree_ID) {
  t <- tree_build(traits = sample_traits, data = df, id_col = 'star_id', is.Mahalanobis.dist = F)
  print(tree_ID)
  return(t)
}
t <- tree_build_mcmc(1)


require(parallel)
require(pbmcapply)

cores <- 10

#let's generate a starting tree
system.time(t1 <- pbmclapply(paste0("Tree_", seq(1,100,1)), tree_build_mcmc, mc.cores = cores))
class(t1) <- 'multiPhylo'
maj1 <- consensus.edges(t1, consensus.tree = consensus(t1, p = 0.5))
maj1$clade_support <- prop.clades(maj1,  t1)
#maj1$params <- gamma_params 
plot(maj1, show.tip.label = F, type = 'fan')

trace <- sum(maj1$clade_support)
Ntraits <- nrow(sample_traits)
Nnodes <- maj1$Nnode
trees <- list()
traits_considered <- list()
traits_considered[[1]] <- sample_traits
trees[[1]] <- maj1
save(trees, file = 'HARPStrees2.Rdata')

nsim <- 10000
clade_support <- c(sum(maj1$clade_support))
Nnode <- c(maj1$Nnode)
pois <- c(poisson_shape)


for(i in length(trace):nsim) {
  sample_traits_old <- sample_traits
  traits_old <- traits
  
  if(length(traits) <= 2){
    
    # Ntraits_sampled <- round(rgamma(1, shape = gamma_params[[1]], 
    #                                 rate = gamma_params[[2]]))
    
    Ntraits_sampled <- rpois(1, poisson_shape) + 1
    
    while(Ntraits_sampled > nrow(sample_traits) | Ntraits_sampled == 0){
      # Ntraits_sampled <- round(rgamma(1, shape = gamma_params[[1]], 
      #                                 rate = gamma_params[[2]]))
      Ntraits_sampled <- rpois(1, poisson_shape) + 1
      
    }
    
    sample <- sample(seq(1, nrow(sample_traits), 1), Ntraits_sampled)
    traits <- rbind(traits, sample_traits[sample,])
    sample_traits <- sample_traits[-sample,]
  } else {
    
    if(length(sample_traits <= 10)){
      
      # Ntraits_sampled <- round(rgamma(1, shape = gamma_params[[1]], 
      #                                 rate = gamma_params[[2]]))
      Ntraits_sampled <- rpois(1, poisson_shape) + 1
      
      while(Ntraits_sampled > nrow(traits) | Ntraits_sampled == 0){
        # Ntraits_sampled <- round(rgamma(1, shape = gamma_params[[1]], 
        #                                 rate = gamma_params[[2]]))
        Ntraits_sampled <- rpois(1, poisson_shape) + 1 
      }
      
      sample <- sample(seq(1, nrow(traits), 1), Ntraits_sampled)
      sample_traits <- rbind(sample_traits, traits[sample,])
      traits <- traits[-sample,]
      
    } else {
      if(runif(1,0,1) > 0.5){
        
        # Ntraits_sampled <- round(rgamma(1, shape = gamma_params[[1]], 
        #                                 rate = gamma_params[[2]]))
        Ntraits_sampled <- rpois(1, poisson_shape) + 1 
        
        while(Ntraits_sampled > nrow(traits) | Ntraits_sampled == 0){
          # Ntraits_sampled <- round(rgamma(1, shape = gamma_params[[1]], 
          #                                 rate = gamma_params[[2]]))
          Ntraits_sampled <- rpois(1, poisson_shape) + 1 
          
        }
        
        sample <- sample(seq(1, nrow(traits), 1), Ntraits_sampled)
        sample_traits <- rbind(sample_traits, traits[sample,])
        traits <- traits[-sample,]
      } else {
        
        # Ntraits_sampled <- round(rgamma(1, shape = gamma_params[[1]], 
        #                                 rate = gamma_params[[2]]))
        Ntraits_sampled <- rpois(1, poisson_shape) + 1
        
        while(Ntraits_sampled > nrow(sample_traits) | Ntraits_sampled == 0){
          # Ntraits_sampled <- round(rgamma(1, shape = gamma_params[[1]], 
          #                                 rate = gamma_params[[2]]))
          Ntraits_sampled <- rpois(1, poisson_shape) + 1
          
        }
        
        sample <- sample(seq(1, nrow(sample_traits), 1), Ntraits_sampled)
        traits <- rbind(traits, sample_traits[sample,])
        sample_traits <- sample_traits[-sample,]
      }
    }
    
  }
  
  
  t2 <- try(pbmclapply(paste0("Tree_", seq(1,100,1)), tree_build_mcmc, mc.cores = cores))
  while(class(t2) == 'try-error'){
    t2 <- try(pbmclapply(paste0("Tree_", seq(1,100,1)), tree_build_mcmc, mc.cores = cores))
  }
  class(t2) <- 'multiPhylo'
  maj2 <- consensus.edges(t2, consensus.tree = consensus(t2, p = 0.5))
  maj2$clade_support <- prop.clades(maj2,  t2)
  #plot(maj2)
  
  
  # if(sum(maj2$clade_support) > sum(maj1$clade_support)){
  #   t1 <- t2
  #   maj1 <- maj2
  #   trace <- c(trace, sum(maj1$clade_support))
  # } else {
  
  #change it so that it always accepts a better option but will sometimes accept a worse option
  p <- sum(maj2$clade_support)/sum(maj1$clade_support)
  
  if(p > runif(1,0.8,1)) { #you could use a different distribution
    t1 <- t2
    maj1 <- maj2
    trace <- c(trace, sum(maj1$clade_support))
    #gamma_params[[1]] <- gamma_params[[1]]/p
    poisson_shape <- poisson_shape/p
    if(poisson_shape > 20){poisson_shape = 1}
    
  } else {
    sample_traits <- sample_traits_old
    traits <- traits_old
    trace <- c(trace, sum(maj1$clade_support))
    #perhaps here you could alter the shape of gamma to tune it
    #gamma_params[[1]] <- gamma_params[[1]]/p
    poisson_shape <- poisson_shape/p
    if(poisson_shape > 20){poisson_shape = 1}
  }
  #}
  #plot(maj1)
  # maj1$params <- gamma_params 
  maj1$params <- poisson_shape 
  
  Ntraits <- c(Ntraits, nrow(sample_traits))
  #Nnodes <- c(Nnodes, maj1$Nnode)
  trees[[i+1]] <- maj1
  
  clade_support <- c(clade_support, sum(maj1$clade_support))
  Nnode <- c(Nnode, maj1$Nnode)
  pois <- c(pois, poisson_shape)
  
  print(set_names(c(i, sum(maj1$clade_support), maj1$Nnode, nrow(sample_traits), poisson_shape), nm = c("count", "clade_support", "Nnode", "Ntraits", 'lambda')))
  traits_considered[[i+1]] <- sample_traits
  
  # par(mfrow = c(3,1))
  # plot(x = seq(1, i+1, 1), y = clade_support, xlim = c(0, nsim), type = 'l', main = 'Clade Support')
  # plot(x = seq(1, i+1, 1), y = Nnode, xlim = c(0, nsim), type = 'l', main = 'Nnode')
  # plot(x = seq(1, i+1, 1), y = Ntraits, xlim = c(0, nsim), type = 'l', main = 'Ntraits')
  # 
  save(trees, file = 'GALAHtrees2.Rdata')
  save(traits_considered, file = 'traits_considered.Rdata')
  write.csv(tibble(clade_support, Nnode, pois, Ntraits), 'trace2.csv')
  # 
  
}

# class(trees) <- 'multiPhylo'
# dev.off()
# #densiTree(trees)
# plot(maj1)
# traits
# sample_traits
# 
# 
# save(MAP_traits, file = 'MAPtraits.Rdata')
# 
# sample_traits
# clade_support[which(clade_support == max(clade_support))]
# MAP_traits <- traits_considered[which(clade_support == max(clade_support))][[1]]
# MAP_tree <- trees[which(clade_support == max(clade_support))][[1]]
# plot(MAP_tree)
# 
# 
# tree_build_final <- function(tree_ID) {
#   t <- tree_build(traits = MAP_traits, data = df, id_col = 'star_id', is.Mahalanobis.dist = F)
#   print(tree_ID)
#   return(t)
# }
# 
# MAP_trees <- pbmclapply(paste0("Tree_", seq(1,1000,1)), tree_build_final, mc.cores = cores)
# class(MAP_trees) <- 'multiPhylo'
# mccTree <- maxCladeCred(MAP_trees)
# plot(mccTree)
# mrcTree <- consensus.edges(MAP_trees, consensus.tree = consensus(MAP_trees, p = 0.5))
# plot(mrcTree, show.tip.label = F, type = 'fan')
# 
# 
# require(ggtree)
# dataset <- 'GALAH_CANON'
# ages <- read.csv('~/Dropbox/GALAH/GALAH_the_canonII/TheCannon_solar_twin_full_set_predict_14xfe_39548spectra_wVAC.csv')
# ages <- ages[which(ages$star_id %in% mccTree$tip.label), c(2,4)]
# ages <- as.tibble(ages)
# colnames(ages)[1] <- 'label'
# 
# #mrcTree$edge.length <- mrcTree$edge.lengthh[-1]
# 
# MAJTree <- full_join(mxcTree, ages, by = "label") 
# tp2 <- ggtree(mxcTree, aes(color = ages$age), ladderize = T, layout = 'circular',  right = T)+
#   scale_color_continuous(type = 'viridis') +
#   geom_tippoint()+
#   labs(color = "Star ages (Gy)", cex = 20)
# tp2
# 
# MCC_Tree <- full_join(mccTree, ages, by = "label") 
# tp3 <- ggtree(MCC_Tree, aes(color = age), ladderize = T, layout = 'circular',  right = T)+
#   scale_color_continuous(type = 'viridis') +
#   geom_tippoint()+
#   labs(color = "Star ages (Gy)", cex = 20)
# tp3
