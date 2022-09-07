#what I need to do

# I need to write code that combines the trace data and the traits data

#it then identifies that MAP set of traits

#it then builds a tree from this data

#at present we only have complete chains for 2 and 3 

setwd("~/Dropbox/GALAH/") 
source("functions.R")
setwd("GALAH_the_canonII/")
colours <- c('#1B1016', '#5669A6', '#89AFDC', '#98B1D3', '#6C2B16')
gradient <- c("#EE05F2", '#6B1C8C', '#011126', '#143D73', '#0BBFBF')

trace_1 <- read.csv('chain_1/trace2.csv')
head(trace_1)
colnames(trace_1)[1] <- 'generation'
trace_1$chain <- 1


trace_2 <- read.csv('chain_2/trace2.csv')
head(trace_2)
colnames(trace_2)[1] <- 'generation'
trace_2$chain <- 2

trace <- rbind(trace_1, trace_2)
head(trace)
trace$chain <- as.factor(trace$chain)

#next we need to trim the first thousand generations from the trace dataset
trace <- trace[which(trace$generation > 1000),]


#want to import the traits we considered
load('chain_1/traits_considered.Rdata')
traits_considered1 <- traits_considered[-c(1:1000)]
load('chain_2/traits_considered.Rdata')
traits_considered2 <- traits_considered[-c(1:1000)]

#next we need to indentify our MAP traits
MAP_gen <- trace[which(trace$clade_support == max(trace$clade_support)),]
MAP_gen 

#so they are all in chain 1 and next to each other
#lets look at them 

MAP_traits <- traits_considered1[MAP_gen$generation]
MAP_traits[[1]] == MAP_traits[[2]] #and yes they are exactly the same set of traits

MAP_traits <- MAP_traits[[1]]

save(MAP_traits, file = 'MAPtraits.Rdata')

#now that we have identified our MAP traits we can use them to build our tree

#we load in the full datatset 
load('full_dataset_GALAH_CANON.Rdata')
df
#and trim by volume
volume_beddell <- read.csv("~/Dropbox/GALAH/dynamics.csv")
volume_beddell <- tibble(volume_beddell)
head(volume_beddell)
beddell_R <- range(volume_beddell$R)
beddell_R
beddell_z <- range(volume_beddell$z)
beddell_z
df <- trim_by_volume(df[[1]], R = beddell_R, z = beddell_z, error = df[[2]])#next we want to sample by volume
df




tree_build_full <- function(tree_ID) {
  t <- tree_build(traits = MAP_traits, data = df, id_col = 'star_id', is.Mahalanobis.dist = F, scaler = 3)
  print(tree_ID)
  return(t)
}
# tree_build_full(1)

# require(pbmcapply)
# trees <- pbmclapply(paste0("Tree_", seq(1,1000,1)), tree_build_full, mc.cores = 4)
# class(trees) <- 'multiPhylo'
# save(trees, file = 'map.trees')

load('map.trees')

# MAJ_tree <- consensus.edges(trees, consensus.tree = consensus(trees, p = 0.5))
# MAJ_tree$clade_support <- prop.clades(MAJ_tree,  trees)
# write.tree(MAJ_tree, 'MAJ.tree')
MAJ_tree <- read.tree('MAJ.tree')


# MCC_tree <- maxCladeCred(trees, rooted = F)
# MCC_tree$clade_support <- prop.clades(MCC_tree, trees)
# write.tree(MCC_tree, 'MCC.tree')
MCC_tree <- read.tree('MCC.tree')
plot(MCC_tree, show.tip.label = F, type = 'fan')
plot(MAJ_tree, show.tip.label = F, type = 'fan')

#plot(cophylo(MAJ_tree, MCC_tree))
clade_suppot <- prop.clades(MCC_tree, trees)
require(castor)
dist_to_root <- get_all_distances_to_root(MCC_tree)[-(1:length(MCC_tree$tip.label))]
clade_suppot <- tibble(clade_suppot = clade_suppot/1000, node_depth = dist_to_root)

pxx <- ggplot(clade_suppot, aes(x = node_depth, y = clade_suppot)) +
  geom_point() +
  geom_smooth(method = 'lm', color = 'grey') +
  ylab('Clade support %') +
  xlab('Distance from root') +
  theme_classic()
pxx

lm <- lm(node_depth ~ clade_suppot, data = clade_suppot)
summary(lm)

ggsave('MCC_distance_from_root_clade_support.pdf', pxx, scale = 0.8)



og_df <- read.csv("~/Dropbox/GALAH/GALAH_the_canonII/TheCannon_solar_twin_full_set_predict_14xfe_39548spectra_wVAC.csv")
require(ggtree)

#MCC_tree <- read.tree('MCC.tree')

#MAJ_tree$edge.length <-  MAJ_tree$edge.length[-length(MAJ_tree$edge.length)]
#MCC_tree$edge.length <-  MCC_tree$edge.length[-length(MCC_tree$edge.length)]

ages <- tibble(label = og_df$star_id, ages = og_df$age)
ages
agesMAJ <- ages[which(ages$label %in% MAJ_tree$tip.label == T),]
agesMCC <- ages[which(ages$label %in% MCC_tree$tip.label == T),]


MAJ_tree_join <- full_join(MAJ_tree, agesMAJ, by = "label") 
MCC_tree_join <- full_join(MCC_tree, agesMCC, by = "label") 
colours <- c('#A62139', '#BF3B5E', '#D96C89', '#AD8C80', '#BDBFBF', '#5D5953')

tp1 <- ggtree(MAJ_tree_join, aes(color = ages), ladderize = T)+
  scale_colour_gradientn(
    colours = rev(colours)
  )  +
  geom_tippoint(size = .5)+
  labs(color = "Star ages (Gy)", cex = 20) +
  #layout_inward_circular(xlim=5) +
  theme(legend.key.width= unit(2, 'cm')) +
  scale_x_reverse() +
  theme(legend.position = 'bottom')  + ylim(-10, 150)
tp1  

tp2 <- ggtree(MCC_tree_join, aes(color = ages), ladderize = T)+
  scale_colour_gradientn(name = 'Gy',
    colours = rev(colours), 
  )  +  geom_tippoint(size = .5)+
  labs(color = "Star ages (Gy)", cex = 20) +
  #layout_inward_circular(xlim=5) +
  theme(legend.key.height= unit(2, 'cm')) + ylim(-10, 150)
tp2  

median_error <- median(na.omit(unlist(c(df[[2]]))))
median_error

MCC_trunc <- di2multi(MCC_tree, tol = median_error)
MCC_trunc_join <- full_join(MCC_trunc, agesMCC, by = "label") 
write.tree(MCC_trunc, 'MCCtrunc.tree')

tp3 <- ggtree(MCC_trunc_join, aes(color = ages), ladderize = T)+
  scale_colour_gradientn(
    colours = rev(colours)
  )  +
  geom_tippoint(size = .5)+
  labs(color = "Star ages (Gy)", cex = 20) +
 # layout_inward_circular(xlim=5) +
  theme(legend.key.height= unit(2, 'cm')) +
  theme(legend.position = 'bottom') + ylim(-10, 150)
 

tp3


MCC_nodes <- tibble(node = unique(MCC_tree$edge[,1]), support = as.numeric(MCC_tree$node.label),
                    age = unique(nodeHeights(MCC_tree)[,1]))

require(ggdist)
p6 <- ggplot(data = MCC_nodes, aes(x = support)) +
  geom_histogram(bins = 60, fill = gradient[3]) +
  xlab("Clade Support") +
  ylab('')+
  xlim(-0.01,1.1) +
  ylim(0, 35) +
  theme_classic()
p6

p7 <- ggplot() +
  geom_histogram(aes(x = as.numeric(MAJ_tree$node.label)), bins = 60, fill = gradient[3]) +
  xlab("Clade Support") +
  ylab('')+
  xlim(-0.01,1.1) +
  ylim(0, 35) +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank() )
p7


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(tp2)

#require(egg)

#ggsave(g2, file = 'Figures/MCC_and_MRC_trees.pdf', scale = 0.7)

require(ggpubr)

g3 <- ggarrange(tp3, tp1, nrow = 1, common.legend = T, legend = 'none', labels = c('MCC', 'MRC'))
g3
 

g4 <- ggarrange(g3, mylegend, ncol = 2, widths = c(7,1))
g4
ggsave(g4, file = 'Figures/MCC_and_MRC_trees.pdf', scale = 0.7)


tp3 + layout_inward_circular(xlim = 10)
tp1 + layout_inward_circular(xlim = 10)
g5 <- ggarrange(tp3 + layout_inward_circular(xlim = 10), 
                tp1 + layout_inward_circular(xlim = 10), 
                common.legend = T, legend = 'none', labels = c('a', 'b'))
g5

require(RPANDA)
MCCspect <- spectR(MCC_trunc)
MRCspect <- spectR(MAJ_tree)

SP_density <- tibble(Eigenvalues = c(MCCspect$eigenvalues, MRCspect$eigenvalues),
                     Tree = c(rep('MCC', length(MCCspect$eigenvalues)), 
                              rep('MRC', length(MRCspect$eigenvalues))))

SP_density

p1 <- ggplot(SP_density, aes(x = Eigenvalues)) +
  geom_density(fill = 'darkgrey') +
  theme_classic() +
  ylab("")
p1 <- p1 +facet_wrap(~Tree) + theme(plot.margin = margin(l = 40, r = 50))
p1



g6 <- ggarrange(p6, p7, nrow = 1) + theme(plot.margin = margin(l = 46, r = 50))
g6

g7 <- ggarrange(g3, p1, g6, nrow = 3, heights = c(5,1, 1), labels = c('a', 'b', 'c'))
g7

p0 <- ggplot() + theme_void()
p0
g0 <- ggarrange(mylegend, p0, p0, nrow = 3, heights = c(5,1,1))
g0
g8 <- ggarrange(g7, g0, ncol = 2, widths = c(8,1))
g8

ggsave('Figures/summary_trees_w_LPG_and_clade_support.pdf', g8, scale = 1)

bplot <- ggarrange(gplot1, gplot5, ncol = 1, heights = c(1.8,2), labels = c("d", 'e'), common.legend = T, legend = 'none')
bplot

bigplot <- ggarrange(g8, bplot, ncol = 2)
bigplot

ggsave('Figures/bigplot.pdf', bigplot, scale = 1)



### Okay some other stuff ###

og_df
ltt.plot(MCC_tree)
ltt.plot(MAJ_tree)



#metalicity
parameters <- tibble(label = og_df$star_id, ages = og_df$age, metalicity = og_df$fe_h)
parametersMAJ <- parameters[which(parameters$label %in% MAJ_tree$tip.label == T),]
parametersMCC <- parameters[which(parameters$label %in% MCC_tree$tip.label == T),]


MAJ_tree_parameters <- full_join(MAJ_tree, parametersMAJ, by = "label") 
MCC_tree_parameters <- full_join(MCC_tree, parametersMCC, by = "label") 


gradient <- c('#BD21BF', '#B527F2', '#B080F2', '#80C7F2', '#1BF2CB')



tp4 <- ggtree(MAJ_tree_parameters, aes(color = metalicity), ladderize = T)+
  scale_colour_gradientn(name = '[Fe/H]',
                         colours = rev(gradient), 
  )  +  geom_tippoint(size = .5)+
  labs(color = "Star ages (Gy)", cex = 20) +
  #layout_inward_circular(xlim=5) +
  scale_x_reverse() +
  theme(legend.key.width= unit(2, 'cm')) + ylim(-10, 150)
tp4

tp5 <- ggtree(MCC_tree_parameters, aes(color = metalicity), ladderize = T)+
  scale_colour_gradientn(name = '[Fe/H]',
                         colours = rev(gradient), 
  )  +  geom_tippoint(size = .5)+
  labs(color = "Star ages (Gy)", cex = 20) +

  #layout_inward_circular(xlim=5) +
  theme(legend.key.width = unit(2, 'cm')) + 
  ylim(-10, 150)
tp5


g4 <- ggarrange(tp5, tp4, common.legend = T, legend = 'none')
g4

p5 <- ggplot(MAJ_tree_parameters@data, aes(y = metalicity, x = ages, colour = metalicity)) +
  geom_point() +
  scale_colour_gradientn(name = '[Fe/H]',
                         colours = rev(gradient), 
  ) +
  xlab("Ages (Gy)") +
  ylab('[FE/H]') +
  theme_classic() +
  theme(legend.key.height = unit(2, 'cm')) 
p5
g5 <- ggarrange(g4, p5, nrow = 2, common.legend = T, legend = 'none', heights = c(4,1), labels = c('a', 'b'))
g5

mylegend<-g_legend(p5)

g6 <- ggarrange(g5, mylegend, ncol = 2, widths = c(8,1))
g6
ggsave('Figures/emp_trees_w_metaliticity.pdf', g6, scale = 0.8)


#see if we can relate the length of our tips to the error in the dataset
n<-length(MAJ_tree$tip.label)
ee<-setNames(MAJ_tree$edge.length[sapply(1:n,function(x,y)   which(y==x),y=MAJ_tree$edge[,2])],MAJ_tree$tip.label)
n<-length(MCC_tree$tip.label)
ef<-setNames(MCC_tree$edge.length[sapply(1:n,function(x,y)   which(y==x),y=MCC_tree$edge[,2])],MCC_tree$tip.label)
length(ee)
length(ee)
require(stringr)
error <- df[[2]]
error <- error[,which(str_remove(colnames(error), '_e') %in% MAP_traits == T)]
error_sums <- rowSums(error)

error_sums <- tibble(label = df[[1]]$star_id, sum_error = error_sums, MRC = ee[match(df[[1]]$star_id, names(ee))], 
                     MCC = ef[match(df[[1]]$star_id, names(ef))])
error_sums$Ages <- parameters$ages[match(error_sums$label, parameters$label)]

error_sums <-error_sums %>% 
  gather('MRC', 'MCC', key = Tree, value = tip_length)

error_sums



p8 <- ggplot(error_sums, aes(x = log(sum_error), y = log(tip_length), color = Ages)) +
  geom_point(alpha = 1) +
  geom_smooth(method = 'lm', color = 'black', alpha = .7) +
  scale_colour_gradientn(name = 'Gy',
                         colours = rev(colours[2:5]), 
  ) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  xlab('ln(Sum Error)') +
  ylab('ln(Tip Length)') +
  theme(legend.key.width= unit(2, 'cm')) 
p8 <- p8 +facet_wrap(~Tree)
p8
ggsave('Figures/tip_length_v_sum_error_scatter.pdf', p8, scale = .8)

#now lest get the residuals 
error_sums
lm <- lm(log(tip_length) ~ log(sum_error), data = error_sums)
summary(lm)
error_sums$residual <- exp(resid(lm))
error_sums
errorsMAJ <- subset(error_sums[which(error_sums$label %in% MAJ_tree$tip.label == T),], Tree == 'MRC')
errorsMCC <- subset(error_sums[which(error_sums$label %in% MCC_tree$tip.label == T),], Tree == 'MCC')
require(phytools)
MAJ_tree_error <- full_join(MAJ_tree, errorsMAJ, by = "label") 
MCC_tree_error <- full_join(MCC_tree, errorsMCC, by = "label") 


new_gradient <- c( '#103778', '#0593A2', '#E3371E')

tp6 <- ggtree(MCC_tree_error, aes(color = log(residual)), ladderize = T)+
  scale_colour_gradientn(name = 'Resid.',
                         colours = (new_gradient), 
  )  +  geom_tippoint(size = .5)+
  labs(color = "Star ages (Gy)", cex = 20) +
  #layout_inward_circular(xlim=5) +
  #scale_x_reverse() +
  theme(legend.key.width= unit(2, 'cm')) + ylim(-10, 150)
tp6

tp7 <- ggtree(MAJ_tree_error, aes(color = log(residual)), ladderize = T)+
  scale_colour_gradientn(name = 'ln(Resid.)',
                         colours = (new_gradient), 
  )  +  geom_tippoint(size = .5)+
  labs(color = "Star ages (Gy)", cex = 20) +
  #layout_inward_circular(xlim=5) +
  scale_x_reverse() +
  theme(legend.key.height= unit(2, 'cm')) + ylim(-10, 150)
tp7


mylegend<-g_legend(tp7)


p8 <- ggplot(error_sums, aes(x = log(sum_error), y = log(tip_length), color = log(residual))) +
  geom_point(alpha = .7) +
  geom_smooth(method = 'lm', color = 'black', alpha = .7) +
  scale_colour_gradientn(name = 'ln(Resid.)',
                         colours = (new_gradient),
  ) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('ln(Sum Error)') +
  ylab('ln(Tip Length)') +
  theme(legend.key.width= unit(2, 'cm')) 
p8 <- p8 +facet_wrap(~Tree)
p8
ggsave('Figures/tip_length_v_sum_error_scatter.pdf', p8, scale = .8)

g9 <- ggarrange(tp6, tp7, common.legend = T, legend = 'none')
g9

g10 <- ggarrange(g9, p8, nrow = 2, heights = c(4,2), labels = c("a", 'b'))
g10

g11 <- ggarrange(g10, mylegend, ncol = 2, widths = c(8,1))
g11

ggsave('Figures/residual_log_error_x_log_edge_length.pdf', g11, scale = 0.8)

error_sums$label[order(error_sums$residual, decreasing = T)][1:10]


#### Chemical Evolution in Distinct Clades #### 

MCC_trunc_play <- MCC_trunc #we need to make a version of our tree where we can actually see the nodelables
MCC_trunc_play$edge.length <- 1/MCC_trunc_play$edge.length
plot(ladderize(MCC_trunc_play), show.tip.label = F)
nodelabels()

#looks like we can spit int clades at nodes 149, 250, 227, and 243
require(geiger)
cladeA <- node.leaves(MCC_trunc, 251)
cladeB <- node.leaves(MCC_trunc, 264)
cladeC <- node.leaves(MCC_trunc, 241)
cladeD <- node.leaves(MCC_trunc, 226)
cladeE <- node.leaves(MCC_trunc, 206)
cladeF <- node.leaves(MCC_trunc, 152)



#this is the one star that isn't included in those clades
MCC_trunc$tip.label[which(MCC_trunc$tip.label %in% c(cladeA, cladeB, cladeC, cladeD, cladeE, cladeF) == F)]
parameters <- parameters[which(parameters$label %in% MCC_tree$tip.label == T),]

parameters$clade[which(parameters$label %in% cladeA == T)] <- 'F'
parameters$clade[which(parameters$label %in% cladeB == T)] <- 'E'
parameters$clade[which(parameters$label %in% cladeC == T)] <- 'D'
parameters$clade[which(parameters$label %in% cladeD == T)] <- 'C'
parameters$clade[which(parameters$label %in% cladeE == T)] <- 'B'
parameters$clade[which(parameters$label %in% cladeF == T)] <- 'A'

parameters$clade <- as.factor(parameters$clade)


MCC_clades <- full_join(MCC_trunc, parameters, by = "label") 


tp8 <- ggtree(MCC_clades, aes(color = clade), ladderize = T)+
  scale_colour_manual(name = 'Clade',
                         values = c(colours[-2], 'black')
  )  +
  geom_tippoint(size = .5)+
  #labs(color = "Star ages (Gy)", cex = 20) +
  #layout_inward_circular(xlim=5) +
  #scale_x_reverse() +
  theme(legend.position = 'bottom') + ylim(-10, 150)
tp8

parameters$clade

p9 <- ggplot(na.omit(parameters), aes(x = metalicity, y = ages, color = clade)) +
  geom_point() +
  scale_colour_manual(name = 'Clade',
                      values = c(colours[-2], 'black')
  )  +
  geom_smooth(method = 'lm') +
  theme_classic() +
  ylab('Ages(Gy)')+
  xlab('[Fe/H]') +
  theme(legend.position = 'bottom')
p9 <- p9 + facet_wrap(~ clade, nrow = 6)
p9

clades <- levels(parameters$clade)

coefficients <- tibble(intercept = vector(), coefficient = vector(),
                       intercept_se = vector(), coefficient_se = vector())
for(i in 1:length(clades)) {
  lm <- lm(ages ~ metalicity, subset(parameters, clade == clades[i]))
  coefficients <- rbind(coefficients, c(coefficients(lm)[[1]], coefficients(lm)[[2]],
                                        sqrt(diag(vcov(lm)))[1], sqrt(diag(vcov(lm)))[2]))
  print(clades[i])
}

coefficients$clades <- clades

library(sjPlot)
library(sjmisc)
library(sjlabelled)

lm <- lm(ages ~ metalicity, parameters)
lmA <- lm(ages ~ metalicity, subset(parameters, clade == 'A'))
lmB <- lm(ages ~ metalicity, subset(parameters, clade == 'B'))
lmC <- lm(ages ~ metalicity, subset(parameters, clade == 'C'))
lmD <- lm(ages ~ metalicity, subset(parameters, clade == 'D'))
lmE <- lm(ages ~ metalicity, subset(parameters, clade == 'E'))
lmF <- lm(ages ~ metalicity, subset(parameters, clade == 'F'))

tab_model(lm, file = 'model_clades.html')
tab_model(lmA, file = 'Amodel_clades.html')
tab_model(lmB, file = 'Bmodel_clades.html')
tab_model(lmC, file = 'Cmodel_clades.html')
tab_model(lmD, file = 'Dmodel_clades.html')
tab_model(lmE, file = 'Emodel_clades.html')
tab_model(lmF, file = 'Fmodel_clades.html')



colnames(coefficients) <- c('intercept', 'coefficient', 'intercent_se', 'coefficient_se', 'clade')

p10 <- p9 + geom_text(
  data    = coefficients,
  mapping = aes(x = -0.2, y = 6, label = paste0(expression('beta =='), round(coefficient, 2))), parse = T
) +
  geom_text(
    data    = coefficients,
    mapping = aes(x = -0.2, y = 4.5, label = paste0(expression('sigma =='), round(coefficient_se, 2))), parse = T,
  ) +
  theme(legend.position = 'none')
p10


g12 <- ggarrange(tp8, p10, common.legend = T, legend = 'bottom', ncol = 2, widths = c(2,1))
g12

ggsave('Figures/clades_age_metalicity.pdf', g12, scale = .8)

# 
# pgls_parameters <- drop_na(parameters)
# pgls_parameters <- pgls_parameters[which(pgls_parameters$label %in% MCC_tree$tip.label == T),]
# 
# pgls_model <- gls(ages ~ metalicity, data = pgls_parameters, correlation = corBrownian(1, phy = keep.tip(MCC_tree, pgls_parameters$label), form = ~label))
# summary(pgls_model)
# coef(pgls_model)
# 
# plot(pgls_parameters$ages ~ pgls_parameters$metalicity)
# abline(a = coef(pgls_model)[1], b = coef(pgls_model)[2])
# 
# pgls_pagel <- gls(ages ~ metalicity, data = pgls_parameters, correlation = corPagel(1, phy = keep.tip(MCC_tree, pgls_parameters$label), form = ~label))
# summary(pgls_pagel)

b1 <- tp4 + geom_tiplab(hjust = 1) +scale_x_reverse() +theme(legend.position = 'bottom') + theme(legend.key.width= unit(2, 'cm'))
b1
ggsave('Figures/MRC_w_metalicity_large.pdf', b1, width = 15, height = 40, scale = 1.3, limitsize = F)

b2 <- tp5 + geom_tiplab(hjust = ) +theme(legend.position = 'bottom')  + theme(legend.key.width= unit(2, 'cm'))
b2
ggsave('Figures/MCC_w_metalicity_large.pdf', b2, width = 15, height = 40, scale = 1.3, limitsize = F)


b3 <- tp1 + geom_tiplab() +theme(legend.position = 'bottom') + theme(legend.key.width= unit(2, 'cm'))
b3
ggsave('Figures/MRC_w_age_large.pdf', b3, width = 15, height = 40, scale = 1.3, limitsize = F)

b4 <- tp2 + geom_tiplab(hjust = ) +theme(legend.position = 'bottom') + theme(legend.key.width= unit(2, 'cm')) 
b4
ggsave('Figures/MCC_w_age_large.pdf', b4, width = 15, height = 40, scale = 1.3, limitsize = F)

b5 <- tp6 + geom_tiplab() +theme(legend.position = 'bottom') + theme(legend.key.width= unit(2, 'cm'))
b5
ggsave('Figures/MCC_w_resid_large.pdf', b5, width = 15, height = 40, scale = 1.3, limitsize = F)

b6 <- tp7 + geom_tiplab(hjust = 1) +theme(legend.position = 'bottom') + theme(legend.key.width= unit(2, 'cm'))
b6
ggsave('Figures/MRC_w_resid_large.pdf', b6, width = 15, height = 40, scale = 1.3, limitsize = F)
 
