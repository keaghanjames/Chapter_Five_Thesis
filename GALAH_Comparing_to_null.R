setwd("~/Dropbox/GALAH/") 
source("functions.R")
setwd("GALAH_the_canonII")
colours <- c('#26151F', '#5872A6', '#89B3D9', '#732A19', '#D99282')

gradient <- c("#EE05F2", '#6B1C8C', '#011126', '#143D73', '#0BBFBF')

require(phytools)
MRCtree <- read.tree('MAJ.tree')
MCCtree <- read.tree('MCCtrunc.tree')

require(RPANDA)
spectMRC <- spectR(MRCtree)
spectMCC <- spectR(MCCtree)


require(foreach)
require(doParallel)

summary <- tibble(principal_eigenvalue = vector(),
                  asymmetry = vector(),
                  peakedness = vector(),
                  eigengap.modalities = vector(),
                  Tree = vector(),
                  Model = vector(),
                  Run = vector())

for(i in 1:100){

  null_MRC <- read.tree(file = paste0('NULL/',i,'consensus.tree'))
  spect_null_MRC <- spectR(null_MRC)
  spect_null_MRC <- unlist(spect_null_MRC[2:5])
  spect_null_MRC <- c(spect_null_MRC, setNames(c('MRC', 'NULL'), nm = c('Tree', 'Model')))

  null_MCC <- read.tree(file = paste0('NULL/',i,'mcc.tree'))
  spect_null_MCC <- spectR(null_MCC)
  spect_null_MCC <- unlist(spect_null_MCC[2:5])
  spect_null_MCC <- c(spect_null_MCC, setNames(c('MCC', 'NULL'), nm = c('Tree', 'Model')))

  pop_MRC <- read.tree(file = paste0('POP_NULL/',i,'consensus.tree'))
  spect_pop_MRC <- spectR(pop_MRC)
  spect_pop_MRC <- unlist(spect_pop_MRC[2:5])
  spect_pop_MRC <- c(spect_pop_MRC, setNames(c('MRC', 'POP'), nm = c('Tree', 'Model')))

  pop_MCC <- read.tree(file = paste0('POP_NULL/',i,'mcc.tree'))
  spect_pop_MCC <- spectR(pop_MCC)
  spect_pop_MCC <- unlist(spect_pop_MCC[2:5])
  spect_pop_MCC <- c(spect_pop_MCC, setNames(c('MCC', 'POP'), nm = c('Tree', 'Model')))

  spect_sum <-  bind_rows(spect_null_MRC,
                spect_null_MCC,
                spect_pop_MRC,
                spect_pop_MCC)
  spect_sum$Run <- i
  summary <- rbind(summary, spect_sum)
  print(i)

}






write.csv(summary, 'spectral_null_summary.csv', row.names = F)
summary <- read.csv('spectral_null_summary.csv')
summary

#make and example of all our nulls


null_MRC <- read.tree(file = paste0('NULL/1consensus.tree'))
null_MCC <- read.tree(file = paste0('NULL/1mcc.tree'))
pop_MRC <- read.tree(file = paste0('POP_NULL/1consensus.tree'))
pop_MCC <- read.tree(file = paste0('POP_NULL/1mcc.tree'))

#to keep it readable we will trim to 30 tips

tips_to_keep <- sample(null_MCC$tip.label, 70)
null_MRC <- keep.tip(null_MRC, tips_to_keep)
null_MCC <- keep.tip(null_MCC, tips_to_keep)
pop_MRC <- keep.tip(pop_MRC, tips_to_keep)
pop_MCC <- keep.tip(pop_MCC, tips_to_keep)

require(ggtree)
nMRC <- ggtree(null_MRC, ladderize = T, color = colours[2])
nMRC
nMCC <- ggtree(null_MCC, ladderize = T, color = colours[3]) +
  scale_x_reverse()
nMCC

pMRC <- ggtree(pop_MRC, ladderize = T, color = colours[4])
pMRC
pMCC <- ggtree(pop_MCC, ladderize = T, color = colours[5]) +
  scale_x_reverse()
pMCC

require(ggpubr)

g1 <- ggarrange(nMRC, nMCC, pMRC, pMCC, nrow = 2, ncol = 2, labels = c('a', 'b', 'c', 'd'))
g1

require(RPANDA)
Spnull_MRC <- spectR(null_MRC)
Spnull_MCC <- spectR(null_MCC)
Sppop_MRC <- spectR(pop_MRC)
Sppop_MCC <- spectR(pop_MCC)

Sp_density <- tibble(eigenvalues = c(Spnull_MRC$eigenvalues,Spnull_MCC$eigenvalues,
                    Sppop_MRC$eigenvalues, Sppop_MCC$eigenvalues), 
                    model = c(rep('NULL', length(c(Spnull_MRC$eigenvalues,Spnull_MCC$eigenvalues))),
                              rep('POP', length(c(Sppop_MRC$eigenvalues,Sppop_MCC$eigenvalues))
                              )
                    ), tree = c(rep('MRC', length(Spnull_MRC$eigenvalues)),
                                rep('MCC', length(Spnull_MCC$eigenvalues)),
                                rep('MRC', length(Sppop_MRC$eigenvalues)),
                                rep('MCC', length(Sppop_MCC$eigenvalues))
                    )
)

Sp_density$Treatment <- paste0(Sp_density$model, '/', Sp_density$tree)

Sp_density$Treatment <- factor(Sp_density$Treatment, ordered = T, 
                               levels = c("NULL/MRC", "NULL/MCC", "POP/MRC", "POP/MCC"))

p1 <- ggplot(Sp_density, aes(x = eigenvalues, colour = Treatment, fill = Treatment)) +
  geom_density(alpha = 0.1) +
  theme_classic() +
  scale_color_manual(values = colours[2:5]) +
  scale_fill_manual(values = colours[2:5]) +
  ylab("") +
  xlab("Eigenvalues") +
  theme(legend.position = 'bottom', axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), plot.margin = margin(50,0,50,0))
p1 <- p1 + facet_wrap(~ Treatment)
p1

g2.1 <- ggarrange(nMRC, pMRC,  nrow = 2)
g2.1
g2.2 <- ggarrange(nMCC, pMCC,  nrow = 2)
g2.2

g2 <- ggarrange(g2.1, p1, g2.2, nrow = 1, common.legend = T, legend = 'none')
g2

ggsave('Figures/null_model_trees_with_SD_plots.pdf', g2, scale = 0.8)


### okay now can systematically compare our true trees to nulls

require(ggdist)

summary$principal_eigenvalue <- as.numeric(summary$principal_eigenvalue)
summary$asymmetry <- as.numeric(summary$asymmetry)
summary$peakedness <- as.numeric(summary$peakedness)
summary$eigengap.modalities <- as.numeric(summary$eigengap.modalities)
summary$Tree <- as.factor(summary$Tree)
summary$Model <- as.factor(summary$Model)
summary$Run <- as.factor(summary$Run)

require(ggridges)

p1 <- ggplot(summary, aes(x = peakedness, y = Model, fill = Tree, color = Tree)) +
 # stat_eye() +
  geom_density_ridges(quantile_lines = F, alpha = 0.7,
                      jittered_points = TRUE,
                      point_size = 0.1, point_alpha = 0.7,
                      position = position_raincloud())+ 
  theme_classic() +
  scale_fill_manual(values = c(colours[3], colours[5])) +
  scale_color_manual(values = c(colours[3], colours[5])) +
  geom_vline(xintercept = spectMRC$peakedness, color = colours[5], 
             linetype = 'dashed', size = 1.2) +
  geom_vline(xintercept = spectMCC$peakedness, color = colours[3], 
             linetype = 'dashed', size = 1.2)  +
  ylab("") +
  xlab("Peakedness")
p1

p2 <- ggplot(summary, aes(x = asymmetry, y = Model, fill = Tree, color = Tree)) +
  geom_density_ridges(quantile_lines = F, alpha = 0.7,
                      jittered_points = TRUE,
                      point_size = 0.1, point_alpha = 0.7,
                      position = position_raincloud())+ 
  theme_classic() +
  scale_fill_manual(values = c(colours[3], colours[5])) +
  scale_color_manual(values = c(colours[3], colours[5])) +
  geom_vline(xintercept = spectMRC$asymmetry, color = colours[5], 
             linetype = 'dashed', size = 1.2) +
  geom_vline(xintercept = spectMCC$asymmetry, color = colours[3], 
             linetype = 'dashed', size = 1.2)  +
  ylab("") +
  xlab("Asymmetry")
p2


p3 <- ggplot(summary, aes(x = principal_eigenvalue, y = Model, fill = Tree, color = Tree)) +
  geom_density_ridges(quantile_lines = F, alpha = 0.7,
                      jittered_points = TRUE,
                      point_size = 0.1, point_alpha = 0.7,
                      position = position_raincloud())+ 
  theme_classic() +
  scale_fill_manual(values = c(colours[3], colours[5])) +
  scale_color_manual(values = c(colours[3], colours[5])) +
  geom_vline(xintercept = spectMRC$principal_eigenvalue, color = colours[5], 
             linetype = 'dashed', size = 1.2) +
  geom_vline(xintercept = spectMCC$principal_eigenvalue, color = colours[3], 
             linetype = 'dashed', size = 1.2)  +
  ylab("") +
  xlab("Principal Eigenvalue")
p3

require(ggpubr)

g1 <- ggarrange(p1,p2,p3, common.legend = T, legend = 'bottom', nrow = 3)
g1

ggsave('Figures/spectral_density_parameters_null_v_empirical.pdf', g1, scale = .8)

### We want a table that summarises means and 95% CIs for each parameter under each combination of tree and model ###


summary <- as.tibble(summary)
summary
sum_stats <- summary %>%
  group_by(Model, Tree) %>%
  summarise(median_eigenvalue = median(principal_eigenvalue), eigenvalue0.05 = quantile(principal_eigenvalue, 0.05),
             eigenvalue0.95 = quantile(principal_eigenvalue, 0.95),
            median_asymmetry = median(asymmetry), asymmetry0.05 = quantile(asymmetry, 0.05),
            asymmetry0.95 = quantile(asymmetry, 0.95),
            peakedness = median(peakedness), peakedness0.05 = quantile(peakedness, 0.05),
            peakedness0.95 = quantile(peakedness, 0.95),
            modalities = median(eigengap.modalities), modalities0.05 = quantile(eigengap.modalities, 0.05),
            modalities0.95 = quantile(eigengap.modalities, 0.95)
            )
sum_stats

write.csv(sum_stats, 'summary_LSDP.csv', row.names = F)


  