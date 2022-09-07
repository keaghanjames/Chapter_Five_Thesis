#need to combine tracefiles I guesds
require(tidyr)
require(dplyr)
setwd('~/Dropbox/GALAH/GALAH_the_canonII/')
df1 <- read.csv('chain_1/trace2.csv', row.names = NULL)
df1$chain = 1
# df2 <- read.csv('chain_2/trace2.csv', row.names = NULL)
# df2$chain <- 2

#colours <- c('#26151F', '#5872A6', '#89B3D9', '#732A19', '#D99282')
colours <- c('#A62139', '#BF3B5E', '#D96C89', '#AD8C80', '#BDBFBF', '#5D5953')


df <- rbind(df1)
colnames(df)[1] <- 'Generations'
df$chain <- as.factor(df$chain)

require(ggplot2)


#we inspect our chains

g1 <- ggplot(data = df, aes(x = Generations, y = clade_support, color = chain)) +
  geom_line() +
  scale_color_manual(values = colours[-2], name = 'Chain') +
  ylab('Sum(Clade Support)') +
  theme_classic()+
  annotate("rect", xmin = -10, xmax = 1000, ymin = 0, ymax = max(df$clade_support),
           alpha = .6,fill = "white")
g1

g2 <- ggplot(data = df, aes(x = Generations, y = Ntraits, color = chain), alpha = 0.5) +
  geom_line() +
  scale_color_manual(values = colours[-2], name = 'Chain') +
  theme_classic()+
  annotate("rect", xmin = -10, xmax = 1000, ymin = 0, ymax = max(df$Ntraits),
           alpha = .6,fill = "white")
g2

g3 <- ggplot(data = df, aes(x = Generations, y = Nnode, color = chain)) +
  geom_line() +
  scale_color_manual(values = colours[-2], name = 'Chain') +
  ylab('Nnodes') +
  theme_classic()+
  annotate("rect", xmin = -10, xmax = 1000, ymin = 0, ymax = max(df$Nnode),
           alpha = .6,fill = "white")
g3

g4 <- ggplot(data = df, aes(x = Generations, y = pois, color = chain)) +
  geom_line() +
  scale_color_manual(values = colours[-2], name = 'Chain') +
  ylab(expression(paste('Poisson (', lambda, ')'))) +
  theme_classic() +
  annotate("rect", xmin = -10, xmax = 1000, ymin = 0, ymax = max(df$pois),
           alpha = .6,fill = "white")
g4

require(ggdist)
d1 <- ggplot(data = df[-c(which(df$Generations <= 1000)),], aes(x = clade_support, y = chain, fill = chain)) +
  stat_eye() +
  scale_fill_manual(values = colours[-2]) +
  theme_classic() +
  ylab('Chain') +
  xlab('Sum(Clade support)') +
  theme(legend.position = 'none')
d1

d2 <- ggplot(data = df[-c(which(df$Generations <= 1000)),], aes(x = Ntraits, y = chain, fill = chain)) +
  stat_eye() +
  scale_fill_manual(values = colours[-2]) +
  theme_classic() +
  ylab('Chain') +
  theme(legend.position = 'none')
d2

d3 <- ggplot(data = df[-c(which(df$Generations <= 1000)),], aes(x = Nnode, y = chain, fill = chain)) +
  stat_eye() +
  scale_fill_manual(values = colours[-2]) +
  theme_classic() +
  ylab('Chain') +
  theme(legend.position = 'none')
d3


pp1 <-ggplot(data = df[-c(which(df$Generations <= 1000)),], aes(x = Nnode, y = Ntraits, color = chain))+
  geom_point(alpha = .3, position = position_jitter()) +
  scale_color_manual(values = colours) +
  geom_smooth(method = 'lm') +
  theme_classic() +
  theme(legend.position = 'none')
pp1

pp2 <-ggplot(data = df[-c(which(df$Generations <= 1000)),], aes(x = clade_support, y = Ntraits, color = chain))+
  geom_point(alpha = .3) +
  scale_color_manual(values = colours) +
  geom_smooth(method = 'lm') +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('Sum(Clade Support)')
pp2

pp2 <-ggplot(data = df[-c(which(df$Generations <= 1000)),], aes(x = clade_support, y = Ntraits))+
  geom_point(alpha = .1, size = .2) +
  scale_color_manual(values = colours[2:5], name = 'Chain') +
  geom_smooth(method = 'lm', color = 'black') +
  theme_classic() +
  theme(legend.position = 'bottom') +
  xlab('Sum(Clade Support)') 
pp2

coefficients <- tibble(intercept = vector(), coef = vector())
  lm <- lm(clade_support ~ Ntraits, df[c(which(df$Generations > 1000)),])
  coefficients <- rbind(coefficients, c(coefficients(lm)[[1]], coefficients(lm)[[2]]))
  #print(clade[i])


colnames(coefficients) <- c('intercept', 'coef')
coefficients$chain <- as.factor(c(1:4))
summary(lm)

.75ppp22 <- pp2 + geom_text(
  data    = coefficients,
  mapping = aes(x = -Inf, y = -Inf, label = paste0(expression('beta =='), round(coef, 2))), parse = T,
  hjust   = -0.1,
  vjust   = -0.1 
) + theme(legend.position = 'none') 
ppp22
ggsave


pp3 <-ggplot(data = df[-c(which(df$Generations <= 1000)),], aes(x = clade_support, y = Nnode, color = chain))+
  geom_point(alpha = .3) +
  scale_color_manual(values = colours) +
  geom_smooth(method = 'lm') +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('Sum(Clade Support)')
pp3

# ### CALCULATE NPD SCORES ###
# 
# # # we inspect our trees
# require(phytools)
# load('chain_1/GALAHtrees2.Rdata')
# 
# trees1 <- trees[-c(1:1000)]
# class(trees1) <- 'multiPhylo'
# trees1
# 
# load('chain_2/GALAHtrees2.Rdata')
# trees2 <- trees[-c(1:1000)]
# class(trees2) <- 'multiPhylo'
# trees2
# 
# require(phangorn)
# source('~/Dropbox/GALAH/NPD_Naser-Khador.R')
# 
# normalised.path.dist(trees2[[1]], trees2[[2]])
# 
# # #we want to get NPD score for all pairs of trees
# #
# # #we are creating a handler function so we can run our pair-wise NPD analysis in paralele
# NPD_pair <- function(pair, trees){
#             NPD <- normalised.path.dist(trees[[pair[1]]], trees[[pair[2]]])
#             return(NPD)
# }
# 
# 
# # #need to set which chain we're working with
# trees <- trees2
# #
# # #now we need to create a list of all possible pairs
# pair_list <- combn(seq(1, length(trees)), 2, simplify = F)
# #and its insanely huge so we will take a sample of 1000
# pair_list <- pair_list[sample(seq(1, length(pair_list)), 1000)]
# pair_list <- do.call(rbind, pair_list)
# head(pair_list)
# # #now we make another handler where we can specify the chain
# NPD_pair_handler <- function(pair){
#   NPD <- NPD_pair(pair = pair, trees = trees)
#   return(NPD)
# }
# 
# NPD_pair_handler <- function(run){
#   NPD <- NPD_pair(pair = pair_list[run,], trees = trees)
#   return(NPD)
# }
# 
# NPD_pair_handler(5)
# 
# require(pbmcapply)
# NPD_scores <- pbmclapply(seq(1,1000,1), NPD_pair_handler, mc.cores = 4)
# 
# get_sum_clade <- function(pair_list){
#   clade_support1 <- vector()
#   clade_support2 <- vector()
#   Nnode1 <- vector()
#   Nnode2 <- vector()
#     for(i in 1:nrow(pair_list)){
#     clade_support1 <- c(clade_support1, sum(trees[[pair_list[i,1]]]$clade_support))
#     clade_support2 <- c(clade_support2, sum(trees[[pair_list[i,2]]]$clade_support))
#     Nnode1 <- c(Nnode1, trees[[pair_list[i,1]]]$Nnode)
#     Nnode2 <- c(Nnode2, trees[[pair_list[i,2]]]$Nnode)
#         }
#   return(tibble(clade_support1, clade_support2, Nnode1, Nnode2))
# }
# 
# 
# NPD_scores <- tibble(chain = 1, NPD_scores = unlist(NPD_scores),
#                       sum_clade_support = get_sum_clade(pair_list))
# NPD_scores <- bind_cols(NPD_scores1, pair_list)
# NPD_scores
# 
# # #lets repeat that processes with each chain
# # trees <- trees2
# # pair_list <- combn(seq(1, length(trees)), 2, simplify = F)
# # pair_list <- pair_list[sample(seq(1, length(pair_list)), 1000)]
# # pair_list <- do.call(rbind, pair_list)
# # head(pair_list)
# # NPD_scores <- pbmclapply(seq(1,1000,1), NPD_pair_handler, mc.cores = 4)
# # NPD_scores2 <- tibble(chain = 2, NPD_scores = unlist(NPD_scores),
# #                       sum_clade_support = get_sum_clade(pair_list))
# # NPD_scores2 <- bind_cols(NPD_scores2, pair_list)
# # NPD_scores2
# 
# 
# # NPD_scores <- bind_rows(NPD_scores1, NPD_scores2, NPD_scores3, NPD_scores4)
# # NPD_scores <- tibble(NPD_scores$chain, NPD_scores$NPD_scores, NPD_scores$sum_clade_support$clade_support1,
# #                      NPD_scores$sum_clade_support$clade_support2,
# #                      NPD_scores$sum_clade_support$Nnode1,
# #                      NPD_scores$sum_clade_support$Nnode2,
# #                      NPD_scores$...4,
# #                      NPD_scores$...5)
# colnames(NPD_scores) <- c('Chain', 'NPD', 'Clade_support_t1',
#                            'Clade_support_t2', 'Nnodes_t1', 'Nnodes_t2',
#                            't1', 't2')
# 
# # NPD_scores$Chain2 <- 0
# # NPD_scores$Chain2[which(NPD_scores$Chain == 2)] <- 4
# # NPD_scores$Chain2[which(NPD_scores$Chain == 3)] <- 3
# # NPD_scores$Chain2[which(NPD_scores$Chain == 4)] <- 2
# # NPD_scores$Chain <- NPD_scores$Chain2
# # NPD_scores$Chain2 <- NULL
# write.csv(NPD_scores, 'NPD_scores.csv')

#### NPD ANALYSIS ###
NPD_scores <- read.csv('NPD_scores.csv')
NPD_scores <- as_tibble(NPD_scores)
NPD_scores$X <- NULL
NPD_scores$Chain <- as.factor(NPD_scores$Chain)
require(ggridges)
d4 <- ggplot(NPD_scores, aes(x = NPD, y = Chain, fill = Chain)) +
  stat_eye() +
  scale_fill_manual(values = colours) +
  geom_vline(xintercept = 1, linetype = 'dashed', size = 1) +
  geom_vline(xintercept = 0.5, linetype = 'dashed', size = 1, color = 'grey') +
  theme_classic() +
  ylab("") +
  xlab('Normalized Path-Distance') +
  theme(legend.position = 'none', axis.text.y = element_blank())
d4

p1 <- ggplot(NPD_scores, aes(x = NPD, y = abs(Clade_support_t1 - Clade_support_t2), col = Chain)) +
  geom_point(alpha = 0.5, size = .7) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = colours) +
  geom_vline(xintercept = 1, linetype = 'dashed', size = 1) +
  ylab('Absolute difference sum(clade support)') +
  xlab('Normalized Path-Distance') +
  geom_vline(xintercept = 0.5, linetype = 'dashed', size = 1, color = 'grey') +
  theme_classic()  +
  theme(legend.position = 'none')
p1

p2 <- ggplot(NPD_scores, aes(x = NPD, y = abs(Nnodes_t1 - Nnodes_t2), col = Chain)) +
  geom_point(alpha = 0.5, size = .7) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = colours) +
  geom_vline(xintercept = 1, linetype = 'dashed', size = 1) +
  ylab('Absolute difference Nnodes') +
  xlab('Normalized Path-Distance') +
  geom_vline(xintercept = 0.5, linetype = 'dashed', size = 1, color = 'grey') +
  theme_classic()  +
  theme(legend.position = 'none')
p2


require(egg)
require(ggpubr)
gplot1 <- ggarrange(g1,g2,g3,g4, ncol = 2, nrow = 2, common.legend = T, legend = 'none')
gplot1
ggsave('Figures/chains.pdf', gplot1, scale = .8)


gplot2 <- ggarrange(p1, p2, ncol = 2, labels = c('b', 'c'))
gplot2

gplot3 <- ggarrange(d4, gplot2, ncol = 1, heights = c(1.5,2), labels = c('a', '', ''))
gplot3
ggsave('Figures/NPD_distribtions.pdf', gplot3, scale = .8)


ggplot6 <- ggarrange(d1, d2, d3, ncol = 1)
ggplot6

ggsave('Figures/chain_parameter_dists.pdf', ggplot6, scale = .8)


ggplot7 <- ggarrange(pp1, pp2, pp3, ncol = 1)
ggplot7

ggplot8 <- ggarrange(ggplot6, ggplot7, ncol = 2, widths = c(1.8,1), labels = c('a', 'b'), common.legend = T)
ggplot8
ggsave('Figures/chain_parameter_dists_cors.pdf', ggplot8, scale = .8)


### Traits Considered ###

load('chain_1/traits_considered.Rdata')
traits_considered1 <- traits_considered[-c(1:1000)]
# load('chain_2/traits_considered.Rdata')
# traits_considered2 <- traits_considered[-c(1:1000)]



#lets add the number of traits to our NPD scores
NPDchain_1 <- subset(NPD_scores, Chain == 1)
NPDchain_1$Ntraits_t1 <- unlist(lapply(traits_considered1[NPDchain_1$t1], nrow))
NPDchain_1$Ntraits_t2 <- unlist(lapply(traits_considered1[NPDchain_1$t2], nrow))

# NPDchain_2 <- subset(NPD_scores, Chain == 2)
# NPDchain_2$Ntraits_t1 <- unlist(lapply(traits_considered2[NPDchain_2$t1], nrow))
# NPDchain_2$Ntraits_t2 <- unlist(lapply(traits_considered2[NPDchain_2$t2], nrow))


NPD_scores_w_traits <- bind_rows(NPDchain_1)
                                 
                      
#okay now we can plot the absolute difference between the two
p3 <- ggplot(NPD_scores_w_traits, aes(x = NPD, y = abs(Ntraits_t1 - Ntraits_t2), col = Chain)) +
  geom_point(alpha = 0.5, size = .7) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = colours) +
  geom_vline(xintercept = 1, linetype = 'dashed', size = 1) +
  ylab('Absolute difference Ntraits') +
  xlab('Normalized Path-Distance') +
  geom_vline(xintercept = 0.5, linetype = 'dashed', size = 1, color = 'grey') +
  theme_classic()  +
  theme(legend.position = 'none')
p3

gplot4 <- ggarrange(p1, p2, p3, ncol = 3)
gplot4

gplot5 <- ggarrange(d4, gplot4, ncol = 1, heights = c(2.5,2))
gplot5
ggsave('Figures/NPD_distribtions_w_traits.pdf', gplot5, scale = .8)


# Okay back to looking at the elements and traits ###

elements1 <- unique(unlist(traits_considered1))
elements2 <- unique(unlist(traits_considered2))

count_elements1 <- summary(as.factor(unlist(traits_considered1)))
count_elements2 <- summary(as.factor(unlist(traits_considered2)))

count_elements1 <- tibble(element = names(count_elements1),
                          counts = count_elements1,
                          percentage = count_elements1/9000,
                          chain = 1)
# count_elements2 <- tibble(element = names(count_elements2),
#                           counts = count_elements2, 
#                           percentage = count_elements2/9000,
#                           chain = 2)

count_elements <- bind_rows(count_elements1)

count_elements$chain <- as.factor(count_elements$chain)
require(stringr)
count_elements$element <- str_to_title(count_elements$element)
p4 <- ggplot(count_elements, aes(x = element, y = counts, fill = chain))+
  geom_bar(stat = 'identity', position = position_dodge()) +
  scale_fill_manual(values = colours) +
  xlab('Element') +
  ylab('Times included') +
  theme_classic()
p4

ggsave('Figures/times_element_included.pdf', p4, scale = .8)

p5 <- ggplot(count_elements, aes(x = element, y = counts/9000, fill = chain))+
  geom_bar(stat = 'identity',  position = position_dodge()) +
  scale_fill_manual(values = colours) +
  xlab('Element') +
  ylab('Percentage visited') +
  #ylim(c(0,100)) +
  theme_classic() 
p5

factor_levels <- unique(c(elements1, elements2))

trait_list1 <- as_tibble(do.call(rbind,traits_considered1))
colnames(trait_list1) <- c('Element1', 'Element2')
trait_list1$Element1 <- factor(trait_list1$Element1, levels = factor_levels)                         
trait_list1$Element2 <- factor(trait_list1$Element2, levels = factor_levels)
trait_list1$Chain <- 1

# trait_list2 <- as_tibble(do.call(rbind,traits_considered2))
# colnames(trait_list2) <- c('Element1', 'Element2')
# trait_list2$Element1 <- factor(trait_list2$Element1, levels = factor_levels)                         
# trait_list2$Element2 <- factor(trait_list2$Element2, levels = factor_levels)
# trait_list2$Chain <- 2


trait_list <- bind_rows(trait_list1)



trait_counts <- trait_list %>% count(Element1, Element2, .drop = F)

# traits_counts <- bind_rows(trait_counts, tibble(Element1 = c('Ci', 'Euii'), Element2 = c('Ci', 'Euii'), n = c(0,0)))
# traits_counts$Element1 <- as.factor(trait_counts$Element1)


trait_counts$Element1 <- str_to_title(trait_counts$Element1)
trait_counts$Element2 <- str_to_title(trait_counts$Element2)
h1 <- ggplot(trait_counts, aes(x = Element1, y = Element2, fill = log(n))) +
  geom_tile() +
  scale_fill_gradient(low = 'grey', high = colours[1], na.value = 'white') +
  theme_classic() +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme(legend.position = 'bottom')
h1


trait_counts1 <- trait_list1 %>% count(Element1, Element2, .drop = F)
trait_counts1$Element1 <- str_to_title(trait_counts1$Element1)
trait_counts1$Element2 <- str_to_title(trait_counts1$Element2)
h1.1 <- ggplot(trait_counts1, aes(x = Element1, y = Element2, fill = log(n))) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = colours[2], na.value = 'white') +
  theme_classic() +
  xlab("") +
  ylab("")+
  theme(legend.position = 'bottom', axis.text.x = element_blank(), axis.text.y = element_blank())+
  coord_fixed()
h1.1

trait_counts2 <- trait_list2 %>% count(Element1, Element2, .drop = F)
trait_counts2$Element1 <- str_to_title(trait_counts2$Element1)
trait_counts2$Element2 <- str_to_title(trait_counts2$Element2)
h1.2 <- ggplot(trait_counts2, aes(x = Element1, y = Element2, fill = log(n))) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = colours[3], na.value = 'white') +
  theme_classic() +
  xlab("") +
  ylab("")+
  theme(legend.position = 'bottom', axis.text.x = element_blank(), axis.text.y = element_blank())+
  coord_fixed()
h1.2


gplot9 <- ggarrange(h1.1,h1.2, common.legend = T, legend = 'none')
gplot9

gplot10 <- ggarrange(h1, gplot9,  ncol = 2, widths = c(1.8,2), labels = c('a', 'b'))
gplot10

ggsave('Figures/heatmaps_traits.pdf', gplot10, scale = .8)

gplot11 <- ggarrange(gplot10, p4, nrow = 2, heights = c(2,1), labels = c("", "c"))
gplot11
ggsave('Figures/heatmaps_traits_bar_chart.pdf', gplot11, scale = 1)

# 
# melt <- melt(traits_counts, id.vars = c('Element1', 'Element2'), measure.vars = 'n', na.rm = T)
# 
# h2 <- ggplot(melt, aes(x = Element1, y = Element2, fill = log(value))) +
#   geom_tile() +
#   #scale_fill_gradient(low = 'white', high = colours[5], na.value = 'white') +
#   theme_classic() +
#   xlab("") +
#   ylab("")+
#   theme(legend.position = 'bottom', axis.text.x = element_blank(), axis.text.y = element_blank()) +
#   coord_fixed()
# h2

count(Pair = str_c(pmin(rowA, rowB), ' - ',
                   pmax(rowA, rowB)), name = "Count")

trait_list$Element1 <- as.character(trait_list$Element1)
trait_list$Element2 <- as.character(trait_list$Element2)

trait_countsX <- trait_list %>% 
  count(Pair = str_c(pmin(Element1, Element2), '-',
                     pmax(Element1, Element2)),
        name = 'count', .drop = F)
trait_countsX


trait_list1$Element1 <- as.character(trait_list1$Element1)
trait_list1$Element2 <- as.character(trait_list1$Element2)

trait_countsX1 <- trait_list1 %>% 
  count(Pair = str_c(pmin(Element1, Element2), '-',
                     pmax(Element1, Element2)),
        name = 'count', .drop = F)
trait_countsX1

trait_list2$Element1 <- as.character(trait_list2$Element1)
trait_list2$Element2 <- as.character(trait_list2$Element2)

trait_countsX2 <- trait_list2 %>% 
  count(Pair = str_c(pmin(Element1, Element2), '-',
                     pmax(Element1, Element2)),
        name = 'count', .drop = F)
trait_countsX2


trait_list3$Element1 <- as.character(trait_list3$Element1)
trait_list3$Element2 <- as.character(trait_list3$Element2)

trait_countsX3 <- trait_list3 %>% 
  count(Pair = str_c(pmin(Element1, Element2), '-',
                     pmax(Element1, Element2)),
        name = 'count', .drop = F)
trait_countsX3

trait_list4$Element1 <- as.character(trait_list4$Element1)
trait_list4$Element2 <- as.character(trait_list4$Element2)

trait_countsX4 <- trait_list4 %>% 
  count(Pair = str_c(pmin(Element1, Element2), '-',
                     pmax(Element1, Element2)),
        name = 'count', .drop = F)
trait_countsX4

trait_countsX$Chain <- 'ALL'
trait_countsX1$Chain <- '1'
trait_countsX2$Chain <- '2'
trait_countsX3$Chain <- '3'
trait_countsX4$Chain <- '4'

trait_countsX$prop <- trait_countsX$count/(9000*4)
trait_countsX1$prop <- trait_countsX1$count/9000
trait_countsX2$prop <- trait_countsX2$count/9000
trait_countsX3$prop <- trait_countsX3$count/9000
trait_countsX4$prop <- trait_countsX4$count/9000


trait_pair_counts <- bind_rows(trait_countsX,
                               trait_countsX1,
                               trait_countsX2,
                               trait_countsX3,
                               trait_countsX4)
trait_pair_counts$Pair <- str_to_title(trait_pair_counts$Pair)
trait_pair_counts$Chain <- as.factor(trait_pair_counts$Chain)
trait_pair_counts$Pair <- as.factor(trait_pair_counts$Pair)
write.csv(trait_pair_counts, 'traits_counted_as_pairs.csv')
h3 <- ggplot(subset(trait_pair_counts, Chain != 'ALL'), aes(x = Chain, y = Pair, fill = prop)) +
  geom_tile() +
  scale_fill_viridis_b(direction = -1) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 4)) +
  ylab("")
h3

ggsave('Figures/pairs_proportion_visited.pdf', h3)

