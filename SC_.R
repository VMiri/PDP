#************************ Structural covariance networks ************************
rm(list = ls())
#************************ 00 ** Libraries  ************************

library(brainGraph)
library(igraph)
library(data.table)
library(RcppEigen)
library(ggplot2)
library(scales)
library(ggrepel)
library(Matrix)
library(abind)
library(expm)
library(boot)
library(permute)
library(Hmisc)
library(MASS)
library(ade4)
library(mediation)
library(tidyverse)
library(dplyr)
library(foreach)
library(reshape)
library(gridExtra)

#************************ 1 ** IMPORT project data  ************************

datadir <- paste0('~/Documents/R_analisi/PDPmeganalysis/covariancearea/data')
covars <- fread(paste0(datadir, '/covars.csv'), stringsAsFactors=TRUE)
covars[, Study.ID := as.character(Study.ID)]
groups <- covars[, levels(Group)]
head(covars)
str(covars)
# exc subjects based on id s*
exclude.subs <- NULL #c(s2, s25, ...)
#modality = 'area' or 'thickness'
raw_data <- import_scn(datadir=datadir, atlas= 'dk', modality='area',
                       exclude.subs=exclude.subs)
head(raw_data)
# place list elements in the global environment
lapply(seq_along(raw_data), function(x)
  assign(names(raw_data)[x], eval(raw_data[[x]]), envir=.GlobalEnv))

# set density of interest 

densities <- seq(0.05, 0.20, 0.01)#10% densityincreasing .01 every cycle
N <- which(abs(densities - 0.10) < 0.001)

#************************ 2 ** RESIDUALS AND CORRS MATRIX *********************
all.dat.resids <- get.resid(lhrh, covars=covars, exclude.cov='Group')
residPlots <- plot(all.dat.resids)
ml <- gridExtra::marrangeGrob(residPlots, nrow=3, ncol=3)
ggsave('residuals.pdf', ml) #resid plot

#***corr matrix of residuals*** 

corrs <- corr.matrix(all.dat.resids, densities=densities)
#explore 
vhcov <- corrs[["vh"]][["P"]]
novhcov <- corrs[["no"]][["P"]]
vhcovR <- corrs[["vh"]][["R"]]
novhcovR <- corrs[["no"]][["R"]]

#visualise
m3 <- vhcovR - novhcovR
library(corrplot)
library(RColorBrewer)
corrplot(vhcovR, method="square", tl.cex=0.6, addCoefasPercent = FALSE,
         )
corrplot(novhcovR, method="square", tl.cex=0.6, addCoefasPercent = FALSE)
corrplot(m3, method="square", tl.cex=0.6, addCoefasPercent = FALSE,
         col=brewer.pal(n=8, name="PuOr"))
write.table(vhcovR, file='vhcovcorr.csv', sep=',')
write.table(novhcovR, file='novhcovcorr.csv', sep=',')

#correlation thresholds 
corrths <- t(sapply(corrs, with, thresholds))
write.table(corrths, file='corrths.csv', sep=',')

#verify matrix difference
library(psych)
cortest(R1=vhcov,R2=novhcov,n1=118L,n2 = 349, fisher = TRUE,cor=TRUE)

# ************ 3compare weighted corrs coefficients with fisher z**************
library(cocor)

for (i in  (1:length(vhcovR[,1]))) {
    print(cocor.indep.groups(r1.jk=vhcovR[i,26],r2.hm =novhcovR[i,26], n1=118, n2=349, 
                             alternative=c("greater") , test="fisher1925"))
  }

#saving each result as vector, then correct for multiple comparisons
x<-unlist(x)
th <- p.adjust(c(x),method="fdr")


#****************************************** 4 graph metrics ******************************************** 

# Create simple, undirected graphs for each group
g <- lapply(corrs, function(x)
  apply(x$r.thresh, 3, graph_from_adjacency_matrix, mode='undirected', diag=F))

    # **IF :graphs weighed by the correlation coefficient

        g <- lapply(corrs, function(x)
             apply(x$r.thresh, 3, function(y)
             graph_from_adjacency_matrix(x$R * y, mode='undirected', diag=F, weighted=T)))


#set graph attributes - takes a while

g <- Map(function(x, y) lapply(x, set_brainGraph_attr, atlas=atlas,
                              modality=modality, group=y), 
         g, as.list(groups))

    # **IF :hemisphere subgraph only
        g.lh <- lapply(g, lapply, function(x) induced.subgraph(x, V(x)$hemi == 'L'))
        g.rh <- lapply(g, lapply, function(x) induced.subgraph(x, V(x)$hemi == 'R'))

#******* measures of interest - data table of the graph metrics for further analyses and plots ******************

dt.G <- rbindlist(lapply(g, graph_attr_dt))
dt.V <- rbindlist(lapply(g, function(x) rbindlist(lapply(x, vertex_attr_dt))))
setkey(dt.V, density, Group)
# Maximum degree for each Group & density
dt.V[density < 0.1, .SD[which.max(degree), .(region, degree)],
     by=.(Group, density)]
# List the regions with the highest participation coefficient at one density
dt.V[density == densities[N], .SD[order(-PC, -degree)[1:5],
                                  .(density, region, lobe, hemi, degree, PC)], by=Group]
# Count vertices with betweenness > mean + sd (i.e. hub regions)
dt.V[density == round(densities[N], 2) & hubs > 0, .N, by=Group] 
# Mean degree by Group and lobe, for one density
dt.V[density == densities[N], .(mean.deg=mean(degree)), by=.(Group, lobe)]

#************************************ 5 tidy and save datasets  **********************************************
# For a given density, vertex-wise network measures 
charCols <- names(which(sapply(dt.V, is.character)))
dt.V.tidy <- melt(dt.V, id.vars=c("density", charCols)) 

# Number of rows = (number of vertices) * (number of variables)
dt.V.tidy[seq(1, nrow(dt.V.tidy), nrow(get(atlas))),]
write.table(dt.V.tidy, file='graph_metrics_area.csv',sep=',')

# or: 
# Number of rows = (number of densities) * (number of global measures)
(dt.G.tidy <- melt(dt.G, c("density",names(which(sapply(dt.G, is.character))))))
write.table(dt.V.tidy, file='graph_metrics_area_dtG.csv',sep=',')



#****************************************** 6 group analysis: permutation ***************************
library(permute)

#**********************bootstrapping for within group on graph measures - estimate variability ****************
#*for MODULARITY 
kNumBoot <- 1e3
bootmod <- brainGraph_boot(densities, all.dat.resids, R=kNumBoot, measure = "mod",.progress=FALSE)
summary(bootmod)
bootmod.p <- plot(bootmod)
p1 <- bootmod.p$se + theme(legend.position=c(1, 1), legend.justification=c(1, 1))
p2 <- bootmod.p$ci + theme(legend.position=c(1, 1), legend.justification=c(1, 1))
gridExtra::grid.arrange(p1, p2)
write.table(bootmod, file='_boot.txt')


#*# PERMUTATIONS 

resids.vh <- all.dat.resids
myPerms <- shuffleSet(n=nrow(resids.vh$resids.all), nset=5e3) 
perms.vh <- brainGraph_permute(densities, resids.vh, perms=myPerms, level="vertex", measure = "btwn.cent",
                               atlas='dk') 
permsvh.res <- summary(perms.vh)
permsvh.res
write.table(permsvh.res, file='premresults.csv',sep=',')
permplot <- plot(perms.vh, measure='degree')
print(permplot[[1]])

#checks on permutation objects - in case of problems with previous lines
as.numeric(resids$resids.all$Group)
groups <- as.numeric(resids$resids.all$Group) # If #1 above worked fine
dim(resids$resids.all[which(groups[myPerms[1, ]] == 1)])   
dim(resids$resids.all[which(groups[myPerms[1, ]] == 2)])
# For as many Groups as you have data for

#permutation for AUC       
auc.perms <- brainGraph_permute(densities, resids.vh, perms=myPerms, level='vertex', atlas='dk', auc=TRUE) 
summary(auc.perms)
write.table(permsvh.res, file='perm__auc.csv', sep=',')
permplot2 <- plot(auc.perms, measure='Betweenness centrality')
print(permplot2[[1]])


# *************************************** 7  plotting ******************************************
# communities
library(RGtk2)
library(cairoDevice)
#**** adjancency matrices 
#adjacency matrix by communalities - change density depending on analyse, fixed 10
matplot.comm1 <- plot_corr_mat(corrs[[1]]$r.thresh[, , 10], type='comm',
                               g=g[[1]][[10]], group=groups[1])
matplot.comm2 <- plot_corr_mat(corrs[[2]]$r.thresh[, , 10], type='comm',
                               g=g[[2]][[10]], group=groups[2])
print(matplot.comm1)
print(matplot.comm2)


#adjacency matrix by lobe 
matplot.lobe1 <- plot_corr_mat(corrs[[1]]$r.thresh[, , 10], type='lobe',
                               g=g[[1]][[10]], group=groups[1])
matplot.lobe2 <- plot_corr_mat(corrs[[2]]$r.thresh[, , 10], type='lobe',
                               g=g[[2]][[10]], group=groups[2])

print(matplot.lobe1)
print(matplot.lobe2)




#* Network plot, requires igrpah which is arelady loaded  *********************

#network representation of the matrix 

vhcov <- corrs[["vh"]][["R"]]
novhcov <- corrs[["no"]][["R"]]

# Make an Igraph object from this matrix, keep only cells that
# have a high corr coeff:

vhcov[vhcov < 0.6] <- 0

network.v <- graph_from_adjacency_matrix(
  m3, weighted=T, mode="undirected", diag=F)


plot(network.v, 
     vertex.label.color = "black", 
     vertex.color = rgb(0.1,0.2,0.6),
     vertex.shape="circle", 
     vertex.size = 6,
     vertex.label.cex=0.4,
     vertex.label.dist=1,
     vertex.label.font=2.5,  
     edge.arrow.size = 0.2,
     edge.color="grey77",
     edge.lty="solid",
     edge.width=2, 
     repel=FALSE,
     layout = layout.circle(network.v))

