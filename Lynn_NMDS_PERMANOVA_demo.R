# Hello! This is a demo of how I created NMDS and PERMANOVA plots for my stream ecological research project. 

full_df <- read.csv(file="Lynn_Research_Data.csv", header=TRUE)
bio_df <- read.csv(file="biology_no_EColi.csv", header=TRUE)

bio_df <- decostand(bio_df[,2:5], method="max")

#Create an NMDS plot of biological variables with 21 sample stream sites
ord <- metaMDS(bio_df, distance = "bray", k=2, trymax=100)
stressplot(ord, "metaMDS")
plot(ord)

ord$stress
#add labels to NMDS plot
sites <- read.csv(file="sites_numbering.csv", header=TRUE)

sites$full

cols = c('red')
points(ord, cex=1.2, pch=16, col=cols[1])

orditorp(ord,display="species", col="blue", cex=1.0, air=.01)
orditorp(ord, display="sites", label=sites$full, air=.01, cex=1.0, pcex=.5, pch="+")

#NMDS plot of biological data completed

#NMDS and PERMANOVA of stream chemical data
chem_df <- read.csv(file="chemical_16_sample_sites_grouped.csv", header=TRUE)
env=chem_df[,1:2]

#create NMDS plot
chem_df <- decostand(chem_df[,3:10], method="max")
ord <- metaMDS(chem_df, k=2, trymax=100, distance='bray')

ord$stress
stressplot(ord,"metaMDS") 

plot(ord)
ordiplot(ord,type="text")

#label plot
cols = c ('red', 'green')
points(ord, cex = 1, pch = 16, col = cols[env$Group])

# decoration
ordispider(ord, groups = env$Group, label = TRUE)
ordihull(ord, groups = env$Group, lty = 'dotted')
orditorp(ord,display="species",col="red",air=.5)

#PERMANOVA
pmv <- adonis2(com ~ Group, data=env, permutations=999, method='bray')
pmv

#PERMANOVA calculated with p-value=0.001 

# An NMDS plot highlights the similarities between samples of complex multidimensional data. 
# The closer two points (samples) are on the plot the more similar those samples are in terms of the underlying data. 

# PERMANOVA (non-parametric) tests the null hypothesis that there are no differences in centroids and dispersions of different groups. 





