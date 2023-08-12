full_df <- read.csv(file = "Lynn_Research_Data.csv", header = TRUE)

bio_df <- read.csv(file = "biology_no_EColi.csv", header = TRUE)

appbio_df <- bio_df[,2:5]

bio_df <- decostand(appbio_df, method="max")

# An NMDS plot highlights the similarities between samples of complex multidimensional data. The closer two points (samples) are on the plot the more similar those samples are in terms of the underlying data. 
ord <- metaMDS(bio_df, distance="bray", k=2, trymax=100)
stressplot(ord,"metaMDS")
plot(ord)

sites <- read.csv(file = "sites_numbering.csv", header = TRUE)
sites
sites$full
sites$short

cols = c ('red', 'green', 'blue','cyan')
points(ord, cex = 1.2, pch = 16, col = cols[2]) 

orditorp(ord,display="species",col="blue",cex=1.0, air=0.01)
orditorp(ord,display="sites",label = sites$full, air=0.01, cex=1.0, pcex=0.5, pch="+")

###

chem_df <- read.csv(file = "chemical_16_sample_sites_grouped.csv", header = TRUE)

ord <- metaMDS(chem_df,k=2,trymax=100, distance = 'bray')
ord$stress
ordiplot(ord,type="text")

cols = c ('red', 'green')
points(ord, cex = 1, pch = 16, col = cols[env$Group])

ordispider(ord, groups = env$Group, label = TRUE)
ordihull(ord, groups = env$Group, lty = 'dotted')
orditorp(ord,display="species",col="red",air=.5)

# PERMANOVA (non-parametric) tests the null hypothesis that there are no differences in centroids and dispersions of different groups
pmv <- adonis2(com ~ Group, data =  env, permutations = 999, method = 'bray')
pmv



