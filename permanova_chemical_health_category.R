# load data and package
library(vegan)

df <- read.csv(file = "chemical_16_sample_sites_health_category.csv", header = TRUE)

df

#######
# Next, subset your data so that you have a data frame with only environmental variables (env), 
# and a data frame with only species abundance data (com). In the code below I'm naming the range 
# of columns in my data containing the respective information.

com = df[,4:11]
##com = df[,3:9]

env = df[,1:3]

com
env

#-- checking data
str(com)
str(env)

#-- more checking
range(com)
#range(env)

## -------- Distance Matrix -------------
##############################################
# plot overall graph and arrows
#############################################3

# standardization
dat <- decostand(com, method="max")

head(dat)

# default use bray curtis distance
ord <- metaMDS(dat,k=2,trymax=100, distance = 'bray')

# distance = "euclidean"
#ord <- metaMDS(dat,k=2,trymax=100,distance = "euclidean")


ord$stress 

stressplot(ord,"metaMDS") 
plot(ord)

# plot NDMS
ordiplot(ord, type = 'n')

# be careful the numbering is not site numbers
#ordiplot(ord,type="text")




# drawing with points
cols = c ('red', 'green')
points(ord, cex = 1.5, pch = 16, col = cols[env$Group])

env$Group
env$Health
env$�..sites

# decoration
ordispider(ord, groups = env$Health, label = TRUE, cex=1.1)
#ordihull(ord, groups = env$Group, lty = 'dotted')
orditorp(ord,display="species",col="blue",cex=1.0,air=0.5)     # add the variable names
orditorp(ord,display="sites",label = env$�..sites, pch="", cex=1.0, col="black")  # add sites numbering


#----------------- PERMANOVA ------------------
pmv <- adonis(com ~ Group, data =  env, permutations = 999, method = 'bray')

pmv

# communities differ statistically significant between two groups
# the Group explains 50% of variance

#---------------- plot permuted F-values
densityplot(permustats(pmv))

