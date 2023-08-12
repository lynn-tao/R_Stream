# load data and package
library(vegan)

df <- read.csv(file = "chemical_biology_envfit.csv", header = TRUE)

df

#######
# Next, subset your data so that you have a data frame with only environmental variables (env), 
# and a data frame with only species abundance data (com). In the code below I'm naming the range 
# of columns in my data containing the respective information.

com = df[,9:12]
env = df[,13:16]

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

# be careful the numbering is not site numbers
ordiplot(ord,type="text")


# plot NDMS
ordiplot(ord, type = 'n')

# drawing with points
cols = c ('red', 'green', 'blue')
points(ord, cex = 1, pch = 16, col = cols[env$Group])

env$Group

# decoration
ordispider(ord, groups = env$Group, label = TRUE)
ordihull(ord, groups = env$Group, lty = 'dotted')
orditorp(ord,display="species",col="red",air=0.5)
#orditorp(ord,display="sites",labels='S', cex=0.75,air=0.01)
#orditorp(ord,display="sites",col="yellow",air=0.01)


# add labels
sites <- read.csv(file = "sites_numbering.csv", header = TRUE)
sites$short

#orditorp(ord,display="species",col="red",air=0.01)
orditorp(ord,display="sites",label = sites$short, cex=0.50,air=0.01, pch="+", pcol="grey")

# not really sure what pch and pcol are
orditorp(ord,display="sites",label = sites$short, pch="+", pcol="grey")

#----------------- PERMANOVA ------------------
pmv <- adonis(com ~ Group, data =  env, permutations = 999, method = 'bray')

pmv

# communities differ statistically significant between two groups
# the Group explains 50% of variance

#---------------- plot permuted F-values
densityplot(permustats(pmv))

