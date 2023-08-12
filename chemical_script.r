
library(vegan)

##############################################
# plot overall graph and arrows
#############################################3
dat2 <- read.csv(file = "chemical_16_sample_sites.csv", header = TRUE, row.names = 1)

dat2


# standardization
dat <- decostand(dat2, method="max")

head(dat)

# default use bray curtis distance
ord <- metaMDS(dat,k=2,trymax=100)

# distance = "euclidean"
ord <- metaMDS(dat,k=2,trymax=100,distance = "euclidean")


ord$stress 

stressplot(ord,"metaMDS") 
plot(ord)

ordiplot(ord,type="text")

orditorp(ord,display="species",col="red",air=0.01)
orditorp(ord,display="sites",labels='S', cex=0.75,air=0.01)

### Add Arrows
# Because we we will be making a custom plot, it is useful to extract the sample scores and variable scores.
variableScores <- ord$species
sampleScores <- ord$points

# Add vectors for the variables. 
# The scaling factor textNudge ensures that the labels are plotted slightly before the tips of the arrows.
arrows(0, 0, variableScores[, 1], variableScores[, 2], length=0.1, angle=20)

#textNudge <- 1.2
#text(variableScores[, 1]*textNudge, variableScores[, 2]*textNudge, rownames(variableScores), cex=0.7)


##########################################
###############################################
# use below two lines to add lables to sites
ordiplot(ord,type="n") # clear the screen

sites <- read.csv(file = "sites_numbering.csv", header = TRUE)
sites$short

orditorp(ord,display="species",col="red",air=0.01)
orditorp(ord,display="sites",label = sites$short, cex=0.50,air=0.01, pch="+", pcol="grey")

# not really sure what pch and pcol are
orditorp(ord,display="sites",label = sites$short, pch="+", pcol="grey")


ordiplot(ord,type="n") # clear the screen



###########################scree test and scree diagram
####################
# function to check for number of axis in NMDS
# NMDS is repeaed 10 times for each dimension to
# check for the influence of the random initial configurations
# x must be a data-boject which can be coerced into a matrix;
# max is the maximum number of dimensions
# dist_meas is the distance measure
scree_value <- function(x, max, dist_meas)
  {
  xx <- as.matrix(x)
  scree.points = NULL
  scree.temp = 1
  for (i in 1:max) 
    {
    sol.mds = NULL
    sol.mds <- replicate(10, metaMDS(xx, k=i, trymax=30, distance = dist_meas)$stress,
                         simplify = TRUE)
    scree.points <- append(scree.points, sol.mds)
    }
  return(scree.points)
  }

############# Now we use IBET data
library(vegan)
dat <- read.csv(file = "chemical_biology_12_sample_sites_with_ecoli.csv", header = TRUE)

#dat <- read.csv(file = "chemical_16_sample_sites.csv", header = TRUE)


#data(dat)
head(dat)
summary(dat)

# if data in similar order of magnitude  (compare the max values)
# no standardisation required
# otherwise perform standardization

dune <- decostand(dat, method="max")
summary(dune)

##############################################
# calculate relationship between dimensions and stress
scree.res <- scree_value(dune, max=5, dist_meas = "bray")

#################################
# plot dimensions vs. stress
par(cex = 1.5, las =1)
plot_vec <- as.vector(sapply(1:(length(scree.res)/10), function(x) rep(x,10), simplify="vector"))

# prerequisite for plotting replicated runs of nmds

plot(plot_vec, scree.res, ylab = "Stress", xla = "Dimensions")
lines(1:(length(scree.res)/10), as.vector((by(scree.res, plot_vec, mean))))
abline(h=0.10, col="red")

# red line indicates threshold for good Stress1 value

# Stress-dimension plot: Two dimensions almost yield as good description of representation of data as you use 
# three dimentions.  This representation of data is quite reliable.

# We select two dimensions as Stress level is around 10




################## MDS 

###############
# if distance measure is not specified bray curtis is used
specnmds <- metaMDS(dune, k=2)

##-----------------------
#specnmds <- metaMDS(dune, k=2, distance = "euclidean" )
################

################### plot Shepard diagram
stressplot(specnmds)
# Non-metric fit is based on stress and calculated as R^2 = sqrt(1 - S^2)
# Linear R^2 = explained variation of regression of observed distances (raw data) against ordination distances (the mapped data)
# There is a good relationship between the raw data and the mapped data

# Prepare for plotting
sumcols <- colSums(dune)
# calculate sum of columns. If two objects are plotted above each other in the ordination, select the object that has higher column sums
par(cex = 1.5)
plot(specnmds, dis = "sp", type = "n")
nmdsplot <- orditorp(specnmds, display = "sp", priority = sumcols, pcol = "gray20", pch = "+", cex = 0.8)

ordiplot(ord,type="n")
orditorp(ord,display="species",col="red",air=0.01)
orditorp(ord,display="sites",cex=0.75,air=0.01)
# species scores are obtained by weighted averaging as for RDA (using WA scores)


#####################################################
#  fitting of enviornmental variables on the NMDS ordination results
#  not used
data(dune.env)
head(dune.env)

# we transform some factor variables into continues numeric variables
dune.env$Moisture <- as.numeric(dune.env$Moisture)
dune.env$Manure <- as.numeric(dune.env$Manure)






############################################### non ecological variables, learning from article  #######33333333
# First, read the data, cull it, and apply data transformation, ending with one that insures that all data values are positive. 
# In this case, the author remove column 1 because it is stratigraphic position (not a geochemical measurement).
# but not applicable in IBET

# row.names
nashville <- read.table('sites_chemical_16_sample_sites.csv', header=TRUE, row.names=1, sep=',')

nashville

#geochem <- nashville[ , -1]
#geochem

# standardization
#geochemdat <- decostand(dat2, method="max")
#geochemdat

#  Run the NMS as for non-ecological data, as described above, and display the results.
chemNMS <- metaMDS(geochemdat, distance='euclidean', k=3, trymax=50, autotransform=FALSE, wasscores=TRUE, noshare=FALSE)

chemNMS

#  experiment with NMDS plot
ordiplot(chemNMS,type="n")

orditorp(chemNMS,display="species",col="red",air=0.01)
#orditorp(ord,display="sites",label = geochemdat$sites, cex=0.50,air=0.01, pch="+", pcol="grey")
orditorp(chemNMS,display="sites",label = 't', cex=1.25,air=0.01, pch="+", pcol="grey")

############################
# ALTERNATIVE PLOTTING, NOT USED
###################
plot(nms, type='n')

points(nms, display=c('sites'), choices=c(1, 2), pch=3, col='red')

text(nms, display=c('species'), choices=c(1, 2), col='blue', cex=0.7)
######################################################

# Because we we will be making a custom plot, it is useful to extract the sample scores and variable scores.
variableScores <- chemNMS$species
sampleScores <- chemNMS$points
