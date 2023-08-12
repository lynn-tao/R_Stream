# load different dataset 
dat2 <- read.csv(file = "chemical_biology_12_sample_sites_with_ecoli.csv", header = TRUE)

# standardization
#dat <- decostand(dat2, method="max")

# calculate Bray Curtis distance
d <- vegdist(dat2)

# use below two lines to add lables to sites
sites <- read.csv(file = "sites_numbering.csv", header = TRUE)
sites$short

# The single linkage clustering can be found with:
csin <- hclust(d,  method = "single")
csin

# The dendrogram can be plotted with:
plot(csin)

