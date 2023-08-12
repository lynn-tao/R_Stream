library(vegan)
library(ggplot2)

df = read.csv("chemical_biology_envfit.csv", header = TRUE)

df

#######
# Next, subset your data so that you have a data frame with only environmental variables (env), 
# and a data frame with only species abundance data (com). In the code below I'm naming the range 
# of columns in my data containing the respective information.

com = df[,9:12]
env = df[,1:8]

com
env

# Perform the NMDS ordination, as discussed more in this tutorial

dat <- decostand(com, method="max")

head(dat)
nmds <- metaMDS(dat,k=2,trymax=100)
stressplot(nmds,"metaMDS") 
plot(nmds)

################
# Now we run the envfit function with our environmental data frame, env.
# The first parameter is the metaMDS object from the NMDS ordination we just performed. Next is env, our environmental data frame. Then we state we want 999 permutations, and to remove any rows with missing data.

en = envfit(nmds, env, permutations = 999, na.rm = TRUE)


# Let's see what the envfit function returns:
en

#plot(nmds)
plot(en)




