library(vegan)

dat <- read.csv(file = "biollogy_no_EColi.csv", header = TRUE)

dat2 <- dat[,2:5]

dat <- decostand(dat2, method="max")

  head(dat)
  ord <- metaMDS(dat, distance="bray", k=2, trymax=100)
  stressplot(ord,"metaMDS")
  plot(ord)

  # use below two lines to add labels to sites

  sites <- read.csv(file = "sites_numbering.csv", header = TRUE)
  sites
  sites$full
  sites$short

  # ordiplot(ord,type="n") # clear the screen
  # plot(ord)
  # drawing with points coloring based on PoIS (Impervious surface)
  cols = c ('red', 'green', 'blue','cyan')
  points(ord, cex = 1.2, pch = 16, col = cols[2])  # cex determines the size of dots

  orditorp(ord,display="species",col="blue",cex=1.0, air=0.01)
  orditorp(ord,display="sites",label = sites$full, air=0.01, cex=1.0, pcex=0.5, pch="+")

  #orditorp(ord,display="sites",label = sites$short, pch="+", pcex=0.5, pcol="grey", air=0.1,color="black")
  #ordiplot(ord,type="n") # clear the screen


  #-- Start here from the article with standization
  #head(dat)

  #ord <- metaMDS(dat,k=2,trymax=100)
  #stressplot(ord,"metaMDS")
  #plot(ord)

  #ordiplot(ord,type="n")

  #orditorp(ord,display="species",col="blue",cex=1.2, air=0.08,pch="+", pcol="grey")
  #orditorp(ord,display="sites",cex=1.0,air=0.01)

#   ====  envfit()

  