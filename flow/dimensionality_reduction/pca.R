### PRINCIPLE COMPONENT ANALYSIS ###
require(stats)
require(ggplot2)
require(lattice)
require(dplyr)

#Read in data
data <- data.frame(read.csv("SEPSIS.csv", header=TRUE, sep=","))
sample <- sample_n(data, 100000)

channels <- c("CD57", "CD161", "CD3", "CCR7", "VA7.2", "CD8",
              "Vdelta2", "CD45RA", "PanGD", "CD4", "CD27")

#Perform PCA with prcomp
pca <- prcomp(sample %>% select(channels), scale=TRUE)

#Calculate the amount of variation contributed from each component
pca_var <- pca$sdev^2
pca_var_per <- round(pca_var/sum(pca_var)*100, 2)
barplot(pca_var_per, main="Scree Plot", 
        xlab="Principle Component", 
        ylab="Percentage Variation")

plot(pca$x[,1], pca$x[,2])
#3D Plot with overlay
sample$PC1 <- pca$x[,1]
sample$PC2 <- pca$x[,2]

#plot_ly(data = pca_data, x=~X, y=~Y, z=~Z, type="scatter3d", color=~label) %>% add_markers()
ggplot(pca_data, aes(x=X, y=Y, color=label)) + geom_point() + xlab("PC1") + ylab("PC2")

ggplot(sample, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=cell_type,fill=CCR7),stroke=0.5, shape=21) + 
  xlab("PC1") + ylab("PC2")

