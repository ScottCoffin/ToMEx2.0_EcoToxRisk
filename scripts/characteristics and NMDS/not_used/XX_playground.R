#################################
#### NOT USED IN MANUSCRIPT ####
#################################


install.packages("nomclust")
library(nomclust)
# sample data
data(data20)
# creating an object with results of hierarchical clustering of
hca.object <- nomclust(data20, measure = "lin", method = "average",
                       clu.high = 5, prox = TRUE)
# assigning variable weights
hca.weights <- nomclust(data20, measure = "lin", method = "average",
                        clu.high = 5, prox = TRUE, var.weights = c(0.7, 1, 0.9, 0.5, 0))
# quick clustering summary
summary(hca.object)
# quick cluster quality evaluation
print(hca.object)
# visualization of the evaluation criteria
eval.plot(hca.object)
# a quick dendrogram
plot(hca.object)
# a dendrogram with three designated clusters
dend.plot(hca.object, clusters = 3)


## NMDS on medians for different measurements including length, width, ratio width/length, density, ... more?
## Samples: Data from Kooi et al. 2021 and other sources?