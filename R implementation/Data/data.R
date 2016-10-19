rm(list=ls())


################################################################################
byFollowed.data <- read.table("NetworkData/byFollowed.txt",header = T, sep='\t',stringsAsFactors = F)
head(byFollowed.data)

byFollower.data <- read.table("NetworkData/byFollowers.txt",header = T, sep=',',stringsAsFactors = F)
head(byFollower.data)


persons.data <- read.csv("NetworkData/persons.csv",header = T,stringsAsFactors = F)
typeof(persons.data)
head(persons.data)

persons.products.data <- read.csv("NetworkData/personProducts.csv",header = T,stringsAsFactors = F)
persons.products.data$id <- as.character(persons.products.data$id)
persons.products.data$personId <- as.character(persons.products.data$personId)
persons.products.data$productId <- as.character(persons.products.data$productId)
persons.products.data$hivedSince <- as.character(persons.products.data$hivedSince)
typeof(persons.products.data$personId)
length(unique((persons.products.data$personId)))
length(intersect(unioned.individuals,unique((persons.products.data$personId))))
length(intersect(persons.data$id,unique((persons.products.data$personId))))
xtable(head(persons.products.data))
xtable(tail(persons.products.data))
################################################################################
##all unique inidividuals in the follower and followed data
unioned.individuals <- union(byFollowed.data$followed, byFollowed.data$follower)
length(unioned.individuals)
################################################################################
install.packages("xtable")
library(xtable)
xtable(head(persons.data))
xtable(tail(persons.data))

xtable(head(byFollower.data))
xtable(tail(byFollower.data))

# install.packages("ff")
# library(ff)
# Whole.data <- read.csv.ffdf(file="NetworkData/list_whole.csv", header=T)
# Whole.data <- as.data.frame(Whole.data)
# typeof(Whole.data$originator)
# Whole.data$originator <- as.character(Whole.data$originator)
# typeof(Whole.data$originator)
# length(unique(Whole.data$originator))
# length(intersect(persons.data$id,unique(Whole.data$originator)))
# length(intersect(unioned.individuals,unique(Whole.data$originator)))
# length(intersect(persons.products.data$personId,unique(Whole.data$originator)))
# 
# head(Whole.data$originator)
# 
# length(unique(tolower(unioned.individuals)))
# length(unique(tolower(persons.data$id)))
# length(unique(tolower(Whole.data$originator)))
# 
# length(unique(tolower(persons.products.data$personId)))
# tolower
################################################################################
length(intersect(unioned.individuals, persons.data$id))
################################################################################

################################################################################
nrow(persons.data)
length(intersect(unioned.individuals, persons.data$id))
################################################################################
# one way: make network of unioned.individuals
# might be useless when dealing with the persons
####
N= as.numeric(length(unioned.individuals))

byFollower.data$source_id <- NA
byFollower.data$sink_id <- NA
source_count = 1
uniq_followers <- unique(byFollower.data$follower)
uniq_followed <- unique(byFollower.data$followed)

uniq_unioned <- unique(unioned.individuals)
indiv.indexed <- data.frame(matrix(0, nrow=length(uniq_unioned), ncol=2))
names(indiv.indexed) <- c('uniq_unioned', 'index')
indiv.indexed$uniq_unioned <- uniq_unioned

indiv.indexed$index <- 1:N
uniq_unioned <- 
names(byFollower.data)
count = 1
for(f in indiv.indexed$uniq_unioned){
    byFollower.data[byFollower.data$follower == f,4] = indiv.indexed[indiv.indexed$uniq_unioned == f, 2]
    byFollower.data[byFollower.data$followed == f,5] = indiv.indexed[indiv.indexed$uniq_unioned == f, 2]
    count = count+1
    print(count)
}
write.csv(x = byFollower.data, file="NetworkData/Network.csv")
length(intersect(byFollower.data$follower, byFollower.data$followed))

byFollower.data <- byFollower.data[order(c(byFollower.data$source_id, byFollower.data$sink_id)),]


byFollower.data <- read.csv("NetworkData/Network.csv", header = T, stringsAsFactors = F)
Network <- byFollower.data
Network <- Network[with(Network, order(source_id, sink_id)), ]
