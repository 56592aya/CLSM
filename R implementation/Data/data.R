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

library(plyr)
sources.count <- table(Network$source_id)
sources_count <- count(Network, "source_id")
sinks.count <- table(Network$sink_id)
sinks_count <- count(Network, "sink_id")

head(sort(sources.count, decreasing = T))
head(sort(sinks.count, decreasing = T))
max(Network$source_id)
max(Network$sink_id)
unsourced <- setdiff(1:79280, Network$source_id)
unsinked <- setdiff(1:79280, Network$sink_id)
tail(intersect(persons.data$id, Network$follower))
persons.data[persons.data$id == "00001020", 'followersNbr']
persons.data[persons.data$id == "00001020", 'followsNbr']
Network[Network$follower == "00001020","sink_id"]
Network[Network$followed == "00001020","source_id"]



persons.data[persons.data$id == "_Haley_", 'followersNbr']
persons.data[persons.data$id == "_Haley_", 'followsNbr']
Network[Network$follower == "_Haley_","sink_id"]
Network[Network$followed == "_Haley_","source_id"]


persons.data[persons.data$id == "00337766", 'followersNbr']
persons.data[persons.data$id == "00337766", 'followsNbr']
Network[Network$follower == "00337766","sink_id"]
Network[Network$followed == "00337766","source_id"]


##################
persons.data$memberSince[1]
persons.data$lastLogin[1]

persons.data$lastLoginYear <- NA
persons.data$lastLoginMonth <- NA
persons.data$lastLoginDay <- NA

persons.data$memberSinceYear <- NA
persons.data$memberSinceMonth <- NA
persons.data$memberSinceDay <- NA
for(i in 1:nrow(persons.data)){
    persons.data$lastLoginYear[i] = strsplit(strsplit(persons.data$lastLogin[i], split =  " ")[[1]], "-" )[[1]][1]
    persons.data$lastLoginMonth[i] = strsplit(strsplit(persons.data$lastLogin[i], split =  " ")[[1]], "-" )[[1]][2]
    persons.data$lastLoginDay[i] = strsplit(strsplit(persons.data$lastLogin[i], split =  " ")[[1]], "-" )[[1]][3]
    persons.data$memberSinceYear[i] = strsplit(strsplit(persons.data$memberSince[i], split =  " ")[[1]], "-" )[[1]][1]
    persons.data$memberSinceMonth[i] = strsplit(strsplit(persons.data$memberSince[i], split =  " ")[[1]], "-" )[[1]][2]
    persons.data$memberSinceDay[i] = strsplit(strsplit(persons.data$memberSince[i], split =  " ")[[1]], "-" )[[1]][3]
    
}
persons.data$sameDayOnly <- 0
for(i in 1:nrow(persons.data)){
    if(persons.data$memberSince[i] == persons.data$lastLogin[i]){
        persons.data$sameDayOnly[i] <- 1
    }
        
}
sum(persons.data$sameDayOnly)/nrow(persons.data)

sameDayGuys <- persons.data[persons.data$sameDayOnly == 1,]
length(intersect(sameDayGuys$id, persons.products.data$personId))
length(intersect(sameDayGuys$id, Whole.data$originator))
length(intersect(persons.data$id, Whole.data$originator))


x
network <- read.csv("NetworkData/Network.csv")
##look if there are a->b and b->a
has.edge <- function(net, a, b){
    has <- F
    for(r in 1:nrow(net)){
        if(net[r,5] == a && net[r,6] == b){
            has=T
        }
    }
    return(has)
}

for(r in 1:nrow(network)){
    print(paste(r, "\n"))
    if(has.edge(network, network$sink_id[r], network$source_id[r])){
        print(paste("also has ", network$sink_id[r], "\t to\t", network$source_id[r], "\n"))
    }
}