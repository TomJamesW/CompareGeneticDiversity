# deterministic.r - compare genetic diversity measurements as differention is increased deterministically. That is each subcommunity starts off with an identical set of equally common alleles and then distinct new alleles are successively added to each subcommunity.

library(rdiversity)
library(vcfR)
library(ggplot2)

#using the vcfR_test object from the vcfR package as a template we will edit the gt section as required. we will just consider the first position.
data("vcfR_test")

vcfR_test@gt <- vcfR_test@gt[1:2,]
vcfR_test@fix <- vcfR_test@fix[1:2,]

vcfR_test@gt <- cbind(vcfR_test@gt,
                      matrix('0', nrow=2, ncol = 5))
colnames(vcfR_test@gt) <- c('FORMAT',paste0('i', 1:8))
vcfR_test@gt[1,-1] <- 1:4
vcfR_test@gt[2,-1] <- 0
vcfR_test@gt[,1] <- 'GT'

d<- gen2dist(vcfR_test)
s<- dist2sim(d, transform = 'l')

partition <- cbind.data.frame(A = c(rep(1, 4), rep(0, 4)), B = c(rep(0, 4), rep(1, 4)))
partitionnorm <- partition/sum(partition)

meta <- metacommunity(partitionnorm, s)

#initialise subcommunities for gendiff calculations
subcoms <- rep(c('a','b'), each = 4)

results <- data.frame(NovelAlleles = NA, Measure = NA, value = NA)



results <- rbind(results, c(0,'Bbar',norm_meta_beta(meta, q=0)$diversity))
results <- rbind(results, c(0,'Gst',0))
results <- rbind(results, c(0,'Gprimest',0))
results <- rbind(results, c(0,'Jost',0))

for (i in 1:20){
  #add individuals with novel alleles
  vcfR_test@gt <- cbind(vcfR_test@gt,
                        matrix(c(2*i+3,0,2*i+4,0), nrow=2, ncol = 2,
                               dimnames = list(NULL, c(paste0('i',(2*i+7):(2*i+8))))))
  #make new partition object
  partition <- rbind(partition,c(1,0),c(0,1))
  partitionnorm <- partition/sum(partition)
  
  #create distance matrix
  d@distance <- rbind(d@distance, rep(1))
  d@distance <- cbind(d@distance, rep(1))
  d@distance <- rbind(d@distance, rep(1))
  d@distance <- cbind(d@distance, rep(1))
  d@distance[2*i+7,2*i+7] <-  d@distance[2*i+8,2*i+8] <- 0
  rownames(d@distance) <- 1:(8+2*i)
  colnames(d@distance) <- 1:(8+2*i)
  #similarity matrix
  s <- dist2sim(d, transform = 'l')
  #create metacommunity object
  m <- metacommunity(partitionnorm, s)
  #calculate Bbar
  results <- rbind(results, c(i,'Bbar',norm_meta_beta(m, q=0)$diversity))
  
  subcoms <- c(subcoms,'a','b')
  
  gd <- genetic_diff(vcfR_test, as.factor(subcoms))
  
  results <- rbind(results, c(i, 'Gst', gd$Gst[1]))
  results <- rbind(results, c(i, 'Gprimest', gd$Gprimest[1]))
  
  Hs <- (gd$Hs_a[1]+gd$Hs_a[1])/2
  jostD <- 2*(gd$Ht[1]-Hs)/(1-Hs)
  results <- rbind(results, c(i, 'Jost', jostD))
}

results <- results[-1,]
results$NovelAlleles <- as.numeric(results$NovelAlleles)
results$value <- as.numeric(results$value)

ggplot(data = results, aes(x=NovelAlleles, y=value)) +
  geom_line(aes(group=Measure, color=Measure))

#########################################################