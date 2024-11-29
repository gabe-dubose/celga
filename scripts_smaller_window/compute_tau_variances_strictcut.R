#!/usr/bin/env Rscript

#load data for homologous group phylogenetic diversities
homologous.group.phylo.data <- read.csv('../data/duplication/homologous_groups_phylogenetic_diversities_striccut.csv')
groups <- c(homologous.group.phylo.data$homologous.group.names)

#load homologous group dictionary
homologous.groups.file <- '../data/duplication/homologous_groups_strict_cutoffs.json'
homologous.groups <- rjson::fromJSON(file=homologous.groups.file)

#load tau data
tau.data <- read.csv('../data/tau_cluster_data.csv')

#initialize vectors to store results
group.tau.mean <- c()
group.tau.variance <- c()
group.stage.bias <- c()
group.stage.bias.frequency <- c()
#iterate through groups
for (group in homologous.group.phylo.data$homologous.group.names){
  genes <- gsub("CELE_", "", homologous.groups[[group]])
  group.tau <- tau.data[tau.data$gene %in% genes,]
  #calculate means and variances
  tau.variance <- var(group.tau$tau, na.rm = TRUE)
  tau.mean <- mean(group.tau$tau, na.rm = TRUE)
  
  #get stage bias information
  stage.frequency.table <- table(group.tau$stage_bias)
  stage.bias <- names(which.max(stage.frequency.table))
  stage.bias.frequency <- max(stage.frequency.table) / length(group.tau$stage_bias)
  
  #add stuff to vectors
  group.tau.mean <- c(group.tau.mean, tau.mean)
  group.tau.variance <- c(group.tau.variance, tau.variance)
  group.stage.bias <- c(group.stage.bias, stage.bias)
  group.stage.bias.frequency <- c(group.stage.bias.frequency, stage.bias.frequency)
}

#add data to dataframe
homologous.group.phylo.data$group.tau.mean <- group.tau.mean
homologous.group.phylo.data$group.tau.variance <- group.tau.variance
homologous.group.phylo.data$group.stage.bias <- group.stage.bias
homologous.group.phylo.data$group.stage.bias.frequence <- group.stage.bias.frequency

#save to dataframe
write.csv(homologous.group.phylo.data, file='../data/homologous_group_phylo_data_tau.csv', row.names=FALSE)
