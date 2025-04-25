#!/usr/bin/env Rscript

#load homologous group dictionary
homologous.groups.file <- '../data/dmel_homologous_groups_orthologs_only.json'
homologous.groups <- rjson::fromJSON(file=homologous.groups.file)

#load tau data
tau.data <- read.csv('../data/dmel_stage_specificity_data.csv')

#initialize vectors to store results
group.tau.mean <- c()
group.tau.variance <- c()
group.stage.bias <- c()
group.stage.bias.frequency <- c()
group.size <- c()

#iterate through groups
for (group in names(homologous.groups)){
  genes <- homologous.groups[[group]]
  group.tau <- tau.data[tau.data$gene %in% genes,]
  n.genes <- length(genes)
  if (n.genes >= 4) {
    #calculate means and variances
    tau.variance <- var(group.tau$tau, na.rm = TRUE)
    tau.mean <- mean(group.tau$tau, na.rm = TRUE)
    
    #get stage bias information
    stage.frequency.table <- table(group.tau$stage_bias)
    stage.bias <- names(which.max(stage.frequency.table))
    if (!is.null(stage.bias)){
      stage.bias.frequency <- max(stage.frequency.table) / length(group.tau$stage_bias)
      
      #add stuff to vectors
      group.tau.mean <- c(group.tau.mean, tau.mean)
      group.tau.variance <- c(group.tau.variance, tau.variance)
      group.stage.bias <- c(group.stage.bias, stage.bias)
      group.stage.bias.frequency <- c(group.stage.bias.frequency, stage.bias.frequency)
      group.size <- c(group.size, n.genes)
    }
  }
}

#make dataframe
homologus.group.data <- data.frame(
  group.tau.mean <- group.tau.mean,
  group.tau.variance <- group.tau.variance,
  group.stage.bias <- group.stage.bias,
  group.stage.bias.frequency <- group.stage.bias.frequency,
  group.size <- group.size
)

#write to file
write.csv(homologus.group.data, file='../data/dmel_homologous_group_tau_var_orthogroup_only.csv', row.names=FALSE)
