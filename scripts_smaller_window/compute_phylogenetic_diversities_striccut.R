#!/usr/bin/env Rscript

#get files
phylogenies.path <- '../data/strictcut_phylogenies'

phylogeny.files <- list.files(path=phylogenies.path)

#initialize vectors to store output
homologous.group.names <- c()
group.size <- c()
phylogenetic.diversity <- c()

#initialize vectors to store individual gene output
gene.name <- c()
distance.from.root <- c()

for (phylogeny.file in phylogeny.files) {
	#define full file path
	full.phylogeny.path <- file.path(phylogenies.path, phylogeny.file)
	#get homologous group name from file name
	homologous.group.name <- paste(strsplit(phylogeny.file, "_")[[1]][1:2], collapse = "_")
	#load tree in as ape object
	tree <- ape::read.tree(full.phylogeny.path)
	#calculate phylogenetic diversity as sun of branch lengths
	branch.lengths <- tree$edge.length
	sum.of.lengths <- sum(branch.lengths)
	#get number of tips as well
	n.tips <- ape::Ntip(tree)
	
	#add values to their list
	homologous.group.names <- c(homologous.group.names, homologous.group.name)
	group.size <- c(group.size, n.tips)
	phylogenetic.diversity <- c(phylogenetic.diversity, sum.of.lengths)
	
	#calculate distance from each tip to root node
	#root tree using midpoint
	#tree <- ape::root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
	tree <- phytools::midpoint.root(tree)
	node.distances <- ape::dist.nodes(tree)
	root.to.tips.distances <- node.distances[1, 1:n.tips]
	names(root.to.tips.distances) <- tree$tip.label
	#add to data
	distance.from.root <- c(distance.from.root, root.to.tips.distances)
	gene.name <- c(gene.name, names(root.to.tips.distances))
}

#save to dataframe and write to csv
data <- data.frame(homologous.group.names, group.size, phylogenetic.diversity)
write.csv(data, file="../data/homologous_groups_phylogenetic_diversities_striccut.csv", row.names=FALSE)

#save other dataframe and write to csv
tree.distances.data <- data.frame(gene.name, distance.from.root)
tree.distances.data$gene.name <- gsub("^CELE_", "", tree.distances.data$gene.name)

#load tau data
tau.data <- read.csv('../data/tau_cluster_data.csv')

#merge data
dist.tau.data <- merge(tree.distances.data, tau.data, by.x = "gene.name", by.y = "gene")

plot(dist.tau.data$distance.from.root, dist.tau.data$tau)
summary(lm(tau~distance.from.root, data=dist.tau.data))

