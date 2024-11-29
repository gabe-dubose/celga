
#load data
data <- read.csv('../data/stage_biased_regions_duplication.csv')


#make linear model
model <- lm(non_homolog_ratio~n_genes, data=data)
summary(model)
plot(model)

hist(model$residuals)
shapiro.test(model$residuals)

#data violated assumptions of linear model, so I'm 
#just going check with correlations
cor.test(data$n_genes, data$non_homolog_ratio, method="spearman")

