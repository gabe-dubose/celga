#load data

data <- read.csv('../data/tau_cluster_data.csv')
#partition data
E_data <- subset(data, stage_bias == "Egg")
L1_data <- subset(data, stage_bias == "L1")
D_data <- subset(data, stage_bias == "Dauer")
L4_data <- subset(data, stage_bias == "L4")
A_data <- subset(data, stage_bias == "Adult")

#embryo model
embryo.model <- glm(embryo_clustered ~ tau, data = E_data, family = binomial)
summary(embryo.model)
#pseudo r2
1 - (embryo.model$deviance / embryo.model$null.deviance)

#L1 model
L1.model <- glm(L1_clustered ~ tau, data = L1_data, family = binomial)
summary(L1.model)
1 - (L1.model$deviance / L1.model$null.deviance)

#Dauer model
D.model <- glm(Dauer_clustered ~ tau, data = D_data, family = binomial)
summary(D.model)
1 - (D.model$deviance / D.model$null.deviance)

#L4 model
L4.model <- glm(L4_clustered ~ tau, data = L4_data, family = binomial)
summary(L4.model)

#Adult model
A.model <- glm(Adult_clustered ~ tau, data = A_data, family = binomial)
summary(A.model)


