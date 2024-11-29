#load data

data <- read.csv('../data/tau_cluster_data.csv')
#partition data
E_data <- subset(data, stage_bias == "Egg")
L1_data <- subset(data, stage_bias == "L1")
D_data <- subset(data, stage_bias == "Dauer")
L4_data <- subset(data, stage_bias == "L4")
A_data <- subset(data, stage_bias == "Adult")

#embryo
E.model <- lm(tau ~ max_expression * embryo_clustered, data = E_data)
summary(E.model)
par(mfrow = c(2, 2))
plot(E.model)

#check model
E.model.reduced <- lm(tau ~ max_expression + embryo_clustered, data = E_data)
anova(E.model.reduced, E.model)

#L1
L1.model <- lm(tau ~ max_expression * L1_clustered, data = L1_data)
summary(L1.model)
par(mfrow = c(2, 2))
plot(L1.model)
#reduce
L1.model.reduced <- lm(tau ~ max_expression + L1_clustered, data = L1_data)
anova(L1.model.reduced, L1.model)

#Dauer
D.model <- lm(tau ~ max_expression * Dauer_clustered, data = D_data)
summary(D.model)
par(mfrow = c(2, 2))
plot(D.model)
#reduce
D.model.reduced <- lm(tau ~ max_expression + Dauer_clustered, data = D_data)
anova(D.model.reduced, D.model)

#L4
L4.model <- lm(tau ~ max_expression * L4_clustered, data = L4_data)
summary(L4.model)
par(mfrow = c(2, 2))
plot(L4.model)
#reduce
L4.model.reduced <- lm(tau ~ max_expression + L4_clustered, data = L4_data)
anova(L4.model.reduced, L4.model)

#Adult
A.model <- lm(tau ~ max_expression * Adult_clustered, data = A_data)
summary(A.model)
par(mfrow = c(2, 2))
plot(A.model)
#reduce
A.model.reduced <- lm(tau ~ max_expression + Adult_clustered, data = A_data)
anova(A.model.reduced, A.model)
