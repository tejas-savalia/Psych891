setwd('C:\\Users\\cemnlcmaplab\\Documents\\Psych 891 Modelling cog dynamics\\Codes\\Midterm')
gsqs = read.csv(file = 'Tables/Q3/gsqs.csv')
gsq_wald = gsqs[4]
gsq_exg = gsqs[3]

softmax_probs = matrix(, nrow = 10, ncol = 2)

logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}

#Punish for more free parameters by adding a lower value to low free parameters model's gsquared error
gsq_wald = gsq_wald + 4
gsq_exg = gsq_exg + 6

for (i in 1:length(gsq_gam[[1]])){
  softmax_probs[i,] = softmax(c(-gsq_wald[i, ], -gsq_exg[i, ]))
}
softmax_probs
gsq_wald[1]  
gsq_exg[1]
report_matrix = matrix(c(0.27, 0.84, 0.41, 0.8, 0.78, 0.15, 0.85, 0.15, 0.59, 0.15, 
                         0.73, 0.16, 0.59, 0.2, 0.22, 0.85, 0.15, 0.85, 0.41, 0.85), nrow = 10)

##########################################################################################

gsq_1 = read.csv(file = 'Tables/Q4/single_v/gsq.csv')
gsq_2 = read.csv(file = 'Tables/Q4/two_v/gsq.csv')
gsq_vfix = gsq_1[2]
gsq_vfree = gsq_2[2]

#punish
gsq_vfree = gsq_vfree + 6
gsq_vfix = gsq_vfix + 4

softmax_probs_2 = matrix(, nrow = 20, ncol = 2)

for (i in 1:length(gsq_vfix[[1]])){
  softmax_probs_2[i, ] = softmax(c(-gsq_vfix[i, ], -gsq_vfree[i, ]))
}

report_matrix_2 = matrix(c(0.71, 0.1, 0.15, 0.81, 0.1, 0.72, 0.81, 0.1, 0.12, 0.8, 0.13, 0.51, 0.15, 0.62, 0.8, 0.1, 0.69, 0.23, 0.15, 0.76,
                           0.29, 0.9, 0.85, 0.19, 0.9, 0.28, 0.19, 0.9, 0.88, 0.2, 0.87, 0.49, 0.85, 0.38, 0.2, 0.1, 0.31, 0.77, 0.85, 0.24), nrow = 20)
scale(softmax(c(-gsq_vfix[1, ], -gsq_vfree[1, ])))
scale = function(x){
  x = (0.85-0.15)*(x - 0)/(1) + 0.15
  return (x)
}

report_matrix = scale(softmax_probs)
report_matrix_2 = scale(softmax_probs_2)
#load('Exp2.RData')

write.csv(report_matrix, file = "exp1_prob.csv", col.names = c("Wald, Exg"))
write.csv(report_matrix_2, file = "exp2_prob.csv")

max(mean(report_matrix[,2]), mean(report_matrix[,1]))
mean_11 = mean(report_matrix[,1])
sd_11 = sd(report_matrix[,1])
mean_12 = mean(report_matrix[,2])
sd_12 = sd(report_matrix[,2])


########################################################################################


#AIC_w = c()
AIC_w = gsq_wald + 2*2 + 2*2*(3)/(200 - 2 - 2)
AIC_e = gsq_exg + 2*3 + 2*3*(4)/(200 - 3 - 2)
