rm(list = ls())
#Starting transition probability. Three possible locations where the stimulus can occur.
p = 1/3
n_trials = 100
data = read.csv(file = file.choose())

#Converts the observed data to transition probabilities. Siple scaling to match the first
#1/3 probability to match the first data point as observed (simulated in this case)
a_to_transition <- function(data){
  a = rev(data$Dist)
  a_trans = a/sum(a)
  a_trans = a_trans * (1/3) / a_trans[n_trials]
  a_trans[a_trans>1]=1
  return (a_trans)
}

observed_transition_prob = a_to_transition(data)

#The RL model. Increasing the transition probability by a constant factor which is a function
#of the trial number. This alpha is the parameter that is to be learned. It is assumed that
#the transition probability toowards other actions are decreased equally.
#This step occurs on correct predictions. For now, it's all modeled on correct predicitons.
#Function outputs predictions of transition probabilities. 
pred_trans = function(alpha){
  p = c(1/3)
  for (i in (1:n_trials-1)){
    p = c(p, p[i] + ((n_trials/100) - i/100) * alpha)
  }
  return(p)
}

#Chi-squared error function to minimize
chi_sq = function(o, p){
  return (sum((o-p)^2/p))
}
predicted_transition_prob = pred_trans(0.01)


chi_sq(observed_transition_prob, predicted_transition_prob)

#Search over a space for the optimal alpha. 
linear_search = function(range){
  alpha = 0.001
  min_chi_sq_val = chi_sq(observed_transition_prob, pred_trans(alpha))
  for (i in seq(from = min(range), to = max(range), by = 0.0001)){
    new_chi_sq_val = chi_sq(observed_transition_prob, pred_trans(i))
    if (new_chi_sq_val < min_chi_sq_val){
      min_chi_sq_val = new_chi_sq_val
      alpha = i
    }
  }
  return (alpha)
}

alpha = linear_search(c(0.0001, 0.1))
p = pred_trans(alpha)

#convert it back to distance on left to travel on the screen. This predicted distance is what is used
#as the starting point of the DDM model.
predicted_dist = min(data$Dist)*p/(1/3)
write.csv(predicted_dist, file = "Predicted_start_points.csv")

