p = c(1/3, 1/3, 1/3)
n_trials = 100
alpha = 0.5
data = read.csv(file = file.choose())
a_to_transition <- function(data){
  a = rev(data$Dist)
  a_trans = a/sum(a)
  a_trans = a_trans * (1/3) / a_trans[n_trials]
  a_trans[a_trans>1]=1
  return (a_trans)
}

observed_transition_prob = a_to_transition(data)
pred_trans = function(alpha){
  p = c(1/3)
  for (i in (1:n_trials)){
    p = c(p, p[i] + ((n_trials/100) - i/100) * alpha)
  }
  return(p)
}
predicted_transition_prob = pred_trans(0.01)

