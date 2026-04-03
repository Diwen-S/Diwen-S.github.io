rm(list = ls())
library(lpSolve)
library(gtools)
n = 3
State <- c(0, 1)
mu <- 0.5 #uniform prior over binary states
actions <- c(0, 0.5, 1)
action_profiles = expand.grid(a1 = actions, 
                              a2 = actions,
                              a3 = actions)

action_grid <- expand.grid(a1 = actions,
                           a2 = actions,
                           a3 = actions,
                           state = State)

action_grid <- action_grid[order(action_grid$state, action_grid$a1,
                                 action_grid$a2, action_grid$a3), ]

action_grid$idx <- 1:nrow(action_grid)
# objective of the organisation wrt actions (v)
obj_org <- rowSums(action_grid[,1:3])

gammaii <- c(0.5, 4/3, 0.25)
gammaij <- 0.1
c <- 0.75


#objective for the workers wrt (actions, state)
obj_worker <- t(apply(action_grid, 1, function(row){
  current_action <- c(row["a1"], row["a2"], row["a3"])
  w <- row["state"]
  
  u <- sapply(1:3, function(i){
    return(current_action[i] * (gammaii[i] * w + gammaij * sum(current_action[-i]))
           - c * current_action[i]^2)
  })
}))

utility_i <- function(action, w, i){
  action[i] * (gammaii[i] * w + gammaij * sum(action[-i])) - c * action[i]^2
}

get_coeffs <- function(active_action_grid, i, ai_p){
  coeffs = rep(0, n_vars)
  for (idx in active_action_grid){
    u_current <- obj_worker[idx, i]
    
    action_profile_deviated <- action_grid[idx, 1:3]
    action_profile_deviated[i] <- ai_p
    w <- action_grid$state[idx]
    u_deviated <- utility_i(action = action_profile_deviated, w = w, i = i)
    
    coeffs[idx] = u_current - u_deviated
  }
  return(coeffs)
}
#Constraints
n_vars <- nrow(action_grid)

## Two plausibility constraints
b_plau <- c(0.5, 0.5)
dir_plau <- c("=", "=")
A_plau <- matrix(0, nrow = 0, ncol = n_vars)

for (s in State){
  rows<- as.numeric(action_grid$state == s)
  A_plau <- rbind(A_plau, rows)
}

## BCE: obedience constraint with full uncertainty over a_{-i}
A_BCE <- matrix(0, nrow = 0, ncol = n_vars)
for (i in 1:3){
  for (ai in actions){
    for (ai_p in actions){
      if (ai == ai_p) next
      
      active_action_grid<- which(action_grid[ ,i]  == ai)
      coeffs <- get_coeffs(active_action_grid = active_action_grid, i = i, ai_p = ai_p)
      A_BCE <- rbind(A_BCE, coeffs)
    }
  }
}
b_BCE <- c(b_plau, rep(0, nrow(A_BCE)))
dir_BCE <- c(dir_plau, rep(">=", nrow(A_BCE)))
A_BCE <- rbind(A_plau, A_BCE)
sol_BCE <- lp("max", obj_org, A_BCE, dir_BCE, b_BCE)
action_grid$bce_prob <- sol_BCE$solution
max_obj = c("BCE" = sol_BCE$objval)

## SMS: in the meeting (a_i != \tilde{a_i}) => constraint only has two non-zero 
##                                             coeffs (a, w = 0) & (a, w = 1)
##      outside meeting (a_i == \tilde{a_i}) => BCE style constraint
## optimal SMS is the maximum of |A| LP problems

optimal_effort = -Inf
for (a_tilde in 1:nrow(action_profiles)){
  A_sms = matrix(0, nrow = 0, ncol = n_vars)
  a_tilde <- action_profiles[a_tilde, ]
  for (i in 1:3){
    for (ai in actions){
      for (ai_p in actions){
        if (ai == ai_p) next
        #get the column with non-zero coefficients
        if (ai == a_tilde[i]){
          active_action_grid <- which(action_grid[ ,i] == ai)
          coeffs <- get_coeffs(active_action_grid = active_action_grid, i = i, ai_p = ai_p)
          A_sms = rbind(A_sms, coeffs)
        }else{
          possible_augments <- unique(action_profiles[ , -i])
          for (aug in 1:nrow(possible_augments)){
            if (i == 1){
              active_action_grid <- which(action_grid[ ,1] == ai &
                                          action_grid[ ,2] == possible_augments[aug, 1]&
                                          action_grid[ ,3] == possible_augments[aug, 2])
            }else if(i == 2){
              active_action_grid <- which(action_grid[ ,2] == ai &
                                          action_grid[ ,1] == possible_augments[aug, 1]&
                                          action_grid[ ,3] == possible_augments[aug, 2])
            }else {
              active_action_grid <- which(action_grid[ ,3] == ai &
                                          action_grid[ ,1] == possible_augments[aug, 1]&
                                          action_grid[ ,2] == possible_augments[aug, 2])
            }
            coeffs <- get_coeffs(active_action_grid = active_action_grid, i = i, ai_p = ai_p)
            A_sms = rbind(A_sms, coeffs)
          }
        }
      }
    }
  }
  b_sms = c(b_plau, rep(0, nrow(A_sms)))
  dir_sms = c(dir_plau, rep(">=", nrow(A_sms)))
  A_sms = rbind(A_plau, A_sms)
  sol_temp <- lp("max", obj_org, A_sms, dir_sms, b_sms)
  if (sol_temp$status == 0 && sol_temp$objval > optimal_effort){
    optimal_effort <- sol_temp$objval
    action_grid$sms_prob <- sol_temp$solution
    optimal_a_tilde <- a_tilde
    sol_sms <- sol_temp
  }
}

max_obj <- c(max_obj, "SMS" = optimal_effort)


# Public Information: A special form of SMS in which all constraints are SMS pattern

optimal_effort = -Inf
A_pi = matrix(0, nrow = 0, ncol = n_vars)
for (i in 1:3){
  for (ai in actions){
    for (ai_p in actions){
      if (ai == ai_p) next
      #get the column with non-zero coefficients
      possible_augments <- unique(action_profiles[ , -i])
      for (aug in 1:nrow(possible_augments)){
        if (i == 1){
          active_action_grid <- which(action_grid[ ,1] == ai &
                                      action_grid[ ,2] == possible_augments[aug, 1]&
                                      action_grid[ ,3] == possible_augments[aug, 2])
        }else if(i == 2){
          active_action_grid <- which(action_grid[ ,2] == ai &
                                      action_grid[ ,1] == possible_augments[aug, 1]&
                                      action_grid[ ,3] == possible_augments[aug, 2])
        }else {
          active_action_grid <- which(action_grid[ ,3] == ai &
                                      action_grid[ ,1] == possible_augments[aug, 1]&
                                      action_grid[ ,2] == possible_augments[aug, 2])
        }
          coeffs <- get_coeffs(active_action_grid = active_action_grid, i = i, ai_p = ai_p)
          A_pi = rbind(A_pi, coeffs)
      }
    }
  }
}
b_pi = c(b_plau, rep(0, nrow(A_pi)))
dir_pi = c(dir_plau, rep(">=", nrow(A_pi)))
A_pi = rbind(A_plau, A_pi)
sol_temp <- lp("max", obj_org, A_pi, dir_pi, b_pi)
if (sol_temp$status == 0 && sol_temp$objval > optimal_effort){
  optimal_effort <- sol_temp$objval
  action_grid$pi_prob <- sol_temp$solution
  optimal_a_tilde <- a_tilde
  sol_pi <- sol_temp
}

max_obj <- c(max_obj, "PI" = optimal_effort)
max_obj


## Delegated Hierarchies (DH)
## optimal DH is the maximum of (n!) * |A| LP problems
## For each order, we loop over all a_bar (serves as the best deviation in recommendation to successors)

# get_coeffs funtion for DH to capture the recommendation deviation to a_bar (due to the str monotonicity construction)

get_coeffs_dh <- function(active_action_grid, i, ai_p, less_informed, a_bar){
  coeffs = rep(0, n_vars)
  for (idx in active_action_grid){
    u_current <- obj_worker[idx, i]
    
    action_profile_deviated <- action_grid[idx, 1:3]
    action_profile_deviated[i] <- ai_p
    
    if (length(less_informed) > 0){
      action_profile_deviated[, less_informed] <- a_bar[less_informed]
    }
    
    w <- action_grid$state[idx]
    u_deviated <- utility_i(action = as.numeric(action_profile_deviated), w = w, i = i)
    
    coeffs[idx] = u_current - u_deviated
  }
  return(coeffs)
}

epsilon <- 1e-4
tot_orders <- permutations(n = 3, r = 3, v = 1:3)
optimal_effort_dh <- -Inf
best_ord <- NULL
best_a_bar <- NULL

for (o_idx in 1:nrow(tot_orders)){
  ord <- tot_orders[o_idx, ]
  
  for (a_bar_idx in 1:nrow(action_profiles)){
    a_bar <- unlist(action_profiles[a_bar_idx, ])
    
    A_dh <- matrix(0, nrow = 0, ncol = n_vars)
    
    # p(a_bar) >= epsilon
    idx_a_bar <- which(action_grid$a1 == a_bar[1] & 
                       action_grid$a2 == a_bar[2] & 
                       action_grid$a3 == a_bar[3])
    row_eps <- rep(0, n_vars)
    row_eps[idx_a_bar] <- 1
    A_dh <- rbind(A_dh, row_eps)
    b_dh <- c(epsilon)
    dir_dh <- c(">=")
    
    # p(a_hat) = 0 for a_hat > a_bar
    for (idx in 1:n_vars){
      a_hat <- unlist(action_grid[idx, 1:3])
      if (any(a_hat > a_bar)){
        row_zero <- rep(0, n_vars)
        row_zero[idx] <- 1
        A_dh <- rbind(A_dh, row_zero)
        b_dh <- c(b_dh, 0)
        dir_dh <- c(dir_dh, "=")
      }
    }
    
    # Incentive constraints
    for (i in 1:3){
      rank <- which(ord == i)
      less_informed <- if(rank < 3) ord[(rank + 1):3] else c()
      
      for (ai in actions){
        for (ai_p in actions){
          if (ai == ai_p && length(less_informed) == 0) next
          
          if (length(less_informed) > 0){
            possible_augments <- unique(action_profiles[, less_informed, drop = FALSE])
            
            for (aug in 1:nrow(possible_augments)){
              match_vec <- (action_grid[, i] == ai)
              for (c_idx in 1:length(less_informed)){
                col_agent <- less_informed[c_idx]
                match_vec <- match_vec & (action_grid[, col_agent] == possible_augments[aug, c_idx])
              }
              active_action_grid <- which(match_vec)
              
              coeffs <- get_coeffs_dh(active_action_grid, i, ai_p, less_informed, a_bar)
              A_dh <- rbind(A_dh, coeffs)
            }
          } else {
            active_action_grid <- which(action_grid[, i] == ai)
            coeffs <- get_coeffs_dh(active_action_grid, i, ai_p, less_informed, a_bar)
            A_dh <- rbind(A_dh, coeffs)
          }
        }
      }
    }
    
    b_dh <- c(b_dh, rep(0, nrow(A_dh) - length(b_dh)))
    dir_dh <- c(dir_dh, rep(">=", nrow(A_dh) - length(dir_dh)))
    
    A_dh <- rbind(A_plau, A_dh)
    dir_dh <- c(dir_plau, dir_dh)
    b_dh <- c(b_plau, b_dh)
    
    sol_temp <- lp("max", obj_org, A_dh, dir_dh, b_dh)
    if (sol_temp$status == 0 && sol_temp$objval > optimal_effort_dh){
      optimal_effort_dh <- sol_temp$objval
      action_grid$dh_prob <- sol_temp$solution
      best_ord <- ord
      best_a_bar <- a_bar
      sol_dh <- sol_temp
    }
  }
}

max_obj <- c(max_obj, "DH" = optimal_effort_dh)
cat("Optimal DH Order:", best_ord, "\nOptimal a_bar:\n")
print(best_a_bar)
max_obj
print(action_grid[action_grid$dh_prob > 1e-6, ])

