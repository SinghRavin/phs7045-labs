
# Note: Control arm is indexed 1 and treatment arms are indexed 2, 3 and 4.

##### ----------- Design 1 -----------#####

design1_func <- function(k=1000, N=228, seed=100, alpha=0.35, beta=0.65,
                         p=0.35, delta=0.9912){
  set.seed(seed)
  y_arm0 <- rbinom(N/4, 1, p=p)
  y_arm1 <- rbinom(N/4, 1, p=p)
  y_arm2 <- rbinom(N/4, 1, p=p)
  y_arm3 <- rbinom(N/4, 1, p=p)

  data_mat <- as.matrix(cbind(y_arm0, y_arm1, y_arm2, y_arm3))
  # Drawing from posterior distribution.
  p_arm0 <- rbeta(k, alpha + sum(y_arm0), beta + N/4 - sum(y_arm0))
  p_arm1 <- rbeta(k, alpha + sum(y_arm1), beta + N/4 - sum(y_arm1))
  p_arm2 <- rbeta(k, alpha + sum(y_arm2), beta + N/4 - sum(y_arm2))
  p_arm3 <- rbeta(k, alpha + sum(y_arm3), beta + N/4 - sum(y_arm3))
  p_vec <- c(sum(as.numeric(p_arm1>p_arm0))/k, sum(as.numeric(p_arm2>p_arm0))/k, 
           sum(as.numeric(p_arm3>p_arm0))/k)
  char <- paste("Arm",which.max(p_vec)+1,"is best") ## +1 because indexing of arms starts with 1.
  return(list(char,p_vec))
}

##### ----------- Design 2 -----------#####

# Non-modular version of design 2.

set.seed(100)
total_patients <- 228
k=1000
n_allocation <- c(0,0,0,0)
arm_weights <- c(0.25, 0.25, 0.25, 0.25)
patients_num <- 40
keep_y_arm1 <- numeric()
keep_y_arm2 <- numeric()
keep_y_arm3 <- numeric()
keep_y_arm4 <- numeric()
p <- 0.35
alpha <- 0.35
beta <- 0.65
patients_assigned <- 0

while (patients_assigned < total_patients){
  
  interim <- min(patients_num, total_patients - patients_assigned)
  
  n_allocation[1] <- interim - (round(arm_weights[2]*interim,0)+
                                  round(arm_weights[3]*interim,0)+
                                  round(arm_weights[4]*interim,0))
  n_allocation[2] <- round(arm_weights[2]*interim,0)
  n_allocation[3] <- round(arm_weights[3]*interim,0)
  n_allocation[4] <- round(arm_weights[4]*interim,0)
  
  y_arm1 <- rbinom(n_allocation[1], 1, p=p)
  y_arm2 <- rbinom(n_allocation[2], 1, p=p)
  y_arm3 <- rbinom(n_allocation[3], 1, p=p)
  y_arm4 <- rbinom(n_allocation[4], 1, p=p)
  
  keep_y_arm1 <- c(keep_y_arm1,y_arm1)
  keep_y_arm2 <- c(keep_y_arm2,y_arm2)
  keep_y_arm3 <- c(keep_y_arm3,y_arm3)
  keep_y_arm4 <- c(keep_y_arm4,y_arm4)
  
  p_arm1 <- rbeta(k, alpha + sum(y_arm1), beta + n_allocation[1] - sum(y_arm1))
  p_arm2 <- rbeta(k, alpha + sum(y_arm2), beta + n_allocation[2] - sum(y_arm2))
  p_arm3 <- rbeta(k, alpha + sum(y_arm3), beta + n_allocation[3] - sum(y_arm3))
  p_arm4 <- rbeta(k, alpha + sum(y_arm4), beta + n_allocation[4] - sum(y_arm4))
  
  d <- as.data.frame(cbind(p_arm1, p_arm2, p_arm3, p_arm4))
  d$max_index <- apply(d, 1, FUN = function(x) which.max(x))
  
  a2 <- sum(as.numeric(d$max_index==2))/k
  a3 <- sum(as.numeric(d$max_index==3))/k
  a4 <- sum(as.numeric(d$max_index==4))/k
  
  a1 <- min(a2*((n_allocation[2]+1)/(n_allocation[1]+1)) + 
                          a3*((n_allocation[3]+1)/(n_allocation[1]+1)) + 
                          a4*((n_allocation[4]+1)/(n_allocation[1]+1)),
            max(a2, a3, a4))
  
  arm_weights[1] <- a1/(a1+a2+a3+a4)
  arm_weights[2] <- a2/(a1+a2+a3+a4)
  arm_weights[3] <- a3/(a1+a2+a3+a4)
  arm_weights[4] <- a4/(a1+a2+a3+a4)
  
  patients_assigned <- patients_assigned + interim
  
}

# Modular version of design 2.

n_allocation_func <- function(n, arm_weights_vec){
  n_allocation_vector <- numeric()
  n_allocation_vector[1] <- n - (round(arm_weights_vec[2]*n,0)+
                                  round(arm_weights_vec[3]*n,0)+
                                  round(arm_weights_vec[4]*n,0))
  n_allocation_vector[2] <- round(arm_weights_vec[2]*n,0)
  n_allocation_vector[3] <- round(arm_weights_vec[3]*n,0)
  n_allocation_vector[4] <- round(arm_weights_vec[4]*n,0)
  return(n_allocation_vector)
}

binom_outcome_sample_func <- function(n_allocation_input_vec, p){
  y_arm1_vec <- rbinom(n_allocation_input_vec[1], 1, p=p)
  y_arm2_vec <- rbinom(n_allocation_input_vec[2], 1, p=p)
  y_arm3_vec <- rbinom(n_allocation_input_vec[3], 1, p=p)
  y_arm4_vec <- rbinom(n_allocation_input_vec[4], 1, p=p)
  return(list(y_arm1_vec, y_arm2_vec, y_arm3_vec, y_arm4_vec))
}

storing_n_allocation_func <- function(existing_vector,new_vector){
  existing_vector <- existing_vector + new_vector
  return(existing_vector)
}

storing_outcome_sample_func <- function(existing_list,new_list){
  existing_list[[1]] <- c(existing_list[[1]],new_list[[1]])
  existing_list[[2]] <- c(existing_list[[2]],new_list[[2]])
  existing_list[[3]] <- c(existing_list[[3]],new_list[[3]])
  existing_list[[4]] <- c(existing_list[[4]],new_list[[4]])
  return(existing_list)
}

posterior_func <- function(k, alpha, beta, y_arm_input_list, n_allocation_vec){
  p_arm1 <- rbeta(k, alpha + sum(y_arm_input_list[[1]]), beta +
                    n_allocation_vec[1] - sum(y_arm_input_list[[1]]))
  p_arm2 <- rbeta(k, alpha + sum(y_arm_input_list[[2]]), beta +
                    n_allocation_vec[2] - sum(y_arm_input_list[[2]]))
  p_arm3 <- rbeta(k, alpha + sum(y_arm_input_list[[3]]), beta +
                    n_allocation_vec[3] - sum(y_arm_input_list[[3]]))
  p_arm4 <- rbeta(k, alpha + sum(y_arm_input_list[[4]]), beta +
                    n_allocation_vec[4] - sum(y_arm_input_list[[4]]))
  return(list(p_arm1, p_arm2, p_arm3, p_arm4))
}

rar_weights_update_func <- function(p_arm_input_list, n_allocation_inputvec){
  
  arm_weights_vector <- numeric()
  
  d <- as.data.frame(cbind(p_arm_input_list[[1]], p_arm_input_list[[2]],
                           p_arm_input_list[[3]], p_arm_input_list[[4]]))
  d$max_index <- apply(d, 1, FUN = function(x) which.max(x))
  
  len <- length(p_arm_input_list[[1]])
  
  a2 <- sum(as.numeric(d$max_index==2))/len
  a3 <- sum(as.numeric(d$max_index==3))/len
  a4 <- sum(as.numeric(d$max_index==4))/len
  
  a1 <- min(a2*((n_allocation_inputvec[2]+1)/(n_allocation_inputvec[1]+1)) + 
              a3*((n_allocation_inputvec[3]+1)/(n_allocation_inputvec[1]+1)) + 
              a4*((n_allocation_inputvec[4]+1)/(n_allocation_inputvec[1]+1)),
            max(a2, a3, a4))
  
  arm_weights_vector[1] <- a1/(a1+a2+a3+a4)
  arm_weights_vector[2] <- a2/(a1+a2+a3+a4)
  arm_weights_vector[3] <- a3/(a1+a2+a3+a4)
  arm_weights_vector[4] <- a4/(a1+a2+a3+a4)
  
  return(arm_weights_vector)
}

best_arm <- function(p_arm_inputlist){
  len <- length(p_arm_inputlist[[1]])
  p_vec <- c(sum(as.numeric(p_arm_inputlist[[2]]>p_arm_inputlist[[1]]))/len,
             sum(as.numeric(p_arm_inputlist[[3]]>p_arm_inputlist[[1]]))/len, 
             sum(as.numeric(p_arm_inputlist[[4]]>p_arm_inputlist[[1]]))/len)
  char <- paste("Arm",which.max(p_vec)+1,"is best") ## +1 because indexing of arms starts with 1.
  return(list(char,p_vec))
}

design2_func <- function(seed=100, total_patients=228, k=1000, patients_num=40, p=0.35,
                         alpha=0.35, beta=0.65, delta=0.9892){
  set.seed(seed)
  keep_y_arm_list <- list(numeric(0), numeric(0), numeric(0), numeric(0))
  keep_n_allocation_vector <- c(0,0,0,0)
  n_allocation <- c(0,0,0,0)
  arm_weights <- c(0.25, 0.25, 0.25, 0.25)
  patients_assigned <- 0
  
  while (patients_assigned < total_patients){
    interim <- min(patients_num, total_patients - patients_assigned)
    
    n_allocation <- n_allocation_func(interim, arm_weights)
    keep_n_allocation_vector <- storing_n_allocation_func(keep_n_allocation_vector,
                                                        n_allocation)
    y_arm_list <- binom_outcome_sample_func(n_allocation, p)
    keep_y_arm_list <- storing_outcome_sample_func(keep_y_arm_list, y_arm_list)
    posterior <- posterior_func(k, alpha, beta, keep_y_arm_list,
                              keep_n_allocation_vector)
    arm_weights <- rar_weights_update_func(posterior, n_allocation)
    patients_assigned <- patients_assigned + interim
  }
  return_list <- best_arm(posterior)
  return_list[[3]] <- keep_n_allocation_vector
  return(return_list)
  }
