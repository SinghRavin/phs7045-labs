
##### ----------- Design 1 -----------#####

# Generating data.

fun_1 <- function(k=1000, N = 228, seed = 100, alpha = 0.35, beta=0.65, p=0.35, delta=0.9912){
  set.seed(seed)
  y_arm0 <- rbinom(N/4, 1, p=p)
  y_arm1 <- rbinom(N/4, 1, p=p)
  y_arm2 <- rbinom(N/4, 1, p=p)
  y_arm3 <- rbinom(N/4, 1, p=p) # some changes..

  data_mat <- as.matrix(cbind(y_arm0, y_arm1, y_arm2, y_arm3))
  # Drawing from posterior distribution.
  p_arm0 <- rbeta(k, alpha + sum(y_arm0), beta + N/4 - sum(y_arm0))
  p_arm1 <- rbeta(k, alpha + sum(y_arm1), beta + N/4 - sum(y_arm1))
  p_arm2 <- rbeta(k, alpha + sum(y_arm2), beta + N/4 - sum(y_arm2))
  p_arm3 <- rbeta(k, alpha + sum(y_arm3), beta + N/4 - sum(y_arm3))
  p_vec <- c(sum(as.numeric(p_arm1>p_arm0))/k, sum(as.numeric(p_arm2>p_arm0))/k, 
           sum(as.numeric(p_arm3>p_arm0))/k)
  return(c(which.max(p_vec),p_vec))
}

##### ----------- Design 2 -----------#####

# Patients 1-40

k=1000
set.seed(100)

n01 <- 10
n11 <- 10
n21 <- 10
n31 <- 10

y_arm01 <- rbinom(n01, 1, p=p)
y_arm11 <- rbinom(n11, 1, p=p)
y_arm21 <- rbinom(n21, 1, p=p)
y_arm31 <- rbinom(n31, 1, p=p)

p_arm01 <- rbeta(k, alpha + sum(y_arm01), beta + n01 - sum(y_arm01))
p_arm11 <- rbeta(k, alpha + sum(y_arm11), beta + n11 - sum(y_arm11))
p_arm21 <- rbeta(k, alpha + sum(y_arm21), beta + n21 - sum(y_arm21))
p_arm31 <- rbeta(k, alpha + sum(y_arm31), beta + n31 - sum(y_arm31))

d <- as.data.frame(cbind(p_arm01, p_arm11, p_arm21, p_arm31))

d$max_index <- apply(d, 1, FUN = function(x) which.max(x))

#v01 <- sum(as.numeric(d$max_index==1))/k
v11 <- sum(as.numeric(d$max_index==2))/k
v21 <- sum(as.numeric(d$max_index==3))/k
v31 <- sum(as.numeric(d$max_index==4))/k

v01 <- min(v11*((n11+1)/(n01+1)) + v21*((n21+1)/(n01+1)) + v31*((n31+1)/(n01+1)),
           max(v11, v21, v31))

v01 <- v01/(v01+v11+v21+v31)
v11 <- v11/(v01+v11+v21+v31)
v21 <- v21/(v01+v11+v21+v31)
v31 <- v31/(v01+v11+v21+v31)

# Patients 41-80 

n02 <- round(v01*40,0)
n12 <- round(v11*40,0)
n22 <- round(v21*40,0)
n32 <- round(v31*40,0)

y_arm02 <- rbinom(n02, 1, p=p)
y_arm12 <- rbinom(n12, 1, p=p)
y_arm22 <- rbinom(n22, 1, p=p)
y_arm32 <- rbinom(n32, 1, p=p)

p_arm02 <- rbeta(k, alpha + sum(y_arm01), beta + n02 - sum(y_arm01))
p_arm12 <- rbeta(k, alpha + sum(y_arm11), beta + n12 - sum(y_arm11))
p_arm22 <- rbeta(k, alpha + sum(y_arm21), beta + n22 - sum(y_arm21))
p_arm32 <- rbeta(k, alpha + sum(y_arm31), beta + n32 - sum(y_arm31))

#----#

total_patient <- 228
k=1000
total_interim <- 6
n_allocation <- numeric(10,10,10,10)
patient_num <- 40
keep_y_data <- data.frame()

for (interim in 1:total_interim){
  y_arm1 <- rbinom(n_allocation[1], 1, p=p)
  y_arm2 <- rbinom(n_allocation[2], 1, p=p)
  y_arm3 <- rbinom(n_allocation[3], 1, p=p)
  y_arm4 <- rbinom(n_allocation[4], 1, p=p)
  
  keep_y_data <- rbind(keep_y_data, cbind(y_arm1, y_arm2, y_arm3, y_arm4))
  
  p_arm1 <- rbeta(k, alpha + sum(y_arm1), beta + n_allocation[1] - sum(y_arm1))
  p_arm2 <- rbeta(k, alpha + sum(y_arm2), beta + n_allocation[2] - sum(y_arm2))
  p_arm3 <- rbeta(k, alpha + sum(y_arm3), beta + n_allocation[3] - sum(y_arm3))
  p_arm4 <- rbeta(k, alpha + sum(y_arm4), beta + n_allocation[4] - sum(y_arm4))
  
  d <- as.data.frame(cbind(p_arm1, p_arm2, p_arm3, p_arm4))
  d$max_index <- apply(d, 1, FUN = function(x) which.max(x))
  
  v2 <- sum(as.numeric(d$max_index==2))/k
  v3 <- sum(as.numeric(d$max_index==3))/k
  v4 <- sum(as.numeric(d$max_index==4))/k
  
  v1 <- min(v2*((n_allocation[2]+1)/(n_allocation[1]+1)) + 
              v3*((n_allocation[3]+1)/(n_allocation[1]+1)) + 
              v4*((n_allocation[4]+1)/(n_allocation[1]+1)),
            max(v2, v3, v4))
  
  v1 <- v1/(v1+v2+v3+v4)
  v2 <- v2/(v1+v2+v3+v4)
  v3 <- v3/(v1+v2+v3+v4)
  v4 <- v4/(v1+v2+v3+v4)
  
  n_allocation[1] <- round(v1*patient_num,0)
  n_allocation[2] <- round(v2*patient_num,0)
  n_allocation[3] <- round(v3*patient_num,0)
  n_allocation[4] <- round(v4*patient_num,0)
  
  if (interim == floor(total_patient/patient_num)){
    patient_num <- total_patient/patient_num
  }
  
}


