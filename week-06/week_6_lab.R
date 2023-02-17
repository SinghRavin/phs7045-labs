
## Jonathan's solution

N <- 228
arms <- 4
h0 <- 0.35
h1 <- 0.35
h2 <- 0.35
h3 <- 0.35

y0 <- rbinom(N,size=1,prob=h0)
y1 <- rbinom(N,size=1,prob=h1)
y2 <- rbinom(N,size=1,prob=h2)
y3 <- rbinom(N,size=1,prob=h3)

y <- cbind("0"=y0,"1"=y1,"2"=y2,"3"=y3)

# Generate draws from posterior and calculate pmax
postDraws <- function(y,nMcmc,h0,h1,h2,h3,n0,n1,n2,n3){
  
  #------------------------------
  # Generate draws from posterior
  #------------------------------
  postDraws0 <- rbeta(n=nMcmc, shape1 = h0 + sum(y[1:n0,"0"]==1), shape2 = (1-h0) + n0 - sum(y[1:n0,"0"]==0))
  postDraws1 <- rbeta(n=nMcmc, shape1 = h1 + sum(y[1:n1,"1"]==1), shape2 = (1-h1) + n1 - sum(y[1:n1,"1"]==0))
  postDraws2 <- rbeta(n=nMcmc, shape1 = h2 + sum(y[1:n2,"2"]==1), shape2 = (1-h2) + n2 - sum(y[1:n2,"2"]==0))
  postDraws3 <- rbeta(n=nMcmc, shape1 = h3 + sum(y[1:n3,"3"]==1), shape2 = (1-h3) + n3 - sum(y[1:n3,"3"]==0))
  
  #-----------------------------------
  # Calculate allocation probabilities
  #-----------------------------------
  # v1-v3 probability each arm is best
  # v0 see RMatch Viele et al. 2020
  pBest <- pmax(postDraws0,postDraws1,postDraws2,postDraws3)
  
  v1 <- mean(pBest==postDraws1)
  v2 <- mean(pBest==postDraws2)
  v3 <- mean(pBest==postDraws3)
  
  v0 <- min(sum( c(v1,v2,v3) * (c( n1, n2, n3) + 1) / (n0 + 1), max(v1, v2, v3)) )
  
  # Standardize
  V0 <- v0 / (sum(v0,v1,v2,v3))
  V1 <- v1 / (sum(v0,v1,v2,v3))
  V2 <- v2 / (sum(v0,v1,v2,v3))
  V3 <- v3 / (sum(v0,v1,v2,v3))
  
  # Calculate probability each arm is greater than control
  p1 <- mean(postDraws1 > postDraws0)
  p2 <- mean(postDraws2 > postDraws0)
  p3 <- mean(postDraws3 > postDraws0)
  
  # Report maximum probablity an arm is greater than control
  pMax <- max(p1,p2,p3)
  # unname n0 objects for consistent object names in output
  n0 <- unname(n0)
  n1 <- unname(n1)
  n2 <- unname(n2)
  n3 <- unname(n3)
  out <- c(V0=V0,V1=V1,V2=V2,V3=V3,p1=p1,p2=p2,p3=p3,pMax=pMax,n0=n0,n1=n1,n2=n2,n3=n3)
  return(out)
  
}

design1 <- function(y, nMcmc=10000, h0, h1, h2, h3){
  
  # By the end of the study, each arm will have equal allocation
  postDraws(y=y,nMcmc=nMcmc,
            n0=nrow(y)/ncol(y),
            n1=nrow(y)/ncol(y),
            n2=nrow(y)/ncol(y),
            n3=nrow(y)/ncol(y),
            h0=h0, h1=h1, h2=h2, h3=h3)
  
}


## My solution of Design 1.

# Generating a function to return a rbinom sample with a size n and probability p.
rbinom_sampler <- function(n, p){
  return(rbinom(n, 1, p=p))
}

# Generating a function to return a rbeta sample with size n and outcome vector y_arm.
rbeta_sampler <- function(n, y_arm){
  return(rbeta(k, alpha + sum(y_arm), beta + n - sum(y_arm)))
}

design1_func <- function(k=1000, N=228, seed=100, alpha=0.35,
                         beta=0.65,trt_effect=c(0.35, 0.35, 0.35,
                                                0.35)){
  
  set.seed(seed)
  
  # Creating outcome data.
  y_arm_list <- mapply(rbinom_sampler, rep(N/4,4), trt_effect, SIMPLIFY = FALSE)
  
  # Drawing from posterior distribution.
  p_arm_list <- mapply(rbeta_sampler, rep(N/4,4), y_arm_list, SIMPLIFY = FALSE)
  
  p_vec <- c(sum(as.numeric(p_arm_list[[2]]>p_arm_list[[1]]))/k,
             sum(as.numeric(p_arm_list[[3]]>p_arm_list[[1]]))/k,
             sum(as.numeric(p_arm_list[[4]]>p_arm_list[[1]]))/k)
  
  char <- paste("Arm",which.max(p_vec)+1,"is best") ## +1 because indexing of arms starts with 1.
  return(list(char,p_vec,rep(N/4,4),p_arm_list))
  #return(c(p_vec,max(p_vec),rep(N/4,4)))
}

# Comparing.

trt_effect <- c(h0,h1,h2,h3)
#y_mat <- as.matrix(mapply(rbinom_sampler, rep(N/4,4), trt_effect, SIMPLIFY = FALSE))

bench::mark("Jonathan's solution"=design1(y,h0=h0,h1=h1,h2=h2,h3=h3),
            "My submission"=design1_func(),
            relative = TRUE, check = FALSE)

## With parallelization

k=1000
N=228
seed=100
alpha=0.35
beta=0.65
trt_effect=c(0.35, 0.35, 0.35,0.35)

library(parallel)
cores <- parallel::detectCores()
cl <- parallel::makePSOCKcluster(cores)  
parallel::clusterSetRNGStream(cl, 123)
parallel::clusterExport(cl, c("k", "N", "seed", "alpha", "beta", "trt_effect", "rbinom_sampler",
                              "rbeta_sampler"))
result_parallel <- parSapply(cl,1:1000,design1_func)

## Comparing parallel and non-parallel version

bench::mark("normal"=replicate(1000, design1_func),
            "parallel"=parSapply(cl,1:1000,design1_func),
            relative = TRUE, check = FALSE)

# Surprisingly, the parallelization doesn't help in my case.