# # Copyright (c) 2019-2020 Carlo Manzo & Tina Košuta 
#
# If you use this code, please cite: 
# (http://dx.doi.org/10.1039/C9CP05616E)
# (https://arxiv.org/abs/1909.13133)
#
# Permission is hereby granted, free of charge, to any person 
# obtaining a copy of this software and associated documentation files 
# (the "Software"), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, 
# publish, distribute, sublicense, and/or sell copies of the Software, 
# and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
  
# The above copyright notice and this permission notice shall be 
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#Upload data
data <- read.table("sample_data.txt")
data <- as.numeric(data[,]) #optional; number of localizations is needed as vector
nb <- table(data)

#mean of logN
p1_1 = 3.3491 

#sd of logN
p1_2 = 0.8462  

############### FUNCTIONS

library(gtools) #rdirichlet
library(pracma) #for erf and fmincom

#calibration curve
pp <- function(x,mu,sigma)  { 0.5*(erf((mu-log(x-1))/(sqrt(2)*sigma))-erf((mu-log(x))/(sqrt(2)*sigma))) }

#log likelihood
log_likelihood <- function(data, theta){
  lambda=ft_u_mat%*%theta
  logL=data%*%log(lambda) 
  return(logL)
}

############### MODEL PARAMETERS

# Number of particles
N <- 30

# MCMC steps per NS iteration
mcmc_steps <- 40

# delta - shape of Dirichlet distribution
delta <- 1.5

# maximum number of parameters in the model we expect
num_par <- 8

############### NESTED SAMPLING ALGORITHM

results <- list()

for(Nmax in 2:num_par){
  
    #calibration
    ft_u_mat <- NULL
    pp_ll <- seq(1, 3*max(data) ,1)
    pp_ft_mat <- matrix(0, ncol=Nmax, nrow=pp_ll)
    pp_ft <- pp(pp_ll, p1_1, p1_2)
    pp_ft_mat <- sapply(1:Nmax, function(x) Re(fft((fft(pp_ft)^x)/length(pp_ft), inverse=TRUE)))
    
    ll_u <- sort(unique(data)) # subset of our values (for our data)
    ft_u_mat <- pp_ft_mat[ll_u,]
    
    #deltas for each parameter in Dirichlet
    delta_dir <- rep(delta, Nmax) #for uniform prior
    
    points <- rdirichlet(N, delta_dir) #each column is new particle
    logL <- sapply(1:N, function(x) log_likelihood(data=nb, theta=points[x,]))
    logP <- log(ddirichlet(points, delta_dir))
    
    logwt <- 0 #Outermost interval of prior mass
    logZnew <- 0
    keep <- list()
    keep2 <- NULL
    Zrat <- Inf
    j <- 0
    
    while(Zrat > log(1e-5)){    #turn the inequallity to create condition
      
      j <- j+1
      
      min <- which.min(logL) #lowest likelihood
      worst_logl <- logL[min]
      logwt <- (-j*log(1+1/N) -log(N)) +worst_logl #Calculate weight of worst walker
      
      #data of first walker, Lmin, logZ, w
      keep[[j]] <- c(points[min,], logL[min])
      keep2[j] <- logwt
      
      threshold <- logL[min]
      
      Wnew <- sample(as.numeric(1:N)[-min],1) #remove walker with logL min and
      points[min,] <- points[Wnew,] #copy new walker; then walk the coppied walker on the place of walker min
      logL[min] <- logL[Wnew]
      logP[min] <- logP[Wnew]
      
      sig <- 0.1
      accept <- 0
      
      for (i in 1:mcmc_steps){
        
        k <- sample(1:Nmax,1)
        
        new <- points[min,]
        step_size <- min(rnorm(1,0,sig), 1-1e-6)
        new[k] <- ifelse(new[k]+step_size>=1, 1-(new[k]+step_size-1), abs(new[k]+step_size))
        new <- new/sum(new)
  
        lambda_new <- ft_u_mat%*%new
        logl_new <- nb%*%log(lambda_new) #the only changing part of likelihood; LOG LIIKELIHOOD
        logP_new <- log(ddirichlet(matrix(new, nrow=1), delta_dir))
        logPdif <- logP_new-logP[min]
        
        if(logl_new >= threshold && runif(1) <= exp(logPdif)){
          points[min,] <- new
          logL[min] <- logl_new
          logP[min] <- logP_new
          accept  <-  accept+1
        }
        
        sig <- ifelse(accept>(i-accept), sig*exp(1/accept), sig/exp(1/(i-accept)))
        
        }
      
      logZnew <- log(sum(exp(keep2-max(keep2)))) + max(keep2)  #new logZ
      logZrem <- -j*log(1+1/N)-log(N)+(log(sum(exp(logL-max(logL)))) + max(logL))  #remaining evidence
      Zrat <- (logZrem-logZnew)
      
    }
    
      logws <- -(1:j)*log(1+1/N) #posterior weights
      logwt2 <- (-j*log(1+1/N)) + logL #calculate the weight of each walker with new voulmes
      
      keep2 <- c(keep2, logwt2) 
      
      keep <- do.call(rbind, keep)
      keep <- rbind(keep, cbind(points,logL))  #add remaining walkers
      
      # Calculate marginal likelihood
      logZ <- log(sum(exp(keep2-max(keep2)))) + max(keep2) #final logZ
      
      post_weights <- exp(keep2 - logZ) #normalised posterior weights
      
      #H calculation
      H <- sum(post_weights*(keep[,ncol(keep)]-logZ))
  
      #Error
      error <- sqrt(H/N)
      
      #N samples from posterior
      ent <- -sum(post_weights*log(post_weights + 1E-300))
      Npost <- floor(exp(ent)) #effective size for posterior samples
      
      # Create posterior samples by resampling; importance weigth
      posterior_samples_ <- list()
      posterior_weights <- NULL
      
      # Counter
      k <- 1
      top <- max(post_weights)
      while(k<=Npost){
        # Choose one of the samples
        which <- sample(1:nrow(keep), 1)
        # Acceptance probability
        prob <- post_weights[which]/top
        if(runif(1) <= prob){
          posterior_samples_[[k]] <- keep[which, ]
          posterior_weights[k] <- post_weights[which]
          k <- k + 1
        }
      }
      
      posterior_samples <- do.call(rbind, posterior_samples_)
      
      #weighted mean from posterior sample
      Wm <- apply(posterior_samples[,-ncol(posterior_samples)]*posterior_weights,2,sum)/sum(posterior_weights)

    results[[Nmax]] <- list("N"=Nmax,"logZ"=logZ,"Wm"=Wm, "ent"= ent, "Npost"=Npost, "H"=H, "error"=error)      
}

#compare logZ for all models
logZ_results <- data.frame("Number of parameters"=sapply(2:num_par, function(x) results[[x]][["N"]]),
                            "logZ"=sapply(2:num_par, function(x) results[[x]][["logZ"]]))


#print results for model with largest logZ
results[[which.max(logZ_results[,"logZ"])]]

