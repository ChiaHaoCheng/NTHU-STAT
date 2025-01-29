#' @title GMM by EM algorithm
#' @description The estimation of parameters for Gaussian mixture model by EM algorithm.
#' @param X X:n x p data matrix or data.frame with numeric elements
#' @param k k:the cluster you specified(nonzero integer)
#' @param max.iter the maximum iteration number(default is 100)
#' @details After the run, if the number of variants belonging to a group is too low,
#'  the message "The s.d.s are unstable." is displayed and the calculation is stopped,
#'  and the current result is recorded.
#'   In addition, if the specified argument "max.iter" fails to complete the convergence,
#'   the message "The iteration reaches maximum iterations" is displayed.
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats dnorm
#' @importFrom stats kmeans
#' @return A list of the result:
#' \item{data.matrix}{The input X}
#' \item{weight}{The vector of weights(k elements)}
#' \item{mu}{The gaussian mean by column for each cluster}
#' \item{sigma}{The variance-covariance matrix for each cluster}
#' \item{density}{The mixture density for all observations}
#' \item{BIC}{The value of BIC}
#' \item{AIC}{The value of AIC}
#' \item{count}{The number of iteration}
#' \item{logL.iter}{The log likelihood obtained by each iteration}
#' \item{data.dim}{The dimension of X,n x p}
#' \item{logL}{The last log likelihood}
#' \item{cluster}{The vector of n clusters assigned to each observation}
#' \item{success}{0:converge;1:No convergence}
#' @examples
#'   data(faithful)
#'   em_GMM(X=faithful,k=2,max.iter=10000)
#' @export
em_GMM <- function(X, k,max.iter=100){
  # initial parameters
  para.init <- function(X, k){
    x = as.matrix(X)
    d <- dim(x)[2]
    w <- rep(1/k, k)

    init_mu <- matrix(0, d, k)
    init_mu <- t(stats::kmeans(x, centers = k)$centers)

    init_sigma <- list()
    if (d!=1){
      for (i in 1:k){
        init_sigma[[i]] <- diag(d)
      }
    }
    else{
      for (i in 1:k){
        init_sigma[[i]] <- 1
      }
    }

    par <- list(weight = w, mu = init_mu, sigma = init_sigma)
    return(par)
  }
  # E step
  E.step <- function(X, par, k){
    x=as.matrix(X)
    n <- dim(x)[1];d=dim(x)[2]
    r <-  matrix(NA, nrow=n, ncol= k)
    L <- rep(0,n)
    if (d != 1){
      for (i in 1:k){
        dmv <- mvtnorm::dmvnorm(x, par$mu[, i], par$sigma[[i]])
        comp <- dmv * par$weight[i]
        L <- L + comp
      } #obtain likelihood

      for (i in 1:k){
        dmv <- mvtnorm::dmvnorm(x, par$mu[, i], par$sigma[[i]])
        r.est <- dmv * par$weight[i] / L
        r[, i] <- r.est
      } #obtain  E(r_k|X)
    }
    else if(d==1){
      for (i in 1:k){
        dmv <- stats::dnorm(x, mean=par$mu[i], sd=sqrt(par$sigma[[i]]))
        comp <- dmv * par$weight[i]
        L <- L + comp
      } #obtain likelihood

      for (i in 1:k){
        dmv <- stats::dnorm(x, mean=par$mu[i], sd=sqrt(par$sigma[[i]]))
        r.est <- dmv * par$weight[i] / L
        r[,i] <- r.est
      } #obtain  E(r_k|X)
    }
    log.GMM <- log(L)
    logL <- sum(log.GMM)

    list(r=r, logLikelihood=logL)
  }
  # M step
  M.step <- function(X, r, k){
    x = as.matrix(X)
    n <- dim(x)[1];d <- dim(x)[2]
    w.t <- rep(0, k)
    mu.t <- matrix(0, nrow=d, ncol = k)
    sigma.t <- list()
    if(d!=1){
      for (i in 1:k){
        new.w <- sum(r[, i])/n
        w.t[i] <- new.w

        new.mu <- 1/sum(r[, i])*colSums(r[, i]* x)
        mu.t[, i] <- new.mu

        new.mu <- matrix(rep(new.mu, n), n, d, byrow = T)
        new.sigma <- t(r[, i]*(x - new.mu)) %*% (r[, i]*(x - new.mu)) * 1/sum(r[, i])
        sigma.t[[i]] <- new.sigma
      }
    }
    else if(d==1){
      for (i in 1:k){
        new.w <- sum(r[, i])/n
        w.t[i] <- new.w

        new.mu <- sum(r[,i]* x)/sum(r[,i])
        mu.t[1, i] <- new.mu

        new.sigma <- sum(r[, i]*(x - mu.t[1, i])^2)/sum(r[, i])
        sigma.t[[i]] <- new.sigma
      }
    }
    para.next <- list(weight = w.t, mu = mu.t, sigma = sigma.t)
    return(para.next)
  }
  x <- as.matrix(X)
  para.0 <- para.init(x, k) #initial para.
  Estep.0 <- E.step(x, para.0, k)  #initial E step
  para.t <- M.step(x, Estep.0$r, k) #t=1 para
  current.logL <- Estep.0$logLikelihood
  logL.all <- Estep.0$logLikelihood
  logL.diff = 0.01
  count = 1
  n= dim(x)[1] ; d = dim(x)[2]
  if (d != 1){
    while(logL.diff > 1e-06){
      Estep.t <- E.step(x, para.t, k)
      para.t <- M.step(x, Estep.t$r, k)
      count <- count + 1
      if (count > max.iter){
        cat("The iteration reaches maximum iterations")
        break
      }
      logL.all <- c(logL.all, current.logL)
      logL.diff <- abs((current.logL - Estep.t$logLikelihood))
      current.logL <- Estep.t$logLikelihood
    }
  }
  else if(d==1){
    error <- c(1,1) #mu and sigma 's diff.
    while(error[1]>10^(-6) | error[2]>10^(-6)){
      mu.t = para.t$mu;sigma.t=as.numeric(para.t$sigma)
      Estep.t <- E.step(x, para.t, k)
      para.t <- M.step(x, Estep.t$r, k)
      count <- count + 1
      if (count > max.iter){
        cat("The iteration reaches maximum iterations","\n")
        break
      }
      else if(sum(sigma.t < 0.000001)>=1){
        cat("The s.d.s are unstable.","\n")
        break
      }
      logL.all <- c(logL.all, current.logL)
      error[1] <- sum(abs(mu.t-para.t$mu))
      error[2] <- sum(abs(sigma.t-as.numeric(para.t$sigma)))
      current.logL <- Estep.t$logLikelihood
    }
  }
  cluster <- factor(apply(Estep.t$r, 1, which.max))
  mix.f <- rep(0,dim(x)[1])
  if(d!=1){
    converg.or.not = ifelse(count < max.iter ,0,1)
    for (i in 1:k){
      dmv <- mvtnorm::dmvnorm(x, para.t$mu[,i], para.t$sigma[[i]])
      comp <- dmv * para.t$weight[i]
      mix.f <- mix.f + comp
    }
  }
  else if(d==1){
    converg.or.not = ifelse(sum(sigma.t < 0.000001)>=1 ,1,0)
    for (i in 1:k){
      dmv <- stats::dnorm(x, para.t$mu[1, i], para.t$sigma[[i]])
      comp <- dmv * para.t$weight[i]
      mix.f <- mix.f + comp
    }
  }
  BIC = -2*current.logL+(d*k/2*(d+3)+k-1)*log(n)
  AIC = -2*current.logL+(d*k/2*(d+3)+k-1)*2
  list(data.matrix=x,weight=para.t$weight,mu=para.t$mu,sigma=para.t$sigma,
       density=mix.f,BIC=BIC,AIC=AIC, count=count,
       logL.iter=logL.all, data.dim = dim(x),
       logL=current.logL, cluster=cluster, success=converg.or.not)
}
