#' @title The plot of Gaussian Mixture model
#' @description Show the result of em_GMM.
#' @param gm The object from em_GMM
#' @param type The type of plot,two choice:"hist","pairs".
#' @details For non-universal dataset,Only histogram with density line is displayed,
#' regardless of what argument type is entered.
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats dnorm
#' @return A plot specified type
#' @examples
#'   data(faithful)
#'   test=em_GMM(X=faithful,k=2,max.iter=10000)
#'   plot_GMM(gm=test,type="hist)
#' @export
plot_GMM <- function(gm,type="pairs"){
  d = gm$data.dim[2]
  if(d==1){
    mix.f <- rep(0,gm$data.dim[1])
    x=sort(gm$data.matrix)
    for (j in 1:length(gm$weight)){
      dmv <- stats::dnorm(x, mean=gm$mu[1,j], sd=sqrt(gm$sigma[[j]]))
      comp <- dmv * gm$weight[j]
      mix.f <- mix.f + comp
    }
    hist(gm$data.matrix, probability=TRUE,20,main="",xlab=colnames(gm$data.matrix))
    lines(x,mix.f,col="red")
  }
  else if(type == "hist"){
    ifelse(d < 3,par(mfrow=c(1,d)),par(mfrow=c(d-3+1,3)))
    for (i in 1:d){
      mix.f <- rep(0,gm$data.dim[1])
      x=sort(gm$data.matrix[,i])
      for (j in 1:length(gm$weight)){
        dmv <- stats::dnorm(x, mean=gm$mu[i,j], sd=sqrt(gm$sigma[[j]][i,i]))
        comp <- dmv * gm$weight[j]
        mix.f <- mix.f + comp
      }
      hist(gm$data.matrix[,i], probability=TRUE,20,main="",xlab=colnames(gm$data.matrix)[i])
      lines(x,mix.f,col="red")
    }
  }
  else if(d != 1 & type =="pairs"){
    panel.lower <- function(x,y,gmm=gm){
      points(x,y,pch=16)
      grid1 <- seq(min(x), max(x), length = 50)
      grid2 <- seq(min(y), max(y), length = 50)
      xy <- expand.grid(x=grid1, y=grid2)
      idx = idy = 0
      n = gmm$data.dim[1]
      for (i in 1:d){
        if(sum(x==gmm$data.matrix[,i]) == n) {idx <- i}
        if(sum(y==gmm$data.matrix[,i]) == n) {idy <- i}
      }
      mix.f <- function(grid,gmm) {
        xy=as.matrix(grid)
        mix.f <- rep(0,dim(xy)[1])

        for (i in 1:length(gmm$weight)){
          mu = gmm$mu[c(idx,idy),i]
          sigma = matrix(0,2,2)
          sigma[1,1] = gmm$sigma[[i]][idx,idx]
          sigma[1,2] = sigma[2,1] = gmm$sigma[[i]][idx,idy]
          sigma[2,2] = gmm$sigma[[i]][idy,idy]
          dmv <- mvtnorm::dmvnorm(xy, mean=mu, sigma=sigma)
          comp <- dmv * gmm$weight[i]
          mix.f <- mix.f + comp
        }
        mix.f
      }
      z <- mix.f(xy,gm)
      contour(grid1, grid2, matrix(z, 50, 50), levels = round(quantile(z, 0.05*(1:19)), 3)
              , add = T, col = "red",drawlabels=T)
    }
    panel.upper <- function(x, y) {
      par(usr = c(0, 1, 0, 1))
      r <- round(cor(x, y), digits = 4)
      text(0.5, 0.5, "Correlation:",cex=1.5)
      text(0.5, 0.3, r,cex=1)
    } #show correlation
    pairs(gm$data.matrix,
          lower.panel = panel.lower,
          upper.panel = panel.upper)
  }
}
