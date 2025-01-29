#' @title Choose the appropriate k cluster
#' @description Choose k via AIC or BIC
#' @param data The data with numeric elements
#' @param lower the lower bound k(default is 2)
#' @param upper the upper bound k(default is 5)
#' @param criterion two choice:"AIC","BIC"(default is "BIC")
#' @details The em_GMM within the upper and lower bounds of the given integer k are performed.
#' If there is an inability to converge or if the sigma is too low,
#' the corresponding message will be displayed on the way. see em_GMM.
#' @return A list of the result:
#' \item{k}{The optimal k via chosed criterion}
#' \item{criterion.value}{The criterion value corresponds to the optimal k}
#' @examples
#'   data(faithful)
#'   optimal_k(faithful,lower=2,upper=6,criterion="AIC")
#' @export
optimal_k = function(data,lower=2,upper=5,criterion="BIC",iter=1000){
  id = 0
  value = c()
  for (i in lower:upper){
    id[i-lower+1]=i
    value[i-lower+1] = em_GMM(data,k=i,max.iter = iter)[[as.character(criterion)]]
  }
  return(list(k=id[which.min(value)],criterion.value=min(value)))
}
