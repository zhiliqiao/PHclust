#' Calculate optimal number of clusters.
#'
#' This function estimates the optimal number of clusters for a given dataset.
#'
#' @param data Data matrix with dimension N*P indicating N features and P samples.
#' @param absolute Logical. Whether we should use absolute (TRUE) or relative (FALSE) abundance of features to determine clusters.
#' @param Kstart Positive integer. The number of clusters for starting the hybrid merging algorithm. Should be relatively large to ensure that Kstart > optimal number of clusters. Uses \emph{max(50, sqrt(N))} by default.
#' @param Treatment Vector of length p, indicating replicates of different treatment groups. For example, \emph{Treatment} = c(1,1,2,2,3,3) indicates 3 treatment groups, each with 2 replicates.
#'
#' @return A positive integer indicating the optimal number of clusters
#'
#' @export
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#'
#' @examples ######## Run the following codes in order:
#' @examples ##
#' @examples ## This is a sample data set which has 100 features, and 4 treatment groups with 4 replicates each.
#' @examples data('sample_data')
#' @examples head(sample_data)
#' @examples set.seed(1)
#' @examples ##
#' @examples ## Finding the optimal number of clusters
#' @examples K <- Hybrid(sample_data, Kstart = 4, Treatment = rep(c(1,2,3,4), each = 4))
#' @examples ##
#' @examples ## Clustering result from EM algorithm
#' @examples result <- PHcluster(sample_data, rep(c(1,2,3,4), each = 4), K, method = 'EM', nstart = 1)
#' @examples print(result$cluster)
#' @examples ##
#' @examples ## Plot the feature abundance level for each cluster
#' @examples plot_abundance(result, sample_data, Treatment = rep(c(1,2,3,4), each = 4))
Hybrid = function(data, absolute = FALSE, Kstart = NULL, Treatment){
  if(is.null(Kstart)){Kstart = min(max(floor(sqrt(nrow(data))), 50), nrow(data))}
  # if(length(unique(Treatment)) == 1){absolute = TRUE}

  i0 = c(match(unique(Treatment), Treatment), length(Treatment) + 1)
  mydata = find_norm(data, Treatment)
  s = mydata$Normalizer
  dis = dis_tau(mydata)

  starting = initial(mydata, Kstart, dis)
  q0 = starting$q0
  mu0 = starting$mu0

  fn = hp_cluster(mydata, q0, mu0, method = 'EM', absolute = absolute)

  final = fn$final
  mu = fn$mu
  q = fn$q
  alpha = fn$alpha

  steps = rep(0, Kstart-1)
  r0 = rep(0, Kstart-1)

  nn = c()

  try = function(i, j, final){
    dataij = data[(final %in% c(i,j)), ]
    alphaij0 = alpha[(final %in% c(i,j))]
    a = g(dataij, s, alphaij0, Treatment, i0)
    alphaij = a$alpha
    muij = a$mu
    qij = a$q
    l = a$l

    return(list(alpha = alphaij, mu = muij, q = qij, l = l))
  }

  for(ind in 1:Kstart){
    if(ind == 1){

      l0 = c()
      alpha = rep(0, nrow(data))

      for(i in 1:length(unique(final))){
        datai = data[final == i, ]
        if(is.vector(datai)) datai = matrix(datai, nrow = 1)
        else datai = as.matrix(datai, ncol = length(Treatment))
        if(is.null(dim(mu))) mu = matrix(mu, ncol = 1)
        alphai = (Est_alpha(s, datai, mu, Treatment))[, i]
        a = g(datai, s, alphai, Treatment, i0)
        alpha[final == i] = a$alpha

        l0 = c(l0, a$l)
      }

      l = matrix(1e40, nrow = length(l0), ncol = length(l0))
      for(i in 2:length(l0)){
        for(j in 1:(i-1)){
          l[i,j] = l0[i] + l0[j] - try(i, j, final)$l
        }
      }
    }

    min_l = min(l)
    merg = as.vector(which(l == min_l, arr.ind = TRUE))

    k1 = which(merg[1] == levels(final))
    k2 = which(merg[2] == levels(final))
    r = nrow(data[(final %in% merg), ])

    nn = c(nn, merg[2])

    l[merg[2], ] = 1e40
    l[, merg[2]] = 1e40

    steps[length(steps) - ind + 1] = min_l
    r0[length(r0) - ind + 1] = r

    mu[k1,] = try(merg[1], merg[2], final)$mu
    mu = mu[-k2, ]
    final[final %in% merg] = merg[1]
    final = as.factor(as.numeric(as.character(final)))

    if (length(levels(final)) == 1) break

    l0[merg[2]] = 0
    i = merg[1]
    datai = data[final == i, ]
    if(is.vector(datai)){datai = matrix(datai, nrow = 1)} else{datai = as.matrix(datai, ncol = length(Treatment))}
    if(is.null(dim(mu))) mu = matrix(mu, ncol = 1)
    alphai = (Est_alpha(s, datai, mu, Treatment))[, which(i == as.numeric(levels(final)))]
    a = g(datai, s, alphai, Treatment, i0)
    l0[merg[1]] = a$l

    for(j in 1:(merg[1]-1)){
      if(sum(j == nn) == 0){
        l[merg[1], j] = l0[j] + l0[merg[1]] - try(j, merg[1], final)$l
      }
    }

    if(i != Kstart){
      for(j in (merg[1]+1):nrow(l)){
        if(sum(j == nn) == 0){l[j, merg[1]] = l0[j] + l0[merg[1]] - try(j, merg[1], final)$l}
      }
    }

  }

  I = length(unique(Treatment))

  p_v = c()
  for(i in 1:length(steps)){
    temp = 1 - pchisq(steps[i], df = r0[i] + 2*I-1)
    p_v = c(p_v, temp)
  }

  if(sum(p_v > 0.01) > 0) k = min(which(p_v > 0.01))
  else k = length(p_v) + 1

  return(k)
}




