#' Poisson hurdle clustering
#'
#' This function gives the clustering result based on a Poisson hurdle model.
#'
#' @param data Data matrix with dimension N*P indicating N features and P samples. The cluster analysis is done feature-wised.
#' @param Treatment Vector of length P. Indicating replicates of different treatment groups. For example, \emph{Treatment} = c(1,1,2,2,3,3) indicates 3 treatment groups, each with 2 replicates.
#' @param nK Positive integer. Number of clusters.
#' @param method Method for the algorithm. Can choose between \emph{"EM"} as Expectation Maximization or \emph{"SA"} as Simulated Annealing.
#' @param absolute Logical. Whether we should use absolute (TRUE) or relative (False) abundance of features to determine clusters.
#' @param cool Real number between (0, 1). Cooling rate for the \emph{"SA"} algorithm. Uses 0.9 by default.
#' @param nstart Positive integer. Number of starts for the entire algorithm. Note that as \emph{nstart} increases the computational time also grows linearly. Uses 1 by default.
#'
#' @return
#' \describe{
#' \item{cluster}{Vector of length N consisting of integers from 1 to nK. Indicating final clustering result. For evaluating the clustering result please check \link[aricode]{NMI} for \emph{Normalized Mutual Information}.}
#' \item{prob}{N*nK matrix. The (i, j)th element representing the probability that observation i belongs to cluster j.}
#' \item{log_l}{Scaler. The Poisson hurdle log-likelihood of the final clustering result.}
#' \item{alpha}{Vector of length N. The geometric mean abundance level for each feature, across all treatment groups.}
#' \item{Normalizer}{vector of length P. The normalizing constant of sequencing depth for each sample.}
#' }
#'
#' @importFrom stats cor optimize pchisq quantile rmultinom
#' @export
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
PHcluster = function(data, Treatment, nK, method = c('EM', 'SA'), absolute = FALSE, cool = 0.9, nstart = 1){
  mydata = find_norm(data, Treatment)

  l0 = c()
  f0 = list()
  Z0 = list()
  alpha0 = list()

  dis = dis_tau(mydata)

  for(tr in 1:nstart){
    starting = initial(mydata, nK, dis)
    q0 = starting$q0
    mu0 = starting$mu0
    fn = hp_cluster(mydata, q0, mu0, method = method, absolute = absolute, cool = cool)
    l0 = c(l0, fn$lglk)
    f0[[tr]] = fn$final
    Z0[[tr]] = fn$Z
    alpha0[[tr]] = fn$alpha
  }

  final = f0[[which.max(l0)]]
  Z1 = Z0[[which.max(l0)]]
  l1 = max(l0)
  alpha1 = alpha0[[which.max(l0)]]

  s = mydata$Normalizer
  Treatment = mydata$Treatment
  i0 = c(match(unique(Treatment), Treatment), length(Treatment) + 1)

  alpha = rep(0, nrow(alpha1))
  for(i in 1:length(unique(final))){
    alphai = alpha1[final == i, ]
    if(is.vector(alphai)) alpha[final == i] = alphai[i]
    else alpha[final == i] = as.vector(alphai[, i])
  }

  return(list(prob = Z1, cluster = final, log_l = l1, alpha = alpha, Normalizer = s))
}







