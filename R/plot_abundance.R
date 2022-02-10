#' Plot of feature abundance level
#'
#' This function plots the feature abundance level for each cluster, after extracting the effect of sample-wise normalization factors and feature-wise geometric mean.
#'
#' @param result Clustering result from function PHclust().
#' @param data Data matrix with dimension N*P indicating N features and P samples.
#' @param Treatment Vector of length P. Indicating replicates of different treatment groups. For example, \emph{Treatment} = c(1,1,2,2,3,3) indicates 3 treatment groups, each with 2 replicates.
#'
#' @return A plot for feature abundance level will be shown. No value is returned.
#'
#' @importFrom stats cor optimize pchisq quantile rmultinom
#' @importFrom graphics par lines
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
plot_abundance = function(result, data, Treatment){
  oldpar <- par(no.readonly = TRUE)

  i0 = c(match(unique(Treatment), Treatment), length(Treatment) + 1)
  s = result$Normalizer
  alpha = result$alpha
  final = result$cluster
  nK = length(unique(final))

  datas = sweep(data, 2, exp(s), '/')
  datas = sweep(datas, 1, exp(alpha), '/')

  all_m = matrix(0, nrow = length(unique(final)), ncol = length(unique(Treatment)))
  each_m = matrix(0, nrow = nrow(data), ncol = length(unique(Treatment)))

  for(j in 1:length(unique(Treatment))){
    dataj = datas[, (i0[j]):(i0[j+1] - 1)]
    dataj = data.matrix(dataj)
    for(i in 1:length(unique(final))){
      dataij = dataj[final == i, ]
      dataij = data.matrix(dataij)
      all_m[i, j] = mean(dataij)
    }
    each_m[, j] = apply(dataj, 1, mean)
  }

  k1 = ceiling(sqrt(nK))
  k2 = ceiling(nK/k1)

  par(mfrow = c(k2, k1))
  for(i in 1:length(unique(final))){
    temp = all_m[i, ]

    plot(temp, type = 'l', col = 'black', ylim = c(min(temp) - (max(temp)-min(temp))*0.5, max(temp) + (max(temp)-min(temp))*0.5),
         xlab = 'Treatment', ylab = 'Expression', cex.lab = 1.25)
    for(j in which(final == i)){
      lines(each_m[j, ], col = c('grey', 0.05))
    }
    lines(all_m[i, ], col = 'black')
  }

  on.exit(par(oldpar))
  print('Plot of feature abundance level for each cluster:')
}

