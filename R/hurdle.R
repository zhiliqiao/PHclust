#' Title
#'
#' @param x a
#' @param by a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
sumRow=function(x,by=NULL){
  if(is.vector(x)) x=matrix(x,ncol=length(x))
  if(is.null(by)) by=rep(1,ncol(x))
  by=as.ordered(as.factor(by))
  mt=outer(by,unique(by),"==")
  xm=x%*%mt
  if(nrow(xm)==1 | ncol(xm)==1) xm=c(xm)
  return(xm)
}



#' Title
#'
#' @param n a
#' @param Treatment a
#' @param Z_mat a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
Est_q = function(n, Treatment, Z_mat, s){
  i0 = c(match(unique(Treatment), Treatment), length(Treatment) + 1)

  q = matrix(0, nrow = ncol(Z_mat), ncol = ncol(n))
  for(k in 1:ncol(Z_mat)){
    z_k = Z_mat[, k]
    f_max = function(x){
      a = x[1]
      b = x[2]
      - sum(sum(log(1+exp(a + b*s)))*z_k) +
        sum(outer(z_k, a + b*s, '*') * (n != 0))
    }
    f_max_gr = function(x){
      a = x[1]
      b = x[2]
      temp = exp(a + b*s) / (1 + exp(a + b*s))
      gr_a = sum(rowSums(n != 0) * z_k) - sum(sum(temp) * z_k)
      gr_b = sum(outer(z_k, s, '*') * (n != 0)) - sum(sum(temp * s) * z_k)
      c(gr_a, gr_b)
    }

    opt_para = optim(c(1,1), f_max, f_max_gr, method = 'BFGS', control=list(fnscale=-1))$par
    q[k, ] = 1 / (1 + exp(- opt_para[1] - opt_para[2]*s))
  }

  q
}




#' Title
#'
#' @param s a
#' @param n a
#' @param mu a
#' @param Treatment a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
Est_alpha = function(s, n, mu, Treatment){
  if(is.null(dim(mu))){mu = matrix(mu, ncol = 1)}
  a1 = rowSums(n)
  a1 = matrix(rep(a1, nrow(mu)), ncol = nrow(mu), byrow = FALSE)

  b1 = matrix(0, nrow(a1), ncol(a1))
  t = as.vector((summary(as.factor(Treatment)))[as.factor(unique(Treatment))])
  for(k in 1:ncol(b1)){
    mu_temp = rep(mu[k, ], times = t)
    b1[, k] = rowSums(sweep((n > 0), 2, exp(s + mu_temp), '*'))
  }

  log(a1/b1)
}










#' Title
#'
#' @param Treatment a
#' @param s a
#' @param ng a
#' @param mu a
#' @param alpha a
#' @param q a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
l_hp = function(Treatment, s, ng, mu, alpha, q){
  t = as.vector((summary(as.factor(Treatment)))[as.factor(unique(Treatment))])

  mu_temp = rep(mu, times = t)

  lambda = exp(s + alpha + mu_temp)
  temp0 = lambda
  temp0[abs(log(lambda)) < log(100)] = log(exp(lambda[abs(log(lambda)) < log(100)]) - 1)
  temp0[lambda < 1/100] = log(lambda[lambda < 1/100])

  l0 = log(1 - q)
  l1 = log(q) + ng * (alpha + mu_temp) - temp0
  sum(l0 * as.numeric(ng == 0) + l1 * as.numeric(ng > 0))
}



#' Title
#'
#' @param Treatment a
#' @param sg a
#' @param ng a
#' @param mu a
#' @param alpha_g a
#' @param q a
#' @param p a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
Z_g = function(Treatment, sg, ng, mu, alpha_g, q, p){
  l = c()
  for(k in 1:length(alpha_g)){
    l = c(l, l_hp(Treatment, sg, ng, mu[k, ], alpha_g[k], q[k, ]))
  }
  l = l - max(l)

  f = exp(l) * p
  g = f/sum(f)
  if(sum(is.na(g))>0) g = rep(1/length(g), length(g))
  return(g)
}



#' Title
#'
#' @param Treatment a
#' @param s a
#' @param n a
#' @param mu a
#' @param alpha a
#' @param p a
#' @param q a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
Expect_Z = function(Treatment, s, n, mu, alpha, p, q){
  Z_mat = matrix(0, nrow(n), nrow(q))
  for (g in 1:nrow(n)){
    Z_mat[g, ] = Z_g(Treatment, s, n[g, ], mu, alpha[g, ], q, p)
  }

  for(k in 1:ncol(Z_mat)){
    if((colSums(Z_mat == 0) == nrow(Z_mat))[k] == 1){
      Z_mat[, k] = 1e-100
    }
  }
  Z_mat = sweep(Z_mat, 1, rowSums(Z_mat), '/')
  return(Z_mat)
}



#' Title
#'
#' @param Treatment a
#' @param s a
#' @param n a
#' @param mu_k a
#' @param alpha a
#' @param q a
#' @param k a
#' @param Z_mat a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
l_hp_k = function(Treatment, s, n, mu_k, alpha, q, k, Z_mat){
  l_mu = 0
  for(g in 1:nrow(n)){
    l_mu = l_mu + Z_mat[g, k] * l_hp(Treatment, s, n[g, ], mu_k, alpha[g, k], q[k, ])
  }
  l_mu
}


#' Title
#'
#' @param Treatment a
#' @param s a
#' @param n a
#' @param mu_ki a
#' @param alpha a
#' @param q a
#' @param k a
#' @param i a
#' @param Z_mat a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
l_hp_ki = function(Treatment, s, n, mu_ki, alpha, q, k, i, Z_mat){
  i0 = c(match(unique(Treatment), Treatment), length(Treatment) + 1)
  data_i = n[, (i0[i]):(i0[i+1] - 1)]
  data_i = as.matrix(data_i, ncol = i0[i+1] - i0[i])
  lambda_i = exp(outer(alpha[, k], s[(i0[i]):(i0[i+1] - 1)], '+') + mu_ki)
  temp0 = lambda_i
  temp0[abs(log(lambda_i)) < log(100)] = log(exp(lambda_i[abs(log(lambda_i)) < log(100)]) - 1)
  temp0[lambda_i < 1/100] = log(lambda_i[lambda_i < 1/100])

  temp = rowSums((data_i > 0) * (sweep(data_i, MARGIN = 1, STATS = alpha[, k] + mu_ki, FUN = '*') - temp0))

  return(sum(Z_mat[, k] * temp))
}



#' Title
#'
#' @param data a
#' @param Treatment a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
find_norm = function(data, Treatment){
  n = data.matrix(data)
  n = matrix(n,ncol=length(Treatment))
  n = n[, order(Treatment)]
  n = n[rowSums(n)!=0, ]
  Treatment = Treatment[order(Treatment)]
  n[n<=0] = NA
  log_q3 = log(apply(n, 2, quantile, 0.75, na.rm = TRUE))
  log_q3 = log_q3 - mean(log_q3)
  n[is.na(n)] = 0
  return(list(Count = n, Treatment = Treatment, Normalizer = log_q3))
}



#' Title
#'
#' @param mydata a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
dis_tau = function(mydata){
  return(cor(t(mydata$Count), method = 'kendall'))
}



#' Title
#'
#' @param mydata a
#' @param nK a
#' @param dis a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
initial = function(mydata, nK, dis){
  Treatment = mydata$Treatment
  i0 = c(match(unique(Treatment), Treatment), length(Treatment) + 1)
  if((sum(is.na(dis)) > 0) == TRUE){index0 = sample(1:nrow(mydata$Count), nK, replace = FALSE)}
  else{
    index0 = sample(1:nrow(mydata$Count), 1)
    for(i in 1:(nK - 1)){
      dis0 = dis[index0, ]
      if(is.vector(dis0)){p = dis0}
      else{
        p = colSums(dis0)
      }

      repeat{
        temp_index = which.min(p)
        if(sum(temp_index == index0) == 0) break
        else p = p[-temp_index]
      }

      index0 = c(index0, temp_index)
    }
  }
  data0 = mydata$Count[index0, ]


  q0 = matrix(0, nK, length(unique(mydata$Treatment)))
  mu0 = matrix(0, nK, length(unique(Treatment)))

  for (i in 1:(length(i0)-1)){
    temp = data0[, (i0[i]):(i0[i+1]-1)]
    temp = as.matrix(temp, ncol = i0[i+1] - i0[i])

    q0[, i] = 1 - apply(temp==0, 1, sum)/ncol(temp)
    q0 = 0.5 + (q0 - 0.5) * (1-1e-10)

    s_temp = (mydata$Normalizer)[(i0[i]):(i0[i+1]-1)]
    mu0[, i] = log(rowSums(temp) / rowSums(exp(sweep((temp == 0) , 2, s_temp, '*'))))
  }
  mu0[!is.finite(mu0)]=0
  mu0 = sweep(mu0, 1, rowMeans(mu0), '-')

  q0 = matrix(0.5 + runif(nK * ncol(mydata$Count), -0.1, 0.1),
              nrow = nK, ncol = ncol(mydata$Count))

  return(list(q0 = q0, mu0 = mu0, index = index0))
}



#' Title
#'
#' @param mydata a
#' @param q a
#' @param mu a
#' @param method a
#' @param absolute a
#' @param cool a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
hp_cluster = function(mydata, q, mu, method = c('EM', 'SA'), absolute, cool = 0.9){
  s = mydata$Normalizer
  n = mydata$Count
  Treatment = mydata$Treatment

  p = rep(1/nrow(mu), nrow(mu))
  if(absolute == FALSE){alpha = Est_alpha(s, n, mu, Treatment)}


  iter = 0
  lglk_list = c(-Inf)
  Z_list = list()
  i0 = c(match(unique(Treatment), Treatment), length(Treatment) + 1)

  run_times = c() # computational time. Delete after use

  repeat{
    iter = iter + 1

    start.time = Sys.time() # computational time. Delete after use
    if(absolute == TRUE){alpha = matrix(0, ncol = nrow(mu), nrow = nrow(n))}
    Z_mat = Expect_Z(Treatment, s, n, mu, alpha, p, q)

    if(method == 'SA') {
      tem = 2*cool^(iter)
      Z_temp = Z_mat^(1/tem)
      z = Z_temp
      for(g in 1:nrow(Z_temp)){
        z[g, ] = t(rmultinom(1, 1, Z_temp[g, ]))
        z[g, ] = z[g, ] + 1e-5
        z[g, ] = z[g, ] / sum(z[g, ])
      }
      Z_mat = z
    }

    Z_list[[iter]] = Z_mat

    p_new = colSums(Z_mat)/nrow(n)

    q_new = Est_q(n, Treatment, Z_mat, s)

    mu_temp = matrix(0, nrow(mu), ncol(mu))
    for(k1 in 1:nrow(mu_temp)){
      for(i1 in 1:ncol(mu_temp)){
        Est1_mu = function(mu_ki){
          -l_hp_ki(Treatment, s, n, mu_ki, alpha, q_new, k1, i1, Z_mat)
        }
        mu_temp[k1, i1] = optimize(Est1_mu, lower = -100, upper = 100)$minimum
      }
    }

    if(absolute == FALSE){
      alpha = sweep(alpha, 2, rowMeans(mu_temp), '+')
      mu_new = sweep(mu_temp, 1, rowMeans(mu_temp), '-')
      alpha_new = Est_alpha(s, n, mu_new, Treatment)
    }
    else{
      alpha_new = matrix(0, ncol = nrow(mu), nrow = nrow(n))
      mu_new = mu_temp
    }

    mu = mu_new
    alpha = alpha_new
    p = p_new
    q = q_new

    lglk_temp = 0
    for(k in 1:nrow(mu_new)){
      fac_Treatment = factor(Treatment)
      levels(fac_Treatment) = q_new[k, ]
      qtemp = as.numeric(as.character(fac_Treatment))
      lglk_temp = lglk_temp + l_hp_k(Treatment, s, n, mu_new[k, ], alpha_new, q_new, k, Z_mat)
    }

    lglk_list = c(lglk_list, lglk_temp)


    prop_diff = abs(1 - lglk_list[length(lglk_list)]/lglk_list[length(lglk_list)-1])

    run_times = c(run_times, Sys.time() - start.time) # computational time. Delete after use

    if(iter > 10 & prop_diff<1e-5) break
  }

  Z_mat_final = Expect_Z(Treatment, s, n, mu, alpha, p, q)
  Z_list[[iter+1]] = Z_mat_final

  Z_opt = Z_list[[which.max(lglk_list)]]

  final = as.factor(as.numeric(apply(Z_opt, 1, which.max)))
  if(length(unique(final)) < ncol(Z_opt)) repeat{
    a = 1:ncol(Z_opt)
    b = a[! a %in% unique(final)]
    b = as.numeric(as.character(b))

    final = as.numeric(as.character(final))
    final[sample(final, length(b))] = b
    final = as.factor(final)

    if(length(unique(final)) == ncol(Z_opt)) break
  }

  lglk = max(lglk_list)

  return(list(final = final, lglk = lglk, lglk_l = lglk_list, Z = Z_opt, q = q, mu = mu, alpha = alpha,
              run_times = run_times))
}






#' Title
#'
#' @param data_i a
#' @param s_i a
#' @param Treatment a
#' @param alphaij a
#' @param mu_i a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
l_mu = function(data_i, s_i, Treatment, alphaij, mu_i){
  lambda = exp(outer(alphaij, s_i, '+') + mu_i)
  temp0 = lambda
  temp0[abs(log(lambda)) < log(100)] = log(exp(lambda[abs(log(lambda)) < log(100)]) - 1)
  temp0[lambda < 1/100] = log(lambda[lambda < 1/100])
  l = sum(data_i) * mu_i - sum((data_i > 0) * temp0)
  return(l)
}



#' Title
#'
#' @param n a
#' @param s a
#' @param alpha a
#' @param Treatment a
#' @param i0 a
#'
#' @return a
#' @importFrom stats cor optimize pchisq quantile rmultinom optim runif
#' @noRd
g = function(n, s, alpha, Treatment, i0){
  mu_m = rep(0, length(unique(Treatment)))
  q_m = rep(0, length(unique(Treatment)))
  for(o in 1:length(unique(Treatment))){
    data_i = n[, (i0[o]):(i0[o+1] - 1)]
    s_i = s[(i0[o]):(i0[o+1] - 1)]

    q_m[o] = sum(data_i > 0) / sum(data_i >= 0)
    q_m[o] = 0.5 + (q_m[o] - 0.5) * (1-1e-10)

    Est2_mu = function(mu){
      temp = l_mu(data_i, s_i, Treatment, alpha, mu)
      return(-temp)
    }

    if(sum(data_i>0) == 0) mu_m[o] = 0
    else mu_m[o] = optimize(Est2_mu, lower = -100, upper = 100)$minimum
  }
  alpha_m = as.vector(Est_alpha(s, n, matrix(mu_m, nrow = 1), Treatment))

  l=0
  for(g in 1:nrow(n)){
    l = l+l_hp(Treatment, s, n[g, ], mu_m, alpha_m[g], q_m)
  }

  return(list(alpha = alpha_m, mu = mu_m, q = q_m, l = l))
}








