

library(dplyr)
install.packages("lpSolve")
library(lpSolve)
load("wiener_holder_norms.txt")
load("tight_mres_norms_const.RData")
load("tight_mres_norms_lin.RData")
tmnc<-load("G:/MASTERS/THESIS/DATA/tight_mres_norms_const.RData",header=F)
head(whn)
summary(tmnc)
tmnl<-load("G:/MASTERS/THESIS/DATA/tight_mres_norms_lin.RData")
whn<-read.table("G:/MASTERS/THESIS/DATA/wiener_holder_norms1.txt")

xw <- wave2sect + rnorm(450)/2
ts.plot(xw, col="grey", ylab="")
lines(wave2sect)
draw_rects(xw.n.0, c(-1, 5), 20, "brown")

ddd


###################### nsp####################
nsp <- function(x, constraints, M, thresh, overlap = FALSE, buffer = 0) {
  d <- dim(constraints)
  
  x.c <- cbind(x, constraints)
  
  x.c.ads <- all_dyadic_scans_array(x.c)
  
  res <- iter_random_checks_scan_array(c(1, d[1]), x.c.ads, M, thresh, overlap, buffer)
  
  intervals <- data.frame(t(order_chron(res)))
  colnames(intervals) <- c("starts", "ends", "values")
  intervals$midpoints <- floor((intervals$starts+intervals$ends)/2)
  
  list(intervals=intervals, threshold.used=thresh)
} 

############# nsp_poly ###############################

nsp_poly <- function(x, M = 1000, thresh = NULL, sigma = NULL, alpha = 0.1, deg = 0, overlap = FALSE) {
  n <- length(x)
  
  x.c <- matrix(x, n, deg+1)
  
  for (i in 1:(deg+1)) {
    
    x.c[,i] <- seq(from=0, to=1, length=n)^(i-1)
    
    
  }
  
  if (is.null(thresh)) {
    
    if (is.null(sigma)) sigma <- mad(diff(x/sqrt(2)))
    if (is.null(alpha)) alpha <- 0.1
    
    thresh <- sigma * thresh_kab(n, alpha)
    
  }
  
  nsp(x, x.c, M, thresh, overlap)
  
}

################ nsp_poly_2 ########################

nsp_poly_2 <- function(x, M = 1000, thresh.type = "univ", thresh.val = NULL, sigma = NULL, alpha = 0.1, deg = 0, overlap = FALSE) {
  n <- length(x)
  
  x.c <- matrix(x, n, deg+1)
  
  for (i in 1:(deg+1)) {
    
    x.c[,i] <- seq(from=0, to=1, length=n)^(i-1)
    
    
  }
  
  if (is.null(thresh.val)) {
    
    if (is.null(sigma)) sigma <- mad(diff(x/sqrt(2)))
    if (is.null(alpha)) alpha <- 0.1
    
    if (thresh.type == "univ") base.thresh <- thresh_kab(n, alpha)
    #		if (thresh.type == "sim") base.thresh <- as.numeric(quantile(cov_dep_multi_norm(x.c, N), 1-alpha))
    if (thresh.type == "sim") base.thresh <- as.numeric(quantile(nearest_simulated_norms(n, deg), 1-alpha))
    
    thresh <- sigma * base.thresh
    
  }
  
  nsp(x, x.c, M, thresh, overlap)
  
}



################## NSP_POLY_AR
nsp_poly_ar <- function(x, ord = 1, M = 1000, thresh = NULL, sigma = NULL, alpha = 0.1, deg = 0, power = 1/2, min.size = 20, overlap = FALSE, buffer = ord) {
  n <- length(x)
  
  x.c <- matrix(x, n, deg+1+ord)
  
  for (i in 1:(deg+1)) {
    
    x.c[,i] <- seq(from=0, to=1, length=n)^(i-1)
    
    
  }
  
  if (ord) for (i in 1:ord) x.c[(1+i):n,deg+1+i] <- x[1:(n-i)]
  
  x.c <- x.c[(ord+1):n,]
  
  x <- x[(ord+1):n]
  
  if (is.null(thresh)) {
    
    if (is.null(sigma)) sigma <- est_var(x, x.c, power, min.size)
    if (is.null(alpha)) alpha <- 0.1
    
    thresh <- sigma * thresh_kab(n-ord, alpha)
  }
  
  res <- nsp(x, x.c, M, thresh, overlap, buffer)
  
  res$intervals[, c(1, 2, 4)] <- res$intervals[, c(1, 2, 4)] + ord
  
  res
  
}


########## nsp_reg ##########################

nsp_tvreg <- function(y, x, M = 1000, thresh = NULL, sigma = NULL, alpha = 0.1, power = 1/2, min.size = 20, overlap = FALSE) {
  n <- length(y)
  
  if (is.null(thresh)) {
    
    if (is.null(sigma)) sigma <- est_var(y, x, power, min.size)
    if (is.null(alpha)) alpha <- 0.1
    
    thresh <- sigma * thresh_kab(n, alpha)
  }
  
  nsp(y, x, M, thresh, overlap)
  
}

########## nsp_selfnorm #####################

nsp_selfnorm <- function(x, constraints, M, thresh, power = 1/2, minsize = 20, eps = 0.03, c = exp(1 + 2 * eps), overlap = FALSE) {
  d <- dim(constraints)
  
  x.c <- cbind(x, constraints)
  
  x.c.ads <- all_dyadic_scans_array(x.c)
  
  Vn2est <- est_var(x, constraints, power, minsize, TRUE)
  
  res <- ircs2sas(c(1, d[1]), x.c.ads, M, thresh, Vn2est, eps, c, overlap)
  
  intervals <- data.frame(t(order_chron(res)))
  colnames(intervals) <- c("starts", "ends", "values")
  intervals$midpoints <- floor((intervals$starts+intervals$ends)/2)
  
  
  list(intervals=intervals, threshold.used=thresh)
}


######### nsp_poly_selfnorm#############
nsp_poly_selfnorm <- function(x, M = 1000, thresh = NULL, power = 1/2, minsize = 20, alpha = 0.1, deg = 0, eps = 0.03, c = exp(1 + 2 * eps), overlap = FALSE) {
  n <- length(x)
  
  x.c <- matrix(x, n, deg+1)
  
  for (i in 1:(deg+1)) {
    
    x.c[,i] <- seq(from=0, to=1, length=n)^(i-1)
    
    
  }
  
  if (is.null(thresh)) {
    
    if (is.null(alpha)) alpha <- 0.1
    wh <- get(paste("wiener.holder_", as.character(eps), sep=""))
    thresh <- as.numeric(quantile(wh, 1-alpha))
    
  }
  
  nsp_selfnorm(x, x.c, M, thresh, power, minsize, eps, c, overlap)
  
}


############### sim_max_holder ###################
sim_max_holder <- function(n, N, eps, c = exp(1 + 2 * eps)) {
  
  # Simulate a sample of size N of values of the Holder-like norm of the Wiener process discretised with step 1/n.
  # See the "Example" in the description of the function "nsp_selfnorm".
  
  max.holder.sample <- rep(0, N)
  
  for (i in 1:N) {
    
    e <- rnorm(n)
    
    max.holder.sample[i] <- max_holder(e, eps, c)
    
  }
  
  max.holder.sample
  
}


################# draw_rects #################
draw_rects <- function(nsp.obj, yrange, density = 10, col = "red", x.axis.start = 1) {
  d <- dim(nsp.obj$intervals)
  if (d[1]) for (i in 1:d[1]) {
    
    rect(nsp.obj$intervals[i,1]+x.axis.start-1, yrange[1], nsp.obj$intervals[i,2]+x.axis.start-1, yrange[2], density=density, col=col)
    
    
  }
  
}

#################### draw_rects_advanced ###############

draw_rects_advanced <- function(x, nsp.obj, half.height = NULL, show.middles = TRUE, col.middles = "blue", lwd = 3, density = 10, col.rects = "red", x.axis.start = 1) {
  loc.est <- round((nsp.obj$intervals[,1] + nsp.obj$intervals[,2])/2)
  
  if (is.null(half.height)) half.height <- 2 * mad(diff(x)/sqrt(2))
  
  centres.y <- x[loc.est]
  
  d <- dim(nsp.obj$intervals)
  if (d[1]) for (i in 1:d[1]) {
    
    rect(nsp.obj$intervals[i,1]+x.axis.start-1, centres.y[i]-half.height, nsp.obj$intervals[i,2]+x.axis.start-1, centres.y[i]+half.height, density=density, col=col.rects)
    
    if (show.middles) lines(rep(loc.est[i]+x.axis.start-1, 2), c(centres.y[i]-half.height, centres.y[i]+half.height), col=col.middles, lwd = lwd)
    
  }
  
}

############### cpt_importance ###############
cpt_importance <- function(nsp.obj) {
  
  # Change-point prominence plot as described in Section 4 of the paper.
  # nsp.obj - quantity returned by one of the nsp_* functions.
  
  d <- dim(nsp.obj$intervals)
  if (d[1]) {
    heights <- nsp.obj$intervals[,2] - nsp.obj$intervals[,1]
    h.ord <- order(heights)
    labels <- paste(as.character(round(nsp.obj$intervals[h.ord,1])), "-", as.character(round(nsp.obj$intervals[h.ord,2])), sep = "")
    barplot(heights[h.ord], names.arg=labels)
  }
  else warning("No change-points to arrange in order of importance.")
  
}

############### order_chron ###############
order_chron <- function(nsp.obj) {
  
  # Order intervals of significance chronologically.
  # nsp.obj - quantity returned by one of the nsp_* functions.
  
  d <- dim(nsp.obj)
  if (d[2]) {
    nsp.obj.ord <- order(nsp.obj[1,])
    nsp.obj <- nsp.obj[,nsp.obj.ord]
  }	
  
  nsp.obj
  
}

############# select_narrowest ###############
select_narrowest <- function(nsp.obj, how.many = dim(nsp.obj)[2], order.chron = FALSE) {
  
  # Select how.many narrowest intervals of significance.
  # nsp.obj - quantity returned by one of the nsp_* functions.
  # order.chron - whether to order them chronologically (TRUE) or in increasing order of width (FALSE).
  
  d <- dim(nsp.obj)
  if (d[2]) {
    nsp.obj.ord <- order(nsp.obj[2,] - nsp.obj[1,])
    nsp.obj <- nsp.obj[,nsp.obj.ord]
    
    if (how.many) nsp.obj <- nsp.obj[,1:min(how.many, d[2])]
    
    if (order.chron) nsp.obj <- order_chron(nsp.obj)
    
  }	
  
  nsp.obj
  
}


###################3 INTERNAL FUNCTIONS ###########
#############cov_dep_multinorm #########
cov_dep_multi_norm <- function(constraints, N = 1000) {
  
  # Covariate-dependent multiresolution norm simulation
  
  d <- dim(constraints)
  
  res <- rep(0, N)
  
  for (i in 1:N) {
    
    x <- rnorm(d[1])
    
    x.c <- cbind(x, constraints)
    
    x.c.ads <- all_dyadic_scans_array(x.c)
    
    res[i] <- check_interval_array(c(1, d[1]), x.c.ads, 0)
    
  }
  
  res
  
}

####### cov_dep_multi_norm_poly ###################
cov_dep_multi_norm_poly <- function(n, deg, N = 10000) {
  
  # As above, but for use specifically in the polynomial models.
  
  x.c <- matrix(0, n, deg+1)
  
  for (i in 1:(deg+1)) {
    
    x.c[,i] <- seq(from=0, to=1, length=n)^(i-1)
    
  }
  
  cov_dep_multi_norm(x.c, N)	
  
}

######## check_specific_vector_poly ##############
check_specific_vector_poly <- function(x, deg) {
  
  n <- length(x)
  
  x.c <- matrix(0, n, deg+1)
  
  for (i in 1:(deg+1)) {
    
    x.c[,i] <- seq(from=0, to=1, length=n)^(i-1)
    
  }
  
  x.c <- cbind(x, x.c)
  
  x.c.ads <- all_dyadic_scans_array(x.c)
  
  check_interval_array(c(1, n), x.c.ads, 0)
  
}


######## nearest_simulated_norms ################3
nearest_simulated_norms <- function(n, deg) {
  
  # This extracts the most suitable variable mres.norm.* to be used
  
  if (n > 2150) stop("The largest permitted n is 2150.")
  
  len <- as.character(round(round(min(max(1, n/100), 21))*100))
  
  if (deg == 0) dg <- "c" else if (deg == 1) dg <- "l"
  
  var.name <- paste("mres.norm.", dg, ".", len, sep="")
  
  get(var.name)
  
}

########## all_dyadic_scan_array ###############33
all_dyadic_scans_array <- function(x) {
  
  d <- dim(x)
  n <- d[1]
  
  if (n) {
    
    add.scales <- floor(logb(n, 2))
    shifts <- rep(0, add.scales+1)
    res <- array(x, c(d[1], d[2], add.scales+1))
    if (add.scales) for (j in 1:add.scales) {
      res[1:(n-2^j+1),,(j+1)] <- 2^(-1/2) * (res[1:(n-2^j+1),,j] + res[(2^(j-1)+1):(n-2^j+1+2^(j-1)),,j])
      res[(n-2^j+2):(n),,(j+1)] <- 0
      shifts[j+1] <- 2^j-1
    }		
    
  }
  else {
    res <- array(0, c(d[1], d[2], 0))	
    shifts <- integer(0)
  }
  
  list(res=res, shifts=shifts)
  
}

####### iter_random_checks_scan_array ################
iter_random_checks_scan_array <- function(ind, ads.array, M, thresh, overlap = FALSE, buffer = 0) {
  
  s <- ind[1]
  e <- ind[2]
  n <- e - s + 1
  
  indices <- ((ads.array$shifts+1) <= (n/2))
  
  ads.array$res <- ads.array$res[s:e,,indices,drop=F]
  
  ads.array$shifts <- ads.array$shifts[indices]
  
  if (n > 1) {
    
    next.int <- random_checks_scan_2stage_array(c(1,n), ads.array, M, thresh)
    
    if (!is.na(next.int$selected.ind))  {
      
      if (!overlap) {
        
        if (next.int$selected.val[1,1]-buffer >= 2) left <- iter_random_checks_scan_array(c(1, next.int$selected.val[1,1]-buffer), ads.array, M, thresh, overlap, buffer) else left <- matrix(NA, 3, 0)
        if (n - next.int$selected.val[1,2]-buffer >= 1) {
          
          right <- iter_random_checks_scan_array(c(next.int$selected.val[1,2]+buffer, n), ads.array, M, thresh, overlap, buffer)
          if (dim(right)[2]) right <- right + c(rep(next.int$selected.val[1,2]-1+buffer, 2), 0)
          
        }
        else right <- matrix(NA, 3, 0)
      }
      
      else {
        
        
        if (floor(mean(next.int$selected.val[1,1:2]))-buffer >= 2) left <- iter_random_checks_scan_array(c(1, floor(mean(next.int$selected.val[1,1:2]))-buffer), ads.array, M, thresh, overlap, buffer) else left <- matrix(NA, 3, 0)
        if (n - floor(mean(next.int$selected.val[1,1:2]))-buffer >= 2) {
          right <- iter_random_checks_scan_array(c(floor(mean(next.int$selected.val[1,1:2]))+1+buffer, n), ads.array, M, thresh, overlap, buffer)
          if (dim(right)[2]) right <- right + c(rep(floor(mean(next.int$selected.val[1,1:2]))+buffer, 2), 0)
          
        }
        else right <- matrix(NA, 3, 0)
        
        
      }
      
      
      return(cbind(t(next.int$selected.val), left, right))
      
      
    }
    
    else(return(matrix(NA, 3, 0)))
    
    
  }
  
  else(return(matrix(NA, 3, 0)))
  
}



random_checks_scan_2stage_array <- function(ind, ads.array, M, thresh) {
  
  s1 <- random_checks_scan_array_1by1(ind, ads.array, M, thresh)
  
  if (!is.na(s1$selected.ind)) {
    
    s <- s1$selected.val[1,1] + ind[1] - 1
    e <- s1$selected.val[1,2] + ind[1] - 1
    
    s2 <- random_checks_scan_array_1by1(c(s,e), ads.array, M, thresh)
    
    if (!is.na(s2$selected.ind)) {
      
      replacement <- s2$selected.val
      replacement[1,1:2] <- replacement[1,1:2] + s - ind[1]
      s1$selected.val <- replacement
      
    }
    
  }
  
  s1	
  
}


######## random_checks_scan_array_1by1 #############
random_checks_scan_array_1by1 <- function(ind, ads.array, M, thresh) {
  
  s <- ind[1]
  e <- ind[2]
  n <- e - s + 1
  
  if (n > 1) {
    
    indices <- ((ads.array$shifts+1) <= (n/2))
    
    ads.array$res <- ads.array$res[s:e,,indices,drop=F]
    
    ads.array$shifts <- ads.array$shifts[indices]
    
    M <- min(M, (n-1)*n/2)
    
    ind <- grid_intervals_sorted(n, M)
    
    M <- dim(ind)[2]
    
    res <- matrix(0, M, 3)
    
    res[,1:2] <- t(ind)
    
    zero.check <- TRUE
    j <- 1
    
    while (zero.check && (j <= M)) {
      
      res[j,3] <- check_interval_array(res[j,1:2], ads.array, thresh)
      zero.check <- (res[j,3] == 0)
      j <- j + 1
      
    }
    
    if (zero.check) {
      
      selected.ind <- NA
      selected.val <- matrix(0, 0, 3)
      
    }
    
    else {
      
      selected.ind <- j-1
      selected.val <- res[selected.ind,,drop=FALSE]
      
    }
    
    
  }
  
  else {
    
    selected.val <- matrix(0, 0, 3)
    selected.ind <- NA
    M <- 0
    
  }
  
  
  list(selected.ind = selected.ind, selected.val = selected.val, M.eff=max(1, M))
  
}

########## all_intervals_flat ###############
all_intervals_flat <- function(n) {
  
  if (n == 2) ind <- matrix(1:2, 2, 1) else {
    M <- (n-1)*n/2	
    ind <- matrix(0, 2, M)
    ind[1,] <- rep(1:(n-1), (n-1):1)
    ind[2,] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)
  }
  ind
  
}


######## all_intervals_sorted ############3
all_intervals_sorted <- function(n) {
  
  d <- all_intervals_flat(n)
  d.ord <- order(d[2,] - d[1,])
  d[,d.ord, drop=FALSE]
  
}

########## grid_intervals_sorted ############
grid_intervals_sorted <- function(n, M) {
  
  if ((n==2) || (M == 0)) ind <- matrix(c(1, n), 2, 1)
  
  else if (M >= (n-1)*n/2) ind <- all_intervals_sorted(n)
  
  else {
    k <- 1
    while (k*(k-1)/2 < M) k <- k+1
    ind2 <- all_intervals_sorted(k)
    ind2.mx <- max(ind2)
    ind <- round((ind2 - 1) * ((n-1) / (ind2.mx-1)) + 1)
  }	
  
  ind	
}

########### check_interval_array ##############
check_interval_array <- function(ind, ads.array, thresh) {
  
  s <- ind[1]
  e <- ind[2]
  n <- e - s + 1
  
  indices <- ((ads.array$shifts+1) <= (n/2))
  
  ads.array$res <- ads.array$res[s:e,,indices,drop=F]
  
  ads.array$shifts <- ads.array$shifts[indices]
  
  dm <- dim(ads.array$res)
  
  f.con.rhs.core <- matrix(0, 0, 1+2*(dm[2]-1))
  
  scales <- length(ads.array$shifts)	
  
  for (i in 1:scales) {
    
    f.con.rhs.current <- ads.array$res[1:(n-ads.array$shifts[i]),,i]
    f.con.rhs.current <- cbind(f.con.rhs.current, -f.con.rhs.current[,-1])
    
    f.con.rhs.core <- rbind(f.con.rhs.core, f.con.rhs.current)		
    
  }
  
  f.con.rhs.core <- rbind(f.con.rhs.core, -f.con.rhs.core)
  
  f.obj <- c(1, rep(0, 2*(dm[2]-1)))
  f.rhs <- f.con.rhs.core[,1]
  f.con <- f.con.rhs.core
  f.con[,1] <- 1	
  d <- dim(f.con.rhs.core)
  
  f.dir <- rep(">=", d[1])
  linf <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution[1]
  linf.t <- linf * (linf > thresh)
  linf.t
  
}

####### lp_selfnorm #############
lp_selfnorm <- function(ind, ads.array, selfnorm.array) {
  
  s <- ind[1]
  e <- ind[2]
  n <- e - s + 1
  
  indices <- ((ads.array$shifts+1) <= (n/2))
  
  ads.array$res <- ads.array$res[s:e,,indices,drop=F]
  
  ads.array$shifts <- ads.array$shifts[indices]
  
  dm <- dim(ads.array$res)
  
  f.con.rhs.core <- matrix(0, 0, 1+2*(dm[2]-1))
  
  scales <- length(ads.array$shifts)	
  
  for (i in 1:scales) {
    
    f.con.rhs.current <- ads.array$res[1:(n-ads.array$shifts[i]),,i]  / selfnorm.array[1:(n-ads.array$shifts[i]),i]
    f.con.rhs.current[!is.finite(f.con.rhs.current)] <- 0
    f.con.rhs.current <- cbind(f.con.rhs.current, -f.con.rhs.current[,-1])
    
    f.con.rhs.core <- rbind(f.con.rhs.core, f.con.rhs.current)		
    
  }
  
  f.con.rhs.core <- rbind(f.con.rhs.core, -f.con.rhs.core)
  
  f.obj <- c(1, rep(0, 2*(dm[2]-1)))
  f.rhs <- f.con.rhs.core[,1]
  f.con <- f.con.rhs.core
  f.con[,1] <- 1	
  d <- dim(f.con.rhs.core)
  
  f.dir <- rep(">=", d[1])
  lp.sol <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution
  lp.sol
  
}



############## create_selfnorm_array ################
create_selfnorm_array_Vn2est <- function(resid, Vn2est, eps, c = exp(1 + 2 * eps)) {
  
  m <- length(resid)
  
  zz <- all_dyadic_scans_array(matrix(resid^2, m, 1))
  zz$res <- zz$res[,,]
  
  zz.norm <- all_dyadic_scans_array(matrix(1, m, 1))
  zz.norm$res <- zz.norm$res[,,]
  
  (1 + eps) * sqrt(zz$res / zz.norm$res) * log(c * pmax(1, Vn2est / (zz$res * zz.norm$res)))^(1/2 + eps)
  
}

###### lin_reg_resid ##############
linreg_resid <- function(ind, ads.array) {
  
  s <- ind[1]
  e <- ind[2]
  
  lmmat <- ads.array$res[s:e,,1]
  
  res <- as.numeric(lm(lmmat[,1] ~ lmmat[,-1]-1)$resid)
  
  if (sum(res^2) == 0) res <- (lmmat[,1] - mean(lmmat[,1]))
  
  if (sum(res^2) == 0) res <- lmmat[,1]
  
  res		
  
}

####### check_interval_array_selfnorm ###########
check_interval_array_selfnorm <- function(ind, ads.array, thresh, Vn2est, eps = 0.03, c = exp(1 + 2 * eps)) {
  
  resid <- linreg_resid(ind, ads.array)
  resid.sna <- create_selfnorm_array_Vn2est(resid, Vn2est, eps, c)
  a <- lp_selfnorm(ind, ads.array, resid.sna)
  sol <- a[1] * (a[1] > thresh)
  sol	
  
}

####### random_check_scan_array_selfnorm #############
random_checks_scan_array_selfnorm_1by1 <- function(ind, ads.array, M, thresh, Vn2est, eps = 0.03, c = exp(1 + 2 * eps)) {
  
  s <- ind[1]
  e <- ind[2]
  n <- e - s + 1
  
  if (n > 1) {
    
    indices <- ((ads.array$shifts+1) <= (n/2))
    
    ads.array$res <- ads.array$res[s:e,,indices,drop=F]
    
    ads.array$shifts <- ads.array$shifts[indices]
    
    M <- min(M, (n-1)*n/2)
    
    ind <- grid_intervals_sorted(n, M)
    
    M <- dim(ind)[2]
    
    res <- matrix(0, M, 3)
    
    res[,1:2] <- t(ind)
    
    zero.check <- TRUE
    j <- 1
    
    while (zero.check && (j <= M)) {
      
      res[j,3] <- check_interval_array_selfnorm(res[j,1:2], ads.array, thresh, Vn2est, eps, c)
      zero.check <- (res[j,3] == 0)
      j <- j + 1
      
    }
    
    if (zero.check) {
      
      selected.ind <- NA
      selected.val <- matrix(0, 0, 3)
      
    }
    
    else {
      
      selected.ind <- j-1
      selected.val <- res[selected.ind,,drop=FALSE]
      
    }
    
    
  }
  
  else {
    
    filtered.res <- selected.val <- matrix(0, 0, 3)
    selected.ind <- NA
    M <- 0
    
  }
  
  list(selected.ind = selected.ind, selected.val = selected.val, M.eff=max(1, M))
  
}

##### random_checks_scan_2stage_array_selfnorm ##########
random_checks_scan_2stage_array_selfnorm <- rcs2sas <- function(ind, ads.array, M, thresh, Vn2est, eps = 0.03, c = exp(1 + 2 * eps)) {
  
  s1 <- random_checks_scan_array_selfnorm_1by1(ind, ads.array, M, thresh, Vn2est, eps, c)
  
  if (!is.na(s1$selected.ind)) {
    
    s <- s1$selected.val[1,1] + ind[1] - 1
    e <- s1$selected.val[1,2] + ind[1] - 1
    
    s2 <- random_checks_scan_array_selfnorm_1by1(c(s,e), ads.array, M, thresh, Vn2est, eps, c)
    
    if (!is.na(s2$selected.ind)) {
      
      replacement <- s2$selected.val
      replacement[1,1:2] <- replacement[1,1:2] + s - ind[1]
      s1$selected.val <- replacement
      
    }
    
  }
  
  s1	
  
}

####### ircs2s ####################
ircs2sas <- function(ind, ads.array, M, thresh, Vn2est, eps = 0.03, c = exp(1 + 2 * eps), overlap = FALSE) {
  
  s <- ind[1]
  e <- ind[2]
  n <- e - s + 1
  
  indices <- ((ads.array$shifts+1) <= (n/2))
  
  ads.array$res <- ads.array$res[s:e,,indices,drop=F]
  
  ads.array$shifts <- ads.array$shifts[indices]
  
  if (n > 1) {
    
    next.int <- rcs2sas(c(1,n), ads.array, M, thresh, Vn2est, eps, c)
    
    if (!is.na(next.int$selected.ind))  {
      
      if (!overlap) {
        
        if (next.int$selected.val[1,1] >= 2) left <- ircs2sas(c(1, next.int$selected.val[1,1]), ads.array, M, thresh, Vn2est, eps, c, overlap) else left <- matrix(NA, 3, 0)
        if (n - next.int$selected.val[1,2] >= 1) {
          
          right <- ircs2sas(c(next.int$selected.val[1,2], n), ads.array, M, thresh, Vn2est, eps, c, overlap)
          if (dim(right)[2]) right <- right + c(rep(next.int$selected.val[1,2]-1, 2), 0)
          
        }
        else right <- matrix(NA, 3, 0)
      }
      
      else {
        
        
        if (floor(mean(next.int$selected.val[1,1:2])) >= 2) left <- ircs2sas(c(1, floor(mean(next.int$selected.val[1,1:2]))), ads.array, M, thresh, Vn2est, eps, c, overlap) else left <- matrix(NA, 3, 0)
        if (n - floor(mean(next.int$selected.val[1,1:2])) >= 2) {
          right <- ircs2sas(c(floor(mean(next.int$selected.val[1,1:2]))+1, n), ads.array, M, thresh, Vn2est, eps, c, overlap)
          if (dim(right)[2]) right <- right + c(rep(floor(mean(next.int$selected.val[1,1:2])), 2), 0)
          
        }
        else right <- matrix(NA, 3, 0)
        
        
      }
      
      
      return(cbind(t(next.int$selected.val), left, right))
      
    }
    
    else(return(matrix(NA, 3, 0)))
    
  }
  
  else(return(matrix(NA, 3, 0)))
  
}

############ max_holder ###########3
max_holder <- function(e, eps, c = exp(1 + 2 * eps)) {
  
  n <- length(e)
  
  eps.cum <- c(0, cumsum(e))
  
  max.stat <- 0
  
  for (i in 0:(n-1)) for (j in (i+1):n) {
    
    scan.stat <- abs(eps.cum[j+1] - eps.cum[i+1]) / sqrt(j-i) / log(c * n / (j-i))^(1/2+eps)
    
    if (scan.stat > max.stat) max.stat <- scan.stat
    
  }
  
  max.stat
  
}


######### est_var ###################
est_var <- function(y, x, power = 1/2, min.size = 20, estVn2 = FALSE) {
  
  n <- length(y)
  w.size <- min(n, max(round(n^power), min.size))
  
  how.many <- n - w.size + 1
  
  res <- rep(0, how.many)
  
  for (i in 1:how.many) {
    
    resp <- y[i:(i+w.size-1)]
    covs <- x[i:(i+w.size-1),]
    
    res[i] <- summary(lm(resp ~ covs))$sigma
    
  }	
  
  if (estVn2) est <- n / (n - w.size + 1) * sum(res^2)
  else est <- median(res)
  
  est
  
}


####### thresh_kab #############
thresh_kab <- function(n, alpha = 0.1, method = "asymp") {
  
  an <- sqrt(2 * log(n)) + (1/2*log(log(n)) + log(0.8197466 / 2 /sqrt(pi))) / sqrt(2 * log(n))
  
  bn <- 1 / sqrt(2 * log(n))
  
  if (method == "bound") beta <- alpha/2
  else if (method == "asymp") beta <- 1 - sqrt(1-alpha)
  
  an + bn * log(1/log(1/(1 - beta)))
  
}





########################################
########### simulation & example ########
########################################

#### smuce_coverage ##########
smuce_coverage <- function(truth, est) {
  
  # Returns "TRUE" if the SMUCE estimator "est" is such that each confidence interval covers a true change-point of "truth", and "FALSE" otherwise. See "sim_coverage" for an example of use.
  
  res <- TRUE
  
  j <- jumpint(est)
  d <- dim(j)
  
  cpts <- which(abs(diff(truth)) > 0)
  
  if (d[1]-1) for (i in 1:(d[1]-1)) {
    
    cur_int <- j[i,3]:j[i,4]
    res <- res && any(cpts %in% cur_int)
  }
  
  res	
  
}

nsp_coverage <- function(truth, est) {
  
  # Returns "TRUE" if the NSP estimator "est" is such that each interval of significance covers a true change-point of "truth", and "FALSE" otherwise. See "sim_coverage" for an example of use.
  
  res <- TRUE
  
  d <- dim(est$intervals)
  
  cpts <- which(abs(diff(truth)) > 0)
  
  if (d[1]) for (i in 1:(d[1])) {
    
    cur_int <- est$intervals[i,1]:(est$intervals[i,2]-1)
    res <- res && any(cpts %in% cur_int)
  }
  
  res	
  
}

######### nsp_coverage_extended 333333333333
nsp_coverage_extended <- function(truth, est) {
  
  # Returns as in "nsp_coverage", plus the number of intervals of significance in "est".	
  
  res <- TRUE
  
  d <- dim(est$intervals)
  
  cpts <- which(abs(diff(truth)) > 0)
  
  if (d[1]) for (i in 1:(d[1])) {
    
    cur_int <- est$intervals[i,1]:(est$intervals[i,2]-1)
    res <- res && any(cpts %in% cur_int)
  }
  
  list(coverage = res, how.many = d[1])
  
}

########### sim_coverage ############
sim_coverage <- function(truth, sigma, M = 100, N = 100) {
  
  # Simulates N realisations of truth + sigma * rnorm(length(truth)) and checks SMUCE and NSP coverage for each sample path, via "nsp_coverage" and "smuce_coverage".
  
  res.nsp <- res.smuce <- rep(0, N)
  
  for (i in 1:N) {
    print(i)
    x <- truth + sigma * rnorm(length(truth))
    
    x.n <- nsp_poly(x, M = M)
    print(x.n)
    res.nsp[i] <- nsp_coverage(truth, x.n)
    
    x.s <- stepFit(x, alpha = 0.1, confband = T)
    print(jumpint(x.s))
    res.smuce[i] <- smuce_coverage(truth, x.s)
    
    print(paste("NSP:", res.nsp[i]))
    print(paste("SMUCE:", res.smuce[i]))
    
  }
  
  list(res.nsp=res.nsp, res.smuce=res.smuce)
  
}

sim_coverage_extended <- function(truth, ar.coef, div.factor, M = 100, N = 100) {
  
  # As in "sim_coverage", but only for NSP and also returns the numbers of intervals of significance in each sample path.
  
  coverage <- how.many <- rep(0, N)
  
  for (i in 1:N) {
    print(i)
    x <- truth + arima.sim(list(ar=ar.coef), n=length(truth)) / div.factor
    
    x.n <- nsp_poly_ar(x, 1, M = M)
    cvg <- nsp_coverage_extended(truth, x.n)
    
    coverage[i] <- cvg$coverage
    how.many[i] <- cvg$how.many
    
  }
  
  list(coverage = coverage, how.many = how.many)
  
}







#################
#################
#################
prost1<-read.csv("G:/MASTERS/THESIS/bimj2320-sup-0001-suppmat/Revision_code/Prostate SEER/prostate_example_data_covs.csv")
library(tidyverse)
pr<-prost1%>%
  arrange(time)
  
table(pr$delta)
library(survival)
library(ggplot2)
s=survfit(Surv(time,delta)~1,data=prost1)
s2=summary(s,censored = FALSE)
table(s1$time,s1$n.event)
Dt<-s2$n.event
ts.plot(s2$n.event)
At <- sqrt(2 * Dt + 3/8)
ts.plot(At)
nsp_poly(At,deg=1) -> At.n


draw_rects(At.n, c(0, 300))
draw_rects_advanced(At,At.n)
# 1st cp
d1=nsp_poly(At[1:67],deg=1)
rr<-(At[1:67])
mn.rr<-mean(rr)
rr1<-as.list(rr)
rr1[[1]]
S <- list()
S[[1]] <- 0
for(i in 1:67){
  S[[i+1]]<-S[[i]] + (rr1[[i]]-mn.rr)  
}
rr[[34]]
s1=S %>% unlist()
max.S<-max(s1)
S.diff<-max.S-min.S
c=which(s1==max.S)

## which corresponds to time=34 has our first change point.
#
########## 

# 2nd cp

At_35<-At[35:67]
d2=nsp_poly(At_34)

rr_2<-At_35[4:13]
mn.rr2<-mean(rr_2)
rr_2l<-as.list(rr_2)
rr1[[1]]
S_2 <- list()
S_2[[1]] <- 0
for(i in 1:10){
  S_2[[i+1]]<-S_2[[i]] + (rr_2l[[i]]-mn.rr2)  
}
rr_2[[7]]
s_2=S_2 %>% unlist()
max.S2<-max(s_2)
S.diff<-max.S-min.S
c=which(s_2==max.S2)
rr_2[[6]]
### so second change point has been found at time point 44.


# 3rd cp
At_72<-At[72:150]
d3=nsp_poly(At_72)

rr_3<-At_72[2:20]
mn.rr3<-mean(rr_3)
rr_3l<-as.list(rr_3)
rr1[[1]]
S_3 <- list()
S_3[[1]] <- 0
for(i in 1:19){
  S_3[[i+1]]<-S_3[[i]] + (rr_3l[[i]]-mn.rr3)  
}
rr_2[[7]]
s_3=S_3 %>% unlist()
max.S3<-max(s_3)
S.diff<-max.S-min.S
c=which(s_3==max.S3)
rr_3[[8]]
### so 3rd  change point has been found at time point 80..


# 4th cp
At_81<-At[81:130]

rr_4<-At_81
mn.rr4<-mean(rr_4)
rr_4l<-as.list(rr_4)
rr1[[1]]
S_4 <- list()
S_4[[1]] <- 0
for(i in 1:50){
  S_4[[i+1]]<-S_4[[i]] + (rr_4l[[i]]-mn.rr4)  
}
rr_2[[7]]
s_4=S_4 %>% unlist()
max.S4<-max(s_4)
S.diff<-max.S-min.S
c=which(s_4==max.S4)
rr_4[[26]]
### so 4th  change point has been found at time point 106.

pr <- At[1:34]
time_1alt <- 1:34
sig_1alt <- lm(pr ~ time_1alt)$fitted

pr2 <- c(rep(NA,34),At[35:44])
tim <- c(rep(NA,34),35:44) 
tim2 <- tim^2 / 44
dd <- lm(pr2  ~ tim + tim2)$fitted
dd1<-c(rep(NA,34),dd)

pr3 <- At[45:80]
tim <- 45:80 
tim2 <- tim^2 / 36
dd2 <- lm(pr3  ~ tim + tim2)$fitted
dd2<-c(rep(NA,44),dd2)

pr4 <- At[81:106]
tim <- 81:106
tim2 <- tim^2 / 26
dd3 <- lm(pr4  ~ tim + tim2)$fitted
dd3<-c(rep(NA,80),dd3)

At[107:250]
AT1<-na.omit(At[107:250])
mean(AT1)
length(AT1)
dd4 <- rep(mean(AT1), 110)
dd4<-c(rep(NA,106),dd4)

pr5 <- At[107:250]
tim <- 107:250
tim2 <- tim^2 / 144
dd5 <- lm(pr5  ~ tim + tim2)$fitted
dd5<-c(rep(NA,106),dd5)

ts.plot(At, xlab="Time (quarters)", ylab="")
lines(sig_1alt, col="red", lwd = 3)
lines(dd1, col="blue", lwd = 3,xlim=c(35,44))
lines(dd2, col="skyblue", lwd = 3,xlim=c(35,44))
lines(dd3, col="yellow", lwd = 3,xlim=c(35,44))
lines(dd4, col="grey", lwd = 3,xlim=c(35,44))
lines(dd5, col="blue", lwd = 3,xlim=c(35,44))


ts.plot(At, xlab="Time (quarters)", ylab="")
lines(sig_1alt, col="red", lwd = 3)
lines(dd1, col="pink", lwd = 3,xlim=c(35,44))



sig_2alt <- 
sig_2alt <- rep(mean(real_dat_sc[77:103]), 103-77+1)
sig_linalt <- c(sig_1alt, sig_2alt)

#########

################# parameter estimation #####

# paraeters for 1st cp

library(SurvRegCensCov)
library(fastDummies)
race_dummy <- dummy_cols(prost1$race)
Z <- matrix(c(prost1$age, race_dummy$.data_Other, race_dummy$.data_White), ncol = 3)
Z1<-prost1%>%matrix(c(prost1$age, race_dummy$.data_Other, race_dummy$.data_White), ncol = 3)

  
pr<-prost1%>%
  mutate(Z1=Z[,1],
         Z2=Z[,2],
         Z3=Z[,3])%>%
  arrange(time)



table(cp3_pr$race)

WR3=WeibullReg(Surv(time, delta) ~ Z1+Z2+Z3, data=cp3_pr)
WR3$coefcp1_pr<-pr%>%
  filter(time<35)
table(cp1_pr$time)
(cp1_pr$Z1[1:20])
Z1[1:20]
WR1=WeibullReg(Surv(time, delta) ~ Z1[1:7786]+Z2[1:7786]+Z3[1:7786], data=cp1_pr)
WR1$coef

#parameters for 2nd cp

cp2_pr<-pr%>%
  filter(time>34)%>%
  filter(time<45)
table(cp2_pr$time)
length(cp3_pr$Z1)
WR2=WeibullReg(Surv(time, delta) ~ Z1+Z2+Z3, data=cp2_pr)
WR2$coef
# parameters for 3rd cp
cp3_pr<-pr%>%
  filter(time>44)%>%
  filter(time<81)
SR$coefficients
# parameters for 4th cp
cp4_pr<-pr%>%
  filter(time>80)%>%
  filter(time<107)
table(cp4_pr$time)

WR4=WeibullReg(Surv(time, delta) ~ Z1+Z2+Z3, data=cp4_pr)
WR4$coef

#PARAMETER FOR 5TH CP
cp5_pr<-pr%>%
  filter(time>106)%>%
  filter(time<251)
table(cp4_pr$time)

WR5=WeibullReg(Surv(time, delta) ~ Z1+Z2+Z3, data=cp5_pr)
WR5$coef
