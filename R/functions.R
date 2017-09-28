#' makeData
#' 
#' Make the data used in the simulation study
#' 
#' @param n Sample size
#' @param t0 The final timepoint
#' @param ... Other args, not currently used
#' @export

makeData <- function(n, t0=6, ...){
	W <- rbinom(n, 4, 0.5) + 1
	Z <- rbinom(n, 1, 0.5)

	# simulate data
	C <- rgeom(n,plogis(-3 + 0.2*as.numeric(W==2) - 0.2*as.numeric(W == 3))) + 1
	T <- rgeom(n,plogis(-2 + 0.4*as.numeric(W %in% c(1:3)) - 
	                    	0.2*as.numeric(W == 4) - Z)) + 1
	J <- rbinom(n,4,plogis(0.2*Z)) + 1
	# if everyone is censored at study end, fix up those times
	if(!is.na(t0)){
		C[C>t0] <- t0
	}

	# observed data time and event indicator
	ftime <- pmin(C,T)
	ftype <- J*as.numeric(T <= C)

	# return a data frame,
	rslt <- list(trt = Z, adjustVars = data.frame(site = W), 
	             ftime = ftime, ftype = ftype)
	return(rslt)
}

#' Get the true value of the trend parameter for a given covariance 
#' weight matrix
#' 
#' @param n Sample size (should be large to ensure numerical accuracy)
#' @param t0 Time to compute incidence
#' @param Upsilon The covariance weighting matrix
#' @param returnL Return the estimated log cumulative incidence ratios
#' 
#' @return If \code{returnL = FALSE}, a vector of trend parameter estimates,
#' the first corresponding to the "intercept" parameter, the second to the "slope"
#' or trend parameter. If \code{returnL = TRUE}, then a list with named entries
#' \code{param} (the projection parameters) and \code{L}, a vector of true values 
#' of L_{j,0}. 
#' @export

getTruth <- function(n = 1e6, t0 = 6, Upsilon = diag(5),
                     returnL = FALSE){
	W <- rbinom(n, 4, 0.5) + 1
	Z <- rbinom(n, 1, 0.5)

	# simulate data
	T <- rgeom(n,plogis(-2 + 0.4*as.numeric(W %in% c(1:3)) - 
	                    	0.2*as.numeric(W == 4) - Z)) + 1
	J <- rbinom(n,4,plogis(0.2*Z)) + 1

	# true value of each incidence
	r_vec <- fj1 <- fj0 <- rep(NA, 5)
	for(j in 1:5){
		fj1[j] <- sum(T <= t0 & Z == 1 & J == j) / sum(Z == 1)
		fj0[j] <- sum(T <= t0 & Z == 0 & J == j) / sum(Z == 0)
		r_vec[j] <- log(fj0[j] / fj1[j])
	}
	# project incidence
	j_vec <- (1:5) - 1

	# design matrix
	X <- cbind(rep(1,5), j_vec)
	R <- matrix(r_vec, ncol = 1)

	# parameters
	Upsilon_inv <- solve(Upsilon)
	psi0 <- solve(t(X)%*%Upsilon_inv%*%X)%*%t(X)%*%Upsilon_inv%*%R
	alpha <- psi0[1]; beta <- psi0[2]
	# plot on VE scale
	# plot(y = 1 - 1/exp(r_vec), x = j_vec, type = "b", ylim = c(0,1))
	# curve(1 - 1/exp(alpha - beta * x), 0, 4, col = 2, add = TRUE)
	if(!returnL){
		return(c(alpha, beta))		
	}else{
		return(list(param = c(alpha, beta),
		            L = R))
	}
}

#' Make a formula that will implement an empirical hazard
#' estimator when called in estimateCensoring within \code{survtmle}.
#' See \code{?survtmle} for more details
#' 
#' @param trt A vector of vaccine assignments
#' @param site A vector of geographical site assignments
#' @param ftime A vector of failure times
#' @param ftype A vector of failure types
#' @return A character formula to be input to \code{glm.ctime} in 
#' calls to \code{survtmle}
#' @export
get.ctimeForm <- function(trt, site, ftime, ftype){
	form <- "-1"
	for(i in unique(trt)){
		for(s in unique(site)){
			form <- c(form, 
			  paste0("I(trt==",i,"& site == ",s," & t==",
			         unique(ftime[ftype==0 & trt==i & site == s]),")",
			         collapse="+"))
		}
	}
	return(paste(form,collapse="+"))
}


#' Helper function to got log ratio of cumulative incidences
#' @param F A vector of cumulative incidence estimates with alternating
#' treatment = 0 and treatment = 1 for each type. 
g <- function(F){
  matrix(log(F[seq(1, length(F), 2)]/F[seq(2, length(F), 2)]), ncol = 1)
}


#' Helper function to compute the gradient for the log ratio
#' transformation of the cumulative incidence estimates.
#' @param F A vector of cumulative incidence estimates with alternating
#' treatment = 0 and treatment = 1 for each type.
grad_g <- function(F){
  K <- length(F) / 2
  tmp_vec <- rep(0, 2*K^2)
  ind <- as.numeric(paste0(sort(rep((1:K) - 1, 2)), 1:length(F)))
  ind[length(ind)] <- 2*K^2
  tmp_vec[ind] <- (-1)^(2:(length(F)+1)) * (1 / F)
  tmp <- matrix(tmp_vec, nrow = K, byrow = TRUE)
  tmp
}

#' Helper function to compute estimate of trend parameter (also returns
#' estimate of intercept). 
#' @param F A vector of cumulative incidence estimates with alternating
#' treatment = 0 and treatment = 1 for each type.
#' @param D A matrix of variance estimates corresponding with the vector 
#' \code{F}. 
#' @param j_vec Vector of genetic distances

h <- function(F, D, j_vec){
  R_n <- g(F)
  K <- length(F)/2
  nabla_g <- grad_g(F)
  Upsilon_n <- nabla_g %*% cov(D) %*% t(nabla_g)
  # design matrix
  X <- cbind(rep(1,length(j_vec)), j_vec)
  Upsilon_inv <- solve(Upsilon_n)
  S_n <- solve(t(X)%*%Upsilon_inv%*%X)%*%t(X)%*%Upsilon_inv
  S_n %*% R_n
}

#' Helper function to compute gradient of trend parameter function. 
#' @param F A vector of cumulative incidence estimates with alternating
#' treatment = 0 and treatment = 1 for each type.
#' @param D A matrix of variance estimates corresponding with the vector 
#' \code{F}.  

grad_h <- function(F, D, j_vec){
  R_n <- g(F)
  K <- length(F)/2
  nabla_g <- grad_g(F)
  Upsilon_n <- nabla_g %*% cov(D) %*% t(nabla_g) 
  # design matrix
  X <- cbind(rep(1,nlength(j_vec)), j_vec)
  Upsilon_inv <- solve(Upsilon_n)
  S_n <- solve(t(X)%*%Upsilon_inv%*%X)%*%t(X)%*%Upsilon_inv
  tmp <- (-1)^(2:(length(F)+1)) * (1 / F)
  tmp2 <- S_n[2, sort(rep(1:K, 2))]
  matrix(tmp*tmp2, ncol = 1)
}    

#' Test for a trend in sieve effect across different levels 
#' of an ordinal failure type. 
#' 
#' @param object An object of class \code{survtmle}. 
#' @param level The nominal coverage probability of the confidence interval.
#' @return An object of class \code{"trend_test"}.
#' \describe{
#'  \item{\code{alpha}}{The intercept from the projection onto the working model.}
#' 	\item{\code{beta}}{The slope from the projection onto the working model, i.e., 
#' 		  the "trend" parameter.}
#' 	\item{\code{L_j}}{A data.frame showing the log ratio of cumulative incidences 
#' 		  that were projected onto the working model.}
#' 	\item{\code{se_beta}}{The influence function-based estimate of the standard 
#' 		  error of \code{beta_n}.}
#' 	\item{\code{ci}}{The confidence interval at the requested level for \code{beta_n}.}
#' 	\item{\code{pval}}{The two-sided p-value from the test of the null hypothesis that
#' 		  \code{beta_n} = 0.}
#' 	\item{\code{level}}{The requested level of the confidence interval.}
#' }
#' @export

trend_test <- function(object, level = 0.95){
	stopifnot(class(object) == "survtmle" | class(object) == "getMO")
	n <- length(object$ic[[1]])
	F <- object$est
	D <- Reduce(cbind, object$ic)
	j_vec <- unique(as.numeric(unlist(lapply(strsplit(row.names(F), " "),"[[", 2))))

	# log ratios
	R_n <- g(F)

	# effect estimates
	est <- h(F, D, j_vec)
	alpha_n <- est[1]
	beta_n <- est[2]

	# standard error of effect estimate
	nabla_h <- grad_h(F, D, j_vec)
	se_beta_n <- sqrt(t(nabla_h) %*% cov(D) %*% nabla_h / n)

	# confidence interval
	ci <- beta_n + c(1,-1)*qnorm((1 - level)/2) * rep(se_beta_n, 2)

	# hypothesis test
	pval <- 2 * pnorm(-abs(beta_n / se_beta_n))

	# output
	out <- list(alpha = alpha_n, beta = beta_n, 
	            L_j = data.frame(j = 1:length(R_n), L_jn = R_n), 
	            se_beta = se_beta_n, ci = ci, pval = pval,
	            level = level)
	class(out) <- "trend_test"
	return(out)
}

#' Print the output of trend_test.
#' @param x An object of class \code{"trend_test"}. 
#' @param digits Number of digits to round output. 
#' @param ... Other options (not currently used)
#' @export
print.trend_test <- function(x, digits = 3, ...){
	tmp <- data.frame(beta = round(x$beta, digits),
	                  ci1 = round(x$ci[1], digits), 
	                  ci2 = round(x$ci[2], digits), 
	                  pval = round(x$pval, digits))
	colnames(tmp)[2:3] <- paste0(c("lower","upper"),"_",round(100*x$level), "%CI")
	cat("Trend in efficacy across failure type levels: \n")
	print(tmp)
}


#' Make a plot illustrating the true trend in vaccine efficacy
#' @export
#' @param file Name of file to save. If \code{NULL}
makeIllustrationPlot <- function(file = "~/Dropbox/Emory/hammingSieve/sim1.pdf"){
	grbg <- getTruth(returnL = TRUE)
	if(!is.null(file)){
		pdf(file, height = 3, width = 4.5)
	}
	par(mar = c(3.1, 3.1, 0.5, 0.5), mgp = c(1.5, 0.5, 0))
	plot(y = grbg$R, x = 0:4, bty = "n", xaxt="n", yaxt = "n",
	     ylim = c(0.2, 1.2), xlim = c(-0.5, 4.5),
	     xlab = "Genetic distance, j", ylab = expression(L["j,0"](tau)))
	axis(side = 1); axis(side = 2)
	# the largest and smallest value of beta_{0n}
	alphaVec <- c(1.10229, 1.0946205)
	betaVec <- c(-0.2021134, -0.1978171)
	abline(a = alphaVec[1], b = betaVec[1], lty = 2)
	abline(a = alphaVec[2], b = betaVec[2], lty = 3)
	legend(x = "topright", lty = 2:3, bty = "n",
	       c(expression(beta["0n"]*" = -0.202"), expression(beta["0n"]*" = -0.199")))
	if(!is.null(file)){
		dev.off()		
	}
}

#' Get multiply outputed cumulative incidence estimates and 
#' influence function estimates. 
#' 
#' @param rslt A list of \code{survtmle} objects 
#' 
#' @return A list with named entries \code{est} containing averaged
#' cumulative incidence estimates and \code{ic} containing averaged
#' influence function estimates. 
#' @export
getMO <- function(rslt){
	# number of outputations
    M <- length(rslt)
    est_matrix <- Reduce(cbind, lapply(rslt, "[[", "est"))

    # MO estimates
    est <- rowMeans(est_matrix)

	# MO influence functions
    ic <- Reduce("+", lapply(rslt, "[[", "ic"))/M

    out <- list(est = est, ic = ic)
    class(out) <- "getMO"
    return(out)
}