wbs.sdll.cpt <- function(x, sigma = stats::mad(diff(x)/sqrt(2)), universal = TRUE, M = NULL, th.const = NULL, th.const.min.mult = 0.3, lambda = 0.9, cusums = "systematic") {
	
	n <- length(x)
    if (n <= 1) {
        no.of.cpt <- 0
        cpt <- integer(0)
    }
    else {
		if (sigma == 0) stop("Noise level estimated at zero; therefore no change-points to estimate.")
		if (universal) {
        	u <- universal.M.th.v3(n, lambda)
        	th.const <- u$th.const
        	M <- u$M
    	}
    	else if (is.null(M) || is.null(th.const)) stop("If universal is FALSE, then M and th.const must be specified.")
    	th.const.min <- th.const * th.const.min.mult
    	th <- th.const * sqrt(2 * log(n)) * sigma
    	th.min <- th.const.min * sqrt(2 * log(n)) * sigma

		if (cusums == "random") cusum.sampling <- random.cusums else if (cusums == "systematic") cusum.sampling <- systematic.cusums

 		rc <- t(wbs.K.int(x, M, cusum.sampling))
 		if (max(abs(rc[,4])) < th) {
    	    no.of.cpt <- 0
        	cpt <- integer(0)

 		}
		else {
			indices <- which(abs(rc[,4]) > th.min)
			if (length(indices) == 1) {
				cpt <- rc[indices, 3]
				no.of.cpt <- 1
			}
			else {
				rc.sel <- rc[indices,,drop=F]
				ord <- order(abs(rc.sel[,4]), decreasing=T)
				z <- abs(rc.sel[ord,4])
				z.l <- length(z)
				dif <- -diff(log(z))
				dif.ord <- order(dif, decreasing=T)
				j <- 1
				while ((j < z.l) & (z[dif.ord[j]+1] > th)) j <- j+1
				if (j < z.l) no.of.cpt <- dif.ord[j] else no.of.cpt <- z.l
				cpt <- sort((rc.sel[ord,3])[1:no.of.cpt])			
			}
		} 
    }
    est <- mean.from.cpt(x, cpt)
	list(est=est, no.of.cpt=no.of.cpt, cpt=cpt)
}


wbs.sdll.cpt.rep <- function(x, sigma = stats::mad(diff(x)/sqrt(2)), universal = TRUE, M = NULL, th.const = NULL, th.const.min.mult = 0.3, lambda = 0.9, repeats = 9) {

	res <- vector("list", repeats)
	
	cpt.combined <- integer(0)
	
	nos.of.cpts <- rep(0, repeats)
	
	for (i in 1:repeats) {
		
		res[[i]] <- wbs.sdll.cpt(x, sigma, universal, M, th.const, th.const.min.mult, lambda, "random")
		cpt.combined <- c(cpt.combined, res[[i]]$cpt)
		nos.of.cpts[i] <- res[[i]]$no.of.cpt				
		
	}

	med.no.of.cpt <- median(nos.of.cpts)
	
	med.index <- which.min(abs(nos.of.cpts - med.no.of.cpt))
	
	med.run <- res[[med.index]]
	
	list(med.run = med.run, cpt.combined = sort(cpt.combined))

}



wbs.K.int <- function(x, M, cusum.sampling) {
	
	n <- length(x)
	if (n == 1) return(matrix(NA, 4, 0))
	else {
		cpt <- t(cusum.sampling(x, M)$max.val)
		return(cbind(cpt, wbs.K.int(x[1:cpt[3]], M, cusum.sampling), wbs.K.int(x[(cpt[3]+1):n], M, cusum.sampling) + c(rep(cpt[3], 3), 0)            ))
	}
	
}



random.cusums <- function(x, M) {

	y <- c(0, cumsum(x))

	n <- length(x)
	
	M <- min(M, (n-1)*n/2)
		
	res <- matrix(0, M, 4)
	
	if (n==2) ind <- matrix(c(1, 2), 2, 1)
	else if (M == (n-1)*n/2) {
		ind <- matrix(0, 2, M)
		ind[1,] <- rep(1:(n-1), (n-1):1)
		ind[2,] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)
	}
	else {
		ind <- ind2 <- matrix(floor(runif(2*M) * (n-1)), nrow=2)
		ind2[1,] <- apply(ind, 2, min)
		ind2[2,] <- apply(ind, 2, max)
		ind <- ind2 + c(1, 2)
	}

	res[,1:2] <- t(ind)
	res[,3:4] <- t(apply(ind, 2, max.cusum, y))

	max.ind <- which.max(abs(res[,4]))

	max.val <- res[max.ind,,drop=F]

	list(res=res, max.val=max.val, M.eff=M)

}


systematic.cusums <- function(x, M) {
	
	y <- c(0, cumsum(x))

	n <- length(x)
	
	M <- min(M, (n-1)*n/2)
			
	ind <- grid.intervals(n, M)

	M <- dim(ind)[2]

	res <- matrix(0, M, 4)

	res[,1:2] <- t(ind)
	res[,3:4] <- t(apply(ind, 2, max.cusum, y))

	max.ind <- which.max(abs(res[,4]))

	max.val <- res[max.ind,,drop=F]

	list(res=res, max.val=max.val, M.eff=M)
		
}


max.cusum <- function(ind, y) {
	
		z <- y[(ind[1]+1):(ind[2]+1)] - y[ind[1]]
		m <- ind[2]-ind[1]+1
		ip <- sqrt(((m-1):1) / m / (1:(m-1))) * z[1:(m-1)] - sqrt((1:(m-1)) / m / ((m-1):1)) * (z[m] - z[1:(m-1)])
		ip.max <- which.max(abs(ip))
		
		c(ip.max + ind[1] - 1, ip[ip.max])

}
