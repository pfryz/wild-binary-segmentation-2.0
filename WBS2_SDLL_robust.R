# The code below is a modification of WBS2.SDLL (https://doi.org/10.1007/s42952-020-00060-x)
# to the problem of detecting change-points in the piecewise-constant median of the data.
# The model is:
# X_t = f_t + e_t,
# where e_t is serially independent, has median zero, and does not need to satisfy any moment 
# conditions and can be arbitrarily heterogeneous.
#
# Example of use:
#
# f <- c(rep(10, 50), rep(5, 30), rep(0, 40), rep(3, 100), rep(6, 50))
# x <- exp(f + rnorm(length(f)))
# x.sol <- wbs2.sol(x, max.cusum.fun = max.cusum.signed)
# wbs2.median.modsel.sdll(x.sol)
#
# Alternatively:
# wbs2.median.modsel.th(x.sol)



random.cusums.for.median <- function(x, M, max.cusum.fun = max.cusum.signed, seed = NULL) {
	
	if(is.integer(seed)) set.seed(seed)

	# M=0 means there is a single cusum that gets taken over [1,length(x)] and therefore we are in the standard BS setting


	n <- length(x)
	

	M <- min(M, (n-1)*n/2)

		

	res <- matrix(0, max(1, M), 4)

	
	if ((n==2) || (M == 0)) ind <- matrix(c(1, n), 2, 1)
	
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

	res[,3:4] <- t(apply(ind, 2, max.cusum.fun, x))



	max.ind <- which.max(abs(res[,4]))



	max.val <- res[max.ind,,drop=F]



	list(res=res, max.val=max.val, M.eff=max(1, M))

}




max.cusum.signed <- function(ind, x) {
	
	
	rig <- sign(x - median(x))

	z <- cumsum(rig)	
	

	m <- ind[2]-ind[1]+1

	ip <- sqrt(((m-1):1) / m / (1:(m-1))) * z[1:(m-1)] - sqrt((1:(m-1)) / m / ((m-1):1)) * (z[m] - z[1:(m-1)])

	ip.max <- which.max(abs(ip))
	
	
	c(ip.max + ind[1] - 1, ip[ip.max])

	
}



wbs2.sol <- function(x, M=100, max.cusum.fun = max.cusum.signed) {
	
	# M = 0 means standard Binary Segmentation; M >= 1 means WBS2
	
	solutions.nested <- TRUE
	
	solution.set <- list()
	
	n <- length(x)
	
	sorted.cusums <- matrix(NA, 0, 4)
	
	if (n <= 1) solution.path <- integer()
	
	else {
		
		rc <- t(wbs.K.int(x, M, max.cusum.fun))

		ord <- order(abs(rc[,4]), decreasing=T)

		sorted.cusums <- abs(rc[ord,, drop=F])

		solution.path <- sorted.cusums[,3]
				
	}	
	
	ret = list(solutions.nested = solutions.nested, solution.path = solution.path, solution.set = solution.set, x = x, M = M, cands = sorted.cusums)

	class(ret) <- "cptpath"
	
	ret
	
}



wbs.K.int <- function(x, M, max.cusum.fun = max.cusum.signed) {

	# M = 0 means standard Binary Segmentation; M >= 1 means WBS2	

	n <- length(x)

	if (n == 1) return(matrix(NA, 4, 0))

	else {

		cpt <- t(random.cusums.for.median(x, M, max.cusum.fun)$max.val)

		return(cbind(cpt, wbs.K.int(x[1:cpt[3]], M, max.cusum.fun), wbs.K.int(x[(cpt[3]+1):n], M, max.cusum.fun) + c(rep(cpt[3], 3), 0)            ))

	}

}




mean.from.cpt <- function(x, cpt) {



	n <- length(x)

	len.cpt <- length(cpt)

	if (len.cpt) cpt <- sort(cpt)

	beg <- endd <- rep(0, len.cpt+1)

	beg[1] <- 1

	endd[len.cpt+1] <- n

	if (len.cpt) {

		beg[2:(len.cpt+1)] <- cpt+1

		endd[1:len.cpt] <- cpt

	}

	means <- rep(0, len.cpt+1)

	for (i in 1:(len.cpt+1)) means[i] <- mean(x[beg[i]:endd[i]])

	rep(means, endd-beg+1)

}


wbs2.median.modsel.sdll <- function(wbs2.sol.obj, th.const = sqrt(3/2), th.const.min.mult = 0.5) {
	
	sdll(wbs2.sol.obj, 1, F, 0, th.const, th.const.min.mult)
	
}


wbs2.median.modsel.th <- function(wbs2.sol.obj, th.const = sqrt(3/2)) {
	
	x <- wbs2.sol.obj$x
	n <- length(x)
	
    if (n <= 1) {

        est <- x

        no.of.cpt <- 0

        cpts <- integer(0)

    }

    else {
	
		th <- th.const * sqrt(2 * log(n))
		
		if (wbs2.sol.obj$cands[1,4] < th) {

    	    no.of.cpt <- 0

        	cpts <- integer(0)
        	
        }
        
        else {
        	
        	indices <- which(wbs2.sol.obj$cands[,4] > th)

        	cpts <- sort(wbs2.sol.obj$cands[indices, 3])

			no.of.cpt <- length(cpts)

        }

	est <- mean.from.cpt(x, cpts)

 	}

	list(est=est, no.of.cpt=no.of.cpt, cpts=cpts)
	
	
}


sdll <- function(sols.object, sigma = stats::mad(diff(sols.object$x)/sqrt(2)), universal = TRUE, M = NULL, th.const = NULL, th.const.min.mult = 0.3, lambda = 0.9) {

	x <- sols.object$x

	n <- length(x)

    if (n <= 1) {

        est <- x

        no.of.cpt <- 0

        cpts <- integer(0)

    }

    else {

		if (sigma == 0) {
			
			s0 <- all.shifts.are.cpts(x)
			est <- s0$est
			no.of.cpt <- s0$no.of.cpt
			cpts <- s0$cpts
			
		} else {

		if (universal) {

        	u <- universal.M.th.v3(n, lambda)

        	th.const <- u$th.const

        	M <- u$M

    	}

    	else if (is.null(M) || is.null(th.const)) stop("If universal is FALSE, then M and th.const must be specified.")

    	th.const.min <- th.const * th.const.min.mult

    	th <- th.const * sqrt(2 * log(n)) * sigma

    	th.min <- th.const.min * sqrt(2 * log(n)) * sigma




 		if (sols.object$cands[1,4] < th) {


    	    no.of.cpt <- 0

        	cpts <- integer(0)



 		}

		else {

			indices <- which(sols.object$cands[,4] > th.min)

			if (length(indices) == 1) {

				cpts <- sols.object$cands[indices, 3]

				no.of.cpt <- 1


			}

			else {

				rc.sel <- sols.object$cands[indices,,drop=F]

				z <- sols.object$cands[indices,4]




				z.l <- length(z)

				dif <- -diff(log(z))

				dif.ord <- order(dif, decreasing=T)

				j <- 1

				while ((j < z.l) & (z[dif.ord[j]+1] > th)) j <- j+1

				if (j < z.l) no.of.cpt <- dif.ord[j] else no.of.cpt <- z.l

				cpts <- sort(sols.object$cands[1:no.of.cpt,3])			


			}

		} 

	est <- mean.from.cpt(x, cpts)

    }

	}

	list(est=est, no.of.cpt=no.of.cpt, cpts=cpts)

}






universal.M.th.v3 <- function(n, lambda = 0.9) {

		

	mat.90 <- matrix(0, 24, 3)

	mat.90[,1] <- c(10, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)

	mat.90[,2] <- c(1.420, 1.310, 1.280, 1.270, 1.250, 1.220, 1.205, 1.205, 1.200, 1.200, 1.200, 1.185, 1.185, 1.170, 1.170, 1.160, 1.150, 1.150, 1.150, 1.150, 1.145, 1.145, 1.135, 1.135)

	mat.90[,3] <- rep(100, 24)

	

	mat.95 <- matrix(0, 24, 3)

	mat.95[,1] <- c(10, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)

	mat.95[,2] <- c(1.550, 1.370, 1.340, 1.320, 1.300, 1.290, 1.265, 1.265, 1.247, 1.247, 1.247, 1.225, 1.225, 1.220, 1.210, 1.190, 1.190, 1.190, 1.190, 1.190, 1.190, 1.180, 1.170, 1.170)

	mat.95[,3] <- rep(100, 24)



	if (lambda == 0.9) A <- mat.90 else A <- mat.95



	d <- dim(A)

	if (n < A[1,1]) {

		th <- A[1,2]

		M <- A[1,3]

	}

	else if (n > A[d[1],1]) {

		th <- A[d[1],2]

		M <- A[d[1],3]

	}

	else {

		ind <- order(abs(n - A[,1]))[1:2]

		s <- min(ind)

		e <- max(ind)

		th <- A[s,2] * (A[e,1] - n)/(A[e,1] - A[s,1]) + A[e,2] * (n - A[s,1])/(A[e,1] - A[s,1])

		M <- A[s,3] * (A[e,1] - n)/(A[e,1] - A[s,1]) + A[e,3] * (n - A[s,1])/(A[e,1] - A[s,1])

	}



	list(th.const=th, M=M)

}



