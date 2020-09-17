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
		
		rc <- t(wbs.K.int.for.median(x, M, max.cusum.fun))

		ord <- order(abs(rc[,4]), decreasing=T)

		sorted.cusums <- abs(rc[ord,, drop=F])

		solution.path <- sorted.cusums[,3]
				
	}	
	
	ret = list(solutions.nested = solutions.nested, solution.path = solution.path, solution.set = solution.set, x = x, M = M, cands = sorted.cusums)

	class(ret) <- "cptpath"
	
	ret
	
}



wbs.K.int.for.median <- function(x, M, max.cusum.fun) {

	# M = 0 means standard Binary Segmentation; M >= 1 means WBS2	

	n <- length(x)

	if (n == 1) return(matrix(NA, 4, 0))

	else {

		cpt <- t(random.cusums.for.median(x, M, max.cusum.fun)$max.val)

		return(cbind(cpt, wbs.K.int.for.median(x[1:cpt[3]], M, max.cusum.fun), wbs.K.int.for.median(x[(cpt[3]+1):n], M, max.cusum.fun) + c(rep(cpt[3], 3), 0)            ))

	}

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
