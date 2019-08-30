#################################
#   WBS2.SDLL code
#################################

wbs.sdll.cpt <- function(x, sigma = stats::mad(diff(x)/sqrt(2)), universal = TRUE, M = NULL, th.const = NULL, th.const.min.mult = 0.3, lambda = 0.9) {
	
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

 		rc <- t(wbs.K.int(x, M))
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

wbs.K.int <- function(x, M) {
	
	n <- length(x)
	if (n == 1) return(matrix(NA, 4, 0))
	else {
		cpt <- t(random.cusums(x, M)$max.val)
		return(cbind(cpt, wbs.K.int(x[1:cpt[3]], M), wbs.K.int(x[(cpt[3]+1):n], M) + c(rep(cpt[3], 3), 0)            ))
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

max.cusum <- function(ind, y) {
	
		z <- y[(ind[1]+1):(ind[2]+1)] - y[ind[1]]
		m <- ind[2]-ind[1]+1
		ip <- sqrt(((m-1):1) / m / (1:(m-1))) * z[1:(m-1)] - sqrt((1:(m-1)) / m / ((m-1):1)) * (z[m] - z[1:(m-1)])
		ip.max <- which.max(abs(ip))
		
		c(ip.max + ind[1] - 1, ip[ip.max])

}


###################################
#   Comparison with other methods
###################################
library(changepoint)
library(breakfast) 
library(FDRSeg)
library(Segmentor3IsBack)
library(stepR)
library(wbs)
library(ecp)
library(strucchange)
library(IDetect)
library(cumSeg)
library(cpm)
library(not)
library(fpop)
library(mosum)
library(jointseg)
library(capushe)



cpt.detect <- function(x) {
		
	setClass("cpt.est", representation(cpt="numeric", nocpt="numeric", mean="numeric", exectime="numeric"), prototype(cpt=numeric(0), nocpt=0, mean=numeric(0), exectime=0))

	print("pelt.mbic")
#	uses package "changepoint"
	pelt.mbic <- new("cpt.est")
	pelt.mbic@exectime <- as.numeric(system.time(z <- cpt.mean(x/mad(diff(x)/sqrt(2)), method="PELT"))[1])
	if (length(z@cpts) > 1) pelt.mbic@cpt <- as.numeric(z@cpts[1:(length(z@cpts)-1)]) else pelt.mbic@cpt <- integer(0)
	pelt.mbic@nocpt <- length(pelt.mbic@cpt)
	pelt.mbic@mean <- mean.from.cpt(x, pelt.mbic@cpt)

	print("pelt.bic")
#	uses package "changepoint"
	pelt.bic <- new("cpt.est")
	pelt.bic@exectime <- as.numeric(system.time(z <- cpt.mean(x/mad(diff(x)/sqrt(2)), penalty="BIC", method="PELT"))[1])
	if (length(z@cpts) > 1) pelt.bic@cpt <- as.numeric(z@cpts[1:(length(z@cpts)-1)]) else pelt.bic@cpt <- integer(0)
	pelt.bic@nocpt <- length(pelt.bic@cpt)
	pelt.bic@mean <- mean.from.cpt(x, pelt.bic@cpt)
	
	print("mosum.lop")
#	uses package "mosum"
	mosum.lop <- new("cpt.est")
	mosum.lop@exectime <- as.numeric(system.time(z <- multiscale.localPrune(x))[1])
	if (length(z$cpts)) mosum.lop@cpt <- z$cpts else mosum.lop@cpt <- integer(0)
	mosum.lop@nocpt <- length(mosum.lop@cpt)
	mosum.lop@mean <- mean.from.cpt(x, mosum.lop@cpt)
			
	print("IDetect")
#	uses package "IDetect"
	id <- new("cpt.est")
	id@exectime <- as.numeric(system.time(z <- ID(x))[1])
	if (z$no_cpt == 0) id@cpt <- integer(0) else id@cpt <- z$cpt
	id@nocpt <- length(id@cpt)
	id@mean <- mean.from.cpt(x, id@cpt)
	
	print("FDRSeg")
#	uses package "FDRSeg"
	fdrs <- new("cpt.est")
	fdrs@exectime <- as.numeric(system.time(z <- fdrseg(x))[1])
	if (length(z$left) == 1) fdrs@cpt <- integer(0) else fdrs@cpt <- as.numeric(z$left[2:length(z$left)]) - 1
	fdrs@nocpt <- length(fdrs@cpt)
	fdrs@mean <- mean.from.cpt(x, fdrs@cpt)
	
	print("S3IB")
#	uses package "Segmentor3IsBack"
	S3IB <- new("cpt.est")
	S3IB@exectime <- as.numeric(system.time(z <- Segmentor(x, model=2, Kmax = trunc(length(x)/3)))[1])
	S3IB@nocpt <- SelectModel(z)-1
	if (S3IB@nocpt) S3IB@cpt <- as.numeric(z@breaks[S3IB@nocpt+1, 1:S3IB@nocpt]) else S3IB@cpt <- integer(0)
	S3IB@mean <- mean.from.cpt(x, S3IB@cpt)	
	
	print("Smuce")
#	uses package "stepR"
	smuce <- new("cpt.est")
	smuce@exectime <- as.numeric(system.time(z <- stepFit(x, alpha=0.5))[1])
	smuce@nocpt <- length(z$rightEnd) - 1
	if (smuce@nocpt) smuce@cpt <- as.numeric(z$rightEnd[1:smuce@nocpt]) else smuce@cpt <- integer(0)
	smuce@mean <- mean.from.cpt(x, smuce@cpt)
	
	print("cumSeg")
#	uses package "cumSeg"
	cumS <- new("cpt.est")
	cumS@exectime <- as.numeric(system.time(z <- jumpoints(x, k = trunc(length(x)/3)))[1])
	if (z$n.psi) cumS@cpt <- as.numeric(z$psi) else cumS@cpt <- integer(0)
	cumS@nocpt <- length(cumS@cpt)
	cumS@mean <- mean.from.cpt(x, cumS@cpt)
	
	print("fpop")
#	uses package "fpop"
#	fpop is not available from CRAN and seems to depend on something that isn't always available; so sometimes does not execute
	fp <- new("cpt.est")
	fp@exectime <- as.numeric(system.time(z <- Fpop(x/mad(diff(x)/sqrt(2)), 2 * log(length(x))))[1])
	if (z$K > 1) fp@cpt <- z$t.est[1:(z$K-1)] else fp@cpt <- integer(0)
	fp@nocpt <- length(fp@cpt)
	fp@mean <- mean.from.cpt(x, fp@cpt)
	
	print("jointseg+DDSE")
#	uses packages "jointseg" and "capushe"
 	js.se <- new("cpt.est")
	Kmax <- trunc(length(x)/3)
	ptm <- proc.time()
  	res <- Fpsn(x, Kmax=Kmax)
  	n <- length(x)
  	dataCapushe <- data.frame(name=1:Kmax,
            pen.shape=lchoose(n-1,0:(Kmax-1)), 
            complexity=1:Kmax, contrast=res$J.est)
  	KC <- try(DDSE(dataCapushe))
    KJ <- Djump(dataCapushe)
  	if(class(KC) == "try-error"){ K <- as.integer(KJ@model) } else { K <- as.integer(KC@model) }
 	tmp <- res$t.est[K, 1:K]
	js.se@exectime <- as.numeric((proc.time() - ptm)[1])
	if (length(tmp) == 1) js.se@cpt <- integer(0) else js.se@cpt <- tmp[1:(length(tmp)-1)]
	js.se@nocpt <- length(js.se@cpt)
	js.se@mean <- mean.from.cpt(x, js.se@cpt)

	print("jointseg+Djump")
#	uses packages "jointseg" and "capushe"
 	js.dj <- new("cpt.est")
	Kmax <- trunc(length(x)/3)
  	n <- length(x)
	
	ptm <- proc.time()
  	res <- Fpsn(x, Kmax=Kmax)
  	dataCapushe <- data.frame(name=1:Kmax,
            pen.shape=lchoose(n-1,0:(Kmax-1)), 
            complexity=1:Kmax, contrast=res$J.est)
  	KJ <- Djump(dataCapushe)
	K <- as.integer(KJ@model)
 	tmp <- res$t.est[K, 1:K]
	js.dj@exectime <- as.numeric((proc.time() - ptm)[1])
	if (length(tmp) == 1) js.dj@cpt <- integer(0) else js.dj@cpt <- tmp[1:(length(tmp)-1)]
	js.dj@nocpt <- length(js.dj@cpt)
	js.dj@mean <- mean.from.cpt(x, js.dj@cpt)
#	This often returns a warning: "1: In Djump(dataCapushe) : There are several maximum jump"

	
	print("tguh")
#	uses package "breakfast"
	tg <- new("cpt.est")
	tg@exectime <- as.numeric(system.time(z <- tguh.cpt(x))[1])
	tg@cpt <- z$cpt
	tg@nocpt <- z$no.of.cpt
	tg@mean <- mean.from.cpt(x, tg@cpt)
	
	print("wbs.thresh.1.0")
#	uses package "breakfast"
	wt10 <- new("cpt.est")
	wt10@exectime <- as.numeric(system.time(z <- wbs.thresh.cpt(x, universal=F, M = 5000, th.const = 1.0, adapt = F))[1])
	wt10@cpt <- z$cpt
	wt10@nocpt <- z$no.of.cpt
	wt10@mean <- mean.from.cpt(x, wt10@cpt)

	print("wbs.thresh.1.3")
#	uses package "breakfast"
	wt13 <- new("cpt.est")
	wt13@exectime <- as.numeric(system.time(z <- wbs.thresh.cpt(x, universal=F, M = 5000, th.const = 1.3, adapt = F))[1])
	wt13@cpt <- z$cpt
	wt13@nocpt <- z$no.of.cpt
	wt13@mean <- mean.from.cpt(x, wt13@cpt)

	print("wbs.bic")
#	uses package "breakfast"
	wb <- new("cpt.est")
	wb@exectime <- as.numeric(system.time(z <- wbs.bic.cpt(x, Kmax = trunc(length(x)/3)))[1])
	wb@cpt <- z$cpt
	wb@nocpt <- z$no.of.cpt
	wb@mean <- mean.from.cpt(x, wb@cpt)

	print("wbs2.90")
	wbs2.90 <- new("cpt.est")
	wbs2.90@exectime <- as.numeric(system.time(z <- wbs.sdll.cpt(x, lambda=0.9))[1])
	wbs2.90@cpt <- z$cpt
	wbs2.90@nocpt <- z$no.of.cpt
	wbs2.90@mean <- mean.from.cpt(x, wbs2.90@cpt)	

	print("wbs2.95")
	wbs2.95 <- new("cpt.est")
	wbs2.95@exectime <- as.numeric(system.time(z <- wbs.sdll.cpt(x, lambda=0.95))[1])
	wbs2.95@cpt <- z$cpt
	wbs2.95@nocpt <- z$no.of.cpt
	wbs2.95@mean <- mean.from.cpt(x, wbs2.95@cpt)	

	
	list(pelt.mbic=pelt.mbic, pelt.bic=pelt.bic, mosum.lop=mosum.lop, id=id, fdrs=fdrs, S3IB=S3IB, smuce=smuce, cumS=cumS, fp=fp, js.se=js.se, js.dj=js.dj, tg=tg, wt10=wt10, wt13=wt13, wb=wb, wbs2.90=wbs2.90, wbs2.95=wbs2.95)

}	

