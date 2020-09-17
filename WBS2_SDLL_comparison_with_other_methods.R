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
