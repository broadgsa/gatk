titvFPEst <- function(titvExpected, titvObserved) { max(min(1 - (titvObserved - 0.5) / (titvExpected - 0.5), 1), 0.001) }

titvFPEstV <- function(titvExpected, titvs) {
    sapply(titvs, function(x) titvFPEst(titvExpected, x))
}

calcHet <- function(nknown, knownTiTv, nnovel, novelTiTv, callable) {
  TP <- nknown + (1-titvFPEst(knownTiTv, novelTiTv)) * nnovel
  2 * TP / 3 / callable
}

marginalTiTv <- function( nx, titvx, ny, titvy ) {
    tvx = nx / (titvx + 1)
    tix = nx - tvx
    tvy = ny / (titvy + 1)
    tiy = ny - tvy
    tiz = tix - tiy
    tvz = tvx - tvy
    return(tiz / tvz)
}
marginaldbSNPRate <- function( nx, dbx, ny, dby ) {
    knownx = nx * dbx / 100
    novelx = nx - knownx
    knowny = ny * dby / 100
    novely = ny - knowny
    knownz = knownx - knowny
    novelz = novelx - novely
    return(knownz / ( knownz + novelz ) * 100)
}

numExpectedCalls <- function(L, theta, calledFractionOfRegion, nIndividuals, dbSNPRate) {
    nCalls <- L * theta * calledFractionOfRegion * sum(1 / seq(1, 2 * nIndividuals))
    return(list(nCalls = nCalls, nKnown = dbSNPRate * nCalls, nNovel = (1-dbSNPRate) * nCalls))
}

normalize <- function(x) {
    x / sum(x)
}

normcumsum <- function(x) {
    cumsum(normalize(x))
}

cumhist <- function(d, ...) {
    plot(d[order(d)], type="b", col="orange", lwd=2, ...)
}

revcumsum <- function(x) {
   return(rev(cumsum(rev(x))))
}

phred <- function(x) {
    log10(max(x,10^(-9.9)))*-10
}

pOfB <- function(b, B, Q) {
    #print(paste(b, B, Q))
    p = 1 - 10^(-Q/10)
    if ( b == B )
        return(p) 
    else
        return(1 - p)
}

pOfG <- function(bs, qs, G) {
    a1 = G[1]
    a2 = G[2]
    
    log10p = 0
    for ( i in 1:length(bs) ) {
        b = bs[i]
        q = qs[i]
        p1 = pOfB(b, a1, q) / 2 + pOfB(b, a2, q) / 2
        log10p = log10p + log10(p1)
    }
    
    return(log10p)
}

pOfGs <- function(nAs, nBs, Q) {
    bs = c(rep("a", nAs), rep("t", nBs))
    qs = rep(Q, nAs + nBs)
    G1 = c("a", "a")
    G2 = c("a", "t")
    G3 = c("t", "t")

    log10p1 = pOfG(bs, qs, G1)
    log10p2 = pOfG(bs, qs, G2)
    log10p3 = pOfG(bs, qs, G3)
    Qsample = phred(1 - 10^log10p2 / sum(10^(c(log10p1, log10p2, log10p3))))
    
    return(list(p1=log10p1, p2=log10p2, p3=log10p3, Qsample=Qsample)) 
}

QsampleExpected <- function(depth, Q) {
    weightedAvg = 0
    for ( d in 1:(depth*3) ) {
        Qsample = 0
        pOfD = dpois(d, depth)
        for ( nBs in 0:d ) {
            pOfnB = dbinom(nBs, d, 0.5)
            nAs = d - nBs
            Qsample = pOfGs(nAs, nBs, Q)$Qsample
            #Qsample = 1
            weightedAvg = weightedAvg + Qsample * pOfD * pOfnB
            print(as.data.frame(list(d=d, nBs = nBs, pOfD=pOfD, pOfnB = pOfnB, Qsample=Qsample, weightedAvg = weightedAvg)))
        }
    }
    
    return(weightedAvg)
}

plotQsamples <- function(depths, Qs, Qmax) {
    cols = rainbow(length(Qs))
    plot(depths, rep(Qmax, length(depths)), type="n", ylim=c(0,Qmax), xlab="Average sequencing coverage", ylab="Qsample", main = "Expected Qsample values, including depth and allele sampling")
    
    for ( i in 1:length(Qs) ) {
        Q = Qs[i]
        y = as.numeric(lapply(depths, function(x) QsampleExpected(x, Q)))
        points(depths, y, col=cols[i], type="b")
    }

    legend("topleft", paste("Q", Qs), fill=cols)
}

pCallHetGivenDepth <- function(depth, nallelesToCall) {
  depths = 0:(2*depth)   
  pNoAllelesToCall = apply(as.matrix(depths),1,function(d) sum(dbinom(0:nallelesToCall,d,0.5)))
  dpois(depths,depth)*(1-pNoAllelesToCall)
}

pCallHets <- function(depth, nallelesToCall) {
  sum(pCallHetGivenDepth(depth,nallelesToCall))
}

pCallHetMultiSample <- function(depth, nallelesToCall, nsamples) {
  1-(1-pCallHets(depth,nallelesToCall))^nsamples
}
