titvFPEst <- function(titvExpected, titvObserved) { 1 - (titvObserved - 0.5) / (titvExpected - 0.5) }

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

cumhist <- function(d) {
    h <- hist(d)
    #plot(h$mids, cumsum(h$count), type="b", col="orange", lwd=2)
}
