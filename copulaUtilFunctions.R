## This is collection of utility functions related to copula computations with R
## Mainly with VineCopula, copula and rugarch packages
## Last edit: 1 December 2016
## Athanassios Stavrakoudis
## http://stavrakoudis.econ.uoi.gr
## astavrak@uoi.gr

## optional, in order to set up the library directory (uncomment)
#.libPaths('/usr/local/lib/R/site-library')

library(stats)
library(dplyr)
library(xtable)
library(timeSeries)
library(FinTS)
library(nortest)
library(moments)
library(WeightedPortTest)

library(VineCopula)
library(gofCopula)
library(copula)
library(rugarch)
#library(rmgarch)
library(vars)
library(ggplot2)
library(gridExtra)

call_time <- function(t=0) {
    ct <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    if (t==0) {
        return (ct)
    } else {
        if (t==1) {
            c = "\nStart :"
        } else if (t==2) {
            c = "\nEnd   :"
        } else {
            c = "\nCall  :"
        }
    }
    cat(c, ct)
    cat("\n---------------------------\n\n")
}


re.Window <- function(x, T1, T2) { return (x[T1:T2]) }

se  <- function(x) sqrt( var(x, na.rm=TRUE) / length(na.omit(x)) )

calcReturn <- function (x, lr=TRUE, per=FALSE, d=-1) {
    if (lr==TRUE) {
        r <- diff(log(x), lag=1)
    } else {
        r <- diff(x) / x[1:(length(x)-1)]
    }
    if (per==TRUE) {
        r <- 100 * r
    }
    if(d > -1) {
        r <- round(r, d)
    }
    return (r)
}

resid_std <- function(resid)
{
    return ( resid / sd(resid) )
}

resid_adj <- function(resid)
{
    return ( (resid-mean(resid)) / sd(resid) )
}

calcRank2 <- function (x, N=c(2,1), k=2)
{
    source("bestArmaOrder.R")
    bao <- bestArmaOrder(x, N=N, k=k)
    ag.fit <- bao$ag
    st <- capture.output(show(ag.fit))
    cat ("arma  : ", bao$bao, "\n")
    cat ("garch : ", bao$gao, "\n")
    cat ("distr : ", bao$distr, "\n")
    cat ("Ljung-Box test                      p-value", "\n")
    i <- grep ("Weighted Ljung-Box Test on Standardized Residuals+", st)
    cat (st[i+3], "\n")
    cat (st[i+4], "\n")
    cat (st[i+5], "\n")
    cat ("ARCH LM test                       p-value", "\n")
    i <- grep ("ARCH Lag+", st)[1]
    cat (st[i+0], "\n")
    cat (st[i+1], "\n")
    cat (st[i+2], "\n")
    pv1 <- cvm.test(ag.fit@fit$resid)$p.v
    pv1 <- round(pv1, 4)
    pv2 <- lillie.test(ag.fit@fit$resid)$p.v
    pv2 <- round(pv2, 4)
    cat ("normality,  CM: ", pv1, "KS:", pv2, "\n")
    cat ("\n")
    return (list(rank=bao$rank, resid=bao$resid))
}

calcRank <- function(x)
{
    x.std <- x / sd(x)
    P     <- ecdf(x.std)
    T     <- length(x.std)
    r     <- P(x.std) * (T/(T+1))
    return (r)
}

calc_uG_Rank <- function(mod)
{
    r <- residuals(mod, standardize=T)
    r <- as.numeric(r)
    P <- ecdf(r)
    T <- length(r)
    u <- P(r) * (T/(T+1))
    return(u)
}

BiCop_Estimate <- function (u1, u2, fam=NA, Nb=0)
{
    if (nargs() < 3) { stop("usage: BiCop_Estimate(u1, u2, fam=fam, Nb=0), need at least 3 arguments \n")}
    Len  <- min(length(u1), length(u2))
    u1   <- u1[1:Len]
    u2   <- u2[1:Len]
    est  <- VineCopula::BiCopEst(u1, u2, family=fam, se=TRUE)
    est  <- summary(est)
    if (Nb > 2) {
        pval <- gof.Copula(u1, u2, fam, Nb)
    } else {
        pval <- NA
    }
    cop  <- list(fam=est$family, fname=est$familyname, p.val=pval, par=est$par, par2=est$par2, tau=est$tau,
                 low=est$taildep$lower, upp=est$taildep$upper, se=est$se, se2=est$se2,
                 logLik=est$logLik, AIC=est$AIC, BIC=est$BIC)
    return (cop)
}

gof.Copula <- function (u1, u2, fam, Nb=100)
{
    if (nargs() < 3) {
        stop("usage: gof.Copula (u1, u2, fam, Nb=100) \n")
    }
    set.seed(12345)
    rnorm(1000)
    if (fam == 2) {
        x     <- matrix(c(u1, u2), ncol=2)
        par   <- cor(u1, u2, method="pearson")
        copu  <- copula::tCopula(par, dim=2, dispstr="un", df.fixed=TRUE)
        gof   <- copula::gofCopula(copu, x, N=Nb, simulation="mult", method="Rn",
                          estim.method="itau", optim.method="BFGS",
                          optim.control=list(reltol=1e-6, maxit=100000))
        p.CvM <- gof$p.value
    } else {
        #p.CvM <- CDVine::BiCopGofKendall(u1, u2, fam, B=Nb)$p.value.CvM
        est  <- try( VineCopula::BiCopGofTest(u1, u2, fam, method="kendall", B=Nb) )
        if(class(est) == "try-error") {
            p.CvM  <- NA
        } else {
            p.CvM  <- est$p.value.CvM
        }

    }
    fName     <- VineCopula::BiCopName(fam)
    #cat("family =", fam, fName, "    p-value =", p.CvM, "\n")
    return (p.CvM)
}

selectCopulaClarke <- function(r1, r2, fams=1:10, Nb=100) {
  VCtest      <- try( VineCopula::BiCopVuongClarke(r1, r2, fams, correction="Schwarz") ) # )
  if(class(VCtest) == "try-error") {
      fam.sel <- 0
      VCtest  <- TDep <- Tau <- p.val <- par <- par2 <- tau <- low <- upp <- NA
  } else {
    print(VCtest)
    Clarke <- as.vector(t(VCtest)[,2])
    famSel <- which(Clarke %in% max(Clarke))
    famSel <- fams[famSel]

        if (length(famSel) > 0) {
            p.val    <- -1
            fam.sel  <- -1
            for(f in famSel) {
                fName <- VineCopula::BiCopName(f)
                p.CvM <- gof.Copula(r1, r2, f, Nb=Nb)
                cat ("family =", f, fName, "    p.value =", p.CvM, "\n\n")
                if (p.CvM > p.val) {
                    p.val   <- p.CvM
                    fam.sel <- f
                }
            }
        } else {
            fam.sel <- famSel
        }
    }
    #cat ("selected:", fam.sel, "\n")

    est  <- try( VineCopula::BiCopEst(r1, r2, fam.sel) )
    if(class(est) == "try-error") {
        VCtest <- TDep <- Tau <- p.val <- par <- par2 <- tau <- low <- upp <- NA
    } else {
        Tau  <- VineCopula::BiCopPar2Tau(fam.sel, est$par, est$par2)
        TDep <- VineCopula::BiCopPar2TailDep(fam.sel, est$par, est$par2)
    }
    cop  <- list(VCtest=VCtest, fam=fam.sel, rank1=r1, rank2=r2, p.val=p.val, par=est$par, par2=est$par2, tau=Tau, low=TDep$lower, upp=TDep$upper)
    #print (reportCopulaRow (cop) )
    return (cop)
}

selectCopulaVuong <- function(u1, u2, fams=1:10, Nb=100) {
    Len    <- min(length(u1), length(u2))
    u1     <- u1[1:Len]
    u2     <- u2[1:Len]
    VCtest <- try( VineCopula::BiCopVuongClarke(u1, u2, fams, correction="Schwarz") )
    if(class(VCtest) == "try-error") {
        fam.sel <- 0
        VCtest  <- TDep <- Tau <- p.val <- par <- par2 <- tau <- low <- upp <- NA
    } else {
        #print(VCtest)
        Vuong   <- as.vector(t(VCtest)[,1])
        famSel  <- which(Vuong %in% max(Vuong))
        famSel  <- fams[famSel]
        fam.sel <- NA

        max1    <- ifelse(sum((Vuong==max(Vuong)))==1, TRUE, FALSE)
        if (length(famSel) > 0) {
            p.val    <- -1
            fam.sel  <- -1
            for(f in famSel) {
                fName <- VineCopula::BiCopName(f)
                p.CvM <- gof.Copula(u1, u2, f, Nb=Nb)
                #cat ("family =", f, fName, "    p.value =", p.CvM, "\n")
                if (p.CvM > p.val) {
                    p.val   <- p.CvM
                    fam.sel <- f
                }
            }
        } else {
            fam.sel <- famSel
        }
    }
    #cat ("selected:", fam.sel, "\n")

    est  <- try( VineCopula::BiCopEst(u1, u2, fam.sel) )
    if(class(est) == "try-error") {
        VCtest <- TDep <- Tau <- p.val <- par <- par2 <- tau <- low <- upp <- NA
    } else {
        Tau  <- VineCopula::BiCopPar2Tau(fam.sel, est$par, est$par2)
        TDep <- VineCopula::BiCopPar2TailDep(fam.sel, est$par, est$par2)
    }
    cop  <-  list(rank1=u1, rank2=u2, VCtest=VCtest, max1=max1, fam=fam.sel, p.val=p.val, par=est$par, par2=est$par2, tau=Tau, low=TDep$lower, upp=TDep$upper)
    #print (reportCopulaRow (cop) )
    return (cop)
}


VC_Test_Matrix <- function(u1, u2, fams=c(1:10), correction="Schwarz")
{

    if (nargs() < 3) { stop("usage: VC_Test_Matrix(u1, u2, fams=c(...), correction=...), need at least the first 3 arguments \n")}
    N   <- length(fams)
    VC  <- matrix(0, nrow=(N+1), ncol=(N+1))
    for (i in 1:(N-1))
    {
        f1 <- fams[i]
        for (j in (i+1):N)
        {
            f2 <- fams[j]
            vc <- VineCopula::BiCopVuongClarke(u1, u2, familyset=c(f1,f2), correction=correction, rotations=FALSE)
            if (vc[1,1]==1) {VC[i,j] <- f1} else if (vc[1,2]==1) {VC[i,j] <- f2}
            if (vc[2,1]==1) {VC[j,i] <- f1} else if (vc[2,2]==1) {VC[j,i] <- f2}
            VC[i, (N+1)] <- VC[i, (N+1)] + vc[1, 1]
            VC[j, (N+1)] <- VC[j, (N+1)] + vc[1, 2]
            VC[(N+1), i] <- VC[(N+1), i] + vc[2, 1]
            VC[(N+1), j] <- VC[(N+1), j] + vc[2, 2]
        }
    }
    diag(VC)     <- NA
    rownames(VC) <- c(fams, "Clarke")
    colnames(VC) <- c(fams, "Vuong")
    return (VC)
}

VC_Test <- function(p, u1, u2, fams=c(1:10), correction="Schwarz")
{
    Len  <- min(length(u1), length(u2))
    u1   <- u1[1:Len]
    u2   <- u2[1:Len]
    tcol <- data.frame( pair=rep(p, 2), correction=rep(correction, 2), test=c("Vuong", "Clarke") )
    vc   <- VineCopula::BiCopVuongClarke(u1, u2, familyset=fams, correction=correction, rotations=FALSE)
    vc   <- cbind(tcol, vc)
    rownames(vc) <- NULL
    return(vc)
}

VC_Test3 <- function(p, u1, u2, fams=c(1:10))
{
    vc1 <- VC_Test (p, u1, u2, fams=fams, correction=FALSE)
    vc2 <- VC_Test (p, u1, u2, fams=fams, correction="Akaike")
    vc3 <- VC_Test (p, u1, u2, fams=fams, correction="Schwarz")
    vc  <- rbind(vc1, vc2, vc3)
    return(vc)
}

VCsimCopula <- function (cop, N)
{
    x     <- rnorm(10)
    a     <- VC_Test_Matrix (calcRank(x), calcRank(1*x+rnorm(10)), fams=fams, correction="Akaike")
    L     <- length(fams) + 1
    a[1:L, 1:L] <- 0
    k1    <- 0 # ak clarke
    k2    <- 0 # ak vuong
    k3    <- 0 # sz clarke
    k4    <- 0 # sz vuong
    fam   <- cop[1]
    par   <- cop[3]
    par2  <- cop[4]
    for (i in 1:N)
    {
        u  <- BiCopSim(124, fam, par, par2)
        u1 <- u[,1]
        u2 <- u[,2]
        a1 <- VC_Test_Matrix (u1, u2, fams=fams, correction="Akaike")
        a2 <- VC_Test_Matrix (u1, u2, fams=fams, correction="Schwarz")
        a  <- a + 0.5 * (a1 + a2)
        m  <- which.max(a1[11,1:10])
        if (m != fam) { k1 <- k1 + 1 }
        m  <- which.max(a1[1:10,11])
        if (m != fam) { k2 <- k2 + 1 }
        m  <- which.max(a2[11,1:10])
        if (m != fam) { k3 <- k3 + 1 }
        m  <- which.max(a2[1:10,11])
        if (m != fam) { k4 <- k4 + 1 }
    }
    a <- round(a/N, 0)
    r <- list(fam=fam, par=par, VC=a, AC=k1, AV=k2, SC=k3, SV=k4)
    return (r)
}

descStat5 <- function (x, k=5)
{
    if (is.list(x) == FALSE) { return (FALSE) }
    n     <- length(x)
    i1    <- seq(1, n*k, by=k)
    i2    <- seq(2, n*k, by=k)
    i3    <- seq(3, n*k, by=k)
    i4    <- seq(4, n*k, by=k)
    i5    <- seq(5, n*k, by=k)
    x     <- unlist(x)
    par   <- as.numeric(x[i1])
    par2  <- as.numeric(x[i2])
    Tau   <- as.numeric(x[i3])
    low   <- as.numeric(x[i4])
    upp   <- as.numeric(x[i5])
    y1    <- c( mean(par), mean(par2), mean(Tau), mean(low), mean(upp) )
    y2    <- c ( sd(par), sd(par2), sd(Tau), sd(low), sd(upp) )
    y3    <- c( se(par), se(par2), se(Tau), se(low), se(upp) )
    y2    <- lapply(y2, FUN=value2NA) %>% unlist
    y3    <- lapply(y3, FUN=value2NA) %>% unlist
    y     <- as.data.frame(rbind(y1, y2, y3))
    colnames(y) <- c("par", "par2", "Tau", "low", "upp")
    rownames(y) <- c("mean", "sd", "se")
    return ( round(y, 4) )
}

# LB test
LB.Test <- function (x, L=12)
{
    v1       <- Box.test(x, lag=L, type = "Box-Pierce")$p.v
    v2       <- Box.test(x, lag=L, type = "Ljung-Box")$p.v
    r        <- c(v1, v2)
    names(r) <- c("Box-Pierce", "Ljung-Box")
    r        <- round(r, 4)
    return (r)
}
# ARCH LM
AR.Test <- function (x, L=12)
{
  r        <- ArchTest(x, lags=L)$p.v[[1]]
  names(r) <- "ARCH-LM"
  return (r)
}

# Kolmogorov Smirnov
KS.Test <- function (x)
{
  r        <- lillie.test(x)$p.v
  names(r) <- "KS"
  return (r)
}

# Cramer von Mises
CM.Test <- function (x)
{
  r        <- cvm.test(x)$p.v
  names(r) <- "CM"
  return (r)
}

descStats <- function (x, Name=NA)
{
    x <- na.omit(x)
    r <- rep(NA, 6)
    r[1]  <- mean(x)
    r[2]  <- sd(x)
    r[3]  <- min(x)
    r[4]  <- max(x)
    r[5]  <- skewness(x)
    r[6]  <- kurtosis(x)
    r[7]  <- KS.Test(x)
    r[8]  <- CM.Test(x)
    r[9]  <- LB.Test(x)
    r[10] <- AR.Test(x)
    r     <- round(r, 4)
    r     <- as.data.frame(x=r)
    if (is.na(Name)==FALSE) {colnames(r) <- Name}
    rownames(r) <- c("Mean", "Std.Dev.", "Min", "Max", "Skewness", "Kurtosis", "KS", "CvM", "Q(12)", "ARCH-LM")
    return (r)
}

reportCopulaRow <- function (x=NA, d=3)
{
    if (is.list(x) == TRUE)
    {
        fam  <- as.integer(x$fam)
        if (is.na(x$p.val)==FALSE) { pval <- round (x$p.val, d) } else {pval<-NA}
        par  <- round (x$par, d)
        par2 <- round (x$par2,d)
        Tau  <- round (x$tau, d)
        TLow <- round (x$taildep$lower, d)
        TUpp <- round (x$taildep$upper, d)
        logL <- round (x$logLik, d)
        AIC  <- round (x$AIC, d)
        BIC  <- round (x$BIC, d)
        y    <- c(fam, par, par2, Tau, TLow, TUpp, pval, logL, AIC, BIC)
        names(y) <- c("fam", "par", "par2", "Tau", "TLow", "TUpp", "pval", "logL", "AIC", "BIC")
        return ( y )
    }
}

reportCopulaEstimate <- function (x=NA, d=3)
{
    if (is.list(x) == TRUE)
    {
        fname  <- BiCopName(x$fam)
        if (is.na(x$p.val)==FALSE) { pval <- round (x$p.val, d) } else {pval <- NA}
        par  <- round (x$par, d)
        par2 <- round (x$par2,d)
        Tau  <- round (x$tau, d)
        TLow <- round (x$low, d)
        TUpp <- round (x$upp, d)
        logL <- round (x$logLik, d)
        AIC  <- round (x$AIC, d)
        BIC  <- round (x$BIC, d)
        y    <- c(     fname,      par,     par2,      Tau,   TLow,       TUpp,       pval,   logL,  AIC,   BIC)
        names(y) <- c("family", "theta", "theta_2", "tau", "lambda_L", "lambda_U", "pval", "logL", "AIC", "BIC")
        return ( y )
    }
}

value2NA <- function (val1, val2=0)
{
    if (val1 == val2) {r <- NA} else { r <- val1 }
    return (r)
}

NCopBoot <- function (family, u1, u2)
{
    T    <- length(u1)
    ir   <- sample(1:T, T, replace=T)
    s1   <- u1[ir]
    s2   <- u2[ir]
    pp   <- VineCopula::BiCopEst(s1, s2, family)
    par  <- pp$par
    par2 <- pp$par2
    Tau  <- VineCopula::BiCopPar2Tau(family, par, par2)
    lu   <- VineCopula::BiCopPar2TailDep(family, par, par2)
    low  <- lu[[1]]
    upp  <- lu[[2]]
    res  <- c(par, par2, Tau, low, upp)
    if (is.na(par)==TRUE) { res <-  NCopBoot(family, u1, u2) }
    return (res)
}

NonBoot <- function(n)
{
    cop <- as.list(cop)
    return ( NCopBoot(cop$fam, cop$rank1, cop$rank2) )
}

do_biq <- function (y1, y2, m=NA)
{
    source("Busetti_Harvey.R")
    T      <- min(length(na.omit(y1)), length(na.omit(y2)))
    m      <- floor( 4 * (T/100)^0.25 )
    biq    <- rep(NA, 3)
    biq[1] <- max(compute_biq(y1[1:T], y2[1:T], tau=0.25, m=m)$bi)
    biq[2] <- max(compute_biq(y1[1:T], y2[1:T], tau=0.50, m=m)$bi)
    biq[3] <- max(compute_biq(y1[1:T], y2[1:T], tau=0.75, m=m)$bi)
    biq    <- round(biq, 4)
    names(biq)  <- c ("0.25", "0.50", "0.75")
    return (biq)
}

dynamicTau <- function(x, y, N=6)
{
    T <- length(x)
    #cat("T=", T, "\n")
    c <- rep(NA, N)
    for (i in 1:N)
    {
        K <- T - N + i
        c[i] <- cor(x[1:K], y[1:K], method="kendall")
        #cat("i=", i, "\t", "K=", K, "c=", c[i], "\n")
    }
    i  <- 1
    K1 <- 1
    K2 <- T - N
    while (K2 < T)
    {
        c[i] <- cor(x[K1:K2], y[K1:K2], method="kendall")
        K1   <- K1 + 1
        K2   <- K2 + 1
        i    <- i  + 1
    }
    return (c)
}

bootTau <- function(x, y, N=100)
{
    T <- length(x)
    c <- rep(NA, N)
    for (i in 1:N)
    {
        k <- sample(1:T, replace=T)
        x1 <- x[k]
        y1 <- y[k]
        c[i] <- cor(x1, y1, method="kendall")
    }
    return (c)
}


residAnalStat <- function (pre, suf, L=12)
{
  N      <- length(pre) * length(suf)
  rnames <- rep(NA, N)
  DF     <- matrix(NA, nrow=N, ncol=5)
  DF     <- as.data.frame(DF)
  i      <- 1
  for (p in pre)
  {
    for (s in suf)
    {
      x         <- eval(parse(text=paste0(p, ".resid.", s)))
      x         <- resid_std(x)
      v         <- c(KS.Test(x), CM.Test(x), AR.Test(x, L=L), LB.Test(x, L=L))
      v         <- t(as.data.frame(v))
      v         <- round(v, 4)
      rnames[i] <- paste0(p, ".", s)
      DF[i,]    <- v
      i         <- i + 1
    }
  }
  colnames(DF) <- colnames(v)
  rownames(DF) <- rnames
  return (DF)
}


varResidTest <- function (C.var, L=12)
{
  x  <- arch.test(C.var, lags.single=L, lags.multi=L, multivariate.only=F)
  v1 <- x$arch.uni[[1]]$p.value[[1]]
  v2 <- x$arch.uni[[2]]$p.value[[1]]
  v3 <- x$arch.mul$p.value[[1]]
  v4 <- serial.test(C.var, type="PT.asymptotic", lags.pt=L)$serial$p.value[[1]]
  v5 <- serial.test(C.var, type="PT.adjusted", lags.pt=L)$serial$p.value[[1]]
  resid <- residuals(C.var)
  res1  <- resid[,1]
  res2  <- resid[,2]
  v6    <- LB.Test(res1)
  v7    <- LB.Test(res2)
  r     <- c(v1, v2, v3, v4, v5, v6, v7)
  names(r) <- c("uni1", "uni2", "mul", "ser-as", "ser-ad", "BP1", "LB1", "BP2", "LB2")
  r     <- round(r, 4)
  return (r)
}


arma21Garch11 <- function (x, pq=c(2,1), k=2, vmod="sGARCH", distr="ged")
{
  # N = maximum arma order
  # k = selection criterion 1 Akaike, 2 Bayes, 3 Shibata, 4 Hannan-Quinn
  x <- as.numeric(na.omit(x))
  arma.order <- c(1, 1)
  crit       <- Inf
  go         <- list(c(1,0), c(0,1), c(1,1), c(1,2), c(2,1), c(2,2))
  go         <- list(c(1,1))
  distr.mod  <- distr   ### c("ged") # "norm", "snorm", "std", "sstd", "ged", "sged"
  if(length(pq)==2) {
    ar <- pq[1]
    ma <- pq[2]
  } else {
    ar <- 0:pq[1]
    ma <- 0:pq[1]
  }
  for (i in ar) # AR terms
  {
    for (j in ma) # MA terms
    {
      for (g in 1:length(go))
      {
        for (d in 1:length(distr.mod)) {
          model <- ugarchspec ( mean.model=list(armaOrder=c(i,j), include.mean=TRUE), # , archm=TRUE
                                variance.model=list(model=vmod, garchOrder=go[[g]]), # , submodel="GARCH"
                                distribution.model=distr.mod[d] )
          #Nt <- 0
          #repeat {
          ag.fit <- ugarchfit (x, spec=model, rseed=1, solver="hybrid", solver.options=list(maxeval=200000, xtol_rel=1e-6),
                               numderiv.control = list(grad.eps=1e-6, grad.d=1e-6, hess.eps=1e-6) )
          arma.crit <- infocriteria(ag.fit)[k] # 1 Akaike, 2 Bayes, 3 Shibata, 4 Hannan-Quinn
          #Nt <- Nt + 1
          #if (convergence(ag.fit) == 0 || Nt >= 2) {break}
          #}
          #cat (i, j, g, d, arma.crit, "\n")
          if (arma.crit < crit ) {
            crit        <- arma.crit
            arma.order  <- c(i, j)
            garch.order <- go[[g]]
            distr       <- distr.mod[d]
            best.fit    <- ag.fit
          }
        }
      }
    }
  }
  resid     <- best.fit@fit$residuals
  agmod     <- best.fit
  r         <- list(resid=resid, agmod=agmod)
  return (agmod)
}

