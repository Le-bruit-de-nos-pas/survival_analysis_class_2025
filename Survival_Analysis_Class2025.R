
library(survival)
library(survminer)
library(lubridate)
library(cmprsk)
library(tidyverse)
library(utile.visuals)
library(Hmisc)
library(corrplot)
library(ggcorrplot)
library(pROC)
library(caret)
library(mgcv)


"CumIncidence" <- function(ftime, fstatus, group, t, strata, rho = 0, 
                           cencode = 0, subset, na.action = na.omit, level,
                           xlab = "Time", ylab = "Probability", 
                           col, lty, lwd, digits = 4)
{
  # check for the required package
  if(!require("cmprsk"))
  { stop("Package `cmprsk' is required and must be installed.\n 
           See help(install.packages) or write the following command at prompt
           and then follow the instructions:\n
           > install.packages(\"cmprsk\")") } 
  
  mf  <- match.call(expand.dots = FALSE)
  mf[[1]] <- as.name("list")
  mf$t <- mf$digits <- mf$col <- mf$lty <- mf$lwd <- mf$level <- 
    mf$xlab <- mf$ylab <- NULL
  mf <- eval(mf, parent.frame())
  g <- max(1, length(unique(mf$group)))
  s <- length(unique(mf$fstatus))
  if(missing(t)) 
  { time <- pretty(c(0, max(mf$ftime)), 6)
  ttime <- time <- time[time < max(mf$ftime)] }
  else { ttime <- time <- t }
  
  fit   <- do.call("cuminc", mf)
  tfit <- timepoints(fit, time)
  
  cat("\n+", paste(rep("-", 67), collapse=""), "+", sep ="")
  cat("\n| Cumulative incidence function estimates from competing risks data |")
  cat("\n+", paste(rep("-", 67), collapse=""), "+\n", sep ="")
  tests <- NULL
  if(g > 1)
  { 
    tests <- data.frame(fit$Tests[,c(1,3,2)], check.names = FALSE)
    colnames(tests) <- c("Statistic", "df", "p-value")
    tests$`p-value` <- format.pval(tests$`p-value`)
    cat("Test equality across groups:\n")
    print(tests, digits = digits) 
  }
  cat("\nEstimates at time points:\n")
  print(tfit$est, digits = digits)
  cat("\nStandard errors:\n")
  print(sqrt(tfit$var), digits = digits)
  
  if(missing(level))
  { 
    if(missing(t))
    { time <- sort(unique(c(ftime, time)))
    x <- timepoints(fit, time) }
    else x <- tfit
    col <- if(missing(col)) rep(1:(s-1), rep(g,(s-1))) else col
    lty <- if(missing(lty)) rep(1:g, s-1) else lty
    lwd <- if(missing(lwd)) rep(1, g*(s-1)) else lwd      
    matplot(time, base::t(x$est), type="s", ylim = c(0,1), 
            xlab = xlab, ylab = ylab, xaxs="i", yaxs="i", 
            col = col, lty = lty, lwd = lwd)
    legend("topleft", legend =  rownames(x$est), x.intersp = 2, 
           bty = "n", xjust = 1, col = col, lty = lty, lwd = lwd)
    out <- list(test = tests, est = tfit$est, se = sqrt(tfit$var))
  }
  else
  { if(level < 0 | level > 1) 
    error("level must be a value in the range [0,1]")
    
    oldpar <- par(ask=TRUE)
    on.exit(par(oldpar))
    if(missing(t))
    { time <- sort(unique(c(ftime, time)))
    x <- timepoints(fit, time) }
    else x <- tfit
    z <- qnorm(1-(1-level)/2)
    lower <- x$est ^ exp(-z*sqrt(x$var)/(x$est*log(x$est)))
    upper <- x$est ^ exp(z*sqrt(x$var)/(x$est*log(x$est)))
    col <- if(missing(col)) rep(1:(s-1), rep(g,(s-1))) 
    else             rep(col, g*(s-1))
    lwd <- if(missing(lwd)) rep(1, g*(s-1)) 
    else             rep(lwd, g*(s-1))      
    
    for(j in 1:nrow(x$est))
    { matplot(time, cbind(x$est[j,], lower[j,], upper[j,]), type="s", 
              xlab = xlab, ylab = ylab, xaxs="i", yaxs="i", 
              ylim = c(0,1), col = col[j], lwd = lwd[j], lty = c(1,3,3))
      legend("topleft", legend =  rownames(x$est)[j], bty = "n", xjust = 1) }
    
    i <- match(ttime, time)
    ci <- array(NA, c(2, length(i), nrow(lower)))
    ci[1,,] <- base::t(lower[,i])
    ci[2,,] <- base::t(upper[,i])
    dimnames(ci) <- list(c("lower", "upper"), ttime, rownames(lower))
    cat(paste("\n", level*100, "% pointwise confidence intervals:\n\n", sep=""))
    print(ci, digits = digits)
    out <- list(test = tests, est = x$est, se = sqrt(tfit$var), ci = ci)
  }
  
  invisible(out)
}


BMT
attach(BMT)

dis = factor(dis, levels = c(0,1), labels = c("ALL", "AML"))

table(dis, status)

fit = CumIncidence(ftime, status, dis, cencode =0, xlab = "Months", level = 0.95, t = c(0:70))              



#     dis ftime status
# 1    0    13      2
# 2    0     1      1
# 3    0    72      0
# 4    0     7      2
# 5    0     8      2
# 6    1    67      0
# 7    0     9      2
# 8    0     5      2
# 9    1    70      0
# 10   1     4      0
# 11   1     7      0
# 12   1    68      0
# 13   0     1      2
# 14   1    10      2
# 15   1     7      2
# 16   1     3      1
# 17   1     4      1
# 18   1     4      1
# 19   1     3      1
# 20   1     3      1
# 21   0    22      2
# 22   1     8      1
# 23   1     2      2
# 24   0     0      2
# 25   0     0      1
# 26   0    35      0
# 27   1    35      0
# 28   0     4      2
# 29   0    14      2
# 30   0    26      2
# 31   0     3      2
# 32   1     2      0
# 33   1     8      0
# 34   1    32      0
# 35   0    12      1




if(!require(cmprsk))
{ stop("the package 'cmprsk' is required, please install it. \nSee help(install.packages).") }

factor2ind <- function(x, baseline)
{
  #### dummy variable encoding ####
  xname <- deparse(substitute(x))
  n <- length(x)
  x <- as.factor(x)
  if(!missing(baseline)) x <- relevel(x, baseline)
  X <- matrix(0, n, length(levels(x)))
  X[(1:n) + n*(unclass(x)-1)] <- 1
  X[is.na(x),] <- NA
  dimnames(X) <- list(names(x), paste(xname, levels(x), sep = ":"))
  return(X[,-1,drop=FALSE])
}

modsel.crr <- function (object, ..., d = log(object$n)) 
{
  if(class(object) != "crr") 
    stop("object is not of class 'crr'")
  objects <- list(object, ...)
  nmodels <- length(objects)
  modnames <- paste("Model ", format(1:nmodels), ": ", 
                    lapply(objects, function(x) x$call), 
                    sep = "", collapse = "\n")
  
  mod0 <- object
  mod0$loglik <- mod0$loglik.null
  mod0$coef <- mod0$call$cov1 <- mod0$call$cov2 <- NULL
  objects <- c(list(mod0), objects)
  nmodels <- nmodels + 1
  
  modnames <- c("Model 0: Null model", modnames)
  ns <- sapply(objects, function(x) x$n) 
  dfs <- sapply(objects, function(x) length(x$coef)) 
  if(any(ns != ns[1]))
    stop("models were not all fitted to the same dataset")
  out <- matrix(rep(NA, 5 * nmodels), ncol = 5)
  loglik <- sapply(objects, function(x) x$loglik)
  crit <- sapply(objects, function(x) -2*x$loglik + d*length(x$coef))
  out[,1] <- ns
  out[,2] <- loglik
  out[,3] <- dfs
  out[,4] <- crit
  out[,5] <- crit - min(crit)
  if(d==log(object$n)) critname <- "BIC"
  else if(d == 2) critname <- "AIC"
  else critname <- "Criterion"
  colnames(out) <- c("Num.obs", "logLik", "Df.fit", critname, paste(critname, "diff"))
  rownames(out) <- 0:(nmodels-1)
  title <- "Model selection table\n"
  topnote <- modnames
  structure(as.data.frame(out), heading = c(title, topnote), 
            class = c("anova", "data.frame"))
}


# 
#      Sex   D   Phase Age Status Source  ftime
# 1     M ALL Relapse  48      2  BM+PB   0.67
# 2     F AML     CR2  23      1  BM+PB   9.50
# 3     M ALL     CR3   7      0  BM+PB 131.77
# 4     F ALL     CR2  26      2  BM+PB  24.03
# 5     F ALL     CR2  36      2  BM+PB   1.47
# 6     M ALL Relapse  17      2  BM+PB   2.23
# 7     M ALL     CR1   7      0  BM+PB 124.83
# 8     F ALL     CR1  17      2  BM+PB  10.53
# 9     M ALL     CR1  26      0  BM+PB 123.90
# 10    F ALL Relapse   8      1  BM+PB   2.00


attach(bmtcrr)

x = cbind(Age, factor2ind(Sex, "M"), 
          factor2ind(D, "ALL"), 
          factor2ind(Phase, "Relapse"),
          factor2ind(Source))


mod1 = crr(ftime, Status, x)
summary(mod1)

