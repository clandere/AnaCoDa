#' Plot Parameter 
#' 
#' @param x A parameter object
#' 
#' @param what Which aspect of the parameter to plot. Default value is
#' "Mutation".
#' 
#' @param samples Number of samples to plot using the posterior mean. Default
#' value is 100.
#'
#' @param mixture.name a vector with names/descriptions of the mixture distributions in the parameter object
#'
#' @param with.ci Plot with or without confidence intervals. Default value
#' is TRUE
#' 
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' 
#' @return This function has no return value.
#' 
#' @description \code{plot} graphs the mutation or selection parameter for a ROC or FONSE
#' parameter object for each mixture element.
#' 
#' @details Graphs are based off the last # samples for the posterior mean.
#' 
plot.Rcpp_ROCParameter <- function(x, what = "Mutation", samples = 100, mixture.name = NULL, with.ci = TRUE, ...)
{
  plotParameterObject(x, what = what, samples= samples, mixture.name=mixture.name, with.ci=with.ci, ...)
}


#' Plot Parameter 
#' 
#' @param x A parameter object
#' 
#' @param what Which aspect of the parameter to plot. Default value is
#' "Mutation".
#' 
#' @param samples Number of samples to plot using the posterior mean. Default
#' value is 100.
#'
#' @param mixture.name a vector with names/descriptions of the mixture distributions in the parameter object
#'
#' @param with.ci Plot with or without confidence intervals. Default value
#' is TRUE
#' 
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' 
#' @return This function has no return value.
#' 
#' @description \code{plot} graphs the mutation or selection parameter for a ROC or FONSE
#' parameter object for each mixture element.
#' 
#' @details Graphs are based off the last # samples for the posterior mean.
#' 
plot.Rcpp_FONSEParameter <- function(x, what = "Mutation", samples = 100, mixture.name = NULL, with.ci = TRUE, ...)
{
  plotParameterObject(x, what = what, samples=samples,mixture.name = mixture.name, with.ci=with.ci, ...)
}

### NOT EXPOSED
plotParameterObject <- function(x, what = "Mutation", samples = 100, mixture.name = NULL, with.ci = TRUE, ...){
  
  numMixtures <- x$numMixtures
  means <- data.frame(matrix(0,ncol=numMixtures,nrow=40))
  sd.values <- data.frame(matrix(0,ncol=numMixtures*2,nrow=40))
  names.aa <- aminoAcids()
  paramType <- ifelse(what == "Mutation", 0, 1)
  #cat("ParamType: ", paramType, "\n")
  
  for (mixture in 1:numMixtures) {
    # get codon specific parameter
    count <- 1
    for (aa in names.aa) {
      if (aa == "M" || aa == "W" || aa == "X") next
      codons <- AAToCodon(aa, T)
      for (i in 1:length(codons)){
       means[count,mixture] <- x$getCodonSpecificPosteriorMean(mixture, samples, codons[i], paramType, TRUE)
        tmp <- x$getCodonSpecificQuantile(mixture, samples, codons[i], paramType, c(0.025, 0.975), TRUE)
        
        ## This approach to storing the quantiles may seem unconventional, but I actually found it to be the most straight forward approach
        ## for plotting later.
        sd.values[count,mixture] <- tmp[1]
        sd.values[count,mixture+numMixtures] <- tmp[2]
        count <- count + 1
      }
    }
  }
  ## Begin graphing
  mat <- matrix(rep(0,numMixtures*numMixtures),
                nrow = numMixtures, ncol = numMixtures, byrow = TRUE)
  count <- 1
  for(i in 1:numMixtures){
    for(j in 1:numMixtures){
      if(i<=j){
        mat[i,j] <-count
        count <- count + 1
      }
    }
  }
  nf <- layout(mat,widths=c(rep(5,numMixtures)),heights=c(rep(5,numMixtures)),respect=FALSE)
  par(mar=c(1,1,1,1))
  for(i in 1:numMixtures){
    for(j in 1:numMixtures){
      if(i==j)
      {
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="",xaxt='n',yaxt='n',ann=FALSE)
        if(is.null(mixture.name)){
          text(x = 0.5, y = 0.5, paste0("Mixture\nElement",i), 
               cex = 1.6, col = "black")
        }else{
          text(x = 0.5, y = 0.5, mixture.name[i], 
               cex = 1.6, col = "black")
        }
      } else if(i<j){
          if(with.ci){
            plot(means[,j],means[,i],ann=FALSE,xlim=range(cbind(sd.values[,j],sd.values[,j+numMixtures])),ylim=range(cbind(sd.values[,i],sd.values[,i+numMixtures])))
            upper.panel.plot(means[,j],means[,i],sd.x=cbind(sd.values[,j],sd.values[,j+numMixtures]),sd.y=cbind(sd.values[,i],sd.values[,i+numMixtures]))
          } else{
            plot(means[,j],means[,i],ann=FALSE,xlim=range(means[,j]),ylim=range(means[,i]))
            upper.panel.plot(means[,j],means[,i])
          }
        }
      }
    }
  }




#TODO: should PA's ploting be here as well?

upper.panel.plot <- function(x, y, sd.x=NULL, sd.y=NULL, ...){
  abline(0, 1, col = "blue", lty = 2)
  points(x, y, ...)
  if(!is.null(sd.y)){
    y.up <- sd.y[,2]
    y.low <- sd.y[,1]
    epsilon <- range(x, na.rm = T) * 0.1
    segments(x, y.low, x, y.up, ...)
  }
  if(!is.null(sd.x)){
    x.up <- sd.x[,2]
    x.low <- sd.x[,1]
    epsilon <- range(y, na.rm = T) * 0.1
    segments(x.low, y, x.up, y, ...)
  }  
  
  lm.line <- lm(y~x, na.action = "na.exclude")
  abline(lm.line, col="blue", lwd = 2)
  
  R2 <- summary(lm.line)$r.squared
  
  b <- lm.line$coef[2]
  rho <- ifelse(b > 0, sqrt(R2), -sqrt(R2)) #make sure rho has correct sign
  

  if(!is.null(sd.x))
  {
    xlim <- range(sd.x, na.rm = T)
  }else{
    xlim <- range(x,na.rm = T)
  }
  if(!is.null(sd.y))
  {
    ylim <- range(sd.y, na.rm = T)
  }else{
    ylim <- range(y,na.rm=T)
  }
  
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  
  std.error <- summary(lm.line)$coefficients[4]
  slope <- round(summary(lm.line)$coefficients[2], 3)
  intercept <- round(summary(lm.line)$coefficients[1], 3)
  t <- (slope - 1)/std.error
  
  if(t > qt(1-(0.05/2), lm.line$df.residual - 1)){
    eq <- paste0("y = ", sprintf("%.3f", intercept), " + ", sprintf("%.3f", slope), "x *")
    text(xlim[1] + width * 0.01, ylim[2] - height * 0.2, eq, pos = 4, cex = 1.5)
  }else{
    eq <- paste0("y = ", sprintf("%.3f", intercept), " + ", sprintf("%.3f", slope), "x")
    text(xlim[1] + width * 0.01, ylim[2] - height * 0.2, eq, pos = 4, cex = 1.5)
  } 
  if(b > 0){
    text(xlim[2] - width * 0.04, ylim[1] + height * 0.05,
         parse(text = paste0("rho == ", sprintf("%.4f", rho))),
         pos = 2, cex = 1.5, font = 2)
  }else{
    text(xlim[2] - width * 0.04, ylim[2] - height * 0.05,
         parse(text = paste0("rho == ", sprintf("%.4f", rho))),
         pos = 2, cex = 1.5, font = 2)
  }
}


lower.panel.plot <- function(x, y, ...)
{
  
}


confidenceInterval.plot <- function(x, y, sd.x=NULL, sd.y=NULL, ...){
  points(x, y, ...)
  if(!is.null(sd.y)){
    y.up <- sd.y[,2]
    y.low <- sd.y[,1]
    epsilon <- range(x, na.rm = T) * 0.1
    segments(x, y.low, x, y.up, ...)
  }
  if(!is.null(sd.x)){
    x.up <- sd.x[,2]
    x.low <- sd.x[,1]
    epsilon <- range(y, na.rm = T) * 0.1
    segments(x.low, y, x.up, y, ...)
  }  
  
  lm.line <- lm(y~x, na.action = "na.exclude") 
  
  b <- lm.line$coef[2]
  
  xlim <- range(x, na.rm = T)
  ylim <- range(y, na.rm = T)
  
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  
  std.error <- summary(lm.line)$coefficients[4]
  slope <- round(summary(lm.line)$coefficients[2], 3)
  intercept <- round(summary(lm.line)$coefficients[1], 3)
  t <- (slope - 1)/std.error
}


#' Plots ACF for codon specific parameter traces
#' 
#' @param parameter object of class Parameter
#' @param csp "Selection" or "Mutation", defaults to "Mutation"
#' @param numMixtures indicates the number of CSP mixtures used
#' @param samples number of samples at the end of the trace used to calculate the acf
#' @param lag.max Maximum amount of lag to calculate acf. Default is 10*log10(N), where N i the number of observations.
#' @param plot logical. If TRUE (default) a plot of the acf is created
#' 
#' @description The function calculates and by defaults plots the acf and estimates the autocorrelation in the trace.
#' 
#' @seealso \code{\link{acfMCMC}}

acfCSP <- function(parameter, csp = "Mutation", numMixtures = 1, samples = NULL, lag.max = 40, plot=TRUE)
{
  paramType <- 0
  if (csp == "Selection" )
  {
    paramType <- 1
  }
  if(is.null(samples)){ samples <- round(10*log10(length(trace))) }
  
  
  names.aa <- aminoAcids()
  trace <- parameter$getTraceObject()
  for (aa in names.aa)
  {
    if (aa == "M" || aa == "W" || aa == "X")
      next
    codons <- AAToCodon(aa, TRUE)
    for (i in 1:length(codons))
    {
      for (j in 1:numMixtures)
      {
        csp.trace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(j, codons[i], paramType, TRUE)
        csp.trace <- csp.trace[(length(csp.trace)-samples):length(csp.trace)]
        csp.acf <- acf(x = csp.trace, lag.max = lag.max, plot = FALSE)
        header <- paste(csp, aa, codons[i], "Mixture:", j, sep = " ")
        plot(x = csp.acf, xlab = "Lag time", ylab = "Autocorrelation", main = header)
      }
    }
  }
}

