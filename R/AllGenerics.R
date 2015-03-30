#' Summary.panda
#'
#' summarizes the results of a PANDA analysis
#'
#' @param object an object of class "panda"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Summary description of panda S4 object
#' @examples
#' data(pandaToyData)
#' panda.res <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001,progress=TRUE)
#' summary(panda.res)
summary.panda <- function(object, ...){
  l <- list(coregNet=dim(object@coregNet),regNet=dim(object@regNet),coopNet=dim(object@coopNet))
  cat("PANDA network for", nrow(object@coregNet),"genes and",nrow(object@coopNet)," transcription factors.")
}
#' print.panda
#'
#' summarizes the results of a PANDA analysis
#'
#' @param x an object of class "panda"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Summary description of panda S4 object
#' @examples
#' data(pandaToyData)
#' panda.res <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001,progress=TRUE)
#' print(panda.res)
print.panda <- function(x, ...){
  l <- list(coregNet=dim(x@coregNet),regNet=dim(x@regNet),coopNet=dim(x@coopNet))
  cat("PANDA network for", nrow(x@coregNet),"genes and",nrow(x@coopNet),"transcription factors.")
  cat("\n\nSlots:\n")
  cat(slotNames(x)[1],"\t: Regulatory network of",nrow(x@coopNet)," transcription factors to", nrow(x@coregNet),"genes.\n")
  cat(slotNames(x)[2],": Co-regulation network of", nrow(x@coregNet),"genes.\n")
  cat(slotNames(x)[3],"\t: Cooperative network of", nrow(x@coopNet),"transcription factors.\n\n")
  numEdges <- sum(x@regNet!=0)
  cat("Regulatory graph contains ",numEdges,"edges.\n")
  if (numEdges==nrow(x@regNet)*ncol(x@regNet)){
    cat("Regulatory graph is complete.")
  } else {
    cat("Regulatory graph is not complete.")
  }
}
#' Plot.panda
#'
#' summarizes the results of a PANDA analysis
#'
#' @param x an object of class "panda"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Plot of the distribution of edge weights in the regulatory network.
#' @examples
#' data(pandaToyData)
#' panda.res <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001,progress=TRUE)
#' plot(panda.res)
plot.panda <- function(x, ...){
  cat("PANDA network for", nrow(x@coregNet),"genes and",nrow(x@coopNet)," transcription factors.\n")
  cat("Mean edge weight = ", mean(x@regNet),"\n")
  cat("Min edge weight = ", min(x@regNet),"\n")
  cat("Max edge weight = ", max(x@regNet),"\n")
  hist(x@regNet, main="Distribution of edge weights")
}
