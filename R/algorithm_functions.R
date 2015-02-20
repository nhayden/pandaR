#' Passing Messages between Biological Networks to Refine Predicted Interactions
#'
#' This function runs the PANDA algorithm 
#' 
#' @param motif A motif dataset, a data.frame, matrix or exprSet containing 3 columns. Each row describes an motif associated with a transcription factor (column 1) a gene (column 2) and a score (column 3) for the motif.
#' @param expr An expression dataset, as a genes (rows) by samples (columns) data.frame
#' @param ppi A Protein-Protein interaction dataset, a data.frame containing 3 columns. Each row describes a protein-protein interaction between transcription factor 1(column 1), transcription factor 2 (column 2) and a score (column 3) for the interaction.
#' @param update.method function for iterative updating.  This must be one of "tanimoto" and "danmethod"(in development)
#' @param alpha value to be used for update variable, alpha (default=0.1)
#' @param hamming value at which to terminate the process based on hamming distance (default 10^-5)
#' @param k sets the maximum number of iterations PANDA can run before exiting.
#' @param progress Boolean to indicate printing of output for algorithm progress.
#' @param output a vector containing which networks to return.  Options include "regulatory", "coregulatory", "cooperative".
#' @param z.scale Boolean to indicate use of z-scores in output.  False will use [0,1] scale.  
#' @param randomize method by which to randomize gene expression matrix.  Default "None".  Must be one of "None", "within.gene", "by.genes".  "within.gene" randomization scrambles each row of the gene expression matrix, "by.gene" scrambles gene labels.
#' @keywords keywords
#' @export
#' @return An object of class "panda" containing matrices describing networks achieved by convergence with PANDA algorithm.  
#' "reg.net" is the regulatory network
#' "coreg.net" is the coregulatory network
#' "coop.net" is the cooperative network
#' (Should wrap this up in a special class)
#' @examples
#' data(yeast)
#' panda.res.cc <- panda(yeast$motif,yeast$exp.cc,yeast$ppi,hamming=.001,progress=T)
#' panda.res.sr <- panda(yeast$motif,yeast$exp.sr,yeast$ppi,hamming=.001)
#' @references
#' Glass K, Huttenhower C, Quackenbush J, Yuan GC. Passing Messages Between Biological Networks to Refine Predicted Interactions. PLoS One. 2013 May 31;8(5):e64832. 
panda <- function( motif, 
                   expr=NULL, 
                   ppi=NULL, 
                   update.method='tanimoto', 
                   alpha=0.1, 
                   hamming=0.00001, 
                   k=NA, 
                   output=c('regulatory','coexpression','cooperative'), 
                   z.scale=T,
                   progress=FALSE,
                   randomize="None"){
  if(progress)
    print('Initializing and validating')
  expr.data  <- expr
  motif.data <- motif
  ppi.data   <- ppi
  
  if(class(expr)=="ExpressionSet")
    expr.data <- expr@assayData
  
  # Create vectors for TF names and Gene names from Motif dataset
  tf.names   <- sort(unique(motif.data[,1]))
  num.TFs    <- length(tf.names)
  if (is.null(expr.data)){
    # Use only the motif data here for the gene list
    gene.names <- sort(unique(motif.data[,2]))
    num.genes  <- length(gene.names)
    num.conditions <- 0
    if (randomize!="None"){
      warning("Randomization ignored because gene expression is not used.")
      randomize <- "None"
    }
  } else {
    # Use the motif data AND the expr data (if provided) for the gene list
    gene.names <- sort(intersect(motif.data[,2],rownames(expr.data)))
    num.genes  <- length(gene.names)
    
    # Filter out the expr genes without motif data
    expr.data <- expr.data[rownames(expr.data) %in% gene.names,]
    
    # Keep everything sorted alphabetically
    expr.data      <- expr.data[order(rownames(expr.data)),]
    num.conditions <- ncol(expr.data);
    if (randomize=='within.gene'){
      expr.data <- t(apply(expr.data, 1, sample))
      if(progress)
        print("Randomizing by reordering each gene's expression")
    } else if (randomize=='by.genes'){
      rownames(expr.data) <- sample(rownames(expr.data))
      expr.data           <- expr.data[order(rownames(expr.data)),]
      if(progress)
        print("Randomizing by reordering each gene labels")
    }
  }

  # Bad data checking
  if (num.genes==0){
    stop("Error validating data.  No matched genes.\n  Please ensure that gene names in expression file match gene names in motif file.")
  }
  
  if(num.conditions==0) {
    warning('No expression data given.  PANDA will run based on an identity co-regulation matrix')
    gene.coreg <- diag(num.genes)
  } else if(num.conditions<3) {
    warning('Not enough expression conditions detected to calculate correlation. Co-regulation network will be initialized to an identity matrix.');
    gene.coreg <- diag(num.genes)
  } else {
    gene.coreg <- cor(t(expr.data), method="pearson", use="pairwise.complete.obs");
    if(progress)
      print('Verified adequate samples')
  }
  
  

  # If no ppi data is given, we use the identity matrix
  if (is.null(ppi.data)){
    ppi.data <- diag(num.TFs)
  }
  
  # Convert 3 column format to matrix format
  regulatory.network <- spread.net(motif.data)

  # sort the genes (columns)
  regulatory.network <- as.matrix(regulatory.network[,order(colnames(regulatory.network))])

  # Filter out any motifs that are not in expr dataset (if given)
  if (!is.null(expr.data)){
    regulatory.network <- regulatory.network[,colnames(regulatory.network) %in% gene.names]
  }

  # store initial motif network (alphabetized for rows and columns)
  starting.motifs <- regulatory.network
  
  # ppi.data Data
  tf.coop.network=diag(num.TFs);
  Idx1=match(ppi.data[,1], tf.names);
  Idx2=match(ppi.data[,2], tf.names);
  Idx=(Idx2-1)*num.TFs+Idx1;
  tf.coop.network[Idx[!is.na(Idx)]]=as.numeric(ppi.data[,3])[!is.na(Idx)];
  Idx=(Idx1-1)*num.TFs+Idx2;
  tf.coop.network[Idx[!is.na(Idx)]]=as.numeric(ppi.data[,3])[!is.na(Idx)];
  colnames(tf.coop.network) <- tf.names
  rownames(tf.coop.network) <- tf.names
  
  ## Run PANDA ##
  tic=proc.time()[3];
  
  if(progress)
    print('Normalizing networks...');
  regulatory.network=normalize.network(regulatory.network,update.method);
  tf.coop.network=normalize.network(tf.coop.network,update.method);
  gene.coreg=normalize.network(gene.coreg,update.method);
  
  if(progress)
    print('Leaning Network!')
  step=0;
  hamming_cur=1;
  if (update.method=='danmethod'){
    if(progress)
      print("using danmethod similarity")
    #   while((hamming_cur>0.00001) && (step<200))
    #   {
    Responsibility=d.function(tf.coop.network, regulatory.network);
    Availability=d.function(t(regulatory.network), gene.coreg);
    hamming_cur=sum(abs(regulatory.network-0.5*(Responsibility+Availability)))/(num.TFs*num.genes);
    regulatory.network=(1-alpha)*regulatory.network+alpha*0.5*(Responsibility+Availability);
    
    ppi.data=d.function(t(regulatory.network), t(regulatory.network));
    #   	ppi.data=update.diagonal(ppi.data, num.TFs, alpha, step);
    tf.coop.network=(1-alpha)*tf.coop.network+alpha*ppi.data;
    
    CoReg2=cor(regulatory.network, regulatory.network);
    CoReg2=update.diagonal(CoReg2, num.genes, alpha, step);
    gene.coreg=(1-alpha)*gene.coreg+alpha*CoReg2;
    print(paste("Step #", step, ", hamming =", round(hamming_cur,5)), sep="");
    step=step+1;
    #   }
  } else if (update.method=='tanimoto'){
    if(progress)
      print("Using tanimoto similarity")
    while(hamming_cur>hamming)
    {
      if ((!is.na(k))&&step>=k){
        stop(paste("Reached maximum iterations, k =",k),sep="")
      }
      Responsibility=t.function(tf.coop.network, regulatory.network);
      Availability=t.function(regulatory.network, gene.coreg);
      hamming_cur=sum(abs(regulatory.network-0.5*(Responsibility+Availability)))/(num.TFs*num.genes);
      regulatory.network=(1-alpha)*regulatory.network+alpha*0.5*(Responsibility+Availability);
      
      ppi.data=t.function(regulatory.network, t(regulatory.network));
      ppi.data=update.diagonal(ppi.data, num.TFs, alpha, step);
      tf.coop.network=(1-alpha)*tf.coop.network+alpha*ppi.data;
      
      CoReg2=t.function(t(regulatory.network), regulatory.network);
      CoReg2=update.diagonal(CoReg2, num.genes, alpha, step);
      gene.coreg=(1-alpha)*gene.coreg+alpha*CoReg2;
      if(progress)
        print(paste("Step #", step, ", hamming =", round(hamming_cur,5)), sep="");
      step=step+1;
    }
  }
  toc=proc.time()[3] - tic;
  if(progress)
    print(paste("Running PANDA on ", num.genes, " Genes and ", num.TFs, " TFs took ", round(toc,2), " seconds!", sep=""));
  
  res.list <- list()
  if (!z.scale){
    regulatory.network <- pnorm(regulatory.network)
    gene.coreg         <- pnorm(gene.coreg)
    tf.coop.network    <- pnorm(tf.coop.network)
  }
  if("regulatory"%in%output){
    res.list$reg.net <- regulatory.network
  }
  if("coregulatory"%in%output){
    res.list$coreg.net <- gene.coreg
  }
  if("cooperative"%in%output){
    res.list$coop.net <- tf.coop.network
  }
  res <- panda.obj(reg.net=regulatory.network,
                   coreg.net=gene.coreg,
                   coop.net=tf.coop.network)
  return(res)  
}

normalize.network<-function(X,update.method)
{
  X <- as.matrix(X)
  if (update.method=='danmethod'){
    X[X<0]<-0
    return(X);
  }
  if (update.method=='tanimoto'){
    
    # overall values
    mu0=mean(X);
    std0=sd(X);
    
    # operations on rows
    mu1=apply(X,1,mean); # operations on rows
    std1=apply(X,1,sd)*sqrt((dim(X)[2]-1)/dim(X)[2]);
    mu1=matrix(rep(mu1, dim(X)[2]), dim(X));
    std1=matrix(rep(std1, dim(X)[2]), dim(X));
    Z1=(X-mu1)/std1;
    
    # operations on columns
    mu2=apply(X,2,mean); # operations on columns
    std2=apply(X,2,sd)*sqrt((dim(X)[1]-1)/dim(X)[1]);
    mu2=matrix(rep(mu2, each=dim(X)[1]), dim(X));
    std2=matrix(rep(std2, each=dim(X)[1]), dim(X));
    Z2=(X-mu2)/std2;
    
    # combine and return
    normMat=Z1/sqrt(2)+Z2/sqrt(2);
    
    # Dan fix to NaN
    normMat[is.na(normMat)]<-0
    return(normMat);
    
  }
}

t.function<-function(X,Y)
{
  Amat=(X %*% Y);
  Bmat=apply(Y*Y,2,sum);
  Bmat=matrix(rep(Bmat, each=dim(X)[1]), dim(Amat));
  Cmat=apply(X*X,1,sum);
  Cmat=matrix(rep(Cmat, dim(Y)[2]), dim(Amat));
  
  Amat=Amat/sqrt(Bmat+Cmat-abs(Amat));
  
  return(Amat);
}

d.function<-function(X,Y)
{
  A <- cor(X,Y)
  A[A<0]<-0
  return(A);
}


update.diagonal<-function(diagMat, num, alpha, step)
{
  diagMat[seq(1, num*num, num+1)]=NaN;
  diagstd=apply(diagMat,2,sd,na.rm=T)*sqrt((num-2)/(num-1));
  diagMat[seq(1, num*num, num+1)]=diagstd*num*exp(2*alpha*step);
  return(diagMat);
}

getNullTansMatrix <- function(expr1, expr2, motif, ppi){
  total.samples <- ncol(expr1) + ncol(expr2)
  combined.data <- cbind(expr1,expr2)
  null.indices.1 <- sample(1:total.samples, ncol(expr1))
  null.expr1 <- combined.data[,null.indices.1]
  null.expr2 <- combined.data[,-null.indices.1]
  reg.net.A <- panda(motif, null.expr1, ppi)$reg.net
  reg.net.B <- panda(motif, null.expr2, ppi)$reg.net
  transformation.matrix(reg.net.A, reg.net.B)
}

spread.net <- function(df){
  df[,3]<- as.numeric(df[,3])
  row_names <- unique(df[,1])
  col_names <- unique(df[,2])
  spread.df <-  data.frame(matrix(0,nrow=length(row_names),ncol=length(col_names)),row.names=row_names)
  colnames(spread.df) <- col_names
  for(i in 1:nrow(df)){
    spread.df[as.character(df[i,1]),as.character(df[i,2])] <- df[i,3]
  }
  spread.df
}

#' Top edges
#'
#' topedges gets a network from a panda obj with a specified cutoff based on magnitude of edgeweight.
#' 
#' @param x an object of class "panda"
#' @param count an optional integer indicating number of top edges to be included in regulatory network.
#' @param cutoff an optional numeric indicating the z-score edge weight cutoff to be used to identify edges. Default is 3.0.  Not used if count is not NA.
#' @param networks an optional vector specifying which networks to be included in output.  May be any combination of c("coregulation","cooperation","regulatory").
#' @keywords keywords
#' @export
#' @return An object of class "panda" containing binary matrices indicating the existence of an edge between two nodes.  For regulatory network the matrix indicates an edge between a transcription factor (row) and a gene (column)
#' @examples
#' data(yeast)
#' panda.res.cc <- panda(yeast$motif,yeast$exp.cc,yeast$ppi,hamming=.001,progress=T)
#' top.panda.res.cc <- topedges(panda.res.cc,1000)
topedges <- function(x, count=NA, cutoff=2.0, networks=c("coregulation","cooperation","regulatory")){
  if(class(x)!="panda"){
    warning(paste(sep="","Cannot run topedges on object of class '",class(x),"'.  Must be of class 'panda'"))
    stop
  }
  if (!is.na(count)){
    cutoff <- sort(x@reg.net)[length(x@reg.net)-(count-1)]
  }
  regulatory.network <- apply(x@reg.net>cutoff, 2,as.numeric)
  rownames(regulatory.network)<-rownames(x@reg.net)
  gene.coreg <- apply(x@coreg.net>cutoff, 2,as.numeric)
  rownames(gene.coreg)<-rownames(x@coreg.net)
  tf.coop.network <- apply(x@coop.net>cutoff, 2,as.numeric)
  rownames(tf.coop.network)<-rownames(x@coop.net)
  
  res <- panda.obj(reg.net=regulatory.network,
                   coreg.net=gene.coreg,
                   coop.net=tf.coop.network)
  return(res)
}

#' Subnetwork
#'
#' subnetwork gets a bipartite network containing only the transcription factors or genes and their respective connections
#' 
#' @param x an object of class "panda"
#' @param nodes character vector containing the transcription factor or gene labels to subset
#' @param sub.tf an optional logical indicating whether to subset by transcription factor.  Default is TRUE.
#' @keywords keywords
#' @export
#' @return An matrix describing the subsetted bipartite network.
#' @examples
#' data(yeast)
#' panda.res.cc <- panda(yeast$motif,yeast$exp.cc,yeast$ppi,hamming=.001,progress=T)
#' top.panda.res.cc <- topedges(panda.res.cc,1000)
#' subnet.panda.res.cc <- subnetwork(top.panda.res.cc,c("YKR099W","YDR423C","YDL056W","YLR182W"))
subnetwork <- function(x, nodes, sub.tf=T){
  if(class(x)!="panda"){
    warning(paste(sep="","Cannot run subnetwork on object of class '",class(x),"'.  Must be of class 'panda'"))
    stop
  }
  if (sub.tf){
    subnet <- x@reg.net[nodes,]
    edgeexists <- apply(subnet,2,sum)>0
    subnet <- subnet[,edgeexists]
  } else {
    subnet <- x@reg.net[,nodes]
    edgeexists <- apply(subnet,1,sum)>0
    subnet <- subnet[edgeexists,]    
  }
  return(subnet)
}


panda.obj <- setClass("panda", slots=c("reg.net","coreg.net","coop.net"))
summary.panda <- function(x){
  l <- list(coreg.net=dim(x@coreg.net),reg.net=dim(x@reg.net),coop.net=dim(x@coop.net))
  print(l)
}
plot.panda <- function(x){
  print("Plotting function has not been implemented for panda... yet.")
}
