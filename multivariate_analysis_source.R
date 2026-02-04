
PCA <- function(input_matrix, scale=FALSE, center=FALSE){
  input_rownames <- rownames(input_matrix)
  INPUT_NROW     <- nrow(input_matrix) # sample size
  INPUT_NCOL     <- ncol(input_matrix) # No of variables
  if (is.null(input_rownames)) input_rownames <- paste("#", 1:INPUT_NROW, sep="")  
  input_matrix   <- subset(input_matrix, complete.cases(input_matrix)) # remove cases with NAs

  if (is.null(colnames(input_matrix))) {
    colnames(input_matrix) <- paste("X", 1:INPUT_NCOL, sep="")
  }

  variable_names <- colnames(input_matrix)
  mean_values    <- colMeans(input_matrix)
  variances      <- apply(input_matrix, 2, var)
  SDs            <- sqrt(variances)
  r              <- cor(input_matrix)

  PCA_result     <- prcomp(input_matrix, scale=scale, center=center)
  eigenvalues    <- sqrt(PCA_result$sdev)
  eigenvector    <- PCA_result$rotation # this is rotation matrix
  contribution   <- eigenvalues/INPUT_NCOL*100
  cumulative_contribution <- cumsum(contribution)
  PCA_loadings <- t(sqrt(eigenvalues)*t(eigenvector))
  PCA_scores   <- 
    scale(input_matrix) %*% eigenvector * sqrt(INPUT_NROW/(INPUT_NROW-1))
  names(mean_values) <- 
    names(variances) <- names(SDs)  <-
    rownames(r)      <- colnames(r) <- rownames(PCA_loadings) <- 
    colnames(input_matrix)
  names(eigenvalues)        <- 
    names(contribution)    <- names(cumulative_contribution) <-
    colnames(PCA_loadings) <- colnames(PCA_scores) <- 
    paste("PC", 1:INPUT_NCOL, sep="")

  return(structure(
    list(
      mean=mean_values, variance=variances,
      standard.deviation=SDs, r=r,
      factor.loadings=PCA_loadings, eval=eigenvalues, 
      evec=eigenvector, nr=INPUT_NROW,  # added for subsequent PCA projection
      contribution=contribution,
      cum.contribution=cumulative_contribution, fs=PCA_scores), 
    class="pca"))
}

# print メソッド
print.PCA <- function(
  pca_returned_object,
  N_PCs_to_visualize=NULL,
  digits_to_visualize=3)
{
  eval <- pca_returned_object$eval
  nv <- length(eval)
  
  if (is.null(N_PCs_to_visualize)) {
    N_PCs_to_visualize <- sum(eval >= 1)
  }

  eval <- eval[1:N_PCs_to_visualize]
  cont <- eval/nv
  cumc <- cumsum(cont)
  fl <- pca_returned_object$factor.loadings[, 1:N_PCs_to_visualize, drop=FALSE]
  rcum <- rowSums(fl^2)
  vname <- rownames(fl)
  max.char <- max(nchar(vname), 12)
  fmt1 <- sprintf("%%%is", max.char)
  fmt2 <- sprintf("%%%is", digits_to_visualize+5)
  fmt3 <- sprintf("%%%i.%if", digits_to_visualize+5, digits_to_visualize)
  cat("\nResult of PCA\n\n")
  cat(sprintf(fmt1, ""),
      sprintf(fmt2, c(sprintf("PC%i", 1:N_PCs_to_visualize), "  Contribution")), "\n", sep="", collapse="")
  
  for (i in 1:nv) {
    cat(sprintf(fmt1, vname[i]),
	sprintf(fmt3, c(fl[i, 1:N_PCs_to_visualize], rcum[i])),
        "\n", sep="", collapse="")
  }

  cat(sprintf(fmt1, "Eigenvalue"),   sprintf(fmt3, eval[1:npca]), "\n", sep="", collapse="")
  cat(sprintf(fmt1, "Contribution"), sprintf(fmt3, cont[1:npca]), "\n", sep="", collapse="")
  cat(sprintf(fmt1, "Cum.contrib."), sprintf(fmt3, cumc[1:npca]), "\n", sep="", collapse="")
	
}

# summary メソッド
summary.PCA <- function(pca_returned_object, output_digits=5){
  print.default(pca_returned_object, digits=output_digits)
}

# plot メソッド
plot.PCA <- function(
  pca_returned_object,
  which=c("loadings", "scores"),
  pc.no=c(1,2),
  draw_axes=TRUE,			        # 座標軸を描き込むかどうか
  label.cex=0.6,			# 主成分負荷量のプロットのラベルのフォントサイズ
  markers=NULL,
  draw.names=FALSE,
  col="black",
  ...)				# plot に引き渡す引数
{
  which <- match.arg(which)

  if (which == "loadings") {
    d <- pca_returned_object$factor.loadings
  } else {
    d <- pca_returned_object$fs
  }
    
  label <- sprintf("PC%i", pc.no)
  plot(d[, pc.no[1]], d[, pc.no[2]], 
       xlab=label[1], ylab=label[2], col=col, ...)
    
  if (which == "loadings") {
    labelPosition <- ifelse(d[, pc.no[1]] < 0, 4, 2)
        
     if (is.null(markers) == FALSE){
     
       for (marker in markers){
         points(x=d[marker, pc.no[1]], y=d[marker, pc.no[2]],
		col="magenta", pch=16)
       }
     }

     if (draw.names == TRUE){
         text(x=d[, pc.no[1]], y=d[, pc.no[2]], 
	      labels=rownames(pca_returned_object$factor.loadings), 
              pos=labelPosition, cex=label.cex
         )
     }
  } else {
    if (which == "scores" && is.null(markers) == FALSE){
      for (marker in markers){
        points(x=d[marker, pc.no[1]], y=d[marker, pc.no[2]], col="magenta", pch=16)
        #text(x=d[marker, pc.no[1]], y=d[marker, pc.no[2]], 
        #     labels=marker, col="red")# At2g14610
      }
    }
  }
}

