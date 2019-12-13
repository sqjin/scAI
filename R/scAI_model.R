#' The scAI Class
#'
#' The scAI object is created from a paired single-cell transcriptomic and epigenomic data.
#' It takes a list of two digital data matrices as input. Genes/loci should be in rows and cells in columns. rownames and colnames should be included.
#' The class provides functions for data preprocessing, integrative analysis, and visualization.
#'
#' The key slots used in the scAI object are described below.
#'
#' @slot raw.data List of raw data matrices, one per dataset (Genes/loci should be in rows and cells in columns)
#' @slot norm.data List of normalized matrices (genes/loci by cells)
#' @slot agg.data Aggregated epigenomic data within similar cells
#' @slot scale.data List of scaled matrices
#' @slot pData data frame storing the information associated with each cell
#' @slot var.features List of informative features to be used, one giving informative genes and the other giving informative loci
#' @slot fit List of inferred low-rank matrices, including W1, W2, H, Z, R
#' @slot fit.variedK List of inferred low-rank matrices when varying the rank K
#' @slot embed List of the reduced 2D coordinates, one per method, e.g., t-SNE/FIt-SNE/umap
#' @slot identity a factor defining the cell identity
#' @slot cluster List of consensus clustering results
#' @slot options List of parameters used throughout analysis
#'
#' @exportClass scAI
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass
#' @useDynLib scAI
scAI <- methods::setClass("scAI",
                          slots = c(raw.data = "list",
                                    norm.data = "list",
                                    agg.data = "matrix",
                                    scale.data = "list",
                                    pData = "data.frame",
                                    var.features = "list",
                                    fit = "list",
                                    fit.variedK = "list",
                                    embed = "list",
                                    identity = "factor",
                                    cluster = "list",
                                    options = "list")
)
#' show method for scAI
#'
#' @param scAI object
#' @param show show the object
#' @param object object
#' @docType methods
#'
setMethod(f = "show", signature = "scAI", definition = function(object) {
  cat("An object of class", class(object), "\n", length(object@raw.data), "datasets.\n")
  invisible(x = NULL)
})



#' creat a new scAI object
#'
#' @param raw.data List of raw data matrices, a paired single-cell transcriptomic and epigenomic data
#' @param do.sparse whether use sparse format
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom methods as new
create_scAIobject <- function(raw.data, do.sparse = T) {
  object <- methods::new(Class = "scAI",
                         raw.data = raw.data)
  if (do.sparse) {
    raw.data <- lapply(raw.data, function(x) {
      as(as.matrix(x), "dgCMatrix")
    })
  }
  object@raw.data <- raw.data
  return(object)
}



#' preprocess the raw.data including quality control and normalization
#'
#' @param object scAI object
#' @param assay List of assay names to be normalized
#' @param minFeatures Filter out cells with expressed features < minimum features
#' @param minCells Filter out features expressing in less than minCells
#' @param minCounts Filter out cells with expressed count < minCounts
#' @param maxCounts Filter out cells with expressed count > minCounts
#' @param libararyflag Whether do library size normalization
#' @param logNormalize whether do log transformation
#'
#' @return
#' @export
#'
#' @examples
preprocessing <- function(object, assay = list("RNA", "ATAC"), minFeatures = 200, minCells = 3, minCounts = NULL, maxCounts = NULL, libararyflag = TRUE, logNormalize = TRUE) {
  if (is.null(assay)) {
    for (i in 1:length(object@raw.data)) {
      object@norm.data[[i]] <- object@raw.data[[i]]
    }

  } else {

    for (i in 1:length(assay)) {
      iniData <- object@raw.data[[assay[[i]]]]
      print(dim(iniData))

      if (class(iniData) == "data.frame") {
        iniData <- as.matrix(iniData)
      }
      # filter cells that have features less than #minFeatures
      msum <- colSums(iniData != 0)
      proData <- iniData[, msum > minFeatures]
      # filter cells that have UMI counts less than #minCounts
      if (!is.null(minCounts)) {
        proData <- proData[, colSums(proData) > minCounts]
      }

      # filter cells that have expressed genes high than #maxGenes
      if (!is.null(maxCounts)) {
        proData <- proData[, colSums(proData) < maxCounts]
      }

      # filter genes that only express less than #minCells cells
      proData <- proData[rowSums(proData != 0) > minCells, ]
      # normalization:we employ a global-scaling normalization method that normalizes the gene expression measurements for each cell by the total expression multiplies this by a scale factor (10,000 by default)
      if (libararyflag) {
        proData <- sweep(proData, 2, colSums(proData), FUN = `/`) * 10000
      }
      if (logNormalize) {
        proData = log(proData + 1)
      }
      object@norm.data[[assay[[i]]]] <- proData

    }
  }
if (length(assay) == 2) {
        X1 <- object@norm.data[[assay[[1]]]]
        X2 <- object@norm.data[[assay[[2]]]]
        cell.keep <- intersect(colnames(X1), colnames(X2))
        object@norm.data[[assay[[1]]]] <- X1[, cell.keep]
        object@norm.data[[assay[[2]]]] <- X2[, cell.keep]
    } else if (length(assay) == 1) {
        X1 <- object@norm.data[[assay[[1]]]]
        assay2 <- setdiff(names(object@raw.data), assay[[1]])
        X2 <- object@raw.data[[assay2]]
        cell.keep <- intersect(colnames(X1), colnames(X2))
        object@norm.data[[assay[[1]]]] <- X1[, cell.keep]
        object@norm.data[[assay2]] <- X2[, cell.keep]
    }
  names(object@norm.data) <- names(object@raw.data)
  return(object)
}



#' add the cell information into pData slot
#'
#' @param object scAi object
#' @param pdata cell information to be added
#' @param pdata.name the name of column to be assigned
#'
#' @return
#' @export
#'
#' @examples
addpData <- function(object, pdata, pdata.name = NULL) {
  if (is.null(x = pdata.name) && is.atomic(x = pdata)) {
    stop("'pdata.name' must be provided for atomic pdata types (eg. vectors)")
  }
  if (inherits(x = pdata, what = c("matrix", "Matrix"))) {
    pdata <- as.data.frame(x = pdata)
  }

  if (is.null(x = pdata.name)) {
    pdata.name <- names(pdata)
  } else {
    names(pdata) <- pdata.name
  }
  object@pData <- pdata
  return(object)
}



#' run scAI model
#'
#' @param object scAI object
#' @param K Rank of the inferred factor
#' @param nrun Number of times to repreat the running
#' @param hvg.use1 whether use high variable genes for RNA-seq data
#' @param hvg.use2 whether use high variable genes for ATAC-seq data
#' @param keep_all Whether keep all the results from multiple runs
#' @param s Probability of Bernoulli distribution
#' @param alpha model parameter
#' @param lambda model parameter
#' @param gamma model parameter
#' @param maxIter Maximum number of iteration
#' @param stop_rule Stop rule to be used
#' @param init List of the initialized low-rank matrices
#' @param rand.seed seed
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach foreach "%dopar%"
#' @importFrom parallel makeForkCluster makeCluster detectCores
#' @importFrom doParallel registerDoParallel
run_scAI <- function(object, K, nrun = 5, hvg.use1 = TRUE, hvg.use2 = FALSE, keep_all = F, s = 0.25, alpha = 1, lambda = 100000, gamma = 1, maxIter = 200, stop_rule = 1, init = NULL, rand.seed = 1) {
  if (!is.null(init)) {
    W1.init = init$W1.init
    W2.init = init$W2.init
    H.init = init$H.init
    Z.init = init$Z.init
    R.init = init$R.init
  } else {
    R.init = NULL
    W1.init = NULL
    W2.init = NULL
    H.init = NULL
    Z.init = NULL
  }
  options(warn = -1)
  # Calculate the number of cores
  numCores <- min(parallel::detectCores(), nrun)
  cl <- tryCatch({
    parallel::makeForkCluster(numCores)
  }, error = function(e) {
    parallel::makeCluster(numCores)
  })
  doParallel::registerDoParallel(cl)
  
  if (hvg.use1) {
        X1 <- as.matrix(object@norm.data[[1]][object@var.features[[1]], ])
    } else {
        X1 <- as.matrix(object@norm.data[[1]])
    }
    if (hvg.use2) {
        X2 <- as.matrix(object@norm.data[[2]][object@var.features[[2]], ])
    } else {
        X2 <- as.matrix(object@norm.data[[2]])
    }

  X1 <- as.matrix(object@norm.data[[1]])
  X2 <- as.matrix(object@norm.data[[2]])
  outs <- foreach(i = 1:nrun, .packages = c("Matrix")) %dopar% {
    set.seed(rand.seed + i - 1)
    scAImodel(X1, X2, K = K, s = s, alpha = alpha, lambda = lambda, gamma = gamma, maxIter = maxIter, stop_rule = stop_rule,
         R.init = R.init, W1.init = W1.init, W2.init = W2.init, H.init = H.init, Z.init = Z.init)
  }

  objs <- foreach(i = 1:nrun, .combine = c) %dopar% {
    W1 <- outs[[i]]$W1
    W2 <- outs[[i]]$W2
    sum(cor(as.matrix(W1))) + sum(cor(as.matrix(W2)))
  }
  parallel::stopCluster(cl)
  N <- ncol(X1)
  C <- matrix(0, N, N)
  for (i in seq_len(nrun)) {
    H <- outs[[i]]$H
    H <- sweep(H, 2, colSums(H), FUN = `/`)
    clusIndex <- apply(H, 2, which.max)
    # compute the consensus matrix
    adjMat <- clust2Mat(clusIndex)
    C <- C + adjMat
  }
  CM <- C/nrun

  if (!keep_all) {
    sprintf("The best seed is %d", which.min(objs))
    outs_final <- outs[[which.min(objs)]]
    object@agg.data <- outs_final$agg.data
    W = list(W1 = outs_final$W1, W2 = outs_final$W2)
    names(W) <- names(object@norm.data)
    object@fit <- list(W = W, H = outs_final$H, Z = outs_final$Z, R = outs_final$R)
    object@cluster$consensus <- CM
    object@options$cost <- objs
    object@options$paras <- outs_final$options
    object@options$paras$nrun <- nrun
    object@options$best.seed <- which.min(objs)
    return(object)
  } else {
    outs_final <- list()
    outs_final$best <- outs[[which.min(objs)]]
    outs_final$best$consensus <- CM
    outs_final$nruns <- outs
    outs_final$options$cost <- objs
    return(outs_final)
  }
}



#' Solving the optimization problem in scAI
#'
#' @param X1 Single-cell transcriptomic data matrix (norm.data)
#' @param X2 Single-cell epigenomic data matrix (norm.data)
#' @param K Rank of inferred factors
#' @param s Probability of Bernoulli distribution
#' @param alpha model parameter
#' @param lambda model parameter
#' @param gamma model parameter
#' @param maxIter Maximum number of iteration
#' @param stop_rule Stop rule to be used
#' @param R.init initialization of R
#' @param W1.init initialization of W1
#' @param W2.init initialization of W2
#' @param H.init initialization of H
#' @param Z.init initialization of Z
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom stats rbinom runif
#' @importFrom rfunctions crossprodcpp
scAImodel <- function(X1, X2, K, s = 0.25, alpha = 1, lambda = 100000, gamma = 1, maxIter = 500, stop_rule = 1,
                 R.init = NULL, W1.init = NULL, W2.init = NULL, H.init = NULL, Z.init = NULL) {
  # Initialization W1,W2,H and Z
  p <- nrow(X1)
  n <- ncol(X1)
  q = nrow(X2)
  if (is.null(W1.init)) {
    W1.init = matrix(runif(p * K), p, K)
  }
  if (is.null(W2.init)) {
    W2.init = matrix(runif(q * K), q, K)
  }
  if (is.null(H.init)) {
    H.init = matrix(runif(K * n), K, n)
  }
  if (is.null(Z.init)) {
    Z.init = matrix(runif(n), n, n)
  }
  if (is.null(R.init)) {
    R.init = matrix(rbinom(n * n, 1, s), n, n)
  }
  W1 = W1.init
  W2 = W2.init
  H = H.init
  Z = Z.init
  R = R.init

  # start the clock to measure the execution time
  ptm = proc.time()
  eps <- .Machine$double.eps
  onesM_K <- matrix(1, K, K)

  XtX2 <- crossprodcpp(X2)
  index <- which(R == 0)
  for (iter in 1:maxIter) {
    # normalize H
    H = H/rowSums(H)
    # update W1
    HHt <- tcrossprod(H)
    X1Ht <- eigenMapMattcrossprod(X1, H)
    W1HHt <- eigenMapMatMult(W1, HHt)
    W1 <- W1 * X1Ht/(W1HHt + eps)

    # update W2
    ZR <- Z
    ZR[index] <- 0
    ZRHt <- eigenMapMattcrossprod(ZR, H)
    X2ZRHt <- eigenMapMatMult(X2, ZRHt)
    W2HHt <- eigenMapMatMult(W2, HHt)
    W2 = W2 * X2ZRHt/(W2HHt + eps)

    # update H
    W1tX1 <- eigenMapMatcrossprod(W1, X1)
    W2tX2 <- eigenMapMatcrossprod(W2, X2)
    W2tX2ZR <- eigenMapMatMult(W2tX2, ZR)
    HZZt <- eigenMapMatMult(H, Z + t(Z))
    W1tW1 <- crossprodcpp(W1)
    W2tW2 <- crossprodcpp(W2)
    temp1 <- H * (alpha * W1tX1 + W2tX2ZR + lambda * HZZt)
    temp2 <- eigenMapMatMult(alpha * W1tW1 + W2tW2 + 2 * lambda * HHt + gamma * onesM_K, H)
    H <- temp1/(temp2 + eps)

    # update Z
    HtH <- crossprodcpp(H)
    X2tW2H <- eigenMapMatcrossprod(W2tX2, H)
    RX2tW2H = X2tW2H
    RX2tW2H[index] = 0
    XtX2ZR <- eigenMapMatMult(XtX2, ZR)
    XtX2ZRR = XtX2ZR
    XtX2ZRR[index] = 0
    Z = Z * (RX2tW2H + lambda * HtH)/(XtX2ZRR + lambda * Z + eps)

    if (stop_rule == 2) {
      obj = alpha * norm(X1 - W1 %*% H, "F")^2 + norm(X2 %*% (Z * R) - W2 %*% H, "F")^2 + lambda * norm(Z - HtH, "F") + gamma * sum(colSums(H) * colSums(H))
      if (iter > 1 && ((obj_old - obj)/obj_old < 10^(-6)) || iter == maxIter) {
        break
      }
      obj_old = obj
    }
  }
  # compute the execution time
  execution.time = proc.time() - ptm

  ZR <- Z
  ZR[index] <- 0
  ZR <- sweep(ZR, 2, colSums(ZR), FUN = `/`)
  X2agg <- eigenMapMatMult(X2, ZR)
  X2agg <- sweep(X2agg, 2, colSums(X2agg), FUN = `/`) * 10000
  X2agg <- log(1 + X2agg)

  barcodes <- colnames(X2)
  feature1 <- rownames(X1)
  feature2 <- rownames(X2)
  names_com <- paste0("factor", seq_len(K))
  attr(X2agg, "dimnames") <- list(feature2, barcodes)
  attr(W1, "dimnames") <- list(feature1, names_com)
  attr(W2, "dimnames") <- list(feature2, names_com)
  attr(H, "dimnames") <- list(names_com, barcodes)
  attr(Z, "dimnames") <- list(barcodes, barcodes)

  outs <- list(agg.data = X2agg, W1 = W1, W2 = W2, H = H, Z = Z, R = R,
               options = list(s = s, alpha = alpha, lambda = lambda, gamma = gamma, maxIter = maxIter, stop_rule = stop_rule, run.time = execution.time))
  return(outs)

}



#' Perform dimensional reduction
#'
#' Dimension reduction by PCA, t-SNE or UMAP
#' @param object scAI object
#' @param return.object whether return scAI object
#' @param data.use input data
#' @param do.scale whether scale the data
#' @param do.center whether scale and center the data
#' @param method Method of dimensional reduction, one of tsne, FItsne and umap
#' @param rand.seed Set a random seed. By default, sets the seed to 42.
#' @param perplexity perplexity parameter in tsne
#' @param theta parameter in tsne
#' @param check_duplicates parameter in tsne

#' @param FItsne.path File path of FIt-SNE
#' @param dim.embed dimensions of t-SNE embedding
#' @param dim.use num of PCs used for t-SNE
#'
#' @param do.fast whether do fast PCA
#' @param dimPC the number of components to keep in PCA
#' @param weight.by.var whether use weighted pc.scores
#'
#' @param n.neighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
#' @param n.components The dimension of the space to embed into.
#' @param distance This determines the choice of metric used to measure distance in the input space.
#' @param n.epochs the number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. If NULL is specified, a value will be selected based on the size of the input dataset (200 for large datasets, 500 for small).
#' @param learning.rate The initial learning rate for the embedding optimization.
#' @param min.dist This controls how tightly the embedding is allowed compress points together.
#' Larger values ensure embedded points are moreevenly distributed, while smaller values allow the
#' algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
#' @param spread he effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.
#' @param set.op.mix.ratio Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets.
#' @param local.connectivity The local connectivity required - i.e. the number of nearest neighbors
#' that should be assumed to be connected at a local level. The higher this value the more connected
#' the manifold becomes locally. In practice this should be not more than the local intrinsic dimension of the manifold.
#' @param repulsion.strength Weighting applied to negative samples in low dimensional embedding
#' optimization. Values higher than one will result in greater weight being given to negative samples.
#' @param negative.sample.rate The number of negative samples to select per positive sample in the
#' optimization process. Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
#' @param a More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param b More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom Rtsne Rtsne
#' @importFrom reticulate py_module_available py_set_seed import
#'

reducedDims <- function(object, data.use = object@fit$H, do.scale = TRUE, do.center = TRUE, return.object = TRUE, method = "umap",
                        dim.embed = 2, dim.use = NULL, perplexity = 30, theta = 0.5, check_duplicates = F, rand.seed = 42L,
                        FItsne.path = NULL,
                        dimPC = 40,do.fast = TRUE, weight.by.var = TRUE,
                        n.neighbors = 30L, n.components = 2L, distance = "correlation",n.epochs = NULL,learning.rate = 1.0,min.dist = 0.3,spread = 1.0,set.op.mix.ratio = 1.0,local.connectivity = 1L,
                        repulsion.strength = 1,negative.sample.rate = 5,a = NULL,b = NULL
                        ) {

  data.use <- as.matrix(data.use)
  if (do.scale) {
    data.use <- t(scale(t(data.use), center = do.center, scale = TRUE))
    data.use[is.na(data.use)] <- 0
  }
  if (!is.null(dim.use)) {
    data.use = object@embed$pca[, dim.use]
  }
  if (method == "pca") {
    cell_coords <- runPCA(data.use, do.fast = do.fast, dimPC = dimPC, seed.use = rand.seed, weight.by.var = weight.by.var)
    object@embed$pca <- cell_coords

  } else if (method == "tsne") {
    set.seed(rand.seed)
    cell_coords <- Rtsne::Rtsne(t(data.use), pca = FALSE, dims = dim.embed, theta = theta, perplexity = perplexity, check_duplicates = F)$Y
    rownames(cell_coords) <- rownames(t(data.use))
    object@embed$tsne <- cell_coords

  } else if (method == "FItsne") {
    if (!exists("FItsne")) {
      if (is.null(fitsne.path)) {
        stop("Please pass in path to FIt-SNE directory as FItsne.path.")
      }
      source(paste0(fitsne.path, "/fast_tsne.R"), chdir = T)
    }
    cell_coords <- fftRtsne(t(data.use), pca = FALSE, dims = dim.embed, theta = theta, perplexity = perplexity, check_duplicates = F, rand_seed = rand.seed)
    rownames(cell_coords) <- rownames(t(data.use))
    object@embed$FItsne <- cell_coords

  } else if (method == "umap") {
    cell_coords <- runUMAP(data.use,
                           n.neighbors = n.neighbors,
                           n.components = n.components,
                           distance = distance,
                           n.epochs = n.epochs,
                           learning.rate = learning.rate,
                           min.dist = min.dist,
                           spread = spread,
                           set.op.mix.ratio = set.op.mix.ratio,
                           local.connectivity = local.connectivity,
                           repulsion.strength = repulsion.strength,
                           negative.sample.rate = negative.sample.rate,
                           a = a,
                           b = b,
                           seed.use = rand.seed)
    object@embed$umap <- cell_coords

  } else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  if (return.object) {
    return(object)
  } else {
    return(cell_coords)
  }

}

#' Dimension reduction using PCA
#'
#' @param data.use input data
#' @param do.fast whether do fast PCA
#' @param dimPC the number of components to keep
#' @param seed.use set a seed
#' @param weight.by.var whether use weighted pc.scores
#' @importFrom stats prcomp
#' @importFrom irlba irlba
#' @return
#' @export
#'
#' @examples
runPCA <- function(data.use, do.fast = T, dimPC = 50, seed.use = 42, weight.by.var = T) {
  set.seed(seed = seed.use)
  if (do.fast) {
    dimPC <- min(dimPC, nrow(data.use) - 1)
    pca.res <- irlba::irlba(t(data.use), nv = dimPC)
    sdev <- pca.res$d/sqrt(max(1, ncol(data.use) - 1))
    if (weight.by.var){
      pc.scores <- pca.res$u %*% diag(pca.res$d)
    } else {
      pc.scores <- pca.res$u
    }
  } else {
    dimPC <- min(dimPC, nrow(data.use) - 1)
    pca.res <- stats::prcomp(x = t(data.use), rank. = dimPC)
    sdev <- pca.res$sdev
    if (weight.by.var) {
      pc.scores <- pca.res$x %*% diag(pca.res$sdev[1:dimPC]^2)
    } else {
      pc.scores <- pca.res$x
    }
  }
  rownames(pc.scores) <- colnames(data.use)
  colnames(pc.scores) <- paste0('PC', 1:ncol(pc.scores))
  cell_coords <- pc.scores
  return(cell_coords)
}

#' Perform dimension reduction using UMAP
#'
#' @param data.use input data
#' @param n.neighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
#' @param n.components The dimension of the space to embed into.
#' @param distance This determines the choice of metric used to measure distance in the input space.
#' @param n.epochs the number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. If NULL is specified, a value will be selected based on the size of the input dataset (200 for large datasets, 500 for small).
#' @param learning.rate The initial learning rate for the embedding optimization.
#' @param min.dist This controls how tightly the embedding is allowed compress points together.
#' Larger values ensure embedded points are moreevenly distributed, while smaller values allow the
#' algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
#' @param spread he effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.
#' @param set.op.mix.ratio Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets.
#' @param local.connectivity The local connectivity required - i.e. the number of nearest neighbors
#' that should be assumed to be connected at a local level. The higher this value the more connected
#' the manifold becomes locally. In practice this should be not more than the local intrinsic dimension of the manifold.
#' @param repulsion.strength Weighting applied to negative samples in low dimensional embedding
#' optimization. Values higher than one will result in greater weight being given to negative samples.
#' @param negative.sample.rate The number of negative samples to select per positive sample in the
#' optimization process. Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
#' @param a More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param b More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param seed.use Set a random seed. By default, sets the seed to 42.
#' @param metric.kwds A dictionary of arguments to pass on to the metric, such as the p value for Minkowski distance
#' @param angular.rp.forest Whether to use an angular random projection forest to initialise the
#' approximate nearest neighbor search. This can be faster, but is mostly on useful for metric that
#' use an angular style distance such as cosine, correlation etc. In the case of those metrics angular forests will be chosen automatically.
#' @param verbose Controls verbosity
#' This function is modified from Seurat package
#' @return
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
#' @examples
runUMAP <- function(
  data.use,
  n.neighbors = 30L,
  n.components = 2L,
  distance = "correlation",
  n.epochs = NULL,
  learning.rate = 1.0,
  min.dist = 0.3,
  spread = 1.0,
  set.op.mix.ratio = 1.0,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5,
  a = NULL,
  b = NULL,
  seed.use = 42L,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  verbose = TRUE){
  if (!reticulate::py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn or reticulate::py_install(packages = 'umap-learn')).")
  }
  set.seed(seed.use)
  reticulate::py_set_seed(seed.use)
  umap_import <- reticulate::import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(n.neighbors),
    n_components = as.integer(n.components),
    metric = distance,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    metric_kwds = metric.kwds,
    angular_rp_forest = angular.rp.forest,
    verbose = verbose
  )
  Rumap <- umap$fit_transform
  umap_output <- Rumap(t(data.use))
  colnames(umap_output) <- paste0('UMAP', 1:ncol(umap_output))
  rownames(umap_output) <- colnames(data.use)
  return(umap_output)
}




#' Identify enriched features in each factor
#'
#' Rank features in each factor by Factor loading analysis
#' @param object scAI object
#' @param assay Name of assay to be analyzed
#' @param features a vector of features
#' @param cutoff.W Threshold of feature loading values
#' @param cutoff.H Threshold of cell loading values
#' @param thresh.pc Threshold of the percent of cells enriched in one factor
#' @param thresh.fc Threshold of Fold Change
#' @param thresh.p Threshold of p-values
#' @param n.top Number of top features to be returned
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom stats sd wilcox.test
#' @importFrom dplyr top_n slice
identifyFactorMarkers <- function(object, assay, features = NULL,
                                     cutoff.W = 0.5, cutoff.H = 0.5,
                                     thresh.pc = 0.05, thresh.fc = 0.25, thresh.p = 0.05, n.top = 10) {
  X <- as.matrix(object@norm.data[[assay]])
  H <- object@fit$H
  W <- object@fit$W[[assay]]
  if (is.null(features)) {
    features <- row.names(W)
  }

  H <- sweep(H, 2, colSums(H), FUN = `/`)
  K = ncol(W)
  lib_W <- base::rowSums(W)
  lib_W[lib_W == 0] <- 1
  lib_W[lib_W < mean(lib_W) - 5 * sd(lib_W)] <- 1  #  omit the nearly null rows
  W <- sweep(W, 1, lib_W, FUN = `/`)
  MW <- base::colMeans(W)
  sW <- apply(W, 2, sd)
  # candidate markers for each factor
  IndexW_record <- vector("list", K)
  for (j in 1:K) {
    IndexW_record[[j]] <- which(W[, j] > MW[j] + cutoff.W * sW[j])
  }

  # divided cells into two groups
  mH <- apply(H, 1, mean)
  sH <- apply(H, 1, sd)
  IndexH_record1 <- vector("list", K)
  IndexH_record2 <- vector("list", K)
  for (i in 1:K) {
    IndexH_record1[[i]] <- which(H[i, ] > mH[i] + cutoff.H * sH[i])
    IndexH_record2[[i]] <- base::setdiff(c(1:ncol(H)), IndexH_record1[[i]])
  }
    my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
    )
  
  # identify factor-specific markers
  factor_markers = vector("list", K)
  markersTopn = vector("list", K)
  factors <- c()
  Features <- c()
  Pvalues <- c()
  Log2FC <- c()
  for (i in 1:K) {
    data1 <- X[IndexW_record[[i]], IndexH_record1[[i]]]
    data2 <- X[IndexW_record[[i]], IndexH_record2[[i]]]
    idx1 <- which(base::rowSums(data1 != 0) > thresh.pc * ncol(data1))  # at least expressed in thresh.pc% cells in one group
    FC <- log2(base::rowMeans(data1)/base::rowMeans(data2))
    idx2 <- which(FC > thresh.fc)
    pvalues <- my.sapply(
        X = 1:nrow(x = data1),
        FUN = function(x) {
            return(wilcox.test(data1[x, ], data2[x, ], alternative = "greater")$p.value)
        }
        )

    idx3 <- which(pvalues < thresh.p)
    idx = intersect(intersect(idx1, idx2), idx3)

    # order
    FC <- FC[idx]
    c = sort(FC, decreasing = T, index.return = T)$ix
    ri = IndexW_record[[i]]

    Pvalues <- c(Pvalues, pvalues[ri[idx[c]]])
    Log2FC <- c(Log2FC, FC[c])
    factors <- c(factors, rep(i, length(c)))
    Features <- c(Features, features[ri[idx[c]]])

  }

  markers.all <- cbind(factors = factors, features = Features, pvalues = Pvalues, log2FC = Log2FC)

  rownames(markers.all) <- Features
  markers.all <- as.data.frame(markers.all)
  markers.top <- markers.all %>% dplyr::group_by(factors) %>% top_n(n.top, log2FC) %>% dplyr::slice(1:n.top)

  markers = list(markers.all = markers.all, markers.top = markers.top)

  return(markers)
}




#' compute the 2D coordinates of embeded cells, genes, loci and factors using VscAI visualization
#'
#' @param object scAI object
#' @param genes.embed A vector of genes to be embedded
#' @param loci.embed A vector of loci to be embedded
#' @param alpha.embed Parameter controlling the distance between the cells and the factors
#' @param snn.embed Number of neighbors
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom Matrix Matrix
getEmbeddings <- function(object, genes.embed = NULL, loci.embed = NULL, alpha.embed = 1.9, snn.embed = 5) {
  H <- object@fit$H
  W1 <- object@fit$W[[1]]
  W2 <- object@fit$W[[2]]
  Z <- object@fit$Z

  if (nrow(H) < 3 & length(object@fit.variedK) == 0) {
    print("VscAI needs at least three factors for embedding. Now rerun scAI with rank being 3...")
    outs <- run_scAI(object, K = 3, nrun = object@options$paras$nrun,
                     s = object@options$paras$s, alpha = object@options$paras$alpha, lambda = object@options$paras$lambda, gamma = object@options$paras$gamma,
                     maxIter = object@options$paras$maxIter, stop_rule = object@options$paras$stop_rule)
    W <- outs@fit$W
    H <- outs@fit$H
    Z <- outs@fit$Z
    W1 <- W[[1]]
    W2 <- W[[2]]
    object@fit.variedK$W <- W
    object@fit.variedK$H <- H
    object@fit.variedK$Z <- Z
  } else if (nrow(H) < 3 & length(object@fit.variedK) != 0) {
    print("Using the previous calculated factors for embedding.")
    W1 <- object@fit.variedK$W[[1]]
    W2 <- object@fit.variedK$W[[2]]
    H <- object@fit.variedK$H
    Z <- object@fit.variedK$Z

  }

  nmf.scores <- H

  Zsym <- Z/max(Z)
  Zsym[Zsym < 10^(-6)] <- 0
  Zsym <- (Zsym + t(Zsym))/2
  diag(Zsym) <- 1
  rownames(Zsym) <- colnames(H)
  colnames(Zsym) <- rownames(Zsym)

  alpha.exp <- alpha.embed  # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
  snn.exp <- snn.embed  # Lower this < 1.0 to move similar cells closer to each other
  n_pull <- nrow(H)  # The number of factors pulling on each cell. Must be at least 3.

  Zsym <- Matrix(data = Zsym, sparse = TRUE)
  snn <- Zsym
  swne.embedding <- swne::EmbedSWNE(nmf.scores, snn, alpha.exp = alpha.exp, snn.exp = snn.exp, n_pull = n_pull, proj.method = "sammon", dist.use = "cosine")
  # project cells and factors
  factor_coords <- swne.embedding$H.coords
  cell_coords <- swne.embedding$sample.coords
  rownames(factor_coords) <- rownames(H)
  rownames(cell_coords) <- colnames(H)

  # project genes onto the low dimension space by VscAI
  if (is.null(genes.embed)) {
    genes.embed <- rownames(W1)
  } else {
    genes.embed <- as.character(as.vector(as.matrix((genes.embed))))
  }

  swne.embedding <- swne::EmbedFeatures(swne.embedding, W1, genes.embed, n_pull = n_pull)
  gene_coords <- swne.embedding$feature.coords
  rownames(gene_coords) <- gene_coords$name

  # project loci
  if (is.null(loci.embed)) {
    loci.embed <- rownames(W2)
  } else {
    loci.embed <- as.character(as.vector(as.matrix((loci.embed))))
  }

  swne.embedding <- swne::EmbedFeatures(swne.embedding, W2, loci.embed, n_pull = n_pull)
  loci_coords <- swne.embedding$feature.coords
  rownames(loci_coords) <- loci_coords$name

  object@embed$VscAI$cells <- cell_coords
  object@embed$VscAI$factors <- factor_coords
  object@embed$VscAI$genes <- gene_coords
  object@embed$VscAI$loci <- loci_coords
  return(object)
}



#' Identify cell clusters
#'
#' @param object scAI object
#' @param partition.type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition (see the Leiden python module documentation for more details)
#' @param seed.use set seed
#' @param n.iter number of iteration
#' @param initial.membership arameters to pass to the Python leidenalg function defaults initial_membership=None
#' @param weights defaults weights=None
#' @param node.sizes Parameters to pass to the Python leidenalg function
#' @param resolution A parameter controlling the coarseness of the clusters
#' @param K Number of clusters if performing hierarchical clustering of the consensus matrix
#' @return
#' @export
#'
#' @examples
#' @importFrom stats as.dist cophenetic cor cutree dist hclust
identifyClusters <- function(object, resolution = 1, partition.type = "RBConfigurationVertexPartition",
                             seed.use = 42L,n.iter = 10L,
                             initial.membership = NULL, weights = NULL, node.sizes = NULL,
                             K = NULL) {
  if (is.null(K) & !is.null(resolution)) {
    data.use <- object@fit$H
    data.use <- as.matrix(data.use)
    data.use <- t(scale(t(data.use), center = TRUE, scale = TRUE))
    snn <- swne::CalcSNN(data.use, k = 20, prune.SNN = 1/15)
    idents <- runLeiden(SNN = snn, resolution = resolution, partition_type = partition.type,
    seed.use = seed.use, n.iter = n.iter,
    initial.membership = initial.membership, weights = weights, node.sizes = node.sizes)
  } else {
    CM <- object@cluster$consensus
    d <- as.dist(1 - CM)
    hc <- hclust(d, "ave")
    idents<- hc %>% cutree(k = K)
    coph <- cophenetic(hc)
    cs <- cor(d, coph)
    object@cluster$cluster.stability <- cs
  }
  names(idents) <- rownames(object@pData)
  object@identity <- factor(idents)
  object@pData$cluster <- factor(idents)
  return(object)
}


#' Run Leiden clustering algorithm
#' This code is modified from Tom Kelly (https://github.com/TomKellyGenetics/leiden), where we added more parameters (seed.use and n.iter) to run the Python version. In addition, we also take care of the singleton issue after running leiden algorithm.
#' @description Implements the Leiden clustering algorithm in R using reticulate to run the Python version. Requires the python "leidenalg" and "igraph" modules to be installed. Returns a vector of partition indices.
#' @param SNN An adjacency matrix compatible with \code{\link[igraph]{igraph}} object.
#' @param seed.use set seed
#' @param n.iter number of iteration
#' @param initial.membership arameters to pass to the Python leidenalg function defaults initial_membership=None
#' @param node.sizes Parameters to pass to the Python leidenalg function
#' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition (see the Leiden python module documentation for more details)
#' @param resolution A parameter controlling the coarseness of the clusters
#' @param weights defaults weights=None
#' @return A parition of clusters as a vector of integers
##' @examples
##' #check if python is availble
##' modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
##' if(modules){
##' #generate partitions
##' partition <- leiden(adjacency_matrix)
##' table(partition)
##'
##' #generate partitions at a lower resolution
##' partition <- leiden(adjacency_matrix, resolution_parameter = 0.5)
##' table(partition)
##'
##' #generate example weights
##' weights <- sample(1:10, sum(adjacency_matrix!=0), replace=TRUE)
##' partition <- leiden(adjacency_matrix, weights = weights)
##' table(partition)
##'
##' #generate example weighted matrix
##' adjacency_matrix[adjacency_matrix == 1] <- weights
##' partition <- leiden(adjacency_matrix)
##' table(partition)
##' }
##'
##'
#' @keywords graph network igraph mvtnorm simulation
#' @importFrom reticulate import r_to_py
##' @export

runLeiden <- function(SNN = matrix(), resolution = 1, partition_type = c(
  'RBConfigurationVertexPartition',
  'ModularityVertexPartition',
  'RBERVertexPartition',
  'CPMVertexPartition',
  'MutableVertexPartition',
  'SignificanceVertexPartition',
  'SurpriseVertexPartition'),
seed.use = 42L,
n.iter = 10L,
initial.membership = NULL, weights = NULL, node.sizes = NULL) {
  if (!reticulate::py_module_available(module = 'leidenalg')) {
    stop("Cannot find Leiden algorithm, please install through pip (e.g. pip install leidenalg).")
  }

  #import python modules with reticulate
  leidenalg <- import("leidenalg", delay_load = TRUE)
  ig <- import("igraph", delay_load = TRUE)

  resolution_parameter <- resolution
  initial_membership <- initial.membership
  node_sizes <- node.sizes
  #convert matrix input (corrects for sparse matrix input)
  adj_mat <- as.matrix(SNN)

  #compute weights if non-binary adjacency matrix given
  is_pure_adj <- all(as.logical(adj_mat) == adj_mat)
  if (is.null(weights) && !is_pure_adj) {
    #assign weights to edges (without dependancy on igraph)
    weights <- t(adj_mat)[t(adj_mat)!=0]
    #remove zeroes from rows of matrix and return vector of length edges
  }

  ##convert to python numpy.ndarray, then a list
  adj_mat_py <- r_to_py(adj_mat)
  adj_mat_py <- adj_mat_py$tolist()

  #convert graph structure to a Python compatible object
  GraphClass <- if (!is.null(weights) && !is_pure_adj){
    ig$Graph$Weighted_Adjacency
  } else {
    ig$Graph$Adjacency
  }
  snn_graph <- GraphClass(adj_mat_py)

  #compute partitions
  partition_type <- match.arg(partition_type)
  part <- switch(
    EXPR = partition_type,
    'RBConfigurationVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBConfigurationVertexPartition,
      initial_membership = initial.membership, weights = weights,
      resolution_parameter = resolution,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'ModularityVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$ModularityVertexPartition,
      initial_membership = initial.membership, weights = weights,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'RBERVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBERVertexPartition,
      initial_membership = initial.membership, weights = weights, node_sizes = node.sizes,
      resolution_parameter = resolution_parameter,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'CPMVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$CPMVertexPartition,
      initial_membership = initial.membership, weights = weights, node_sizes = node.sizes,
      resolution_parameter = resolution,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'MutableVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$MutableVertexPartition,
      initial_membership = initial.membership,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'SignificanceVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SignificanceVertexPartition,
      initial_membership = initial.membership, node_sizes = node.sizes,
      resolution_parameter = resolution,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'SurpriseVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SurpriseVertexPartition,
      initial_membership = initial.membership, weights = weights, node_sizes = node.sizes,
      n_iterations = n.iter,
      seed = seed.use
    ),
    stop("please specify a partition type as a string out of those documented")
  )
  partition <- part$membership+1
  idents <- partition

  if (min(table(idents)) == 1) {
    idents <- assignSingletons(idents, SNN)
  }
  idents <- factor(idents)
  names(idents) <- row.names(SNN)
  return(idents)
}

# Group single cells that make up their own cluster in with the cluster they are most connected to.
#
# @param idents  clustering result
# @param SNN     SNN graph
# @return        Returns scAI object with all singletons merged with most connected cluster

assignSingletons <- function(idents, SNN) {
  # identify singletons
  singletons <- c()
  for (cluster in unique(idents)) {
    if (length(which(idents %in% cluster)) == 1) {
      singletons <- append(x = singletons, values = cluster)
    }
  }
  #singletons = names(table(idents))[which(table(idents)==1)]
  # calculate connectivity of singletons to other clusters, add singleton to cluster it is most connected to
  cluster_names <- unique(x = idents)
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode="numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  for (i in singletons) {
    print(i)
    for (j in cluster_names) {
      subSNN = SNN[
        which(idents %in% i),
        which(idents %in% j)
        ]
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
    which(idents %in% i)[which(idents %in% i)] <- closest_cluster
  }
  if (length(x = singletons) > 0) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(idents)),
      "final clusters."
    ))
  }
  return(idents)
}


#' Convert membership into an adjacent matrix
#'
#' @param memb Membership vector
#'
#' @return
#' @export
#'
#' @examples
clust2Mat <- function(memb) {
  N <- length(memb)
  return(as.numeric(outer(memb, memb, FUN = "==")) - outer(1:N, 1:N, "=="))
}


#' Identify enriched features in each cell cluster
#'
#' @param object scAI object
#' @param assay Name of assay to be analyzed
#' @param features a vector of features
#' @param test.use which test to use ("bimod" or "wilcox")
#' @param thresh.pc Threshold of the percent of cells enriched in one cluster
#' @param thresh.fc Threshold of Fold Change
#' @param thresh.p Threshold of p-values
#' @param n.top Number of top features to be returned
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @return
#' @export
#'
#' @examples
identifyClusterMarkers <- function(object, assay, features = NULL, test.use = "bimod",
                                   thresh.pc = 0.05, thresh.fc = 0.25, thresh.p = 0.05, n.top = 10) {
  X <- as.matrix(object@norm.data[[assay]])
  H <- object@fit$H
  W <- object@fit$W[[assay]]
  identity <- object@identity

  if (is.null(features)) {
    features <- row.names(W)
  }

  lib_W <- rowSums(W)
  lib_W[lib_W == 0] <- 1
  lib_W[lib_W < mean(lib_W) - 5 * sd(lib_W)] <- 1  #  omit the nearly null rows
  W <- sweep(W, 1, lib_W, FUN = `/`)
  MW <- colMeans(W)
  sW <- apply(W, 2, sd)
  # candidate markers for each factor
  IndexW_record <- vector("list", ncol(W))
  for (j in 1:ncol(W)) {
    IndexW_record[[j]] <- which(W[, j] > MW[j])
  }
  IndexW_record <- unique(c(IndexW_record, recursive = TRUE))

  # divided cells into two groups
  K <- length(unique(identity))

  IndexH_record1 <- vector("list", K)
  IndexH_record2 <- vector("list", K)
  for (i in 1:K) {
    IndexH_record1[[i]] <- which(identity == i)
    IndexH_record2[[i]] <- base::setdiff(c(1:length(identity)), IndexH_record1[[i]])
  }

  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )

  # identify factor-specific markers
  factor_markers = vector("list", K)
  markersTopn = vector("list", K)
  factors <- c()
  Features <- c()
  Pvalues <- c()
  Log2FC <- c()
  for (i in 1:K) {
    data1 <- X[IndexW_record, IndexH_record1[[i]]]
    data2 <- X[IndexW_record, IndexH_record2[[i]]]
    idx1 <- which(rowSums(data1 != 0) > thresh.pc * ncol(data1))  # at least expressed in thresh.pc cells in one group
    FC <- log2(rowMeans(data1)/rowMeans(data2))
    idx2 <- which(FC > thresh.fc)
    if (test.use == "wilcox") {
      pvalues <- unlist(
        x = my.sapply(
          X = 1:nrow(x = data1),
          FUN = function(x) {
            return(wilcox.test(data1[x, ], data2[x, ], alternative = "greater", correct = T)$p.value)
          }
        )
      )

    } else if (test.use == 'bimod') {
      pvalues <- unlist(
        x = my.sapply(
          X = 1:nrow(x = data1),
          FUN = function(x) {
            return(DifferentialLRT(
              x = as.numeric(x = data1[x,]),
              y = as.numeric(x = data2[x,])
            ))
          }
        )
      )
    }

    idx3 <- which(pvalues < thresh.p)
    idx = intersect(intersect(idx1, idx2), idx3)

    # order
    FC <- FC[idx]
    c = sort(FC, decreasing = T, index.return = T)$ix
    ri = IndexW_record

    Pvalues <- c(Pvalues, pvalues[ri[idx[c]]])
    Log2FC <- c(Log2FC, FC[c])
    factors <- c(factors, rep(i, length(c)))
    Features <- c(Features, features[ri[idx[c]]])

  }

  markers.all <- data.frame(factors = factors, features = as.character(Features),
                            pvalues = as.numeric(as.character(Pvalues)), log2FC = as.numeric(as.character(Log2FC)))

  markers.top <- markers.all %>% group_by(factors) %>% top_n(n.top, log2FC) %>% slice(1:n.top)

  markers = list(markers.all = markers.all, markers.top = markers.top)


  return(markers)
}

# function to run mcdavid et al. DE test
#' likelood ratio test
#'
#' @param x a vector
#' @param y a vector
#' @param xmin threshold for the values in the vector
#'
#' @return
#'
#' @importFrom stats pchisq
DifferentialLRT <- function(x, y, xmin = 0) {
  lrtX <- bimodLikData(x = x)
  lrtY <- bimodLikData(x = y)
  lrtZ <- bimodLikData(x = c(x, y))
  lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  return(pchisq(q = lrt_diff, df = 3, lower.tail = F))
}

#' likelood ratio test
#' @importFrom stats sd dnorm
#' @param x a vector
#' @param xmin threshold for the values in the vector

bimodLikData <- function(x, xmin = 0) {
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  xal <- length(x = x2) / length(x = x)
  xal[xal > 1 - 1e-5] <- 1 - 1e-5
  xal[xal < 1e-5] <- 1e-5
  likA <- length(x = x1) * log(x = 1 - xal)
  if (length(x = x2) < 2) {
    mysd <- 1
  } else {
    mysd <- sd(x = x2)
  }
  likB <- length(x = x2) *
    log(x = xal) +
    sum(dnorm(x = x2, mean = mean(x = x2), sd = mysd, log = TRUE))
  return(likA + likB)
}


#' Infer regulatory relationships
#'
#' @param object scAI object
#' @param gene.use genes to be inferred
#' @param candinate_loci cadinate loci
#' @param cutoff_H threshold of coefficient values
#' @param cutoff the difference of correlation to be considered as significant
#' @param thresh_corr threshold of correlation coefficient
#'
#' @return
#' @export
#'
#' @examples
inferRegulations <- function(object, gene.use, candinate_loci, cutoff_H = 0.5, cutoff = 0.1, thresh_corr = 0.2) {
  H <- as.matrix(object@fit$H)
  H <- sweep(H, 2, colSums(H), FUN = `/`)
  K = nrow(H)
  mH <- apply(H, 1, mean)
  sH <- apply(H, 1, sd)
  IndexH_record <- vector("list", K)
  for (i in 1:K) {
    IndexH_record[[i]] <- which(H[i, ] > mH[i] + cutoff_H * sH[i])
  }
  gene.use <- as.character(gene.use)
  X1 <- object@norm.data[[1]]
  X1 <- X1[gene.use, ]
  X2a <- object@agg.data

  regulations <- vector("list", length(gene.use))
  names(regulations) <- gene.use
  for (i in 1:length(gene.use)) {
    regulations_i = vector("list", K)
    names(regulations_i) <- rownames(H)
    for (j in 1:K) {
      loci_j <- candinate_loci$markers.all$features[candinate_loci$markers.all$factors == j]
      # compute the correlation between
      x1 <- X1[i, ]
      x2a <- X2a[loci_j, ]
      cors1 <- cor(x1, t(x2a))

      # set the values of this gene and its candiate loci to be zero
      X1_new <- X1
      X2a_new <- X2a
      X1_new[gene.use[i], IndexH_record[[j]]] <- 0
      X2a_new[loci_j, IndexH_record[[j]]] <- 0
      x1_new <- X1_new[gene.use[i], ]
      x2a_new <- X2a_new[loci_j, ]
      cors2 <- cor(x1_new, t(x2a))
      cors3 <- cor(x1, t(x2a_new))
      cors1[is.na(cors1)] <- 0
      cors2[is.na(cors2)] <- 0
      cors3[is.na(cors3)] <- 0
      D <- rbind(cors1 - cors2, cors1 - cors3)
      flag <- (rowSums(abs(D) > cutoff) > 0) & abs(cors1) > thresh_corr
      regulations_i[[j]]$link <- as.character(loci_j[flag])
      regulations_i[[j]]$intensity <- cors1[flag]
    }
    regulations[[i]] <- regulations_i
  }
  return(regulations)
}



#' select highly variable features
#'
#' @param object scAI objecy
#' @param assay Name of assay
#' @param do.plot Whether showing plot
#' @param do.text Whether adding feature names
#' @param x.low.cutoff The minimum expression level
#' @param x.high.cutoff The maximum expression level
#' @param y.cutoff The minimum fano factor values
#' @param y.high.cutoff The maximum fano factor values
#' @param num.bin Number of bins
#' @param pch.use Shape of dots in ggplot
#' @param col.use Color of dots
#' @param cex.text.use Size of text
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom graphics smoothScatter text
selectFeatures <- function(object, assay = "RNA", do.plot = TRUE, do.text = TRUE,
                           x.low.cutoff = 0.01, x.high.cutoff = 3.5, y.cutoff = 1, y.high.cutoff = Inf,
                           num.bin = 20, pch.use = 16, col.use = "black", cex.text.use = 0.5) {
  # This function is modified from Seurat Package
  data <- object@norm.data[[assay]]
  genes.use <- rownames(data)

  gene.mean <- rep(x = 0, length(x = genes.use))
  names(x = gene.mean) <- genes.use
  gene.dispersion <- gene.mean
  gene.dispersion.scaled <- gene.mean
  bin.size <- 1000
  max.bin <- floor(x = length(x = genes.use)/bin.size) + 1

  for (i in 1:max.bin) {
    my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1)) + 1
    my.inds <- my.inds[my.inds <= length(x = genes.use)]
    genes.iter <- genes.use[my.inds]
    data.iter <- data[genes.iter, , drop = F]
    gene.mean[genes.iter] <- apply(X = data.iter, MARGIN = 1, FUN = mean.function)
    gene.dispersion[genes.iter] <- apply(X = data.iter, MARGIN = 1, FUN = dispersion.function)
  }
  gene.dispersion[is.na(x = gene.dispersion)] <- 0
  gene.mean[is.na(x = gene.mean)] <- 0
  data_x_bin <- cut(x = gene.mean, breaks = num.bin)

  names(x = data_x_bin) <- names(x = gene.mean)
  mean_y <- tapply(X = gene.dispersion, INDEX = data_x_bin, FUN = mean)
  sd_y <- tapply(X = gene.dispersion, INDEX = data_x_bin, FUN = sd)
  gene.dispersion.scaled <- (gene.dispersion - mean_y[as.numeric(x = data_x_bin)])/sd_y[as.numeric(x = data_x_bin)]
  gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0
  names(x = gene.dispersion.scaled) <- names(x = gene.mean)
  mv.df <- data.frame(gene.mean, gene.dispersion, gene.dispersion.scaled)
  rownames(x = mv.df) <- rownames(x = data)
  hvg.info <- mv.df


  names(x = gene.mean) <- names(x = gene.dispersion) <- names(x = gene.dispersion.scaled) <- rownames(hvg.info)
  pass.cutoff <- names(x = gene.mean)[which(x = ((gene.mean > x.low.cutoff) & (gene.mean < x.high.cutoff)) & (gene.dispersion.scaled > y.cutoff) & (gene.dispersion.scaled < y.high.cutoff))]
  if (do.plot) {
    smoothScatter(x = gene.mean, y = gene.dispersion.scaled, pch = pch.use, cex = 0.5, col = col.use, xlab = "Average expression", ylab = "Dispersion", nrpoints = Inf)
  }
  if (do.text) {
    text(x = gene.mean[pass.cutoff], y = gene.dispersion.scaled[pass.cutoff], labels = pass.cutoff, cex = cex.text.use)
  }
  hvg.info <- hvg.info[order(hvg.info$gene.dispersion, decreasing = TRUE), ]
  object$var.features <- list(assay = pass.cutoff)
  return(object)

}


#' Select number of the rank
#'
#' @param object scAI object
#' @param rangeK A predefined range of K
#'
#' @return
#' @export
#'
#' @examples
selectK <- function(object, rangeK = c(2:15)) {
  coph <- vector("double", length(rangeK))
  i <- 0
  for (k in rangeK) {
    i <- i+1
    outs <- run_scAI(object, K = k, nrun = object@options$paras$nrun,
                     s = object@options$paras$s, alpha = object@options$paras$alpha, lambda = object@options$paras$lambda, gamma = object@options$paras$gamma,
                     maxIter = object@options$paras$maxIter, stop_rule = object@options$paras$stop_rule)
    CM <- outs@cluster$consensus
    d1 <- dist(CM)
    hc <- hclust(d1, method = "average")
    d2 <- cophenetic(hc)
    coph[i] <- cor(d1, d2)
  }

  df <- data.frame(k = rangeK, Coph = coph)
  gg <- ggplot(df, aes(x = k, y = Coph)) + geom_line(size=1) +
    geom_point() +
    theme_classic() + labs(x = 'K', y='Stability score (Coph)') +
    labs(title = feature.name) + scAI_theme_opts() + theme(text = element_text(size = 10)) + labs(x = 'K', y='Stability score (Coph)') +
    theme(legend.position = "right", legend.title = NULL)
  gg
}



#' Reorder features according to the loading values
#'
#' @param W Basis matrix
#' @param cutoff Threshold of the loading values
#'
#' @return
#' @export
#'
#' @examples
reorderFeatures <- function(W, cutoff = 0.5) {
  M <- nrow(W)
  K = ncol(W)
  MW <- colMeans(W)
  vw <- apply(W, 2, sd)
  # order features
  IndexW_record <- vector("list", K)
  IndexW_record[[1]] <- which(W[, 1] > MW[1] + cutoff * VW[1])
  n <- length(IndexW_record[[1]])
  A <- c(1:M)
  c <- match(A, IndexW_record[[1]])
  A <- A[-c]
  for (j in 2:K) {
    IndexW_record[[j]] <- which(W[, j] > MW[j] + cutoff * VW[j])
    s <- 0
    for (k in 1:j - 1) {
      s <- s + length(IndexW_record[[k]])
    }
    ir2 <- match(IndexW_record[[j]], 1:s)
    IndexW_record[[j]] <- IndexW_record[[j]][-ir2]
    n <- length(IndexW_record[[j]])
    B <- union(1:s, IndexW_record[[j]])
    A <- 1:M
    c <- match(A, B)
    A <- A[-c]
    W0 <- W
    W[s + 1:s + n, ] <- W0[IndexW_record[[j]], ]
    W[s + n + 1:M, ] <- W0[A, ]
  }
  results <- list()
  list[["W"]] <- W
  list[["index_record"]] <- IndexW_record
  return(results)
}

