reduceDimension_old <- function(cds,
                                max_components=2,
                                reduction_method=c("DDRTree", "ICA", 'tSNE', "SimplePPT", 'L1-graph', 'SGL-tree'),
                                norm_method = c("log", "vstExprs", "none"),
                                residualModelFormulaStr=NULL,
                                pseudo_expr=1,
                                relative_expr=TRUE,
                                auto_param_selection = TRUE,
                                verbose=FALSE,
                                scaling = TRUE,
                                ...){
  extra_arguments <- list()
  set.seed(2016) #ensure results from RNG sensitive algorithms are the same on all calls
  
  FM <- normalize_expr_data(cds, norm_method, pseudo_expr)
  
  #FM <- FM[unlist(sparseApply(FM, 1, sd, convert_to_dense=TRUE)) > 0, ]
  xm <- Matrix::rowMeans(FM)
  xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
  FM <- FM[xsd > 0,]
  
  if (is.null(residualModelFormulaStr) == FALSE) {
    if (verbose)
      message("Removing batch effects")
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr),
                                       data = pData(cds), drop.unused.levels = TRUE)
    
    fit <- limma::lmFit(FM, X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
  }else{
    X.model_mat <- NULL
  }
  
  if(scaling){
    FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
    FM <- FM[!is.na(row.names(FM)), ]
  } else FM <- as.matrix(FM)
  
  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }
  
  FM <- FM[apply(FM, 1, function(x) all(is.finite(x))), ] #ensure all the expression values are finite values
  if (is.function(reduction_method)) {
    reducedDim <- reduction_method(FM, ...)
    colnames(reducedDim) <- colnames(FM)
    reducedDimW(cds) <- as.matrix(reducedDim)
    reducedDimA(cds) <- as.matrix(reducedDim)
    reducedDimS(cds) <- as.matrix(reducedDim)
    reducedDimK(cds) <- as.matrix(reducedDim)
    dp <- as.matrix(dist(reducedDim))
    cellPairwiseDistances(cds) <- dp
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- dp_mst
    cds@dim_reduce_type <- "function_passed"
  }else{
    reduction_method <- match.arg(reduction_method)
    if (reduction_method == "tSNE") {
      #first perform PCA
      if (verbose)
        message("Remove noise by PCA ...")
      
      # # Calculate the variance across genes without converting to a dense
      # # matrix:
      # FM_t <- Matrix::t(FM)
      # cell_means <- Matrix::rowMeans(FM_t)
      # cell_vars <- Matrix::rowMeans((FM_t - cell_means)^2)
      # Filter out genes that are constant across all cells:
      #genes_to_keep <- expression_vars > 0
      #FM <- FM[genes_to_keep,]
      #expression_means <- expression_means[genes_to_keep]
      #expression_vars <- expression_vars[genes_to_keep]
      # Here✬s how to take the top PCA loading genes, but using
      # sparseMatrix operations the whole time, using irlba.
      
      
      if("num_dim" %in% names(extra_arguments)){ #when you pass pca_dim to the function, the number of dimension used for tSNE dimension reduction is used
        num_dim <- extra_arguments$num_dim #variance_explained
      }
      else{
        num_dim <- 50
      }
      
      FM <- (FM)
      irlba_res <- prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1),
                                center = TRUE, scale. = TRUE)
      irlba_pca_res <- irlba_res$x
      
      # irlba_res <- irlba(FM,
      #                    nv=min(num_dim, min(dim(FM)) - 1),
      #                        nu=0,
      #                        center=cell_means,
      #                        scale=sqrt(cell_vars),
      #                        right_only=TRUE)
      # irlba_pca_res <- irlba_res$v
      # row.names(irlba_pca_res) <- genes_to_keep
      
      # pca_res <- prcomp(t(FM), center = T, scale = T)
      # std_dev <- pca_res$sdev
      # pr_var <- std_dev^2
      # prop_varex <- pr_var/sum(pr_var)
      # prop_varex <- irlba_res$sdev^2 / sum(irlba_res$sdev^2)
      
      topDim_pca <- irlba_pca_res#[, 1:num_dim]
      
      # #perform the model formula transformation right before tSNE:
      # if (is.null(residualModelFormulaStr) == FALSE) {
      #   if (verbose)
      #     message("Removing batch effects")
      #   X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr),
      #                                      data = pData(cds), drop.unused.levels = TRUE)
      
      #   fit <- limma::lmFit(topDim_pca, X.model_mat, ...)
      #   beta <- fit$coefficients[, -1, drop = FALSE]
      #   beta[is.na(beta)] <- 0
      #   topDim_pca <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
      # }else{
      #   X.model_mat <- NULL
      # }
      
      #then run tSNE
      if (verbose)
        message("Reduce dimension by tSNE ...")
      
      tsne_res <- Rtsne::Rtsne(as.matrix(topDim_pca), dims = max_components, pca = F,...)
      
      tsne_data <- tsne_res$Y[, 1:max_components]
      row.names(tsne_data) <- colnames(tsne_data)
      
      reducedDimA(cds) <- t(tsne_data) #this may move to the auxClusteringData environment
      
      #set the important information from densityClust to certain part of the cds object:
      cds@auxClusteringData[["tSNE"]]$pca_components_used <- num_dim
      cds@auxClusteringData[["tSNE"]]$reduced_dimension <- t(topDim_pca)
      #cds@auxClusteringData[["tSNE"]]$variance_explained <- prop_varex
      
      cds@dim_reduce_type <- "tSNE"
    }
    
    else if (reduction_method == "ICA") {
      # FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
      # FM <- FM[!is.na(row.names(FM)), ]
      
      if (verbose)
        message("Reducing to independent components")
      init_ICA <- ica_helper(Matrix::t(FM), max_components,
                             use_irlba = TRUE, ...)
      x_pca <- Matrix::t(Matrix::t(FM) %*% init_ICA$K)
      W <- Matrix::t(init_ICA$W)
      weights <- W
      A <- Matrix::t(solve(weights) %*% Matrix::t(init_ICA$K))
      colnames(A) <- colnames(weights)
      rownames(A) <- rownames(FM)
      S <- weights %*% x_pca
      rownames(S) <- colnames(weights)
      colnames(S) <- colnames(FM)
      reducedDimW(cds) <- as.matrix(W)
      reducedDimA(cds) <- as.matrix(A)
      reducedDimS(cds) <- as.matrix(S)
      reducedDimK(cds) <- as.matrix(init_ICA$K)
      adjusted_S <- Matrix::t(reducedDimS(cds))
      dp <- as.matrix(dist(adjusted_S))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "ICA"
    }
    else if (reduction_method == "DDRTree") {
      # FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
      # FM <- FM[!is.na(row.names(FM)), ]
      
      if (verbose)
        message("Learning principal graph with DDRTree")
      
      # TODO: DDRTree should really work with sparse matrices.
      if(auto_param_selection & ncol(cds) >= 100){
        if("ncenter" %in% names(extra_arguments)) #avoid overwrite the ncenter parameter
          ncenter <- extra_arguments$ncenter
        else
          ncenter <- cal_ncenter(ncol(FM))
        #add other parameters...
        ddr_args <- c(list(X=FM, dimensions=max_components, ncenter=ncenter, verbose = verbose),
                      extra_arguments[names(extra_arguments) %in% c("initial_method", "maxIter", "sigma", "lambda", "param.gamma", "tol")])
        #browser()
        ddrtree_res <- do.call(DDRTree, ddr_args)
      } else{
        ddrtree_res <- DDRTree(FM, max_components, verbose = verbose, ...)
      }
      if(ncol(ddrtree_res$Y) == ncol(cds))
        colnames(ddrtree_res$Y) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      else
        colnames(ddrtree_res$Y) <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      
      colnames(ddrtree_res$Z) <- colnames(FM)
      reducedDimW(cds) <- ddrtree_res$W
      reducedDimS(cds) <- ddrtree_res$Z
      reducedDimK(cds) <- ddrtree_res$Y
      cds@auxOrderingData[["DDRTree"]]$objective_vals <- ddrtree_res$objective_vals
      
      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "DDRTree"
      cds <- findNearestPointOnMST(cds)
    }else {
      stop("Error: unrecognized dimensionality reduction method")
    }
  }
  cds
}







# Helper function to normalize the expression data prior to dimensionality
# reduction
normalize_expr_data <- function(cds,
                                norm_method = c("log", "vstExprs", "none"),
                                pseudo_expr = 1,
                                relative_expr = TRUE){
  FM <- exprs(cds)
  use_for_ordering <- NULL
  # If the user has selected a subset of genes for use in ordering the cells
  # via setOrderingFilter(), subset the expression matrix.
  # if (is.null(fData(cds)$use_for_ordering) == FALSE &&
  #     nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
  #   FM <- FM[fData(cds)$use_for_ordering, ]
  # }
  # 
  # norm_method <- match.arg(norm_method)
  # if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
  #   
  #   # If we're going to be using log, and the user hasn't given us a pseudocount
  #   # set it to 1 by default.
  #   if (is.null(pseudo_expr)){
  #     if(norm_method == "log")
  #       pseudo_expr = 1
  #     else
  #       pseudo_expr = 0
  #   }
  #   
  #   checkSizeFactors(cds)
  #   
  #   if (norm_method == "vstExprs") {
  #     if (relative_expr == FALSE)
  #       message("Warning: relative_expr is ignored when using norm_method == 'vstExprs'")
  #     
  #     if (is.null(fData(cds)$use_for_ordering) == FALSE &&
  #         nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
  #       VST_FM <- vstExprs(cds[fData(cds)$use_for_ordering,], round_vals = FALSE)
  #     }else{
  #       VST_FM <- vstExprs(cds, round_vals = FALSE)
  #     }
  #     
  #     if (is.null(VST_FM) == FALSE) {
  #       FM <- VST_FM
  #     }
  #     else {
  #       stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
  #     }
  #   }else if (norm_method == "log") {
  #     # If we are using log, normalize by size factor before log-transforming
  #     
  #     if (relative_expr)
  #       FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))
  #     
  #     if(is.null(pseudo_expr))
  #       pseudo_expr <- 1 
  #     FM <- FM + pseudo_expr
  #     FM <- log2(FM)
  #   }else if (norm_method == "none"){
  #     # If we are using log, normalize by size factor before log-transforming
  #     FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))
  #     FM <- FM + pseudo_expr
  #   }
  # }else if (cds@expressionFamily@vfamily == "binomialff") {
  #   if (norm_method == "none"){
  #     #If this is binomial data, transform expression values into TF-IDF scores.
  #     ncounts <- FM > 0
  #     ncounts[ncounts != 0] <- 1
  #     FM <- Matrix::t(Matrix::t(ncounts) * log(1 + ncol(ncounts)/rowSums(ncounts)))
  #   }else{
  #     stop("Error: the only normalization method supported with binomial data is 'none'")
  #   }
  # }else if (cds@expressionFamily@vfamily == "Tobit") {
  #   FM <- FM + pseudo_expr
  #   if (norm_method == "none"){
  #     
  #   }else if (norm_method == "log"){
  #     FM <- log2(FM)
  #   }else{
  #     stop("Error: the only normalization methods supported with Tobit-distributed (e.g. FPKM/TPM) data are 'log' (recommended) or 'none'")
  #   }
  # }else if (cds@expressionFamily@vfamily == "uninormal") {
  #   if (norm_method == "none"){
  #     FM <- FM + pseudo_expr
  #   }else{
  #     stop("Error: the only normalization method supported with gaussian data is 'none'")
  #   }
  # }
  # if(norm_method != "none")
  #normalize_expr_data
  return (FM)
}


cal_ncenter <- function(ncells, ncells_limit = 100){
  round(2 * ncells_limit * log(ncells)/ (log(ncells) + log(ncells_limit)))
}




extract_ddrtree_ordering <- function(cds, root_cell, verbose=T)
{
  
  dp <- cellPairwiseDistances(cds)
  dp_mst <- minSpanningTree(cds)
  
  curr_state <- 1
  
  res <- list(subtree = dp_mst, root = root_cell)
  
  states = rep(1, ncol(dp))
  names(states) <- V(dp_mst)$name
  
  pseudotimes = rep(0, ncol(dp))
  names(pseudotimes) <- V(dp_mst)$name
  
  parents = rep(NA, ncol(dp))
  names(parents) <- V(dp_mst)$name
  
  mst_traversal <- graph.dfs(dp_mst,
                             root=root_cell,
                             neimode = "all",
                             unreachable=FALSE,
                             father=TRUE)
  mst_traversal$father <- as.numeric(mst_traversal$father)
  curr_state <- 1
  
  for (i in 1:length(mst_traversal$order)){
    curr_node = mst_traversal$order[i]
    curr_node_name = V(dp_mst)[curr_node]$name
    
    if (is.na(mst_traversal$father[curr_node]) == FALSE){
      parent_node = mst_traversal$father[curr_node]
      parent_node_name = V(dp_mst)[parent_node]$name
      parent_node_pseudotime = pseudotimes[parent_node_name]
      parent_node_state = states[parent_node_name]
      curr_node_pseudotime = parent_node_pseudotime + dp[curr_node_name, parent_node_name]
      if (degree(dp_mst, v=parent_node_name) > 2){
        curr_state <- curr_state + 1
      }
    }else{
      parent_node = NA
      parent_node_name = NA
      curr_node_pseudotime = 0
    }
    
    curr_node_state = curr_state
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_node_state
    parents[curr_node_name] <- parent_node_name
  }
  
  ordering_df <- data.frame(sample_name = names(states),
                            cell_state = factor(states),
                            pseudo_time = as.vector(pseudotimes),
                            parent = parents)
  row.names(ordering_df) <- ordering_df$sample_name
  # ordering_df <- plyr::arrange(ordering_df, pseudo_time)
  return(ordering_df)
}
