library(Matrix)
library(qlcMatrix)
library(data.table)
library(fastDummies)
library(future.apply)

`%ni%` <- Negate(`%in%`)

# Function for generating combinations' dummies
get.dummies <- function(nck,lst_cte,lst_fixed,lst_comp,prohibs=TRUE){
  
  if (length(nck)==1){
    colnames(nck) <- "X"
  }
  
  for (proh in prohs){
    nck <- nck[!apply(nck,1,function(comb) all(proh %in% comb)),]
  }
  
  nck <- dummy_cols(nck,select_columns=colnames(nck),
                    remove_selected_columns=TRUE)
  
  head <- unique(gsub("X[0-9]+", "", colnames(nck)))
  
  nck <- sapply(head, 
                function(x)rowSums(nck[endsWith(colnames(nck), x)])) %>% data.frame()
  
  colnames(nck) <- sapply(colnames(nck),
                          function(col) str_sub(col,start=nchar(prefix)+1))
  
  for (x in lst_fixed){
    nck[[as.name(x)]] <- 1
  }
  
  for (x in lst_comp){
    nck[[as.name(x)]] <- 1
  }
  
  for (x in setdiff(lst_cte, 
                    colnames(nck) %>% as.numeric())) {
    nck[[as.name(x)]] <- 0
  }
  
  nck <- nck[, order(colnames(nck) %>% as.numeric())]
  
  return(nck)
}

# cbind for data.frames with different numbers of rows
cbind2 <- function (...){
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n  <- max(sapply(nm, nrow))
  do.call(cbind, 
          lapply(nm, function(x) rbind(x, matrix(, n - nrow(x), ncol(x))))) %>% data.frame()
}

# Iteration of TURF Analysis
turf.iter <- function(i,nck_bin,utils,prefix,ttest.data=FALSE){
  
  comb <- nck_bin[i,]
  
  items <- copy(comb)
  items[] <- sapply(1:ncol(comb), 
                    function(i) ifelse(comb[i]==1,
                                       paste0(prefix,colnames(comb)[i]),0))
  
  items <- items %>% select(-which(comb == 0)) %>% 
    select_if(colnames(.) %ni% lst_comp) %>% select(-none)

  colnames(items) <- 1:ncol(items)
  
  case_util <- mapply(`*`, util, 
                      comb*prod_wei) %>% Matrix(data=.,sparse=TRUE)
  
  if (firstChoice) {
    
    case_max  <- rowMax(case_util)
    case_util <- apply(case_util,2,
                       function(col) ifelse(col == case_max,1,0)) %>% Matrix(data=.,
                                                                             sparse=TRUE)
    
  }
  
  case_sop <- wei * (case_util / rowSums(case_util))
  
  item_sop <- colMeans(case_sop) %>% data.frame()
  colnames(item_sop) <- i
  
  cte_sop <- sum(item_sop[lst_cte,]) %>% data.frame()
  colnames(cte_sop) <- 'SoP'
  
  output <- cbind(cte_sop,items)
  
  if (ttest.data) {
    
    ttest_sop <- rowSums(case_sop[,lst_cte]) %>% data.frame()
    colnames(ttest_sop) <- paste0('comb',i)
    
    output <- list(output,ttest_sop)
    
  } 
  
  return(output)
  
} 

# Function that runs TURF Analysis
run.turf <- function(nck_bin,utils,prefix,ttest.data=FALSE){
  
  n_comb <- nrow(nck_bin)
  results <- future_lapply(1:n_comb, 
                           function(i) turf.iter(i,nck_bin,utils,prefix,ttest.data)) 
  
  if (ttest.data) {
    
    cte_SoP   <- do.call(rbind,future_lapply(results,"[[",1))
    cte_SoP   <- cte_SoP[order(cte_SoP$SoP,decreasing=TRUE),]
    
    ttest_SoP <- do.call(cbind2,future_lapply(results,"[[",2))
    
    return(list(cte_SoP,ttest_SoP))
    
  } else {
    
    cte_SoP <- do.call(rbind,results)
    cte_SoP <- cte_SoP[order(cte_SoP$SoP,decreasing=TRUE),]
    
    return(cte_SoP)
    
  }

}

# Function for determining "optimal" k-partition
steps <- function(array,m,k,threshold=60000){
    
  # Comments
  mCk <- function(m,k) prod(setdiff(1:m, 1:(m-k))) / prod(1:k)
  
  # Comments
  if (mCk(m,k) <= threshold) {
    array <- append(array,k)
    return(array[order(array)])
  } else {
    ki <- sapply(1:(k-1), 
                 function(i) mCk(m,i)) %>% ifelse(.<=threshold,.,0) %>% which.max()
    ki <- if (k - ki > 1) ki else ki - 1=
    array <- append(array,ki)
    steps(array,m-ki,k-ki)
  }
  
}

# Function that implements greedy approach for TURF Analysis 
stepwise.turf <- function(k,utils,lst_cte,lst_fixed,lst_comp,prefix,start_from=NULL){
  
  m <- length(setdiff(lst_cte,lst_fixed))
  
  k_iter <- steps(c(),m,k) # rep(1,k)
  stepwise <- if (length(k_iter) > 1) TRUE else FALSE
  
  iter_fixed <- copy(lst_fixed)
  
  start_from <- start_from[[as.name(k_iter[1])]][[2]]
  
  if (stepwise & !is.null(start_from)) {
    
    best <- start_from[1,2:ncol(start_from)]
    iter_fixed <- sapply(best,
                         function(i) as.numeric(str_sub(i, start=nchar(prefix)+1)))
    
    k_iter <- k_iter[2:length(k_iter)]
  
  }
  
  results <- NULL
  
  for (i in seq_along(k_iter)) {
    
    elem <- setdiff(lst_cte,iter_fixed)
    
    nck_bin <- data.frame(t(combn(elem,k_iter[i]))) %>% 
      get.dummies(.,lst_cte,iter_fixed,lst_comp)
    
    nck_bin$none <- if (None) 1 else 0
    
    if (i < length(k_iter)) {
      
      results <- run.turf(nck_bin,utils,prefix)
      
      best <- results[1,2:ncol(results)]
      iter_fixed <- sapply(best,
                           function(i) as.numeric(str_sub(i, start=nchar(prefix)+1)))
      
    } else {
      
      results <- run.turf(nck_bin,utils,prefix)
      
    }
    
  }
  
  if (stepwise) {
    
    results <- results %>% 
      rbind(.,run.swapping(results,utils,lst_cte,lst_fixed,lst_comp,prefix))
    
    results <- results[!duplicated(results),]
    results <- results[order(results$SoP,decreasing=TRUE),]
    
    results <- list("stepwise+swapping",results)
      
  } else {
    
    results <- list("full-search",results)
    
  }
  
  return(results)
}

# Iteration of swapping algorithm
swap <- function(comb,lst_cte,lst_fixed,lst_comp){
  
  comb <- comb[,2:ncol(comb)]
  comb <- mapply(function(x) as.numeric(str_sub(x, start=nchar(prefix)+1) ), comb)
  
  alt <- setdiff(lst_cte,comb)
  base <- setdiff(comb,lst_fixed)
  
  output <- c()
  
  for (i in seq_along(base)) {
    new = copy(base)
    for (j in alt) {
      new[i] <- j
      output <- rbind(output,new[order(new)])
    }
  }
  
  output <- output %>% data.frame()
  rownames(output) <- 1:nrow(output)
  
  return(output)
  
}

# Function that runs swapping
run.swapping <- function(results,utils,lst_cte,lst_fixed,lst_comp,prefix){
  
  top100 <- min(nrow(results),100)
  
  swapped <- future_lapply(1:top100,
                           function(i) swap(results[i,],
                                            lst_cte,lst_fixed,lst_comp))
  
  swapped <- do.call(rbind,swapped)
  swapped <- swapped[!duplicated(swapped),]
  
  swap_bin <- swapped %>% get.dummies(.,lst_cte,lst_fixed,lst_comp)
  swap_bin$none <- if (None) 1 else 0
  
  output <- run.turf(swap_bin,utils,prefix)

  return(output)

}

# Function that performs a paired t-test over the top K combos
paired.ttest <- function(results,utils,top,alpha,prefix) {
  
  results <- results[1:min(nrow(results),top),]
  rn <- rownames(results)
  
  results <- results[,2:ncol(results)] %>%
    mapply(function(x) as.numeric(str_sub(x, start=nchar(prefix)+1)),.)
  
  results_bin <- results %>% 
    data.frame() %>% get.dummies(.,lst_cte,lst_fixed,lst_comp,prohibs=F)
  
  results_bin$none <- if (None) 1 else 0
  
  ttest_data <- run.turf(results_bin,utils,prefix,ttest.data=TRUE)
  
  output <- ttest_data[[1]]
  ttest_data <- ttest_data[[2]]
  
  p_val   <- c(1)
  sig_val <- c('==')
  
  i0 <- 1
  for (i in 2:nrow(output)) {
    
    test <- t.test(ttest_data[[as.name(paste0('comb',i0))]],
                   ttest_data[[as.name(paste0('comb',i))]], paired=TRUE)
    
    p_val[[i]] <- test$p.value
    
    if (test$p.value < alpha) {
      
      i0 <- i
      sig_val[[i]] <- '***'
      
    } else {
      
      sig_val[[i]] <- '=='
      
    }
  }
  
  output$pvalue <- p_val
  output$sigtest <- sig_val
  
  rownames(output) <- rn
  
  return(output)
  
}