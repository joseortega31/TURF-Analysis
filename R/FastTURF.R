# Fast TURF Analysis (based on of e.grant's Python tool)
# Author: j.ortega

# 0. Setup ---------------------------------------------------------------------

#rm(list = ls())

library(tictoc)
library(future)
library(tidyverse)

tic.clearlog()
plan(multisession)

if ( !(file.exists('output')) ){
  dir.create('output')
}

source('functions.R')

# 1. Project inputs ------------------------------------------------------------

proj <- 'Remco' #------------------------------------------------ Project Number
prefix <- 'KV' #-------------------------------------------------- KV / C / SKU
utilities <- 'utilities.csv' #--------------------------- Name of utilities file

n <- 19 #--------------------------------- Total number of KVs / Concepts / SKUs
k <- 3 #----------------------- Number of elements to be drawn (EXCLUDES FIXED)
top <- 100 #-------------------------------------------------- Number of results 
None <- TRUE #------------------------------------------ None option? FALSE = no
firstChoice <- FALSE #--------------------------------- First choice? FALSE = no

sigtest <- 0.05 #----------------------------- Significance (either 0.05 or 0.1)
corrmap <- FALSE #------------------------- Generate correlation map? FALSE = no

# 2. Buckets & restrictions ----------------------------------------------------

lst_fixed <<- c()
len_fixed <- length(lst_fixed)

# TODO: buckets as csv inputs
lst_comp <<- c(17,18,19) # -------------------- Competitor KVs / Concepts / SKUs

n_buks <- 0 # ----------------------------------------- Number of client buckets
#buk_1 <- 1:16 # ------------------------------------------------- Bucket name 1

# Capping of buckets on line per bucket
# Maximum number of concepts per bucket allowed - default = k
#lim_buk1 <- n

# TODO: restrictions as csv inputs
n_rest <- 0 # ------------------------------------- Total number of restrictions 
# Example: rest_1 = list(buk_1, buk_2)

# TODO: prhobitions as csv inputs
n_prohs <- 0
proh_1 <- c(1,9,16)

lst <- 1:n
head <- append(lst,'none')
lst_cte <- setdiff(unique(lst),unique(lst_comp))

buks <- lapply(seq_len(n_buks), function(i) get(paste0('buk_',i)))
lims <- lapply(seq_len(n_buks), function(i) get(paste0('lim_buk',i)))
prohs <- lapply(seq_len(n_prohs), function(i) get(paste0('proh_',i))) 

prod_wei  <- rep(1,n+1) %>% data.frame() %>% t()
colnames(prod_wei) <- head

# 3. Read and clean utilities --------------------------------------------------

util <- read.csv(utilities)

wei  <- util$wei

util <- util %>% select(-wei, -id)
colnames(util) <- head

util <- exp(util)

if (corrmap == TRUE){
  corr <- cor(util)
  write.csv(corr,'output/CorrMatrix.csv')
  heatmap(corr,Colv=NA,Rowv=NA)
}

# 4. Run fast TURF analysis ----------------------------------------------------

if (length(k) == 1) {
  
  tic(k)
  
  results <- stepwise.turf(k,util,lst_cte,
                           lst_fixed,lst_comp,prefix)
  
  toc(log=TRUE, quiet=TRUE)
  
  method <- results[[1]]
  results <- results[[2]] 
  
} else {
  
  results <- lapply(k, function(i) NULL)
  names(results) <- k
  
  for (ki in k) {
    
    tic(ki)
    
    results[[ki]] <- stepwise.turf(ki,util,lst_cte,lst_fixed,
                                   lst_comp,prefix,start_from=results)
    
    toc(log=TRUE, quiet=TRUE)
    
  }
  
  method  <- sapply(results,"[[",1)
  results <- sapply(results,"[[",2)
  
}

t <- tic.log(format=F) %>% 
  lapply(., function(k) k$toc - k$tic) %>% unlist()

print(paste0('TURF for ',(k+len_fixed),' items using ',
             method,' method took: ',format(t,digits=3),' secs'))

# 5. Perform t-test for top combs ---------------------------------------------

t <- Sys.time()

if (class(results) == 'list') {
  
  results <- future_sapply(results, 
                           function(k) paired.ttest(k,util,top,sigtest,prefix))
  for (i in k){
    write.csv(results[[i]],
              file=paste0("./output/",proj,"_TURF",i,".csv"))
  }
  
} else {
  
  results <- paired.ttest(results,util,top,sigtest,prefix)
  write.csv(results,file=paste0("./output/",proj,"_TURF",k,".csv"))
  
}


t <- difftime(Sys.time(),t,units="secs")

print(paste0('Paired ttest for top 100 portfolios took ',format(t,digits=3)))