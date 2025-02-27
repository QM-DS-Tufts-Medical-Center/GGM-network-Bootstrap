library(rlist)
library(rootSolve)

famBootstrap <- function(file.path, n.vertice, boot.size=1000, fam.sample.ratio=1, replace=T, rho=0, alpha=0.05){
  datasets <- list.files(file.path)
  R.mat <- c()
  z.var.mat <- c()
  vec.list <- list()
  for(k in datasets){
    #print(k)
    dat <- read.csv(paste0(file.path, "/", k))
    
    #calculate Ri,j
    X.mat <- dat[,!names(dat)%in%"fam.ID"]
    S.mat <- cov(X.mat)
    K.mat <- solve(S.mat)
    y <- diag(sqrt(1/diag(K.mat)))
    PC.mat <- y%*%K.mat%*%y
    lower.tri.vec <- PC.mat[lower.tri(PC.mat)]
    vec.R <- sqrt(nrow(dat)-n.vertice-1)*0.5*log((1+lower.tri.vec)/(1-lower.tri.vec))
    
    ####sample family ids for each bootstrap sample
    resample_fam_id <- lapply( 1:boot.size, function(i)  sample( unique(dat$fam.ID), length(unique(dat$fam.ID))*fam.sample.ratio, replace = replace ) )
    
    ####generate all bootstrap samples by drawing all subject from family ids
    ##each list element is a bootstrap sample
    
    boot_samples <- list()
    for( i in resample_fam_id ){
      index_list <- rle(sort(i))
      boot_sample <- c()
      while(max(index_list$lengths)>=1){
        index <- index_list$values[index_list$lengths>=1]
        boot_sample <- rbind(boot_sample, dat[which( dat$fam.ID %in% index ),])
        index_list$lengths <- index_list$lengths-1
      }
      boot_samples <- list.append(boot_samples, boot_sample)
    }
    
    z.mat.boot <- c()
    n.edges <- (n.vertice-1)*n.vertice/2
    for(i in boot_samples){
      X.mat <- i[,!names(dat)%in%"fam.ID"]
      S <- cov(X.mat) #S <- pcor(X.mat)$estimate
      pc <- solve(S)
      y <- diag(sqrt(1/diag(pc)))
      PC.boot <- y%*%pc%*%y
      lower.tri.vec <- PC.boot[lower.tri(PC.boot)]
      vec.z <- sqrt(nrow(i)-n.vertice-1)*0.5*log((1+lower.tri.vec)/(1-lower.tri.vec))
      z.mat.boot <- rbind(z.mat.boot, vec.z)
    }
    z.var.vec.boot <- apply(z.mat.boot, 2, sd)
    R.mat <- rbind(R.mat, vec.R)
    z.var.mat <- rbind(z.var.mat, z.var.vec.boot)
  }
  z.rho <- 0.5*log((1+rho)/(1-rho))*sqrt(nrow(dat)-n.vertice-1)
  
  if(rho==0){
    ratio <- apply(abs(R.mat)/z.var.mat > qnorm(1-alpha/2), 2, sum)/length(datasets)
  }else{
    ratio <- apply(pnorm( (z.rho-abs(R.mat))/z.var.mat )+pnorm((-abs(R.mat)-z.rho)/z.var.mat) < alpha, 2, sum)/length(datasets)
  }
}

#example
#file.path <- "/restricted/projectnb/necs/Zeyuan_Analysis/Network/fam_data_simulation/dataset_fam_sim_p_o_pairs_hsq_function_Dec26/fam250_h.sq0.75_4v4e_5/"
#r0=famBootstrap(file.path, n.vertice=4, boot.size=100, fam.sample.ratio=1, rho=0, alpha=0.05)
#r0.1=famBootstrap(file.path, n.vertice=4, boot.size=100, fam.sample.ratio=1, rho=0.1, alpha=0.05)
#r0.3=famBootstrap(file.path, n.vertice=4, boot.size=100, fam.sample.ratio=1, rho=0.316, alpha=0.05)
#r0.5=famBootstrap(file.path, n.vertice=4, boot.size=100, fam.sample.ratio=1, rho=0.5, alpha=0.05)
#r0.78=famBootstrap(file.path, n.vertice=4, boot.size=100, fam.sample.ratio=1, rho=0.78, alpha=0.05)
