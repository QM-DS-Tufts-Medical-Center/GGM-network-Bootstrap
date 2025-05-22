library(rlist)

clusterBootstrap <- function(dat, sample.id = "sample.id", cluster.id = "cluster.id", boot.size=1000, cluster.sample.ratio=1, replace=T, alpha=0.05, rho=0){
  #### data should have only two ID columns and the variables(vertices) of interest
  
  #### sample IDs and cluster IDs are mandatory as "sample.id", "cluster.id"
  if(sum(colnames(dat) %in% c(sample.id, cluster.id))!=2){
    return("error: sample.id or cluster.id missing")
  }else{
    
    n.vertice = ncol(dat)-2
    n.edges <- (n.vertice-1)*n.vertice/2
    
    #calculate Ri,j
    X.mat <- as.data.frame(dat)[,!names(dat)%in%c(sample.id, cluster.id)]
    S.mat <- cov(X.mat)
    K.mat <- solve(S.mat)
    y <- diag(sqrt(1/diag(K.mat)))
    PC.mat <- -y%*%K.mat%*%y + 2*diag(1, nrow(y))
    lower.tri.vec <- PC.mat[lower.tri(PC.mat)]
    vec.R <- sqrt(nrow(dat)-n.vertice-1)*0.5*log((1+lower.tri.vec)/(1-lower.tri.vec))
    
    #### sample family ids for each bootstrap sample
    resample_fam_id <- lapply(1:boot.size, function(i){sample(unique(dat$cluster.id), length(unique(dat$cluster.id))*cluster.sample.ratio, replace = replace)})
    
    #### generate sample set for each Bootstrap iteration
    boot_samples <- list()
    for( i in resample_fam_id ){
      index_list <- rle(sort(i))
      boot_sample <- c()
      while(max(index_list$lengths)>=1){
        index <- index_list$values[index_list$lengths>=1]
        boot_sample <- rbind(boot_sample, dat[which( dat$cluster.id %in% index ),])
        index_list$lengths <- index_list$lengths-1
      }
      boot_samples <- list.append(boot_samples, boot_sample)
    }
    
    z.mat.boot <- c()
    p.mat.boot <- c()
    for(i in boot_samples){
      V.mat <- as.data.frame(i)[,!names(dat)%in%c(sample.id, cluster.id)]
      S <- cov(V.mat) 
      pc <- solve(S)
      y <- diag(sqrt(1/diag(pc)))
      PC.boot <- y%*%pc%*%y
      
      lower.tri.vec <- PC.boot[lower.tri(PC.boot)]
      vec.z <- sqrt(nrow(i)-n.vertice-1)*0.5*log((1+lower.tri.vec)/(1-lower.tri.vec))
      
      z.mat.boot <- rbind(z.mat.boot, vec.z)
    }
    z.var.vec.boot <- apply(z.mat.boot, 2, sd)
    
    p.edges <- c()
    for(r in rho){
      z.rho <- 0.5*log((1+r)/(1-r))*sqrt(nrow(dat)-n.vertice-1)

      if(r==0){
        p.edges <- cbind(p.edges, 2*(pnorm( -abs(vec.R)/z.var.vec.boot) ) )
      }else{
        p.edges <- cbind(p.edges, pnorm( (z.rho-abs(vec.R))/z.var.vec.boot )+pnorm((-abs(vec.R)-z.rho)/z.var.vec.boot) )
      }
    }
    
    #### estimates of the partial correlations
    p.cor.mat <- PC.mat
    
    edges <- apply(p.edges ,2 ,function(x){ifelse(x > alpha, 0, 1)})
    
    num.edges <- apply(edges, 2, sum)
    
    v.names <- names(dat)[!names(dat)%in%c(sample.id, cluster.id)]
    
    adj.mat.low <- matrix(0,n.vertice,n.vertice)  
    adj.mat.rho <- list()

    pval.mat.low <- matrix(0,n.vertice,n.vertice)  
    pval.mat.rho <- list()
    
    for(r in rho){
      adj.mat.low[lower.tri(adj.mat.low, diag=FALSE)] <- as.numeric(edges[,which(rho==r)])
      adj.mat.up <- t(adj.mat.low)
      adj.mat <- adj.mat.low + adj.mat.up + diag(1, n.vertice)
      rownames(adj.mat) <- v.names
      colnames(adj.mat) <- v.names
      adj.mat.rho <- list.append(adj.mat.rho, adj.mat)
      
      pval.mat.low[lower.tri(pval.mat.low, diag=FALSE)] <- as.numeric(p.edges[,which(rho==r)])
      pval.mat.up <- t(pval.mat.low)
      pval.mat <- pval.mat.low + pval.mat.up + diag(NA, n.vertice)
      rownames(pval.mat) <- v.names
      colnames(pval.mat) <- v.names
      pval.mat.rho <- list.append(pval.mat.rho, pval.mat)
    }
    
    rownames(p.cor.mat) <- v.names
    colnames(p.cor.mat) <- v.names
    
    out.list <- list(rho, num.edges, p.cor.mat, adj.mat.rho, pval.mat.rho)
    names(out.list) <- c("rho", "num.edges", "p.cor.mat", "adj.mat.rho", "pval.mat.rho")
    
    return(out.list)
  }
}

#check edges:  plot(rho, out.list$num.edges)
