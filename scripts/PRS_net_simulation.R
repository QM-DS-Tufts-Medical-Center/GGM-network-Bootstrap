
fam_dat_sim_prsnet <- function(sim.size, h.sq, fam.num, out.path){
  pc.mat <- readRDS("/restricted/projectnb/necs/Zeyuan_Analysis/Network/PRS_network_LLFS/prsnet_simulation/prs.pc.mat.rds")
  y.cov.mat <- solve(pc.mat)
  lambda <- eigen(y.cov.mat)$values
  eigenv <- eigen(y.cov.mat)$vectors
  temp.m <- eigenv %*% diag(lambda^(1/2))
  
  sigma.e.mat <- (1-h.sq^2/2)*y.cov.mat
  sigma.g.mat <- h.sq/2*y.cov.mat-h.sq^2/2*y.cov.mat
  
  n.off1 <- 5
  fam1.ID <- rep(1:(fam.num/2), 2+n.off1)
  sigma.e.mat.k1 <- kronecker(diag(1, n.off1, n.off1), sigma.e.mat)
  sigma.g.mat.k1 <- kronecker(matrix(1, n.off1, n.off1)-diag(1, n.off1, n.off1), sigma.g.mat)
  sigma.cond.mat1 <- sigma.e.mat.k1 + sigma.g.mat.k1
  lambda.cond1 <- eigen(sigma.cond.mat1)$values
  eigenv.cond1 <- eigen(sigma.cond.mat1)$vectors
  temp.m.cond1 <- eigenv.cond1 %*% diag(lambda.cond1^(1/2))
  
  n.off2 <- 1
  fam2.ID <- rep(((fam.num/2)+1):fam.num, 2+n.off2)
  sigma.e.mat.k2 <- kronecker(diag(1, n.off2, n.off2), sigma.e.mat)
  sigma.g.mat.k2 <- kronecker(matrix(1, n.off2, n.off2)-diag(1, n.off2, n.off2), sigma.g.mat)
  sigma.cond.mat2 <- sigma.e.mat.k2 + sigma.g.mat.k2
  lambda.cond2 <- eigen(sigma.cond.mat2)$values
  eigenv.cond2 <- eigen(sigma.cond.mat2)$vectors
  temp.m.cond2 <- eigenv.cond2 %*% diag(lambda.cond2^(1/2))
  
  for(i in 1:sim.size){
    y_father <- t(matrix(replicate(fam.num/2, temp.m  %*% rnorm(nrow(y.cov.mat), 0,1)), nrow = nrow(y.cov.mat)))
    y_mother <- t(matrix(replicate(fam.num/2, as.vector(temp.m  %*% rnorm(nrow(y.cov.mat), 0,1))), nrow = nrow(y.cov.mat)))
    y_child_res1 <- t(matrix(replicate(fam.num/2, as.vector(temp.m.cond1  %*% rnorm(n.off1*nrow(y.cov.mat), 0,1))), nrow = n.off1*nrow(y.cov.mat)))
    y_child1 <- h.sq*(0.5*y_father + 0.5*y_mother) + y_child_res1[,1:30]
    y_child2 <- h.sq*(0.5*y_father + 0.5*y_mother) + y_child_res1[,31:60]
    y_child3 <- h.sq*(0.5*y_father + 0.5*y_mother) + y_child_res1[,61:90]
    y_child4 <- h.sq*(0.5*y_father + 0.5*y_mother) + y_child_res1[,91:120]
    y_child5 <- h.sq*(0.5*y_father + 0.5*y_mother) + y_child_res1[,121:150]
    data_out1 <- rbind(y_father, y_mother, y_child1, y_child2, y_child3, y_child4, y_child5)
    data_out_ID1 <- cbind(fam1.ID, data_out1)
    colnames(data_out_ID1) <- c("fam.ID",paste0("Y", 1:30))
    
    y_father <- t(matrix(replicate(fam.num/2, temp.m  %*% rnorm(nrow(y.cov.mat), 0,1)), nrow = nrow(y.cov.mat)))
    y_mother <- t(matrix(replicate(fam.num/2, as.vector(temp.m  %*% rnorm(nrow(y.cov.mat), 0,1))), nrow = nrow(y.cov.mat)))
    y_child_res1 <- t(matrix(replicate(fam.num/2, as.vector(temp.m.cond2  %*% rnorm(n.off2*nrow(y.cov.mat), 0,1))), nrow = n.off2*nrow(y.cov.mat)))
    y_child1 <- h.sq*(0.5*y_father + 0.5*y_mother) + y_child_res1
    data_out2 <- rbind(y_father, y_mother, y_child1)
    data_out_ID2 <- cbind(fam2.ID, data_out2)
    colnames(data_out_ID1) <- c("fam.ID",paste0("Y", 1:30))
    
    data_out_ID <- rbind(data_out_ID1, data_out_ID2)
    write.csv(data_out_ID, paste0(out.path,"/dataset",i,"_prsnet_h.sq",h.sq,".csv"), row.names = F)
  }
}


#check
#cov(y_father)
#cov(y_mother)
#cov(y_child)
#cov(y_father, y_mother)
#cov(y_father, y_child)
#cov(y_mother, y_child)
#cov(data_out)
