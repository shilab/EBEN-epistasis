##Set up the running path where the data and the scipt locate in
setwd('~/Coloncancer/');
library('EBEN');

## Input the miRNA data and pathological data
mi <- read.table("./bc_matrix.txt",header=T);
mi <- as.matrix(mi);

## Remove the samples without pathological data
mi2 <- mi[,which(mi[nrow(mi),]!='NA')];

target <- as.matrix(mi2[nrow(mi2),(2:ncol(mi2))]);
target1 <- log(as.numeric(target),base=exp(1));
x <- mi2[(1:nrow(mi2)-1),];
x11 <- x[,2:ncol(x)];
x11 <- matrix(as.numeric(x11),nrow(x11));

## Filter the miRNA data with more than 20% missing data
x1 <- NULL;
for(i in 1:nrow(x11)){
  if(sum(as.numeric(x11[i,])!=0)){x1 <- rbind(x1,x[i,]);}
}

x2 <- NULL;
criteria <- trunc((ncol(x1)-1) *0.8);  
for(i in 1:nrow(x1)){
  if(sum(as.numeric(x1[i,(2:ncol(x1))])!=0) > criteria){
    x2 <- rbind(x2,x1[i,]);
  }
}
colnames(x2) <- colnames(mi2);
new_matrix_miRNA <- rbind(x2,c("stage",target));

## Quantile normalization:
x3 <- x2[,2:ncol(x2)];
for( sl in 1:nrow(x3) ) {
  mat = matrix(as.numeric(x3[sl,]),1);
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(x3)+1));
  x3[sl,] = mat;
}
rm(sl, mat);

## Main effect estimated using EBEN:
x4 <- matrix(as.numeric(x3),nrow = nrow(x3));
CV = EBelasticNet.GaussianCV(t(x4), target1, nFolds = 5,Epis = "no");
Blup1 = EBelasticNet.Gaussian(t(x4), target1,lambda = CV$Lambda_optimal,alpha = CV$Alpha_optimal, Epis = "no",verbose = 0);
Blup_main_sig = Blup1$weight[which(Blup1$weight[,6] <= 0.05),];

## Substract the main effect:
x5 <- t(x4);
index_main <- Blup_main_sig[,1];
effect_main <- Blup_main_sig[,3];
target_new <- as.matrix(target1) - x5[,index_main] %*% (as.matrix(effect_main));

## Epistatic effect estimated using EBEN:
CV_epis = EBelasticNet.GaussianCV(t(x4), target_new, nFolds = 5,Epis = "yes");
Blup_epis = EBelasticNet.Gaussian(t(x4), target_new,lambda =  CV_epis$Lambda_optimal,alpha = CV_epis$Alpha_optimal, Epis = "yes",verbose = 0)
Blup_epis_sig = Blup_epis$weight[which(Blup_epis$weight[,6] <= 0.05),];

## Final run:
mir <- as.matrix(t(x3));
mir <- matrix(as.numeric(mir),nrow = nrow(mir));
main_epi_miR_id = rbind(Blup_main_sig[,1:2],Blup_epis_sig[,1:2]);

new_x6 <- NULL;
for ( i in 1:nrow(main_epi_miR_id)) {
  if (main_epi_miR_id[i,1]==main_epi_miR_id[i,2]){
    new_x6 <- cbind(new_x6,mir[,main_epi_miR_id[i,1]]);
  }
  if (main_epi_miR_id[i,1]!=main_epi_miR_id[i,2]){
    col <- mir[,main_epi_miR_id[i,1]] * mir[,main_epi_miR_id[i,2]];
    new_x6 <- cbind(new_x6,col);
  }
}
new_x7 <- t(new_x6);
for( sl in 1:nrow(new_x7) ) {
  mat = matrix(as.numeric(new_x7[sl,]),1);
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(new_x7)+1));
  new_x7[sl,] = mat;
}
rm(sl, mat);
new_x8 <- t(new_x7);

CV_full = EBelasticNet.GaussianCV(new_x8, target1, nFolds = 5,Epis = "no");
Blup_full = EBelasticNet.Gaussian(new_x8, target1,lambda =  CV_full$Lambda_optimal,alpha = CV_full$Alpha_optimal, Epis = "no",verbose = 0)
Blup_full_sig =  Blup_full$weight[which(Blup_full$weight[,6] <= 0.05),];

idma <- matrix(NA,nrow = nrow(Blup_full_sig ),6);
for(i in 1:nrow(Blup_full_sig)){
  idma[i,] = c(main_epi_miR_id[Blup_full_sig[i,1],1:2],Blup_full_sig[i,3:6]);
}

## Ouput the final result including main and epistatic effect
write.table(idma,"idma05",quote=F,sep='\t',row.names = F,col.names = F);
