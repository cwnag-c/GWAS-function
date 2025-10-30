library(data.table)
library(magrittr)
library(tidyverse)
cat("*******************************************************************\n")
cat("* CWAS \n")
cat("* Version 1.8\n")
cat("* Update Date:2024.07.17\n")
cat("* (C) 2024 Chao Wang \n")
cat("* Nanjing Medical University\n")
cat("*******************************************************************\n")

#----------------1. Cross-trait:PLACO-------------
GWAS_placo<-function(sumstats1,sumstats2,
                     GWAS_txt_dir,
                     gwas_prefix,
                     out_path,
                     PLACO_path,
                     ncore,
                     rm_MHC=T){
  
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.5\n")
  cat("* Update Date:2024.07.07\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  library(data.table)
  library(magrittr)
  library(tidyverse)
  library(ggplot2)
  library(readxl)
  library(parallel)
  code<-fread(GWAS_txt_dir)
  load(PLACO_path)#####载入PLACO
  
  f1<-paste0(gwas_prefix,code[code$abr==sumstats1,absolute_dir])
  f2<-paste0(gwas_prefix,code[code$abr==sumstats2,absolute_dir])
  f1<-data.table::fread(f1)
  f2<-data.table::fread(f2)
  f12<-inner_join(f1,f2,by=c("SNP","CHR","POS","A1","A2"),suffix=c(".f1", ".f2"))
  rm(f1)
  rm(f2)
  gc()
  
  ###增加Z_col,并剔除Z^2>80
  f12 %<>%
    dplyr::mutate(Z.f1=BETA.f1/SE.f1) %>% 
    dplyr::mutate(Z.f2=BETA.f2/SE.f2) %>%   
    dplyr::mutate(Z2.f1=Z.f1^2) %>% 
    dplyr::mutate(Z2.f2=Z.f2^2) %>%  
    dplyr::filter(Z2.f1<80) %>% 
    dplyr::filter(Z2.f2<80)
  
  ###剔除MHC(chr6 25000000-35000000)
  if(rm_MHC){
    f12 %<>% dplyr::filter(!(CHR==6&POS>25000000&POS<35000000)) 
  }
  f12 %<>% 
    dplyr::select(SNP,CHR,POS,A1,A2,
                  BETA.f1,Z.f1,P.f1,N.f1,EAF.f1,
                  BETA.f2,Z.f2,P.f2,N.f2,EAF.f1)
  fwrite(f12,
         file = paste0(out_path,"/",sumstats1, "-",sumstats2,".txt.gz"),
         compress = "gzip",
         sep = "\t",
         quote = F,
         col.names = T,
         row.names = F)
  
  ####PLACO#####
  k <- 2 
  p <-nrow(f12)
  Z.matrix<-as.matrix(f12[,c("Z.f1","Z.f2")])
  P.matrix<-as.matrix(f12[,c("P.f1","P.f2")])
  colnames(Z.matrix) <- paste("Z",1:k,sep="")
  colnames(P.matrix) <- paste("P",1:k,sep="")
  
  
  #### decorrelating the Z-scores
  R <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)###计算Z的相关度
  "%^%" <- function(x, pow) with(eigen(x), vectors %*% (values^pow * t(vectors)))####创建公式
  Z.matrix.decor <- Z.matrix %*% (R %^% (-0.5)) ###decor
  colnames(Z.matrix.decor) <- paste("Z",1:k,sep="")
  
  VarZ <- var.placo(Z.matrix.decor, P.matrix, p.threshold=1e-4)
  
  ####删除不相关的
  rm(f12)
  rm(Z.matrix)
  rm(P.matrix)
  rm(R)
  rm(k)
  rm("%^%")
  rm(cor.pearson)
  rm(var.placo)
  gc()
  
  #### 保存函数内部环境的对象
  objects_to_save <- list(Z.matrix.decor = Z.matrix.decor, VarZ = VarZ, p = p, sumstats1 = sumstats1, sumstats2 = sumstats2)
  save(objects_to_save, file = paste0(out_path,"/",sumstats1, "-", sumstats2, ".RData"))
  
  print(paste(m, n, sep = "-"))
  
  # 设置多核并行计算
  cl <- makeCluster(ncore)
  clusterExport(cl, c("sumstats1","sumstats2","PLACO_path","out_path"))
  clusterEvalQ(cl, {
    library(data.table)
    library(magrittr)
    library(tidyverse)
    load(paste0(out_path,"/",sumstats1, "-", sumstats2, ".RData"))
    load(PLACO_path)
  })
  out <- parSapply(cl, 1:p, function(i) placo(Z = Z.matrix.decor[i,], VarZ = VarZ))
  stopCluster(cl)
  
  ####整合结果到原始文件
  f12<-fread(paste0(out_path,"/",sumstats1, "-",sumstats2,".txt.gz"))
  out %>% 
    t() %>% 
    as.data.table() %>% 
    mutate(p.placo=as.numeric(p.placo)) %>% 
    mutate(T.placo=as.numeric(T.placo)) %>% 
    cbind(f12)->f12
  f12 %>% filter(p.placo<5*10^-8)->f12_sig###输出显著的点
  fwrite(f12,
         file = paste0(sumstats1, "-",sumstats2,".txt.gz"),
         compress = "gzip",
         sep = "\t",
         quote = F,
         col.names = T,
         row.names = F)
  write.table(f12_sig,
              file = paste0(sumstats1, "-",sumstats2,"_sig.txt"),
              sep = "\t",
              quote = F,
              col.names = T,
              row.names = F)
  
  rm(f12)
  rm(f12_sig)
  rm(out)
  gc()
  print(paste0(sumstats1, "-",sumstats2))  
} 

