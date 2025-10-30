library(data.table)
library(magrittr)
library(tidyverse)
cat("*******************************************************************\n")
cat("* CWAS \n")
cat("* Version 1.2\n")
cat("* Update Date:2025.10.29\n")
cat("* (C) 2025 Chao Wang \n")
cat("* Nanjing Medical University\n")
cat("*******************************************************************\n")

#----------------1. Cross-trait:PLACO-------------
GWAS_placo<-function(sumstats1_name,sumstats2_name,
                     sumstats1_dir=NULL,sumstats2_dir=NULL,
                     GWAS_txt_dir=NULL,
                     gwas_prefix=NULL,
                     out_path="./",
                     rscript_dir,
                     Placo.plus=F,
                     ncore=8,
                     rm_MHC=T){
  
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 2.0\n")
  cat("* Update Date:2025.10.29\n")
  cat("* (C) 2025 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  ################process###################
  source(rscript_dir)
  library(data.table)
  library(magrittr)
  library(tidyverse)
  library(readxl)
  library(parallel)
  if(!is.null(GWAS_txt_dir)){##GWAS_txt_dir模式
    code<-fread(GWAS_txt_dir)
    f1<-paste0(gwas_prefix,code[code$abr==sumstats1_name,absolute_dir])
    f2<-paste0(gwas_prefix,code[code$abr==sumstats2_name,absolute_dir])
    f1<-fread(f1)
    f2<-fread(f2)
  }else{##直接输入dir模式
    f1<-fread(sumstats1_dir)
    f2<-fread(sumstats2_dir)
  }
  ##########看列名##
  l1<-all(c("SNP","CHR","POS","A1","A2","BETA","SE","P","N","EAF") %in% colnames(f1))
  l2<-all(c("SNP","CHR","POS","A1","A2","BETA","SE","P","N","EAF") %in% colnames(f2))
  z1<-all(c("SNP","CHR","POS","A1","A2","Z","P","N","EAF") %in% colnames(f1))
  z2<-all(c("SNP","CHR","POS","A1","A2","Z","P","N","EAF") %in% colnames(f2))
  f12<-inner_join(f1,f2,by=c("SNP","CHR","POS","A1","A2"),suffix=c(".f1", ".f2"))
  rm(f1);rm(f2)
  
  ###剔除MHC(chr6 25000000-35000000)
  if(rm_MHC){
    f12 %<>% dplyr::filter(!(CHR==6&POS>25000000&POS<35000000)) 
  }
  ###增加Z_col,并剔除Z^2>80
  if(l1 & l2){
    f12 %<>%
      dplyr::mutate(Z.f1=BETA.f1/SE.f1) %>% 
      dplyr::mutate(Z.f2=BETA.f2/SE.f2) %>%   
      dplyr::mutate(Z2.f1=Z.f1^2) %>% 
      dplyr::mutate(Z2.f2=Z.f2^2) %>%  
      dplyr::filter(Z2.f1<80) %>% 
      dplyr::filter(Z2.f2<80) %>% 
      dplyr::select(SNP,CHR,POS,A1,A2,
                    BETA.f1,SE.f1,Z.f1,P.f1,N.f1,EAF.f1,
                    BETA.f2,SE.f2,Z.f2,P.f2,N.f2,EAF.f2)
  }else{
    if(z1 & z2){
      f12 %<>%
        dplyr::mutate(Z2.f1=Z.f1^2) %>% 
        dplyr::mutate(Z2.f2=Z.f2^2) %>%  
        dplyr::filter(Z2.f1<80) %>% 
        dplyr::filter(Z2.f2<80) %>% 
        dplyr::select(SNP,CHR,POS,A1,A2,
                      Z.f1,P.f1,N.f1,EAF.f1,
                      Z.f2,P.f2,N.f2,EAF.f2)
    }else{
      stop("缺失列")
    }
  }
  fwrite(f12,
         file = paste0(out_path,"/",sumstats1_name, "-",sumstats2_name,".txt.gz"),
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
  
  VarZ <- var.placo(Z.matrix, P.matrix)
  CorZ <- cor.pearson(Z.matrix, P.matrix,returnMatrix=FALSE)
  
  #### 保存函数内部环境的对象
  objects_to_save <- list(Z.matrix = Z.matrix, VarZ = VarZ, 
                          CorZ = CorZ , p = p, 
                          sumstats1_name = sumstats1_name, 
                          sumstats2_name = sumstats2_name)
  save(objects_to_save, file = paste0(out_path,"/",sumstats1_name, "-", sumstats2_name, ".RData"))
  
  ###run
  if(Placo.plus){
    cl <- makeCluster(ncore)
    clusterExport(cl, c("sumstats1_name","sumstats2_name","rscript_dir","out_path"))
    clusterEvalQ(cl, {
      library(data.table)
      library(magrittr)
      library(tidyverse)
      source(rscript_dir)
      load(paste0(out_path,"/",sumstats1_name, "-", sumstats2_name, ".RData"))
    })
    out <- parSapply(cl, 1:p, function(i) placo.plus(Z = Z.matrix[i,], VarZ = VarZ, CorZ=CorZ))
    stopCluster(cl)
    
    ####整合结果到原始文件
    f12<-fread(paste0(out_path,"/",sumstats1_name, "-",sumstats2_name,".txt.gz"))
    out %>% 
      t() %>% 
      as.data.table() %>% 
      mutate(p.placo.plus=as.numeric(p.placo.plus)) %>% 
      mutate(T.placo.plus=as.numeric(T.placo.plus)) %>% 
      cbind(f12)->f12
    f12 %>% filter(p.placo.plus<5*10^-8)->f12_sig###输出显著的点
  }else{
    # 设置多核并行计算
    cl <- makeCluster(ncore)
    clusterExport(cl, c("sumstats1_name","sumstats2_name",
                        "placo","placo.plus",".p.bessel",".pdfx",".pdfx.cor",".p.bessel.cor",
                        "out_path"))
    clusterEvalQ(cl, {
      library(data.table)
      library(magrittr)
      library(tidyverse)
      load(paste0(out_path,"/",sumstats1_name, "-", sumstats2_name, ".RData"))
    })
    out <- parSapply(cl, 1:p, function(i) placo(Z = Z.matrix[i,], VarZ = VarZ))
    stopCluster(cl)
    
    ####整合结果到原始文件
    f12<-fread(paste0(out_path,"/",sumstats1_name, "-",sumstats2_name,".txt.gz"))
    out %>% 
      t() %>% 
      as.data.table() %>% 
      mutate(p.placo=as.numeric(p.placo)) %>% 
      mutate(T.placo=as.numeric(T.placo)) %>% 
      cbind(f12)->f12
    f12 %>% filter(p.placo<5*10^-8)->f12_sig###输出显著的点
  }
  fwrite(f12,
         file = paste0(out_path,"/",sumstats1_name, "-",sumstats2_name,".txt.gz"),
         compress = "gzip",
         sep = "\t",
         quote = F,
         col.names = T,
         row.names = F)
  write.table(f12_sig,
              file = paste0(out_path,"/",sumstats1_name, "-",sumstats2_name,"_sig.txt"),
              sep = "\t",
              quote = F,
              col.names = T,
              row.names = F)
  rm(f12);rm(f12_sig);rm(out)
  print(paste0(sumstats1_name, "-",sumstats2_name))  
} 

