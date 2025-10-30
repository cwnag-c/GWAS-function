library(data.table)
library(magrittr)
library(tidyverse)
cat("*******************************************************************\n")
cat("* CWAS \n")
cat("* Version 2.0\n")
cat("* Update Date:2024.08.13\n")
cat("* (C) 2024 Chao Wang \n")
cat("* Nanjing Medical University\n")
cat("*******************************************************************\n")

#-----------------1. 7z解压------------
GWAS_7z<-function(file_path,destination_dir,py_dir){
  library(reticulate)
  # 运行Python脚本
  source_python(py_dir)
  
  # 调用Python函数解压7z文件
  py$extract_7z(file_path, destination_dir)
  
}

#---------------2.质控前看看数据情况----------
GWAS_look<-function(gwas_dir,n=5){
  Tmp<-fread(gwas_dir,nrows = n)
  temp1<-NULL
  frq_cols1 <- grep("^FRQ_U_", names(Tmp), value = TRUE)
  frq_cols2 <- grep("^FRQ_A_", names(Tmp), value = TRUE)
  if(length(frq_cols1)!=0){
    if(length(frq_cols2)!=0){
      temp1<-paste0("提示：似乎是PGC的daner格式，建议使用GWAS_clean_daner")
    }
  }
  temp_list<-list(Tmp,temp1)
  return(temp_list)
}
#------------3.1 SNP匹配chr和pos----------------
XWAS_map_chrpos<-function(f1=NULL,MAP_dir,gwas_dir,pre_col_name,post_col_name){
  if(is.null(f1)==T){
    f1<-fread(gwas_dir,select =pre_col_name,col.names = post_col_name,fill=TRUE)
  }
  if(!exists("MAP_file")){
    MAP_file<-fread(MAP_dir)
    MAP_file %<>% rename(SNP=rsmid,CHR=chr,POS=pos)
    }

  f1<-inner_join(f1,MAP_file,by="SNP")
  return(f1)
}
#------------3.2 chr和pos匹配SNP----------------
XWAS_map_snp<-function(f1=NULL,MAP_dir,gwas_dir,pre_col_name,post_col_name,
                          character_chr=F){
  if(is.null(f1)==T){
    f1<-fread(gwas_dir,select =pre_col_name,col.names = post_col_name,fill=TRUE)
  }
  f1 %<>% mutate(CHR=as.integer(CHR)) %>% mutate(POS=as.integer(POS))
  if(!exists("MAP_file")){
    MAP_file<-fread(MAP_dir)
    MAP_file %<>% rename(SNP=rsmid,CHR=chr,POS=pos)
    }
  if(character_chr){
    f1 %<>%mutate(CHR=str_remove(string=CHR,pattern="chr")) ####CHR去除“chr”
  }
  f1<-inner_join(f1,MAP_file,by=c("CHR","POS"))
  return(f1)
}
#-----------------5.1 gwas清洗质控----------
XWAS_clean<-function(f1=NULL,gwas_dir,Samplesize,map_dir,ref_dir,
                     pre_col_name,post_col_name,
                     output_dir,folders_name,sumstat_name,
                     is_character_chr=FALSE, is_OR=F, is_MLOG10P=F ,is_info=F,
                     is_SE=F,is_Ncases_controls=F,is_AFcases_controls=F){
  cat("****************************Strat***************************************\n")
  if(is.null(f1)==T){
    f1<-fread(gwas_dir,select =pre_col_name,col.names = post_col_name,fill=TRUE)
  }
  
  #-------------------清洗------------------#
  if(is_character_chr){
    f1 %<>%mutate(CHR=str_remove(string=CHR,pattern="chr")) ####CHR去除“chr”
  }
  
  if(is_OR){
    f1 %<>%mutate(BETA=log(BETA))###如果是OR，转化为BETA
  } 
  
  ##有些P是chr形式，需要统一成num
  f1 %<>% mutate(P=as.numeric(P))  
  if(is_MLOG10P){
    f1 %<>% 
      rename(MLOG10P=P) %>% 
      mutate(P=10^(-MLOG10P))
  }
  
  #----缺少se，转化计算se-----#
  if(is_SE){
    f1 %<>% mutate(SE= abs(BETA/qnorm(P/2)))
  }
  #------若有病例对照单独的频率，则合成一个EAF-------------#
  if(is_AFcases_controls){
    f1 %<>% mutate(EAF=(AF1+AF2)/2)
  }
  
  #-----若存在单独的Ncanses和Ncontrols列，则合并成N--------#
  if(is_Ncases_controls){
    f1 %<>% mutate(N=Ncases+Ncontrols)
    n_cases=max(f1$Ncases)
    n_controls=max(f1$Ncontrols)
  }
  #-----添加样本量--------#
  if ("N" %in% names(f1)) {
    print(max(f1$N))
  }else{
    f1 %<>% mutate(N=Samplesize) }

  ###创建文件夹和工作目录
  dir1<-paste0(output_dir,folders_name)
  if(!dir.exists(dir1)){
    dir.create(dir1)
  }
  
  if(!dir.exists(paste0(dir1,"/",sumstat_name))){
    dir.create(paste0(dir1,"/",sumstat_name))
  }
  setwd(paste0(dir1,"/",sumstat_name))
  
  ###创建list
  summary.dat<-data.frame("筛选"=c("样本量",
                                  "初始snp数",
                                 "保留MAF>0.01",
                                 "保留无缺失值",
                                 "保留数值有意义",
                                 "保留 A,C,T,G 碱基",
                                 "保留X染色体的snp",
                                 "N>30&MAC>6",
                                 "info>0.8",
                                 "去除重复值",
                                 "map后的snp数",
                                 "chr、pos与ref一致",
                                 "与ref取交集的snp数",
                                 "不是回文碱基的snp数量",
                                 "回文碱基，保留MAF<0.4且与ref差值绝对值小于0.1的snp数量",
                                 "清洗完的snp数量",
                                 "剔除ref与EAF差值大于0.2"))
  
  ###样本量
  summary.dat[summary.dat$筛选=="样本量",sumstat_name]=max(f1$N)
  ###初始snp数量
  summary.dat[summary.dat$筛选=="初始snp数",sumstat_name]=nrow(f1)
  
  ####MAF>0.01
  f1%<>% filter(EAF>0.01&EAF<0.99)
  
  summary.dat[summary.dat$筛选=="保留MAF>0.01",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留MAF>0.01",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="初始snp数",
                                                                                   sumstat_name]*100, digits=3)
  
  ###去除缺失值
  f1 %<>% 
    filter(is.na(CHR)==F) %>% 
    filter(is.na(SNP)==F) %>% 
    filter(is.na(POS)==F) %>% 
    filter(is.na(A1)==F) %>% 
    filter(is.na(A2)==F) %>% 
    filter(is.na(EAF)==F) %>%
    filter(is.na(BETA)==F) %>%  
    filter(is.na(SE)==F) %>% 
    filter(is.na(P)==F) %>% 
    filter(is.na(N)==F)
  summary.dat[summary.dat$筛选=="保留无缺失值",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留无缺失值",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留MAF>0.01",
                                                                                   sumstat_name]*100, digits=3)
  
  ####看P值的形式
  f1 %<>% 
    filter(P>=0&P<=1) %>% 
    filter(SE>0) %>% 
    filter(BETA!=0) %>% 
    filter(N>0)
  
  summary.dat[summary.dat$筛选=="保留数值有意义",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留数值有意义",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留无缺失值",
                                                                                   sumstat_name]*100, digits=3)
  
  
  ####Filter SNPs with ID
  f1 %<>% 
    mutate(A1=toupper(A1)) %>% 
    mutate(A2=toupper(A2)) %>% 
    filter(A1=="A"|A1=="T"|A1=="G"|A1=="C") %>% 
    filter(A2=="A"|A2=="T"|A2=="G"|A2=="C")
  
  summary.dat[summary.dat$筛选=="保留 A,C,T,G 碱基",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留 A,C,T,G 碱基",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留数值有意义",
                                                                                   sumstat_name]*100, digits=3)
  
  ####性染色体如果是23或者小写x，统一转化成X
  if(is.character(f1$CHR)){
    f1 %>% 
      filter(CHR=="X")->f1.1
    f1 %>% 
      filter(CHR=="x") %>%
      mutate(CHR="X")->f1.2
    f1<-rbind(f1.1,f1.2);rm(f1.1);rm(f1.2);gc()
  }else{
    if(is.numeric(f1$CHR)){
      f1 %<>%mutate(CHR="X")
    }else{print("CHR is error")}
  }
  
  summary.dat[summary.dat$筛选=="保留X染色体的snp",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留X染色体的snp",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留 A,C,T,G 碱基",
                                                                                   sumstat_name]*100, digits=3)
  ####N>30&MAC>6
  f1 %>% 
    filter(EAF<=0.5) %>% 
    mutate(MAF=EAF)->f1.1
  
  f1 %>% 
    filter(EAF>0.5) %>% 
    mutate(MAF=1-EAF)->f1.2
  
  f1<-rbind(f1.1,f1.2);rm(f1.1);rm(f1.2);gc()
  
  f1 %<>%
    mutate(MAC=2*MAF*N) %>% 
    filter(N>=30) %>% 
    filter(MAC>6)
  
  summary.dat[summary.dat$筛选=="N>30&MAC>6",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="N>30&MAC>6",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留X染色体的snp",
                                                                                   sumstat_name]*100, digits=3)  
  
  
  ####INFO>0.8
  if(is_info){f1 %<>%filter(INFO>0.8)}
  
  summary.dat[summary.dat$筛选=="info>0.8",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="info>0.8",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="N>30&MAC>6",
                                                                                   sumstat_name]*100, digits=3)
  
  ####去重
  f1 %<>% 
    dplyr::arrange(desc(MAF)) %>% 
    distinct(SNP,.keep_all = T)
  summary.dat[summary.dat$筛选=="去除重复值",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="去除重复值",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="info>0.8",
                                                                                   sumstat_name]*100, digits=3)
  
  ####P值检验
  f1 %<>%mutate(manual.p=2*pnorm(-abs(BETA/SE)))
  
  png(paste0(sumstat_name,"_provided_p_values_vs_manually_calc_p_val.png"), width=2880, height=2880, res=360, type="cairo")
  smoothScatter(f1$manual.p,
                f1$P,
                ylim=c(0,1), 
                xlim=c(0,1), 
                ylab=paste0(sumstat_name," Provided P values"),
                xlab="Beta-calculated P values",
                main=paste0(sumstat_name," provided p values vs manually calculated p values"),nrpoints = Inf)
  dev.off()
  
  ####snp 的chr和pos检验
  map_ref<-fread(map_dir)
  map_ref %<>% rename(SNP=rsmid) 
  map_ref %<>%distinct(SNP,.keep_all = T)
  f1 %<>% inner_join(map_ref,by="SNP",multiple = "all")
  summary.dat[summary.dat$筛选=="map后的snp数",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="map后的snp数",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="去除重复值",
                                                                                   sumstat_name]*100, digits=3)
  f1 %<>%
    filter(CHR==chr) %>% 
    filter(POS==pos)
  summary.dat[summary.dat$筛选=="chr、pos与ref一致",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="chr、pos与ref一致",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="map后的snp数",
                                                                                   sumstat_name]*100, digits=3)
  rm(map_ref);gc()
  
  ####EAF check and A1与ref一致
  ref<-fread(ref_dir)
  f1 %<>% 
    mutate(cptid=paste(CHR,POS,sep = ":")) %>% 
    inner_join(ref,by="cptid",multiple = "all")
  rm(ref);gc()
  summary.dat[summary.dat$筛选=="与ref取交集的snp数",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="与ref取交集的snp数",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="chr、pos与ref一致",
                                                                                   sumstat_name]*100, digits=3)
  ###质控前的ref eaf与EAF散点图
  png(paste0("ref_eaf_vs_",sumstat_name,"_eaf_before_qc.png"), width=2880, height=2880,res=360, type="cairo")
  smoothScatter(f1$eaf,
                f1$EAF,
                ylim=c(0,1), xlim=c(0,1), 
                ylab=paste0(sumstat_name," EAF"),
                xlab="Reference EAF",
                main=paste0("Reference EAF vs ",sumstat_name," SNPs EAF before QC"),nrpoints = Inf)
  dev.off()
  
  #----------回文碱基和不是的分别处理--------#
  f1 %>% filter(((A1=="A" & A2=="T")|(A1=="T" & A2=="A")|
                    (A1=="G" & A2=="C")|(A1=="C" & A2=="G"))) %>% 
    filter(EAF<0.4|EAF>0.6) %>% 
    filter(eaf<0.4|eaf>0.6)->f_amb
  
  f1 %<>% filter(!((A1=="A" & A2=="T")|(A1=="T" & A2=="A")|
                    (A1=="G" & A2=="C")|(A1=="C" & A2=="G")))
  
  
  ####处理不是回文碱基##
  f1 %>% filter(A1==ea&A2==oa)->f1.1 ###A1与ref一致
  f1 %>% filter(A1==oa&A2==ea)->f1.2 ###A2与ref一致
  f1 %>% filter(!((A1==oa&A2==ea)|(A1==ea&A2==oa)))->a ###不一致的
  f1.1 %<>% 
    mutate(newbeta=BETA) %>% 
    mutate(neweaf=EAF)
  f1.2%<>% 
    mutate(newbeta=BETA/-1) %>% 
    mutate(neweaf=1-EAF)
  f1 <-rbind(f1.1,f1.2)
  rm(f1.1);rm(f1.2);gc()
  
  #---正负链翻转---#
  if(nrow(a)!=0){
    a %<>% 
      mutate(a1 = case_when(A1 == "A" ~ "1", A1 == "T" ~ "-1", 
                            A1 == "C" ~ "-2",A1 == "G" ~ "2",
                            TRUE ~ as.character(A1) )) %>% 
      mutate(a2 = case_when(A2 == "A" ~ "1", A2 == "T" ~ "-1", 
                            A2 == "C" ~ "-2",A2 == "G" ~ "2",
                            TRUE ~ as.character(A2) )) %>% 
      mutate(EA = case_when(ea == "A" ~ '1', ea == "T" ~ "-1", 
                            ea == "C" ~ "-2",ea == "G" ~ "2",
                            TRUE ~ as.character(ea) )) %>% 
      mutate(OA = case_when(oa == "A" ~ "1", oa == "T" ~ "-1", 
                            oa == "C" ~ "-2",oa == "G" ~ "2",
                            TRUE ~ as.character(oa) ))
    a$a1<-as.numeric(a$a1); a$a2<-as.numeric(a$a2)
    a$EA<-as.numeric(a$EA); a$OA<-as.numeric(a$OA)
    a %>% filter((a1+OA==0 & a2+EA==0))->a1  ###需要变EAF和beta
    a %>% filter((a1+EA==0 & a2+OA==0))->a2  ###虽然在不同链上，但不需要变EAF和beta
    if(nrow(a1)!=0){
      a1%<>% 
        mutate(newbeta=BETA/-1) %>% 
        mutate(neweaf=1-EAF)
    }
    if(nrow(a2)!=0){
      a2 %<>% 
        mutate(newbeta=BETA) %>% 
        mutate(neweaf=EAF)
    }
    a <-bind_rows(a1,a2);rm(a1);rm(a2)
    a %<>% select(-c("a1","a2","EA","OA")) 
    f1<-bind_rows(f1,a);rm(a)
  }
  
  #------处理回文碱基----#
  f_amb %<>% filter(((ea=="A" & oa=="T")|(ea=="T" & oa=="A")|
                      (ea=="G" & oa=="C")|(ea=="C" & oa=="G")))
  f_amb %>% 
    filter((EAF<0.5 & eaf<0.5)|(EAF>0.5 & eaf>0.5)) %>% 
    mutate(newbeta=BETA) %>% 
    mutate(neweaf=EAF)->f_amb1
  f_amb %>% 
    filter(!((EAF<0.5 & eaf<0.5)|(EAF>0.5 & eaf>0.5))) %>% 
    mutate(newbeta=BETA/-1) %>% 
    mutate(neweaf=1-EAF)->f_amb2
  
  f_amb<-bind_rows(f_amb1,f_amb2);rm(f_amb1);rm(f_amb2)
  #----剔除与ref相差超过0.1的回文碱基----#
  f_amb %>% filter(abs(neweaf-eaf)<=0.1)

  summary.dat[summary.dat$筛选=="不是回文碱基的snp数量",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="不是回文碱基的snp数量",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="与ref取交集的snp数",
                                                                                     sumstat_name]*100, digits=3)
  
  summary.dat[summary.dat$筛选=="回文碱基，保留MAF<0.4且与ref差值绝对值小于0.1的snp数量",sumstat_name]=nrow(f_amb)
  summary.dat[summary.dat$筛选=="回文碱基，保留MAF<0.4且与ref差值绝对值小于0.1的snp数量",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f_amb)/summary.dat[summary.dat$筛选=="与ref取交集的snp数",
                                                                                     sumstat_name]*100, digits=3)
  f1<-bind_rows(f1,f_amb);rm(f_amb)
  summary.dat[summary.dat$筛选=="清洗完的snp数量",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="清洗完的snp数量",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="与ref取交集的snp数",
                                                                                  sumstat_name]*100, digits=3)
  ###剔除REF与EAF差值在0.2以上的snp
  f1 %<>% filter(abs(neweaf-eaf)<=0.2)
  summary.dat[summary.dat$筛选=="剔除ref与EAF差值大于0.2",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="剔除ref与EAF差值大于0.2",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="清洗完的snp数量",
                                                                                   sumstat_name]*100, digits=3)
  ####质控完后的散点图
  png(paste0("ref_eaf_vs_",sumstat_name,"_eaf_after_qc.png"), width=2880, height=2880, res=360, type="cairo")
  smoothScatter(f1$eaf, 
                f1$neweaf,
                ylim=c(0,1), xlim=c(0,1), 
                ylab=paste0(sumstat_name," EAF"),
                xlab="Reference EAF",
                main=paste0("Reference EAF vs ",sumstat_name," SNPs EAF after QC"),nrpoints = Inf)
  dev.off()
  
  ###输出###
  f1 %<>% 
    select(SNP,CHR,POS,ea,oa,neweaf,newbeta,SE,P,N) %>% 
    rename(A1=ea,
           A2=oa,
           BETA=newbeta,
           EAF=neweaf) %>% 
    arrange(CHR,POS)
  
  
  f1 %<>% 
    mutate(CHR=23) %>% 
    mutate(POS=as.integer(POS)) %>%
    mutate(N=as.integer(N)) %>%  
    mutate(EAF=round(EAF,4))
  
  fwrite(f1,
         file = paste0(sumstat_name,"_chrX.txt.gz"),
         compress = "gzip",
         sep = "\t",
         quote = F,
         col.names = T,
         row.names = F)
  
  write.table(summary.dat,
              file = paste0(sumstat_name,"_summary.txt"),
              sep = "\t",
              quote = F,
              col.names = T,
              row.names = F)
  gc()
  #------如果有病例对照列，则返回---#
  if(is_Ncases_controls){
    n_list=c(n_cases,n_controls)
    return(n_list)
    cat("****************************End***************************************\n")
  }
}

#-------------5.2 gwas清洗质控不提供等位基因频率/用了参考基因组的等位基因频率-----
XWAS_clean_noEAF<-function(f1=NULL,gwas_dir,Samplesize,map_dir,ref_dir,
                           pre_col_name,post_col_name,
                           output_dir,folders_name,sumstat_name,
                           is_character_chr=FALSE, is_OR=F, is_MLOG10P=F ,is_info=F,
                           is_SE=F,is_Ncases_controls=F){
  cat("********************************Start*******************************\n")
  if(is.null(f1)==T){
    f1<-fread(gwas_dir,select =pre_col_name,col.names = post_col_name,fill=TRUE)
  }
  
  #-------------------清洗------------------#
  if(is_character_chr){
    f1 %<>%dplyr::mutate(CHR=str_remove(string=CHR,pattern="chr")) ####CHR去除“chr”
  }
  
  if(is_OR){
    f1 %<>%dplyr::mutate(BETA=log(BETA))###如果是OR，转化为BETA
  } 
  
  ##有些P是chr形式，需要统一成num
  f1 %<>% dplyr::mutate(P=as.numeric(P))  
  if(is_MLOG10P){
    f1 %<>% 
      dplyr::rename(MLOG10P=P) %>% 
      dplyr::mutate(P=10^(-MLOG10P))
  }   
  #----缺少se，转化计算se-----#
  if(is_SE){
    f1 %<>% dplyr::mutate(SE= abs(BETA/qnorm(P/2)))
  }
  #-----若存在单独的Ncanses和Ncontrols列，则合并成N--------#
  if(is_Ncases_controls){
    f1 %<>% dplyr::mutate(N=Ncases+Ncontrols)
    n_cases=max(f1$Ncases)
    n_controls=max(f1$Ncontrols)
  }
  #-----添加样本量--------#
  if ("N" %in% names(f1)) {
    print(max(f1$N))
  }else{
    f1 %<>% dplyr::mutate(N=Samplesize) }
  

  ###创建文件夹和工作目录
  dir1<-paste0(output_dir,folders_name)
  if(!dir.exists(dir1)){
    dir.create(dir1)
  }
  
  if(!dir.exists(paste0(dir1,"/",sumstat_name))){
    dir.create(paste0(dir1,"/",sumstat_name))
  }
  setwd(paste0(dir1,"/",sumstat_name))
  
  ###创建list
  summary.dat<-data.frame("筛选"=c("样本量",
                                 "初始snp数",
                                 "保留无缺失值",
                                 "保留数值有意义",
                                 "保留 A,C,T,G 碱基",
                                 "保留X染色体的snp",
                                 "info>0.8",
                                 "去除重复值",
                                 "map后的snp数",
                                 "chr、pos与ref一致",
                                 "与ref取交集的snp数",
                                 "剔除回文碱基后的snp",
                                 "N>30&MAC>6",
                                 "MAF>0.01"))
  
  ###样本量
  summary.dat[summary.dat$筛选=="样本量",sumstat_name]=max(f1$N)
  ###初始snp数量
  summary.dat[summary.dat$筛选=="初始snp数",sumstat_name]=nrow(f1)
  
  
  ###去除缺失值
  f1 %<>% 
    dplyr::filter(is.na(CHR)==F) %>% 
    dplyr::filter(is.na(SNP)==F) %>% 
    dplyr::filter(is.na(POS)==F) %>% 
    dplyr::filter(is.na(A1)==F) %>% 
    dplyr::filter(is.na(A2)==F) %>% 
    dplyr::filter(is.na(BETA)==F) %>%  
    dplyr::filter(is.na(SE)==F) %>% 
    dplyr::filter(is.na(P)==F) %>% 
    dplyr::filter(is.na(N)==F)
  summary.dat[summary.dat$筛选=="保留无缺失值",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留无缺失值",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="初始snp数",
                                                                                   sumstat_name]*100, digits=3)
  
  ####看P值的形式
  f1 %<>% 
    dplyr::filter(P>=0&P<=1) %>% 
    dplyr::filter(SE>0) %>% 
    dplyr::filter(BETA!=0) %>% 
    dplyr::filter(N>0)
  
  summary.dat[summary.dat$筛选=="保留数值有意义",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留数值有意义",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留无缺失值",
                                                                                   sumstat_name]*100, digits=3)
  
  
  ####Filter SNPs with ID
  f1 %<>% 
    dplyr::mutate(A1=toupper(A1)) %>% 
    dplyr::mutate(A2=toupper(A2)) %>% 
    dplyr::filter(A1=="A"|A1=="T"|A1=="G"|A1=="C") %>% 
    dplyr::filter(A2=="A"|A2=="T"|A2=="G"|A2=="C")
  
  summary.dat[summary.dat$筛选=="保留 A,C,T,G 碱基",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留 A,C,T,G 碱基",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留数值有意义",
                                                                                   sumstat_name]*100, digits=3)
  ####性染色体如果是23或者小写x，统一转化成X
  if(is.character(f1$CHR)){
    f1 %>% 
      filter(CHR=="X")->f1.1
    f1 %>% 
      filter(CHR=="x") %>%
      mutate(CHR="X")->f1.2
    f1<-rbind(f1.1,f1.2);rm(f1.1);rm(f1.2);gc()
  }else{
    if(is.numeric(f1$CHR)){
      f1 %<>%mutate(CHR="X")
    }else{print("CHR is error")}
  }

  summary.dat[summary.dat$筛选== "保留X染色体的snp",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选== "保留X染色体的snp",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留 A,C,T,G 碱基",
                                                                                   sumstat_name]*100, digits=3)
  ####INFO>0.8
  if(is_info){f1 %<>% dplyr::filter(INFO>0.8)}
  
  summary.dat[summary.dat$筛选=="info>0.8",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="info>0.8",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留X染色体的snp",
                                                                                   sumstat_name]*100, digits=3)
  
  ####去重
  f1 %<>% 
    #arrange(desc(MAF)) %>% 
    dplyr::distinct(SNP,.keep_all = T)
  summary.dat[summary.dat$筛选=="去除重复值",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="去除重复值",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="info>0.8",
                                                                                   sumstat_name]*100, digits=3)
  
  ####P值检验
  f1 %<>% dplyr::mutate(manual.p=2*pnorm(-abs(BETA/SE)))
  
  png(paste0(sumstat_name,"_provided_p_values_vs_manually_calc_p_val.png"), width=2880, height=2880, res=360, type="cairo")
  smoothScatter(f1$manual.p,
                f1$P,
                ylim=c(0,1), 
                xlim=c(0,1), 
                ylab=paste0(sumstat_name," Provided P values"),
                xlab="Beta-calculated P values",
                main=paste0(sumstat_name," provided p values vs manually calculated p values"),nrpoints = Inf)
  dev.off()
  
  ####snp 的chr和pos检验
  map_ref<-fread(map_dir)
  map_ref %<>% rename(SNP=rsmid)
  map_ref %<>% dplyr::distinct(SNP,.keep_all = T)
  f1 %<>% dplyr::inner_join(map_ref,by="SNP",multiple = "all")
  summary.dat[summary.dat$筛选=="map后的snp数",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="map后的snp数",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="去除重复值",
                                                                                   sumstat_name]*100, digits=3)
  f1 %<>%
    dplyr::filter(CHR==chr) %>% 
    dplyr::filter(POS==pos)
  summary.dat[summary.dat$筛选=="chr、pos与ref一致",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="chr、pos与ref一致",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="map后的snp数",
                                                                                   sumstat_name]*100, digits=3)
  rm(map_ref);gc()
  
  ####用参考面板的EAF填补
  ref<-fread(ref_dir)
  f1 %<>% 
    dplyr::mutate(cptid=paste(CHR,POS,sep = ":")) %>% 
    dplyr::inner_join(ref,by="cptid",multiple = "all")
  rm(ref);gc()
  summary.dat[summary.dat$筛选=="与ref取交集的snp数",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="与ref取交集的snp数",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="chr、pos与ref一致",
                                                                            sumstat_name]*100, digits=3)
  
  
  #------剔除回文碱基------#
  f1 %<>% dplyr::filter(!((A1=="A" & A2=="T")|(A1=="T" & A2=="A")|
                     (A1=="G" & A2=="C")|(A1=="C" & A2=="G")))
  
  summary.dat[summary.dat$筛选=="剔除回文碱基后的snp",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="剔除回文碱基后的snp",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="与ref取交集的snp数",
                                                                                   sumstat_name]*100, digits=3)
  
  
  
  ####N>30&MAC>6
  f1 %>% 
    dplyr::filter(eaf<=0.5) %>% 
    dplyr::mutate(MAF=eaf)->f1.1
  
  f1 %>% 
    dplyr::filter(eaf>0.5) %>% 
    dplyr::mutate(MAF=1-eaf)->f1.2
  f1<-rbind(f1.1,f1.2);rm(f1.1);rm(f1.2);gc()
  f1 %<>%
    dplyr::mutate(MAC=2*MAF*N) %>% 
    dplyr::filter(N>=30) %>% 
    dplyr::filter(MAC>6)
  
  summary.dat[summary.dat$筛选=="N>30&MAC>6",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="N>30&MAC>6",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="剔除回文碱基后的snp",
                                                                                   sumstat_name]*100, digits=3)

  #------MAF要大于0.01---------#
  f1 %<>% dplyr::filter(MAF>0.01)
  summary.dat[summary.dat$筛选=="MAF>0.01",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="MAF>0.01",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="N>30&MAC>6",
                                                                                   sumstat_name]*100, digits=3)
  
  ####处理不是回文碱基##
  f1 %>% dplyr::filter(A1==ea&A2==oa)->f1.1 ###A1与ref一致
  f1 %>% dplyr::filter(A1==oa&A2==ea)->f1.2 ###A2与ref一致
  f1 %>% dplyr::filter(!((A1==oa&A2==ea)|(A1==ea&A2==oa)))->a ###不一致的
  f1.1 %<>% 
    dplyr::mutate(newbeta=BETA)
  f1.2%<>% 
    dplyr::mutate(newbeta=BETA/-1)
  f1 <-rbind(f1.1,f1.2)
  rm(f1.1);rm(f1.2);gc()
  
  #---正负链翻转---#
  if(nrow(a)!=0){
    a %<>% 
      dplyr::mutate(a1 = case_when(A1 == "A" ~ "1", A1 == "T" ~ "-1", 
                            A1 == "C" ~ "-2",A1 == "G" ~ "2",
                            TRUE ~ as.character(A1) )) %>% 
      dplyr::mutate(a2 = case_when(A2 == "A" ~ "1", A2 == "T" ~ "-1", 
                            A2 == "C" ~ "-2",A2 == "G" ~ "2",
                            TRUE ~ as.character(A2) )) %>% 
      dplyr::mutate(EA = case_when(ea == "A" ~ '1', ea == "T" ~ "-1", 
                            ea == "C" ~ "-2",ea == "G" ~ "2",
                            TRUE ~ as.character(ea) )) %>% 
      dplyr::mutate(OA = case_when(oa == "A" ~ "1", oa == "T" ~ "-1", 
                            oa == "C" ~ "-2",oa == "G" ~ "2",
                            TRUE ~ as.character(oa) ))
    a$a1<-as.numeric(a$a1); a$a2<-as.numeric(a$a2)
    a$EA<-as.numeric(a$EA); a$OA<-as.numeric(a$OA)
    a %>% dplyr::filter((a1+OA==0 & a2+EA==0))->a1  ###需要变EAF和beta
    a %>% dplyr::filter((a1+EA==0 & a2+OA==0))->a2  ###虽然在不同链上，但不需要变EAF和beta
    if(nrow(a1)!=0){
      a1%<>% 
        dplyr::mutate(newbeta=BETA/-1)
    }
    if(nrow(a2)!=0){
      a2 %<>% 
        dplyr::mutate(newbeta=BETA)
    }
    a <-bind_rows(a1,a2);rm(a1);rm(a2)
    a %<>% dplyr::select(-c("a1","a2","EA","OA")) 
    f1<-bind_rows(f1,a);rm(a)
  }
  
  
  

  ###输出###
  f1 %<>% 
    dplyr::select(SNP,CHR,POS,ea,oa,eaf,newbeta,SE,P,N) %>% 
    dplyr::rename(A1=ea,
           A2=oa,
           BETA=newbeta,
           EAF=eaf) %>% 
    dplyr::arrange(CHR,POS)
  
  
  f1 %<>% 
    mutate(CHR=23) %>% 
    mutate(POS=as.integer(POS)) %>%
    mutate(N=as.integer(N)) %>%  
    mutate(EAF=round(EAF,4))
  
  fwrite(f1,
         file = paste0(sumstat_name,"_af1000_chrX.txt.gz"),
         compress = "gzip",
         sep = "\t",
         quote = F,
         col.names = T,
         row.names = F)
  
  write.table(summary.dat,
              file = paste0(sumstat_name,"_summary.txt"),
              sep = "\t",
              quote = F,
              col.names = T,
              row.names = F)
  gc()
  #------如果有病例对照列，则返回---#
  if(is_Ncases_controls){
    n_list=c(n_cases,n_controls)
    return(n_list)
    cat("********************************End*******************************\n")
    
  }
}

#----------------5.3 PGC daner格式清洗-----------
#CHR        SNP        BP     A1     A2 A1F40463 FRQ_U_313436  INFO      OR     SE      P   ngt....
XWAS_clean_daner<-function(gwas_dir,map_dir,ref_dir,
                     output_dir,folders_name,sumstat_name,is_seperate_Ncase_control=T){
  cat("*******************************Start********************************\n")
  f1<-fread(gwas_dir,fill=TRUE)
  
  #------从列名中提出病例对照数----#
  frq_cols1 <- grep("^FRQ_A_", names(f1), value = TRUE)
  frq_cols2 <- grep("^FRQ_U_", names(f1), value = TRUE)
  Ncases <- as.numeric(sub("^FRQ_A_", "", frq_cols1))
  Ncontrols <- as.numeric(sub("^FRQ_U_", "", frq_cols2))
  #-----------统一病例和对照的频率的列名，并生成一个总EAF，根据加权样本量-----#
  names(f1)[names(f1) %in% frq_cols1] <- ("A1F")
  names(f1)[names(f1) %in% frq_cols2] <- ("A2F")
  if(is_seperate_Ncase_control){
    f1 %<>%
      mutate(EAF= (Ncases*A1F+Ncontrols*A2F)/(Ncases+Ncontrols)) %>% 
      mutate(N=Nca+Nco) %>% 
      mutate(BETA=log(OR)) %>% 
      rename(POS=BP) %>% 
      select(SNP,CHR,POS,A1,A2,EAF,BETA,SE,P,INFO,N)
  }else{
    f1 %<>%
      mutate(EAF= (Ncases*A1F+Ncontrols*A2F)/(Ncases+Ncontrols)) %>% 
      mutate(N=Ncases+Ncontrols) %>% 
      mutate(BETA=log(OR)) %>% 
      rename(POS=BP) %>% 
      select(SNP,CHR,POS,A1,A2,EAF,BETA,SE,P,INFO,N)
  }
  
  gc()
  #-----清洗-------------#
  
  ###创建文件夹和工作目录
  dir1<-paste0(output_dir,folders_name)
  if(!dir.exists(dir1)){
    dir.create(dir1)
  }
  
  if(!dir.exists(paste0(dir1,"/",sumstat_name))){
    dir.create(paste0(dir1,"/",sumstat_name))
  }
  setwd(paste0(dir1,"/",sumstat_name))
  
  ###创建list
  summary.dat<-data.frame("筛选"=c("样本量",
                                 "初始snp数",
                                 "保留MAF>0.01",
                                 "保留无缺失值",
                                 "保留数值有意义",
                                 "保留 A,C,T,G 碱基",
                                 "保留X染色体的snp",
                                 "N>30&MAC>6",
                                 "info>0.8",
                                 "去除重复值",
                                 "map后的snp数",
                                 "chr、pos与ref一致",
                                 "与ref取交集的snp数",
                                 "不是回文碱基的snp数量",
                                 "回文碱基，保留MAF<0.4且与ref差值绝对值小于0.1的snp数量",
                                 "清洗完的snp数量",
                                 "剔除ref与EAF差值大于0.2"))
  
  ###样本量
  summary.dat[summary.dat$筛选=="样本量",sumstat_name]=max(f1$N)
  ###初始snp数量
  summary.dat[summary.dat$筛选=="初始snp数",sumstat_name]=nrow(f1)
  
  ####MAF>0.01
  f1%<>% filter(EAF>0.01&EAF<0.99)
  
  summary.dat[summary.dat$筛选=="保留MAF>0.01",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留MAF>0.01",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="初始snp数",
                                                                                   sumstat_name]*100, digits=3)
  
  ###去除缺失值
  f1 %<>% 
    filter(is.na(CHR)==F) %>% 
    filter(is.na(SNP)==F) %>% 
    filter(is.na(POS)==F) %>% 
    filter(is.na(A1)==F) %>% 
    filter(is.na(A2)==F) %>% 
    filter(is.na(EAF)==F) %>%
    filter(is.na(BETA)==F) %>%  
    filter(is.na(SE)==F) %>% 
    filter(is.na(P)==F) %>% 
    filter(is.na(N)==F)
  summary.dat[summary.dat$筛选=="保留无缺失值",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留无缺失值",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留MAF>0.01",
                                                                                   sumstat_name]*100, digits=3)
  
  ####看P值的形式
  f1 %<>% 
    filter(P>=0&P<=1) %>% 
    filter(SE>0) %>% 
    filter(BETA!=0) %>% 
    filter(N>0)
  
  summary.dat[summary.dat$筛选=="保留数值有意义",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留数值有意义",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留无缺失值",
                                                                                   sumstat_name]*100, digits=3)
  
  
  ####Filter SNPs with ID
  f1 %<>% 
    mutate(A1=toupper(A1)) %>% 
    mutate(A2=toupper(A2)) %>% 
    filter(A1=="A"|A1=="T"|A1=="G"|A1=="C") %>% 
    filter(A2=="A"|A2=="T"|A2=="G"|A2=="C")
  
  summary.dat[summary.dat$筛选=="保留 A,C,T,G 碱基",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留 A,C,T,G 碱基",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留数值有意义",
                                                                                   sumstat_name]*100, digits=3)
  ####性染色体如果是23或者小写x，统一转化成X
  if(is.character(f1$CHR)){
    f1 %>% 
      filter(CHR=="X")->f1.1
    f1 %>% 
      filter(CHR=="x") %>%
      mutate(CHR="X")->f1.2
    f1<-rbind(f1.1,f1.2);rm(f1.1);rm(f1.2);gc()
  }else{
    if(is.numeric(f1$CHR)){
      f1 %<>%mutate(CHR="X")
    }else{print("CHR is error")}
  }
  

  summary.dat[summary.dat$筛选=="保留X染色体的snp",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留X染色体的snp",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留 A,C,T,G 碱基",
                                                                                   sumstat_name]*100, digits=3)
  ####N>30&MAC>6
  f1 %>% 
    filter(EAF<=0.5) %>% 
    mutate(MAF=EAF)->f1.1
  
  f1 %>% 
    filter(EAF>0.5) %>% 
    mutate(MAF=1-EAF)->f1.2
  
  f1<-rbind(f1.1,f1.2);rm(f1.1);rm(f1.2);gc()
  
  f1 %<>%
    mutate(MAC=2*MAF*N) %>% 
    filter(N>=30) %>% 
    filter(MAC>6)
  
  summary.dat[summary.dat$筛选=="N>30&MAC>6",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="N>30&MAC>6",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留X染色体的snp",
                                                                                   sumstat_name]*100, digits=3)  
  
  ####INFO>0.8
  f1 %<>%filter(INFO>0.8)
  
  summary.dat[summary.dat$筛选=="info>0.8",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="info>0.8",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="N>30&MAC>6",
                                                                                   sumstat_name]*100, digits=3)
  
  ####去重
  f1 %<>% 
    arrange(desc(MAF)) %>% 
    distinct(SNP,.keep_all = T)
  summary.dat[summary.dat$筛选=="去除重复值",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="去除重复值",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="info>0.8",
                                                                                   sumstat_name]*100, digits=3)
  
  ####P值检验
  f1 %<>%mutate(manual.p=2*pnorm(-abs(BETA/SE)))
  
  png(paste0(sumstat_name,"_provided_p_values_vs_manually_calc_p_val.png"), width=2880, height=2880, res=360, type="cairo")
  smoothScatter(f1$manual.p,
                f1$P,
                ylim=c(0,1), 
                xlim=c(0,1), 
                ylab=paste0(sumstat_name," Provided P values"),
                xlab="Beta-calculated P values",
                main=paste0(sumstat_name," provided p values vs manually calculated p values"),nrpoints = Inf)
  dev.off()
  
  ####snp 的chr和pos检验
  map_ref<-fread(map_dir)
  map_ref %<>% rename(SNP=rsmid) 
  map_ref %<>% distinct(SNP,.keep_all = T)
  f1 %<>% inner_join(map_ref,by="SNP",multiple = "all")
  summary.dat[summary.dat$筛选=="map后的snp数",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="map后的snp数",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="去除重复值",
                                                                                   sumstat_name]*100, digits=3)
  f1 %<>%
    filter(CHR==chr) %>% 
    filter(POS==pos)
  summary.dat[summary.dat$筛选=="chr、pos与ref一致",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="chr、pos与ref一致",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="map后的snp数",
                                                                                   sumstat_name]*100, digits=3)
  rm(map_ref);gc()
  
  ####EAF check and A1与ref一致
  ref<-fread(ref_dir)
  f1 %<>% 
    mutate(cptid=paste(CHR,POS,sep = ":")) %>% 
    inner_join(ref,by="cptid",multiple = "all")
  rm(ref);gc()
  summary.dat[summary.dat$筛选=="与ref取交集的snp数",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="与ref取交集的snp数",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="chr、pos与ref一致",
                                                                                   sumstat_name]*100, digits=3)
  ###质控前的ref eaf与EAF散点图
  png(paste0("ref_eaf_vs_",sumstat_name,"_eaf_before_qc.png"), width=2880, height=2880,res=360, type="cairo")
  smoothScatter(f1$eaf,
                f1$EAF,
                ylim=c(0,1), xlim=c(0,1), 
                ylab=paste0(sumstat_name," EAF"),
                xlab="Reference EAF",
                main=paste0("Reference EAF vs ",sumstat_name," SNPs EAF before QC"),nrpoints = Inf)
  dev.off()
  
  #----------回文碱基和不是的分别处理--------#
  f1 %>% filter(((A1=="A" & A2=="T")|(A1=="T" & A2=="A")|
                   (A1=="G" & A2=="C")|(A1=="C" & A2=="G"))) %>% 
    filter(EAF<0.4|EAF>0.6) %>% 
    filter(eaf<0.4|eaf>0.6)->f_amb
  
  f1 %<>% filter(!((A1=="A" & A2=="T")|(A1=="T" & A2=="A")|
                     (A1=="G" & A2=="C")|(A1=="C" & A2=="G")))
  
  
  ####处理不是回文碱基##
  f1 %>% filter(A1==ea&A2==oa)->f1.1 ###A1与ref一致
  f1 %>% filter(A1==oa&A2==ea)->f1.2 ###A2与ref一致
  f1 %>% filter(!((A1==oa&A2==ea)|(A1==ea&A2==oa)))->a ###不一致的
  f1.1 %<>% 
    mutate(newbeta=BETA) %>% 
    mutate(neweaf=EAF)
  f1.2%<>% 
    mutate(newbeta=BETA/-1) %>% 
    mutate(neweaf=1-EAF)
  f1 <-rbind(f1.1,f1.2)
  rm(f1.1);rm(f1.2);gc()
  
  #---正负链翻转---#
  if(nrow(a)!=0){
    a %<>% 
      mutate(a1 = case_when(A1 == "A" ~ "1", A1 == "T" ~ "-1", 
                            A1 == "C" ~ "-2",A1 == "G" ~ "2",
                            TRUE ~ as.character(A1) )) %>% 
      mutate(a2 = case_when(A2 == "A" ~ "1", A2 == "T" ~ "-1", 
                            A2 == "C" ~ "-2",A2 == "G" ~ "2",
                            TRUE ~ as.character(A2) )) %>% 
      mutate(EA = case_when(ea == "A" ~ '1', ea == "T" ~ "-1", 
                            ea == "C" ~ "-2",ea == "G" ~ "2",
                            TRUE ~ as.character(ea) )) %>% 
      mutate(OA = case_when(oa == "A" ~ "1", oa == "T" ~ "-1", 
                            oa == "C" ~ "-2",oa == "G" ~ "2",
                            TRUE ~ as.character(oa) ))
    a$a1<-as.numeric(a$a1); a$a2<-as.numeric(a$a2)
    a$EA<-as.numeric(a$EA); a$OA<-as.numeric(a$OA)
    a %>% filter((a1+OA==0 & a2+EA==0))->a1  ###需要变EAF和beta
    a %>% filter((a1+EA==0 & a2+OA==0))->a2  ###虽然在不同链上，但不需要变EAF和beta
    if(nrow(a1)!=0){
      a1%<>% 
        mutate(newbeta=BETA/-1) %>% 
        mutate(neweaf=1-EAF)
    }
    if(nrow(a2)!=0){
      a2 %<>% 
        mutate(newbeta=BETA) %>% 
        mutate(neweaf=EAF)
    }
    a <-bind_rows(a1,a2);rm(a1);rm(a2)
    a %<>% select(-c("a1","a2","EA","OA")) 
    f1<-bind_rows(f1,a);rm(a)
  }
  
  #------处理回文碱基----#
  f_amb %<>% filter(((ea=="A" & oa=="T")|(ea=="T" & oa=="A")|
                       (ea=="G" & oa=="C")|(ea=="C" & oa=="G")))
  f_amb %>% 
    filter((EAF<0.5 & eaf<0.5)|(EAF>0.5 & eaf>0.5)) %>% 
    mutate(newbeta=BETA) %>% 
    mutate(neweaf=EAF)->f_amb1
  f_amb %>% 
    filter(!((EAF<0.5 & eaf<0.5)|(EAF>0.5 & eaf>0.5))) %>% 
    mutate(newbeta=BETA/-1) %>% 
    mutate(neweaf=1-EAF)->f_amb2
  
  f_amb<-bind_rows(f_amb1,f_amb2);rm(f_amb1);rm(f_amb2)
  #----剔除与ref相差超过0.1的回文碱基----#
  f_amb %>% filter(abs(neweaf-eaf)<=0.1)
  
  summary.dat[summary.dat$筛选=="不是回文碱基的snp数量",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="不是回文碱基的snp数量",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="与ref取交集的snp数",
                                                                                   sumstat_name]*100, digits=3)
  
  summary.dat[summary.dat$筛选=="回文碱基，保留MAF<0.4且与ref差值绝对值小于0.1的snp数量",sumstat_name]=nrow(f_amb)
  summary.dat[summary.dat$筛选=="回文碱基，保留MAF<0.4且与ref差值绝对值小于0.1的snp数量",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f_amb)/summary.dat[summary.dat$筛选=="与ref取交集的snp数",
                                                                                      sumstat_name]*100, digits=3)
  f1<-bind_rows(f1,f_amb);rm(f_amb)
  summary.dat[summary.dat$筛选=="清洗完的snp数量",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="清洗完的snp数量",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="与ref取交集的snp数",
                                                                                   sumstat_name]*100, digits=3)
  ###剔除REF与EAF差值在0.2以上的snp
  f1 %<>% filter(abs(neweaf-eaf)<=0.2)
  summary.dat[summary.dat$筛选=="剔除ref与EAF差值大于0.2",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="剔除ref与EAF差值大于0.2",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="清洗完的snp数量",
                                                                                   sumstat_name]*100, digits=3)
  ####质控完后的散点图
  png(paste0("ref_eaf_vs_",sumstat_name,"_eaf_after_qc.png"), width=2880, height=2880, res=360, type="cairo")
  smoothScatter(f1$eaf, 
                f1$neweaf,
                ylim=c(0,1), xlim=c(0,1), 
                ylab=paste0(sumstat_name," EAF"),
                xlab="Reference EAF",
                main=paste0("Reference EAF vs ",sumstat_name," SNPs EAF after QC"),nrpoints = Inf)
  dev.off()
  
  ###输出###
  f1 %<>% 
    select(SNP,CHR,POS,ea,oa,neweaf,newbeta,SE,P,N) %>% 
    rename(A1=ea,
           A2=oa,
           BETA=newbeta,
           EAF=neweaf) %>% 
    arrange(CHR,POS)
  
  
  f1 %<>% 
    mutate(CHR=23) %>% 
    mutate(POS=as.integer(POS)) %>%
    mutate(N=as.integer(N)) %>%  
    mutate(EAF=round(EAF,4))
  
  fwrite(f1,
         file = paste0(sumstat_name,"_chrX.txt.gz"),
         compress = "gzip",
         sep = "\t",
         quote = F,
         col.names = T,
         row.names = F)
  
  write.table(summary.dat,
              file = paste0(sumstat_name,"_summary.txt"),
              sep = "\t",
              quote = F,
              col.names = T,
              row.names = F)
  gc()
  N_list<-list(Ncases,Ncontrols)
  return(N_list)
  cat("*******************************End********************************\n")
  
}


