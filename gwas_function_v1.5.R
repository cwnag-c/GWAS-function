library(data.table)
library(magrittr)
library(tidyverse)
cat("*******************************************************************\n")
cat("* CWAS \n")
cat("* Version 1.0\n")
cat("* Update Date:2024.04.06\n")
cat("* (C) 2024 Chao Wang \n")
cat("* Nanjing Medical University\n")
cat("*******************************************************************\n")


#---------------1.质控前看看数据情况----------
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
#------------2.1.1 SNP匹配chr和pos----------------
GWAS_map_chrpos<-function(f1=NULL,MAP_dir,gwas_dir,pre_col_name,post_col_name,
                          character_chr=F){
  if(is.null(f1)==T){
    f1<-fread(gwas_dir,select =pre_col_name,col.names = post_col_name,fill=TRUE)
  }
  if(!exists("MAP_file")){
    MAP_file<-fread(MAP_dir)}
  if(character_chr){
    f1 %<>%mutate(CHR=as.integer(str_remove(string=CHR,pattern="chr"))) ####CHR去除“chr”
  }
  f1<-inner_join(f1,MAP_file,by="SNP")
  return(f1)
}
#------------2.1.2 chr和pos匹配SNP----------------
GWAS_map_snp<-function(f1=NULL,MAP_dir,gwas_dir,pre_col_name,post_col_name,
                          character_chr=F){
  if(is.null(f1)==T){
    f1<-fread(gwas_dir,select =pre_col_name,col.names = post_col_name,fill=TRUE)
  }
  
  if(!exists("MAP_file")){
    MAP_file<-fread(MAP_dir)}
  if(character_chr){
    f1 %<>%mutate(CHR=as.integer(str_remove(string=CHR,pattern="chr"))) ####CHR去除“chr”
  }
  f1<-inner_join(f1,MAP_file,by=c("CHR","POS"))
  return(f1)
}
#------------2.2.1 基因组位置转换,无需rsid,hg38to37-----------------
GWAS_liftover_hg38to37<-function(liftover_dir,gwas_dir,pre_col_name,post_col_name,
                        character_chr=F){
  gwas_f1<-fread(gwas_dir,select =pre_col_name,col.names = post_col_name,fill=TRUE) ##读入文件
  if(!exists("liftover_file")){
    liftover_file<-fread(liftover_dir)}
  if(character_chr){
    gwas_f1 %<>%mutate(CHR=as.integer(str_remove(string=CHR,pattern="chr"))) ####CHR去除“chr”
  }
  gwas_f1 %<>% filter(A1 %in% c("A","T","G","C")) %>% filter(A2 %in% c("A","T","G","C")) 
  gwas_f1 %<>% rename(pos38=POS,chr38=CHR) 
  gwas_f1<-inner_join(gwas_f1,liftover_file,by=c("chr38","pos38"),
                      relationship = "many-to-many")
  
  gwas_f1 %<>% select(-c("chr38","pos38")) %>% rename(POS=pos37,CHR=chr37)
  return(gwas_f1)
}
#------------2.2.2 基因组位置转换,无需rsid,hg37to38-----------------
GWAS_liftover_hg37to38<-function(liftover_dir,gwas_dir,pre_col_name,post_col_name,
                                 character_chr=F){
  gwas_f1<-fread(gwas_dir,select =pre_col_name,col.names = post_col_name,fill=TRUE) ##读入文件
  if(!exists("liftover_file")){
    liftover_file<-fread(liftover_dir)}
  if(character_chr){
    gwas_f1 %<>%mutate(CHR=as.integer(str_remove(string=CHR,pattern="chr"))) ####CHR去除“chr”
  }
  gwas_f1 %<>% filter(A1 %in% c("A","T","G","C")) %>% filter(A2 %in% c("A","T","G","C")) 
  gwas_f1 %<>% rename(pos37=POS,chr37=CHR) 
  gwas_f1<-inner_join(gwas_f1,liftover_file,by=c("chr37","pos37"),
                      relationship = "many-to-many")
  
  gwas_f1 %<>% select(-c("chr37","pos37")) %>% rename(POS=pos38,CHR=chr38)
  return(gwas_f1)
}

#-----------------3.1gwas清洗质控----------
GWAS_clean<-function(f1=NULL,gwas_dir,Samplesize,map_dir,ref_dir,
                     pre_col_name,post_col_name,
                     output_dir,folders_name,sumstat_name,
                     is_character_chr=FALSE, is_OR=F, is_MLOG10P=F ,is_info=F,
                     is_SE=F,is_Ncases_controls=F,is_AFcases_controls=F){
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.0.1\n")
  cat("* Update Date:2024.04.06\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  if(is.null(f1)==T){
    f1<-fread(gwas_dir,select =pre_col_name,col.names = post_col_name,fill=TRUE)
  }
  
  #-------------------清洗------------------#
  if(is_character_chr){
    f1 %<>%mutate(CHR=as.integer(str_remove(string=CHR,pattern="chr"))) ####CHR去除“chr”
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
                                 "保留常染色体的snp",
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
  
  ####CHR1-22
  if(is.character(f1$CHR)){
    f1 %<>% 
      filter(CHR!="X") %>%
      filter(CHR!="x") %>% 
      mutate(CHR=as.integer(CHR))
  }else{
    if(is.numeric(f1$CHR)){
      f1 %<>%
        filter(CHR>=1&CHR<=22) %>% 
        mutate(CHR=as.integer(CHR))
    }else{print("CHR is error")}
    
  }
  
  summary.dat[summary.dat$筛选=="保留常染色体的snp",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留常染色体的snp",
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
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留常染色体的snp",
                                                                                   sumstat_name]*100, digits=3)
  
  ####INFO>0.8
  if(is_info){f1 %<>%filter(INFO>0.8)}
  
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
  map_ref %<>%rename(SNP=rsid)
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
    mutate(POS=as.integer(POS)) %>%
    mutate(N=as.integer(N)) %>%  
    mutate(EAF=round(EAF,4))
  
  fwrite(f1,
         file = paste0(sumstat_name,"_v1.txt.gz"),
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
  }
}

#-------------3.2 gwas清洗质控不提供等位基因频率/用了参考基因组的等位基因频率-----
GWAS_clean_noEAF<-function(f1=NULL,gwas_dir,Samplesize,map_dir,ref_dir,
                           pre_col_name,post_col_name,
                           output_dir,folders_name,sumstat_name,
                           is_character_chr=FALSE, is_OR=F, is_MLOG10P=F ,is_info=F,
                           is_SE=F,is_Ncases_controls=F){
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.0.1\n")
  cat("* Update Date:2024.04.06\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  if(is.null(f1)==T){
    f1<-fread(gwas_dir,select =pre_col_name,col.names = post_col_name,fill=TRUE)
  }
  
  #-------------------清洗------------------#
  if(is_character_chr){
    f1 %<>%dplyr::mutate(CHR=as.integer(str_remove(string=CHR,pattern="chr"))) ####CHR去除“chr”
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
                                 "保留常染色体的snp",
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
  
  ####CHR1-22
  if(is.character(f1$CHR)){
    f1 %<>% 
      dplyr::filter(CHR!="X") %>%
      dplyr::filter(CHR!="x") %>% 
      dplyr::mutate(CHR=as.integer(CHR))
  }else{
    if(is.numeric(f1$CHR)){
      f1 %<>%
        dplyr::filter(CHR>=1&CHR<=22) %>% 
        dplyr::mutate(CHR=as.integer(CHR))
    }else{print("CHR is error")}
    
  }
  
  summary.dat[summary.dat$筛选=="保留常染色体的snp",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留常染色体的snp",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留 A,C,T,G 碱基",
                                                                                   sumstat_name]*100, digits=3)
  
  ####INFO>0.8
  if(is_info){f1 %<>% dplyr::filter(INFO>0.8)}
  
  summary.dat[summary.dat$筛选=="info>0.8",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="info>0.8",
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留常染色体的snp",
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
  map_ref %<>% dplyr::rename(SNP=rsid)
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
    dplyr::mutate(POS=as.integer(POS)) %>%
    dplyr::mutate(N=as.integer(N)) %>%  
    dplyr::mutate(EAF=round(EAF,4))
  
  fwrite(f1,
         file = paste0(sumstat_name,"_af1000_v1.txt.gz"),
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
  }
}

#----------------3.3 PGC daner格式清洗-----------
#CHR        SNP        BP     A1     A2 A1F40463 FRQ_U_313436  INFO      OR     SE      P   ngt....
GWAS_clean_daner<-function(gwas_dir,map_dir,ref_dir,
                     output_dir,folders_name,sumstat_name,is_seperate_Ncase_control=T){
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.0.1\n")
  cat("* Update Date:2024.04.10\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  
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
                                 "保留常染色体的snp",
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
  
  ####CHR1-22
  if(is.character(f1$CHR)){
    f1 %<>% 
      filter(CHR!="X") %>%
      filter(CHR!="x") %>% 
      mutate(CHR=as.integer(CHR))
  }else{
    if(is.numeric(f1$CHR)){
      f1 %<>%
        filter(CHR>=1&CHR<=22) %>% 
        mutate(CHR=as.integer(CHR))
    }else{print("CHR is error")}
    
  }
  
  summary.dat[summary.dat$筛选=="保留常染色体的snp",sumstat_name]=nrow(f1)
  summary.dat[summary.dat$筛选=="保留常染色体的snp",
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
              paste0(sumstat_name,"_pct_of_prev_step")]=round(nrow(f1)/summary.dat[summary.dat$筛选=="保留常染色体的snp",
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
  map_ref %<>%rename(SNP=rsid)
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
    mutate(POS=as.integer(POS)) %>%
    mutate(N=as.integer(N)) %>%  
    mutate(EAF=round(EAF,4))
  
  fwrite(f1,
         file = paste0(sumstat_name,"_v1.txt.gz"),
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
}


#-----------------3.4 phe-gwas清洗，如BrainMRI批量清洗----------
GWAS_cohort_look<-function(phenotype_list, filename_list,
                            type,input_dir,
                            output_dir,folders_name,
                            unzip_filename #解压后的文件名
){
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.3\n")
  cat("* Update Date:2024.04.20\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  
  #-----------------1.读入第一个文件，然后输出gwas summary-----#
  if(type=="zip"){
    f1_dir=paste0(input_dir,"/",filename_list[1])
    temp_path <- tempdir()# 创建临时路径
    unzip(zipfile = f1_dir, exdir = temp_path)
    f1_dir=paste0(temp_path,"\\",unzip_filename[1])
    f1<-fread(f1_dir)
    #-----关闭临时文件-------#
    unlink(temp_path, recursive = TRUE)
    
  }else{
    f1<-fread(f1_dir)
  }
  return(head(f1))
  
}

GWAS_cohort_clean<-function(phenotype_list, filename_list,
                            type,input_dir,
                            output_dir,folders_name,
                            unzip_filename, #解压后的文件名
                            pre_col_name,post_col_name,
                            is_N=T,#有无N列
                            seperate_N=T #N样本量是否不一致
                            ){
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.3\n")
  cat("* Update Date:2024.04.20\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  
  #-----------------1.读入第一个文件，然后输出gwas summary-----#
  if(type=="zip"){
    f1_dir=paste0(input_dir,"/",filename_list[1])
    temp_path <- tempdir()# 创建临时路径
    unzip(zipfile = f1_dir, exdir = temp_path)
    f1_dir=paste0(temp_path,"\\",unzip_filename[1])
    f1<-fread(f1_dir)
    #-----关闭临时文件-------#
    unlink(temp_path, recursive = TRUE)
  }else{
    f1<-fread(f1_dir)
  }
  
  colnames(f1)<-post_col_name 
  f1 %<>% 
    filter(EAF>0.01 & EAF<0.99) %>% 
    mutate(A1=toupper(A1)) %>% 
    mutate(A2=toupper(A2)) %>% 
    mutate(POS=as.integer(POS)) %>% 
    mutate(CHR=as.integer(CHR)) %>% 
    filter(A1 %in% c("A","T","G","C")) %>% 
    filter(A1 %in% c("A","T","G","C")) %>% 
    distinct(SNP,.keep_all = T) %>% 
    mutate(BETA=round(BETA,4)) %>% 
    mutate(SE=round(SE,4)) %>% 
    mutate(EAF=round(EAF,4)) %>% 
    mutate(P=format(P, scientific = TRUE, digits = 4)) %>% 
    mutate(P=as.numeric(P))
    
  
  if(seperate_N){
    f1 %>% select(SNP,CHR,POS,A1,A2,EAF)->SNP_annotation
    f1 %<>% select(SNP,N,BETA,SE,P)
  }else{
    f1 %>% select(SNP,CHR,POS,A1,A2,EAF,N)->SNP_annotation
    f1 %<>% select(SNP,BETA,SE,P)
  }

  #----输出----#
  fwrite(SNP_annotation,
         file = paste0(output_dir,"/SNP_annotation.txt.gz"),
         col.names = T,row.names = F,quote = F,sep = "\t",
         compress = "gzip")
  fwrite(f1,
         file = paste0(output_dir,"/",phenotype_list[1],".txt.gz"),
         col.names = T,row.names = F,quote = F,sep = "\t",
         compress = "gzip")
  
  SNP_annotation %<>% select(-c("EAF"))
  #-------循环--------------#
  for (i in 2:length(filename_list)) {
    if(type=="zip"){
      f1_dir=paste0(input_dir,"/",filename_list[i])
      temp_path <- tempdir()# 创建临时路径
      unzip(zipfile = f1_dir, exdir = temp_path)
      f1_dir=paste0(temp_path,"\\",unzip_filename[i])
      f1<-fread(f1_dir)
      #-----关闭临时文件-------#
      unlink(temp_path, recursive = TRUE)
    }else{
      f1<-fread(f1_dir)
    }
    colnames(f1)<-post_col_name 
    f1<-inner_join(f1,SNP_annotation,by=c("SNP","CHR","POS","A1","A2"))
    f1 %<>% 
      select(SNP,N,BETA,SE,P) %>% 
      mutate(BETA=round(BETA,4)) %>% 
      mutate(SE=round(SE,4)) %>% 
      mutate(P=format(P, scientific = TRUE, digits = 4)) %>% 
      mutate(P=as.numeric(P))
    fwrite(f1,
           file = paste0(output_dir,"/",phenotype_list[i],".txt.gz"),
           col.names = T,row.names = F,quote = F,sep = "\t",
           compress = "gzip")
  }
  }
#----------------4.1 MR:twosampleMR----------
GWAS_twosampleMR_iv<-function(
                          exposure_list,outcome_list, #pattern 1 必须输入的
                           gwas_prefix,out_path,
                           pattern,bi=T,
                           gwas_threshold=5*10^-8,
                           bfile_dir,plink_dir,MR_ancestry,
                           Clump_kb = 10000, Clump_r2 = 0.001){
  library(ieugwasr)
  library(TwoSampleMR)
  library(data.table)
  library(magrittr)
  library(tidyverse)
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.1\n")
  cat("* Update Date:2024.04.07\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  
  if(pattern==1){#分别输入暴露和结局，全跑
    out1<-lapply(1:nrow(exposure_list), function(e){
      exposure=paste0(gwas_prefix,exposure_list$absolute_dir[e])
      exposure<-fread(exposure)
      exposure<-as.data.frame(exposure)
      exposure$phenotype=exposure_list$abr[e]
      #----twosampleMR---#
      f1_exposure<-format_data(dat=exposure, type = "exposure",snps = NULL,header = TRUE,
                               phenotype_col="phenotype",
                               snp_col = "SNP",beta_col = "BETA",se_col = "SE",eaf_col = "EAF",
                               effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P",
                               chr_col = "CHR",pos_col = "POS",samplesize_col="N",
                               min_pval = 1e-500,log_pval = FALSE)
      f1_exposure<-f1_exposure[f1_exposure$pval.exposure< gwas_threshold,]
      #-----如果没有显著的点，就returnNULL，有则继续清洗-----#
      if(nrow(f1_exposure)!=0){ #如果没有显著的点，就跳过#
        f1_exposure %<>% mutate(rsid=SNP) %>%  mutate(pval=pval.exposure)
        f1_exposure_clump<-ld_clump(dat = f1_exposure,
                                    clump_kb = Clump_kb,
                                    clump_r2 = Clump_r2,
                                    pop = MR_ancestry,
                                    bfile =bfile_dir,
                                    plink_bin =plink_dir)
        rm(exposure);rm(f1_exposure);gc()
        
        #-------读入out_数据----#
        res<-lapply(1:nrow(outcome_list), function(c){
          outcome=paste0(gwas_prefix,outcome_list$absolute_dir[c])
          outcome<-fread(outcome)
          outcome<-as.data.frame(outcome)
          outcome$phenotype=outcome_list$abr[c]
          f2_outcome<-format_data(outcome,type = "outcome",snps = NULL,header = TRUE,
                                  phenotype_col="phenotype",
                                  snp_col = "SNP",beta_col = "BETA",se_col = "SE",eaf_col = "EAF",
                                  effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P",
                                  chr_col = "CHR",pos_col = "POS",samplesize_col="N",
                                  min_pval = 1e-500,log_pval = FALSE)
          dat<-harmonise_data(exposure_dat = f1_exposure_clump,outcome_dat =f2_outcome,action=1)
          #------------MR steiger---------#
          dat<-steiger_filtering(dat)
          #---------------F值计算--------------#
          dat %<>%
            mutate(R2=(2*(beta.exposure^2)*eaf.exposure*(1 - eaf.exposure))) %>% 
            mutate(Fvalue=((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))
          return(dat)
        })
        List_Name<- paste0(exposure_list$abr[e]," to ", outcome_list$abr)
        names(res)<-List_Name
        return(res)
      }else{
        return(NULL)
      }
    })
    if(bi){
      out2<-lapply(1:nrow(outcome_list), function(e){
        exposure=paste0(gwas_prefix,outcome_list$absolute_dir[e])
        exposure<-fread(exposure)
        exposure<-as.data.frame(exposure)
        exposure$phenotype=outcome_list$abr[e]
        #----twosampleMR---#
        f1_exposure<-format_data(dat=exposure, type = "exposure",snps = NULL,header = TRUE,
                                 phenotype_col="phenotype",
                                 snp_col = "SNP",beta_col = "BETA",se_col = "SE",eaf_col = "EAF",
                                 effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P",
                                 chr_col = "CHR",pos_col = "POS",samplesize_col="N",
                                 min_pval = 1e-500,log_pval = FALSE)
        f1_exposure<-f1_exposure[f1_exposure$pval.exposure< gwas_threshold,]
        #-----如果没有显著的点，就returnNULL，有则继续清洗-----#
        if(nrow(f1_exposure)!=0){ #如果没有显著的点，就跳过#
          f1_exposure %<>% mutate(rsid=SNP) %>%  mutate(pval=pval.exposure)
          f1_exposure_clump<-ld_clump(dat = f1_exposure,
                                      clump_kb = Clump_kb,
                                      clump_r2 = Clump_r2,
                                      pop = MR_ancestry,
                                      bfile =bfile_dir,
                                      plink_bin =plink_dir)
          rm(exposure);rm(f1_exposure);gc()
          
          #-------读入out_数据----#
          res<-lapply(1:nrow(exposure_list), function(c){
            outcome=paste0(gwas_prefix,exposure_list$absolute_dir[c])
            outcome<-fread(outcome)
            outcome<-as.data.frame(outcome)
            outcome$phenotype=exposure_list$abr[c]
            f2_outcome<-format_data(outcome,type = "outcome",snps = NULL,header = TRUE,
                                    phenotype_col="phenotype",
                                    snp_col = "SNP",beta_col = "BETA",se_col = "SE",eaf_col = "EAF",
                                    effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P",
                                    chr_col = "CHR",pos_col = "POS",samplesize_col="N",
                                    min_pval = 1e-500,log_pval = FALSE)
            dat<-harmonise_data(exposure_dat = f1_exposure_clump,outcome_dat =f2_outcome,action=1)
            #------------MR steiger---------#
            dat<-steiger_filtering(dat)
            #---------------F值计算--------------#
            dat %<>%
              mutate(R2=(2*(beta.exposure^2)*eaf.exposure*(1 - eaf.exposure))) %>% 
              mutate(Fvalue=((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))
            return(dat)
          })
          List_Name<- paste0(outcome_list$abr[e]," to ", exposure_list$abr)
          names(res)<-List_Name
          return(res)
        }else{
          return(NULL)
        }
      })
    }else{
      out2<-NULL
    }
    out<-c(out1,out2)
    
    integrated_list <- list() # 创建一个空列表，用于整合元素
    # 将嵌套的列表展平并添加到整合列表中
    for (sublist in out) {
      integrated_list <- c(integrated_list, sublist)
    }
    out<-integrated_list
    tmp_list<-list(out,exposure_list,outcome_list,Clump_kb,Clump_r2,
                             MR_ancestry,gwas_threshold,bfile_dir,plink_dir,
                             iv_rdata_name,out_path)
    tmp_list_name<-list("ivs","exposure_list","outcome_list","Clump_kb","Clump_r2",
                   "MR_ancestry","gwas_threshold","bfile_dir","plink_dir",
                   "iv_rdata_name","out_path")
    names(tmp_list)<-tmp_list_name
    return(tmp_list)
    
  }
  
}

GWAS_twosampleMR<-function(pattern,out_path,mr_name){
  library(TwoSampleMR)
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.1\n")
  cat("* Update Date:2024.04.07\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  rdata_dir<-paste0(out_path,mr_name,".Rdata")
  load(rdata_dir)
  IVs<-ivs$ivs
  if(pattern==1){
    exposure_list<-ivs$exposure_list
    outcome_list<-ivs$outcome_list
    exp_c<-exposure_list$abr
    ouc_c<-outcome_list$abr
    
    tmp1<-lapply(1:length(exp_c), function(e){
      exp_tmp<-exp_c[e]
      tmp_c<-paste0(exp_tmp," to ",ouc_c)
    })
    tmp2<-lapply(1:length(ouc_c), function(e){
      ouc_tmp<-ouc_c[e]
      tmp_c<-paste0(ouc_tmp," to ",exp_c)
    })
    tmp2<-unlist(tmp2);tmp1<-unlist(tmp1)
    list_name<-c(tmp1,tmp2)
    
    #-------twosampleMR-------------#
    MR_res<-lapply(list_name, function(i){
      dat=IVs[[i]]
      if(length(dat)==0){
        restBindsub<-NULL
        leaveone_out<-NULL
        mr_presso<-NULL
      }else{
        if(nrow(dat)==0){
          restBindsub<-NULL
          leaveone_out<-NULL
          mr_presso<-NULL
        }else{
          dat %<>% dplyr::filter(ambiguous=="FALSE")
          if(nrow(dat)>2){
            #----------twosamoleMR------#
            res<-mr(dat)
            het<-mr_heterogeneity(dat)
            ple<-mr_pleiotropy_test(dat)
            restBindsub<-bind_rows(res,het,ple)
            leaveone_out<-mr_leaveoneout(dat)
          }else{if(nrow(dat)==2){
            res<-mr(dat)
            het<-mr_heterogeneity(dat)
            ple<-mr_pleiotropy_test(dat)
            restBindsub<-bind_rows(res,het,ple)
            leaveone_out<-mr_leaveoneout(dat)
          }else{if(nrow(dat)==1){
            restBindsub<-mr(dat)
            leaveone_out<-NULL
          }else{if(nrow(dat)==0){
            restBindsub<-NULL
            leaveone_out<-NULL
          }
          }
          }
            
            
          }
        }

      }
      return(list(mr_res=restBindsub,leave_one_out=leaveone_out))
    })
    mr_res_list <- lapply(MR_res, function(x) x$mr_res)
    leaveone_out_list <- lapply(MR_res, function(x) x$leave_one_out)
    mr_res_list<-bind_rows(mr_res_list)
    leaveone_out_list<-bind_rows(leaveone_out_list)
    write.csv(mr_res_list,file = paste0(out_path,mr_name,".mr_res.csv"))
    write.csv(leaveone_out_list,file = paste0(out_path,mr_name,".leaveone_out_list.csv"))
    
  }

}

#----------------4.2 MR:GSMR2----------------
#----------------4.3 MR:CAUSE---------------
#----------------4.4 MR:MVMR---------------
#----------------5.1 SNP2GENE:SMR--------------
GWAS_SMR<-function(smr_dir,smr_bfile,smr_gwas_file,smr_qtl_dir,smr_qtl_name,
                   smr_out,smr_gwas_name,smr_qtl_type,
                   diff_freq=0.2,diff_freq_prop=0.05,heidi_mtd=1, 
                   specify_probe=F, myprobe_list, 
                   smr_type="all"){ #all即部分chr，seperate_chr则分chr
  
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.1\n")
  cat("* Update Date:2024.04.07\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  
  if(smr_type=="all"){
    #------smr文件制作,读入临时路径----------#
    smr_gwas<-fread(smr_gwas_file)
    smr_gwas %<>%
      dplyr::rename(freq=EAF,b=BETA,se=SE,p=P,n=N) %>% 
      dplyr::select(SNP,A1,A2,freq,b,se,p,n)
    
    temp_path <- tempdir()# 创建临时路径
    fwrite(smr_gwas,
           file = paste0(temp_path,"\\tmp.ma"),
           sep = "\t",col.names = T,row.names = F,quote = F)
    smr_gwas_dir=paste0(temp_path,"\\tmp.ma")
    #-----设置输出路径,创建路径----------#
    smr_out1<-paste0(smr_out,"/",smr_gwas_name)
    smr_out2<-paste0(smr_out,"/",smr_gwas_name,"/",smr_qtl_type)
    smr_out3<-paste0(smr_out,"/",smr_gwas_name,"/",smr_qtl_type,"/",smr_qtl_name)
    if(!dir.exists(smr_out1)){dir.create(smr_out1)}
    if(!dir.exists(smr_out2)){dir.create(smr_out2)}
    if(!dir.exists(smr_out3)){dir.create(smr_out3)}
    smr_out_name=paste0(smr_out3,"/",smr_gwas_name,".",smr_qtl_type,".",smr_qtl_name)
    
    #------生成命令------------#
    if(specify_probe){
      c1<-paste(smr_dir,
                "--bfile",smr_bfile,
                "--gwas-summary",smr_gwas_dir,
                "--beqtl-summary",smr_qtl_dir,
                "--heidi-mtd",heidi_mtd,
                "--diff-freq",diff_freq,
                "--diff-freq-prop",diff_freq_prop,
                "--extract-probe",myprobe_list, ##指定探针跑smr
                "--out",smr_out_name)
    }else{
      c1<-paste(smr_dir,
                "--bfile",smr_bfile,
                "--gwas-summary",smr_gwas_dir,
                "--beqtl-summary",smr_qtl_dir,
                "--heidi-mtd",heidi_mtd,
                "--diff-freq",diff_freq,
                "--diff-freq-prop",diff_freq_prop,
                #"--extract-probe",myprobe_list, ##指定探针跑smr
                "--out",smr_out_name)
    }
    #---------------输出结果---------------#
    system(c1)
    #-----关闭临时文件-------#
    unlink(temp_path, recursive = TRUE)
  }
  
}

#----------------5.2 SNP2GENE:SMR-multi--------------
GWAS_SMR_multi<-function(smr_dir,smr_bfile,smr_gwas_file,smr_qtl_dir,smr_qtl_name,
                   smr_out,smr_gwas_name,smr_qtl_type,
                   diff_freq=0.2,diff_freq_prop=0.05,heidi_mtd=1, 
                   specify_probe=F, myprobe_list, 
                   smr_type="all"){ #all即部分chr，seperate_chr则分chr
  
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.1\n")
  cat("* Update Date:2024.04.07\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  
  if(smr_type=="all"){
    #------smr文件制作,读入临时路径----------#
    smr_gwas<-fread(smr_gwas_file)
    smr_gwas %<>%
      dplyr::rename(freq=EAF,b=BETA,se=SE,p=P,n=N) %>% 
      dplyr::select(SNP,A1,A2,freq,b,se,p,n)
    
    temp_path <- tempdir()# 创建临时路径
    fwrite(smr_gwas,
           file = paste0(temp_path,"\\tmp.ma"),
           sep = "\t",col.names = T,row.names = F,quote = F)
    smr_gwas_dir=paste0(temp_path,"\\tmp.ma")
    #-----设置输出路径,创建路径----------#
    smr_out1<-paste0(smr_out,"/",smr_gwas_name)
    smr_out2<-paste0(smr_out,"/",smr_gwas_name,"/",smr_qtl_type)
    smr_out3<-paste0(smr_out,"/",smr_gwas_name,"/",smr_qtl_type,"/",smr_qtl_name)
    if(!dir.exists(smr_out1)){dir.create(smr_out1)}
    if(!dir.exists(smr_out2)){dir.create(smr_out2)}
    if(!dir.exists(smr_out3)){dir.create(smr_out3)}
    smr_out_name=paste0(smr_out3,"/",smr_gwas_name,".",smr_qtl_type,".",smr_qtl_name)
    
    #------生成命令------------#
    if(specify_probe){
      c1<-paste(smr_dir,
                "--bfile",smr_bfile,
                "--gwas-summary",smr_gwas_dir,
                "--beqtl-summary",smr_qtl_dir,
                "--heidi-mtd",heidi_mtd,
                "--diff-freq",diff_freq,
                "--diff-freq-prop",diff_freq_prop,
                "--extract-probe",myprobe_list, ##指定探针跑smr
                "--out",smr_out_name,
                "--smr-multi")
    }else{
      c1<-paste(smr_dir,
                "--bfile",smr_bfile,
                "--gwas-summary",smr_gwas_dir,
                "--beqtl-summary",smr_qtl_dir,
                "--heidi-mtd",heidi_mtd,
                "--diff-freq",diff_freq,
                "--diff-freq-prop",diff_freq_prop,
                #"--extract-probe",myprobe_list, ##指定探针跑smr
                "--out",smr_out_name,
                "--smr-multi")
    }
    #---------------输出结果---------------#
    system(c1)
    #-----关闭临时文件-------#
    unlink(temp_path, recursive = TRUE)
  }
  
}

#----------------5.3 SNP2GENE:SMR-trans--------------
GWAS_SMR_trans<-function(smr_dir,smr_bfile,smr_gwas_file,smr_qtl_dir,smr_qtl_name,
                         smr_out,smr_gwas_name,smr_qtl_type,
                         diff_freq=0.2,diff_freq_prop=0.05,heidi_mtd=1, 
                         specify_probe=F, myprobe_list, 
                         smr_type="all"){ #all即部分chr，seperate_chr则分chr
  
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.1\n")
  cat("* Update Date:2024.04.07\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  
  if(smr_type=="all"){
    #------smr文件制作,读入临时路径----------#
    smr_gwas<-fread(smr_gwas_file)
    smr_gwas %<>%
      dplyr::rename(freq=EAF,b=BETA,se=SE,p=P,n=N) %>% 
      dplyr::select(SNP,A1,A2,freq,b,se,p,n)
    
    temp_path <- tempdir()# 创建临时路径
    fwrite(smr_gwas,
           file = paste0(temp_path,"\\tmp.ma"),
           sep = "\t",col.names = T,row.names = F,quote = F)
    smr_gwas_dir=paste0(temp_path,"\\tmp.ma")
    #-----设置输出路径,创建路径----------#
    smr_out1<-paste0(smr_out,"/",smr_gwas_name)
    smr_out2<-paste0(smr_out,"/",smr_gwas_name,"/",smr_qtl_type)
    smr_out3<-paste0(smr_out,"/",smr_gwas_name,"/",smr_qtl_type,"/",smr_qtl_name)
    if(!dir.exists(smr_out1)){dir.create(smr_out1)}
    if(!dir.exists(smr_out2)){dir.create(smr_out2)}
    if(!dir.exists(smr_out3)){dir.create(smr_out3)}
    smr_out_name=paste0(smr_out3,"/",smr_gwas_name,".",smr_qtl_type,".",smr_qtl_name)
    
    #------生成命令------------#
    if(specify_probe){
      c1<-paste(smr_dir,
                "--bfile",smr_bfile,
                "--gwas-summary",smr_gwas_dir,
                "--beqtl-summary",smr_qtl_dir,
                "--heidi-mtd",heidi_mtd,
                "--diff-freq",diff_freq,
                "--diff-freq-prop",diff_freq_prop,
                "--extract-probe",myprobe_list, ##指定探针跑smr
                "--out",smr_out_name,
                "--trans")
    }else{
      c1<-paste(smr_dir,
                "--bfile",smr_bfile,
                "--gwas-summary",smr_gwas_dir,
                "--beqtl-summary",smr_qtl_dir,
                "--heidi-mtd",heidi_mtd,
                "--diff-freq",diff_freq,
                "--diff-freq-prop",diff_freq_prop,
                #"--extract-probe",myprobe_list, ##指定探针跑smr
                "--out",smr_out_name,
                "--trans")
    }
    #---------------输出结果---------------#
    system(c1)
    #-----关闭临时文件-------#
    unlink(temp_path, recursive = TRUE)
  }
  
}

#----------------5.4 SNP2GENE: mbat-combo--------------
GWAS_mbat_combo<-function(gcta_dir,bfile_dir,gene_list_dir,
                          gwas_prefix,gwas_dir,
                          out_path,ancestry,file_name){
  
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.2\n")
  cat("* Update Date:2024.04.15\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  #------gcta文件制作,读入临时路径----------#
  gcta_gwas<-fread(paste0(gwas_prefix,gwas_dir))
  gcta_gwas %<>%
    dplyr::select(SNP,A1,A2,EAF,BETA,SE,P,N) %>% 
    dplyr::rename(freq=EAF)
  
  temp_path <- tempdir()# 创建临时路径
  fwrite(gcta_gwas,
         file = paste0(temp_path,"\\tmp.ma"),
         sep = "\t",col.names = T,row.names = F,quote = F)
  gcta_gwas_dir=paste0(temp_path,"\\tmp.ma")
  #-----设置输出路径,创建路径----------#
  
  if(!dir.exists(paste0(out_path,"/",ancestry))){
    dir.create(paste0(out_path,"/",ancestry))}
    
  if(!dir.exists(paste0(out_path,"/",ancestry,"/",file_name))){
    dir.create(paste0(out_path,"/",ancestry,"/",file_name))
    setwd(paste0(out_path,"/",ancestry,"/",file_name))
    out_dir=paste0(out_path,"/",ancestry,"/",file_name,"/",file_name)
    #-----gcta命令-------#
    c1=gcta_dir
    c2=paste0("--bfile ",bfile_dir)
    c3=paste0("--mBAT-combo ",gcta_gwas_dir)
    c4=paste0("--mBAT-gene-list ",gene_list_dir)
    c5="--mBAT-print-all-p "
    c6=paste0("--out ",out_dir)
    command=paste(c1,c2,c3,c4,c5,c6,sep = " ")
    system(command)
  }else{
    print(paste0("文件夹已经存在，可能已经做过分析，已跳过这个"))}
  #-----关闭临时文件-------#
  unlink(temp_path, recursive = TRUE)
} 
#----------------6.1 Cross-trait:PLACO-------------
GWAS_placo<-function(sumstats1,sumstats2,
                     GWAS_txt_dir,
                     gwas_prefix,
                     out_path,
                     PLACO_path,
                     ncore){
  
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
  f12 %<>% 
    dplyr::filter(!(CHR==6&POS>25000000&POS<35000000)) %>% 
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
  save(objects_to_save, file = paste0(out_path, "/", sumstats1, "-", sumstats2, ".RData"))
  
  print(paste(m, n, sep = "-"))
  
  # 设置多核并行计算
  cl <- makeCluster(ncore)
  clusterExport(cl, c("sumstats1","sumstats2","PLACO_path"))
  clusterEvalQ(cl, {
    library(data.table)
    library(magrittr)
    library(tidyverse)
    load(paste0(sumstats1, "-", sumstats2, ".RData"))
    load(PLACO_path)
  })
  out <- parSapply(cl, 1:p, function(i) placo(Z = Z.matrix.decor[i,], VarZ = VarZ))
  stopCluster(cl)
  
  ####整合结果到原始文件
  f12<-fread(paste0(sumstats1, "-",sumstats2,".txt.gz"))
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

