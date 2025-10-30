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
#----------------1. SNP2GENE:SMR--------------
GWAS_SMR<-function(smr_dir,smr_bfile,smr_gwas_file,smr_qtl_dir,smr_qtl_name,
                   smr_out,smr_gwas_name,smr_qtl_type,
                   diff_freq=0.2,diff_freq_prop=0.05,heidi_mtd=1, 
                   specify_probe=F, myprobe_list, 
                   smr_type="all", #all即qtl不分chr，seperate chr则分chr
                   seperate_chr #若smr_type="seperate chr"，则需要填这个
                   ){ 
  
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
    if(!dir.exists(smr_out3)){
      dir.create(smr_out3)
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
    }else{
      print(paste0("文件夹已经存在，可能已经做过分析，已跳过这个"))}
    
    #-----关闭临时文件-------#
    unlink(temp_path, recursive = TRUE)
  }else{
    if(smr_type=="seperate chr"){
      #------smr文件制作,读入临时路径----------#
      smr_gwas<-fread(smr_gwas_file)
      temp_path <- tempdir()# 创建临时路径
      for (chr in 1:22) {
        smr_gwas %>% 
          dplyr::rename(freq=EAF,b=BETA,se=SE,p=P,n=N) %>% 
          dplyr::filter(CHR==chr) %>% 
          dplyr::select(SNP,A1,A2,freq,b,se,p,n)->tmp_chr
        fwrite(tmp_chr,
               file = paste0(temp_path,"\\tmp",chr,".ma"),
               sep = "\t",col.names = T,row.names = F,quote = F)
      }
      smr_gwas_dir=paste0(temp_path,"\\tmp",seperate_chr,".ma")
      #-----设置输出路径,创建路径----------#
      smr_out1<-paste0(smr_out,"/",smr_gwas_name)
      smr_out2<-paste0(smr_out,"/",smr_gwas_name,"/",smr_qtl_type)
      smr_out3<-paste0(smr_out,"/",smr_gwas_name,"/",smr_qtl_type,"/",smr_qtl_name)
      if(!dir.exists(smr_out1)){dir.create(smr_out1)}
      if(!dir.exists(smr_out2)){dir.create(smr_out2)}
      if(!dir.exists(smr_out3)){
        dir.create(smr_out3)
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
      }else{
        print(paste0("文件夹已经存在，可能已经做过分析，已跳过这个"))}
      #-----关闭临时文件-------#
      unlink(temp_path, recursive = TRUE)
    }
  }
  }
  

#----------------2. SNP2GENE:SMR-multi--------------
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
    if(!dir.exists(smr_out3)){
      dir.create(smr_out3)
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
      }else{
      print(paste0("文件夹已经存在，可能已经做过分析，已跳过这个"))
    }
    #-----关闭临时文件-------#
    unlink(temp_path, recursive = TRUE)
  }
  
}

#----------------3. SNP2GENE:SMR-trans--------------
GWAS_SMR_trans<-function(smr_dir,smr_bfile,smr_gwas_file,smr_qtl_dir,smr_qtl_name,
                         smr_out,smr_gwas_name,smr_qtl_type,
                         diff_freq=0.2,diff_freq_prop=0.05,heidi_mtd=1, 
                         specify_probe=F, myprobe_list, 
                         smr_type="all", #all即qtl不分chr，seperate chr则分chr
                         seperate_chr #若smr_type="seperate chr"，则需要填这个
                         ){ 
  
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
    if(!dir.exists(smr_out3)){
      dir.create(smr_out3)
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
    }else{
      print(paste0("文件夹已经存在，可能已经做过分析，已跳过这个"))
      }
    
    #-----关闭临时文件-------#
    unlink(temp_path, recursive = TRUE)
  }else{
    if(smr_type=="seperate chr"){
      #------smr文件制作,读入临时路径----------#
      smr_gwas<-fread(smr_gwas_file)
      temp_path <- tempdir()# 创建临时路径
      for (chr in 1:22) {
        smr_gwas %>% 
          dplyr::rename(freq=EAF,b=BETA,se=SE,p=P,n=N) %>% 
          dplyr::filter(CHR==chr) %>% 
          dplyr::select(SNP,A1,A2,freq,b,se,p,n)->tmp_chr
        fwrite(tmp_chr,
               file = paste0(temp_path,"\\tmp",chr,".ma"),
               sep = "\t",col.names = T,row.names = F,quote = F)
      }
      smr_gwas_dir=paste0(temp_path,"\\tmp",seperate_chr,".ma")
      #-----设置输出路径,创建路径----------#
      smr_out1<-paste0(smr_out,"/",smr_gwas_name)
      smr_out2<-paste0(smr_out,"/",smr_gwas_name,"/",smr_qtl_type)
      smr_out3<-paste0(smr_out,"/",smr_gwas_name,"/",smr_qtl_type,"/",smr_qtl_name)
      if(!dir.exists(smr_out1)){dir.create(smr_out1)}
      if(!dir.exists(smr_out2)){dir.create(smr_out2)}
      if(!dir.exists(smr_out3)){
        dir.create(smr_out3)
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
      }else{
        print(paste0("文件夹已经存在，可能已经做过分析，已跳过这个"))}
      #-----关闭临时文件-------#
      unlink(temp_path, recursive = TRUE)
    }
  }
  
}

#----------------4. SNP2GENE: mbat-combo--------------
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