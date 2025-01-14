library(data.table)
library(magrittr)
library(tidyverse)
library(ieugwasr)
library(TwoSampleMR)
cat("*******************************************************************\n")
cat("* CWAS \n")
cat("* Version 1.8\n")
cat("* Update Date:2024.07.17\n")
cat("* (C) 2024 Chao Wang \n")
cat("* Nanjing Medical University\n")
cat("*******************************************************************\n")
#----------------1.1 twosampleMR:iv----------
GWAS_twosampleMR_iv<-function(
    exposure_list=NULL,outcome_list=NULL, #输入gwas_list
    exposure=NULL,outcome=NULL, #输入gwas summary数据
    gwas_prefix,
    bi_mr=F,
    gwas_threshold=5*10^-8,
    bfile_dir,plink_dir,MR_ancestry,
    Clump_kb = 10000, Clump_r2 = 0.001){
  if(any(is.null(exposure),is.null(outcome))){
    out1<-lapply(1:nrow(exposure_list), function(e){
      exposure=paste0(gwas_prefix,exposure_list$absolute_dir[e])
      exposure<-fread(exposure)
      exposure<-as.data.frame(exposure)
      exposure$phenotype=exposure_list$abr[e]
      #----twosampleMR---#
      exposure<-format_data(dat=exposure, type = "exposure",snps = NULL,header = TRUE,
                            phenotype_col="phenotype",
                            snp_col = "SNP",beta_col = "BETA",se_col = "SE",eaf_col = "EAF",
                            effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P",
                            chr_col = "CHR",pos_col = "POS",samplesize_col="N",
                            min_pval = 1e-500,log_pval = FALSE)
      exposure<-exposure[exposure$pval.exposure< gwas_threshold,]
      gc()
      #-----如果没有显著的点，就returnNULL，有则继续清洗-----#
      if(nrow(exposure)!=0){ #如果没有显著的点，就跳过#
        res<-lapply(1:nrow(outcome_list), function(c){
          #-------------------读入outcome--------------#
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
          rm(outcome);gc()
          #-----------------clump前先合并------------#
          snp_outcome<-f2_outcome$SNP
          f1_exposure<- exposure[exposure$SNP %in% snp_outcome,]
          if(nrow(f1_exposure)!=0){
            #-------------clump----------------------------#
            f1_exposure %<>% mutate(rsid=SNP) %>%  mutate(pval=pval.exposure)
            f1_exposure_clump<-ld_clump(dat = f1_exposure,
                                        clump_kb = Clump_kb,
                                        clump_r2 = Clump_r2,
                                        pop = MR_ancestry,
                                        bfile =bfile_dir,
                                        plink_bin =plink_dir)
            rm(f1_exposure);gc()
            dat<-harmonise_data(exposure_dat = f1_exposure_clump,outcome_dat =f2_outcome,action=1)###默认不翻，因为清洗的时候翻过了
            #------------MR steiger---------#
            dat<-steiger_filtering(dat)
            #---------------F值计算--------------#
            dat %<>%
              mutate(R2=(2*(beta.exposure^2)*eaf.exposure*(1 - eaf.exposure))) %>% 
              mutate(Fvalue=((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))
            return(dat)
          }else{return(NULL)}
        })
        List_Name<- paste0(exposure_list$abr[e]," to ", outcome_list$abr)
        names(res)<-List_Name
        return(res)
      }else{return(NULL)}
    })
    if(bi_mr){
      out2<-lapply(1:nrow(outcome_list), function(e){
        exposure=paste0(gwas_prefix,outcome_list$absolute_dir[e])
        exposure<-fread(exposure)
        exposure<-as.data.frame(exposure)
        exposure$phenotype=outcome_list$abr[e]
        #----twosampleMR---#
        exposure<-format_data(dat=exposure, type = "exposure",snps = NULL,header = TRUE,
                              phenotype_col="phenotype",
                              snp_col = "SNP",beta_col = "BETA",se_col = "SE",eaf_col = "EAF",
                              effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P",
                              chr_col = "CHR",pos_col = "POS",samplesize_col="N",
                              min_pval = 1e-500,log_pval = FALSE)
        exposure<-exposure[exposure$pval.exposure< gwas_threshold,]
        gc()
        #-----如果没有显著的点，就returnNULL，有则继续清洗-----#
        if(nrow(exposure)!=0){ #如果没有显著的点，就跳过#
          res<-lapply(1:nrow(exposure_list), function(c){
            #-------读入out_数据----#
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
            rm(outcome);gc()
            #-----------------clump前先合并------------#
            snp_outcome<-f2_outcome$SNP
            f1_exposure<- exposure[exposure$SNP %in% snp_outcome,]
            if(nrow(f1_exposure)!=0){
              #-------------clump----------------------------#
              f1_exposure %<>% mutate(rsid=SNP) %>%  mutate(pval=pval.exposure)
              f1_exposure_clump<-ld_clump(dat = f1_exposure,
                                          clump_kb = Clump_kb,
                                          clump_r2 = Clump_r2,
                                          pop = MR_ancestry,
                                          bfile =bfile_dir,
                                          plink_bin =plink_dir)
              rm(f1_exposure);gc()
              dat<-harmonise_data(exposure_dat = f1_exposure_clump,outcome_dat =f2_outcome,action=1)
              #------------MR steiger---------#
              dat<-steiger_filtering(dat)
              #---------------F值计算--------------#
              dat %<>%
                mutate(R2=(2*(beta.exposure^2)*eaf.exposure*(1 - eaf.exposure))) %>% 
                mutate(Fvalue=((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))
              return(dat)
            }else(return(NULL))
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
                   iv_rdata_name)
    tmp_list_name<-list("ivs","exposure_list","outcome_list","Clump_kb","Clump_r2",
                        "MR_ancestry","gwas_threshold","bfile_dir","plink_dir",
                        "iv_rdata_name")
    names(tmp_list)<-tmp_list_name
    return(tmp_list)
  }else{
    exposure<-as.data.frame(exposure)
    exposure$phenotype="exposure"
    #----twosampleMR---#
    exposure<-format_data(dat=exposure, type = "exposure",snps = NULL,header = TRUE,
                          phenotype_col="phenotype",
                          snp_col = "SNP",beta_col = "BETA",se_col = "SE",eaf_col = "EAF",
                          effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P",
                          chr_col = "CHR",pos_col = "POS",samplesize_col="N",
                          min_pval = 1e-500,log_pval = FALSE)
    exposure<-exposure[exposure$pval.exposure< gwas_threshold,]
    gc()
    #-----如果没有显著的点，就returnNULL，有则继续清洗-----#
    if(nrow(exposure)!=0){ #如果没有显著的点，就跳过#
      #-------------------读入outcome--------------#
      outcome<-as.data.frame(outcome)
      outcome$phenotype="outcome"
      f2_outcome<-format_data(outcome,type = "outcome",snps = NULL,header = TRUE,
                              phenotype_col="phenotype",
                              snp_col = "SNP",beta_col = "BETA",se_col = "SE",eaf_col = "EAF",
                              effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P",
                              chr_col = "CHR",pos_col = "POS",samplesize_col="N",
                              min_pval = 1e-500,log_pval = FALSE)
      rm(outcome);gc()
      #-----------------clump前先合并------------#
      snp_outcome<-f2_outcome$SNP
      f1_exposure<- exposure[exposure$SNP %in% snp_outcome,]
      if(nrow(f1_exposure)!=0){
        #-------------clump----------------------------#
        f1_exposure %<>% mutate(rsid=SNP) %>%  mutate(pval=pval.exposure)
        f1_exposure_clump<-ld_clump(dat = f1_exposure,
                                    clump_kb = Clump_kb,
                                    clump_r2 = Clump_r2,
                                    pop = MR_ancestry,
                                    bfile =bfile_dir,
                                    plink_bin =plink_dir)
        rm(f1_exposure);gc()
        dat<-harmonise_data(exposure_dat = f1_exposure_clump,outcome_dat =f2_outcome,action=1)###默认不翻，因为清洗的时候翻过了
        #------------MR steiger---------#
        dat<-steiger_filtering(dat)
        #---------------F值计算--------------#
        dat %<>%
          mutate(R2=(2*(beta.exposure^2)*eaf.exposure*(1 - eaf.exposure))) %>% 
          mutate(Fvalue=((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))
        return(dat)
      }else{return(NULL)}
    }else{return(NULL)}
  }
}
#----------------1.2 twosampleMR:pheGWAS_iv----------
pheGWAS_twosampleMR_iv<-function(
    pheGWAS_annotation_dir, #######区别于普通iv的
    exposure_list,outcome_list, 
    gwas_prefix,out_path,
    bi_mr=F,
    gwas_threshold=5*10^-8,
    bfile_dir,plink_dir,MR_ancestry,
    Clump_kb = 10000, Clump_r2 = 0.001,
    Action=2 ##harmonise_data里的action
){
  out1<-lapply(1:nrow(exposure_list), function(e){
    pheGWAS_annotation<-fread(paste0(gwas_prefix,pheGWAS_annotation_dir))####读入注释文件
    exposure=paste0(gwas_prefix,exposure_list$absolute_dir[e])
    exposure<-fread(exposure)
    #注释文件和exposure文件合并----#
    exposure<-inner_join(pheGWAS_annotation,exposure,by="SNP")
    exposure$phenotype=exposure_list$abr[e]
    exposure<-as.data.frame(exposure)
    #----twosampleMR---#
    exposure<-format_data(dat=exposure, type = "exposure",snps = NULL,header = TRUE,
                             phenotype_col="phenotype",
                             snp_col = "SNP",beta_col = "BETA",se_col = "SE",eaf_col = "EAF",
                             effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P",
                             chr_col = "CHR",pos_col = "POS",samplesize_col="N",
                             min_pval = 1e-500,log_pval = FALSE)
    exposure<-exposure[exposure$pval.exposure< gwas_threshold,]
    gc()
    #-----如果没有显著的点，就returnNULL，有则继续清洗-----#
    if(nrow(exposure)!=0){ #如果没有显著的点，就跳过#
      res<-lapply(1:nrow(outcome_list), function(c){
        #-------读入out_数据----#
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
        rm(outcome);gc()
        #------clump前先合并-----------#
        snp_outcome<-f2_outcome$SNP
        f1_exposure<- exposure[exposure$SNP %in% snp_outcome,]
        if(nrow(f1_exposure)!=0){
          #-------------clump----------------------------#
          f1_exposure %<>% mutate(rsid=SNP) %>%  mutate(pval=pval.exposure)
          f1_exposure_clump<-ld_clump(dat = f1_exposure,
                                      clump_kb = Clump_kb,
                                      clump_r2 = Clump_r2,
                                      pop = MR_ancestry,
                                      bfile =bfile_dir,
                                      plink_bin =plink_dir)
          rm(f1_exposure);gc()
          dat<-harmonise_data(exposure_dat = f1_exposure_clump,outcome_dat =f2_outcome,action=1)###默认不翻，因为清洗的时候翻过了
          #------------MR steiger---------#
          dat<-steiger_filtering(dat)
          #---------------F值计算--------------#
          dat %<>%
            mutate(R2=(2*(beta.exposure^2)*eaf.exposure*(1 - eaf.exposure))) %>% 
            mutate(Fvalue=((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))
          return(dat)
        }else{return(NULL)}
      })
      List_Name<- paste0(exposure_list$abr[e]," to ", outcome_list$abr)
      names(res)<-List_Name
      return(res)
    }else{
      return(NULL)
    }
  })
  if(bi_mr){
    out2<-lapply(1:nrow(outcome_list), function(e){
      exposure=paste0(gwas_prefix,outcome_list$absolute_dir[e])
      exposure<-fread(exposure)
      exposure<-as.data.frame(exposure)
      exposure$phenotype=outcome_list$abr[e]
      #----twosampleMR---#
      exposure<-format_data(dat=exposure, type = "exposure",snps = NULL,header = TRUE,
                               phenotype_col="phenotype",
                               snp_col = "SNP",beta_col = "BETA",se_col = "SE",eaf_col = "EAF",
                               effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P",
                               chr_col = "CHR",pos_col = "POS",samplesize_col="N",
                               min_pval = 1e-500,log_pval = FALSE)
      exposure<-exposure[exposure$pval.exposure< gwas_threshold,]
      gc()
      #-----如果没有显著的点，就returnNULL，有则继续清洗-----#
      if(nrow(exposure)!=0){ #如果没有显著的点，就跳过#
        res<-lapply(1:nrow(exposure_list), function(c){
          #-------读入out_数据----#
          pheGWAS_annotation<-fread(paste0(gwas_prefix,pheGWAS_annotation_dir))####读入注释文件
          outcome=paste0(gwas_prefix,exposure_list$absolute_dir[c])
          outcome<-fread(outcome)
          #注释文件合并----#
          outcome<-inner_join(pheGWAS_annotation,outcome,by="SNP")
          outcome$phenotype=exposure_list$abr[c]
          outcome<-as.data.frame(outcome)
          f2_outcome<-format_data(outcome,type = "outcome",snps = NULL,header = TRUE,
                                  phenotype_col="phenotype",
                                  snp_col = "SNP",beta_col = "BETA",se_col = "SE",eaf_col = "EAF",
                                  effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P",
                                  chr_col = "CHR",pos_col = "POS",samplesize_col="N",
                                  min_pval = 1e-500,log_pval = FALSE)
          rm(outcome);gc()
          #-----------------clump前先合并------------#
          snp_outcome<-f2_outcome$SNP
          f1_exposure<- exposure[exposure$SNP %in% snp_outcome,]
          if(nrow(f1_exposure)!=0){
            #-------------clump----------------------------#
            f1_exposure %<>% mutate(rsid=SNP) %>%  mutate(pval=pval.exposure)
            f1_exposure_clump<-ld_clump(dat = f1_exposure,
                                        clump_kb = Clump_kb,
                                        clump_r2 = Clump_r2,
                                        pop = MR_ancestry,
                                        bfile =bfile_dir,
                                        plink_bin =plink_dir)
            rm(f1_exposure);gc()
            dat<-harmonise_data(exposure_dat = f1_exposure_clump,outcome_dat =f2_outcome,action=1)###默认不翻，因为清洗的时候翻过了
            #------------MR steiger---------#
            dat<-steiger_filtering(dat)
            #---------------F值计算--------------#
            dat %<>%
              mutate(R2=(2*(beta.exposure^2)*eaf.exposure*(1 - eaf.exposure))) %>% 
              mutate(Fvalue=((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))
            return(dat)
          }else{return(NULL)}
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
#----------------1.3 twosampleMR_Primary----------
GWAS_twosampleMR_Primary<-function(out_path,mr_name,bi_mr=F){
  rdata_dir<-paste0(out_path,mr_name,".Rdata")
  load(rdata_dir)
  IVs<-ivs$ivs
  exposure_list<-ivs$exposure_list
  outcome_list<-ivs$outcome_list
  exp_c<-exposure_list$abr
  ouc_c<-outcome_list$abr
  
  if(bi_mr){
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
  }else{
    tmp1<-lapply(1:length(exp_c), function(e){
      exp_tmp<-exp_c[e]
      tmp_c<-paste0(exp_tmp," to ",ouc_c)
    })
    list_name<-unlist(tmp1)
  }
  #-------提出结果-------------#
  MR_res<-lapply(list_name, function(i){
    dat=IVs[[i]]
    if(length(dat)==0){
      twosamplemr_res<-NULL
      heterogeneity_res<-NULL
      pleiotropy_res<-NULL
      leaveone_out<-NULL
      mr_singleSNP<-NULL
    }else{
      if(nrow(dat)==0){
        twosamplemr_res<-NULL
        heterogeneity_res<-NULL
        pleiotropy_res<-NULL
        leaveone_out<-NULL
        mr_singleSNP<-NULL
      }else{
        dat %<>% dplyr::filter(ambiguous=="FALSE")
        if(nrow(dat)>=2){
          #----------twosamoleMR------#
          twosamplemr_res<-mr(dat)
          heterogeneity_res<-mr_heterogeneity(dat)
          pleiotropy_res<-mr_pleiotropy_test(dat)
          leaveone_out<-mr_leaveoneout(dat)
          ivs_isq<-Isq(dat$beta.exposure,dat$se.exposure)
          mr_singleSNP<-mr_singlesnp(dat)
          heterogeneity_res[3,"method"]="I2"
          heterogeneity_res[3,"I2"]=ivs_isq
          
        }else{
          if(nrow(dat)==1){
            twosamplemr_res<-mr(dat)
            heterogeneity_res<-NULL
            pleiotropy_res<-NULL
            leaveone_out<-NULL
            mr_singleSNP<-NULL
          }else{
            if(nrow(dat)==0){
              twosamplemr_res<-NULL
              heterogeneity_res<-NULL
              pleiotropy_res<-NULL
              leaveone_out<-NULL
              mr_singleSNP<-NULL
            }
          }
        }
      }
    }
    return(list(mr_res=twosamplemr_res,
                leave_one_out=leaveone_out,
                heterogeneity_res=heterogeneity_res,
                pleiotropy_res=pleiotropy_res,
                mr_singleSNP=mr_singleSNP))
  })
  names(MR_res)<-list_name
  return(MR_res)
}

#----------------1.4 mrpresso_parallel----------
GWAS_mrpresso_parallel<-function(out_path,
                                 mr_name,
                                 ncore,
                                 list_name=NULL #若不设置目标表型对，则把所有表型对都跑了
                                 ){
  library(TwoSampleMR)
  library(parallel)
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.0\n")
  cat("* Update Date:2024.08.08\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  if(is.null(list_name)){
    load(paste0(Out_path,iv_rdata_name,".Rdata"))
    IVs<-ivs$ivs
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
  }
  # 设置多核并行计算
  cl <- parallel::makeCluster(ncore)
  clusterExport(cl, c("Out_path","iv_rdata_name"))
  clusterEvalQ(cl, {
    library(TwoSampleMR)
    library(magrittr)
    library(tidyverse)
    load(paste0(Out_path,iv_rdata_name,".Rdata"))
    IVs<-ivs$ivs
  })
  MR_res <- parSapply(cl, list_name, function(i){
    dat=IVs[[i]]
    if(length(dat)==0){
      mr_presso<-NULL
    }else{
      if(nrow(dat)==0){
        mr_presso<-NULL
      }else{
        dat %<>% dplyr::filter(ambiguous=="FALSE")
        dat %<>% dplyr::filter(steiger_dir=="TRUE")
        if(nrow(dat)>3){
          mr_presso<-run_mr_presso(dat,NbDistribution = 3000,SignifThreshold = 0.05)
        }else{
          mr_presso<-NULL
        }
      }
    }
    return(presso_res=mr_presso)
  })
  stopCluster(cl)
  names(MR_res)<-list_name
  mr_presso<-lapply(list_name, function(i){
    tmp<-MR_res[[i]]$`Main MR results`
    if(is.null(tmp)==F){
      tmp$e2o<-i
      return(tmp)
    }
  })
  mr_presso<-bind_rows(mr_presso)
  write.csv(mr_presso,file = paste0(out_path,"mr_presso.csv"))
  return(mr_presso)
}
#----------------2 MR:GSMR2----------------
GWAS_GSMR2<-function(gcta_dir,
                     gcta_bfile,
                     exposure_dir_dat,
                     outcome_dir_dat,
                     gsmr_direction=2,
                     out_path,
                     file_name,
                     gwas_prefix){ 
  
  cat("*******************************************************************\n")
  cat("* CWAS \n")
  cat("* Version 1.9\n")
  cat("* Update Date:2024.07.25\n")
  cat("* (C) 2024 Chao Wang \n")
  cat("* Nanjing Medical University\n")
  cat("*******************************************************************\n")
  #-----------创建临时路径---------#
  temp_path <- tempdir()
  
  #------smr文件制作,读入临时路径----------#
  for (i in 1:nrow(exposure_dir_dat)) {
    gsmr_gwas<-fread(paste0(gwas_prefix,exposure_dir_dat$absolute_dir[i]))
    gsmr_gwas %<>%
      dplyr::rename(freq=EAF,b=BETA,se=SE,p=P) %>% 
      dplyr::select(SNP,A1,A2,freq,b,se,p,N)
    
    fwrite(gsmr_gwas,
           file = paste0(temp_path,"\\",exposure_dir_dat$abr[i],".raw"),
           sep = "\t",col.names = T,row.names = F,quote = F)
  }
  
  for (i in 1:nrow(outcome_dir_dat)) {
    gsmr_gwas<-fread(paste0(gwas_prefix,"/",outcome_dir_dat$absolute_dir[i]))
    gsmr_gwas %<>%
      dplyr::rename(freq=EAF,b=BETA,se=SE,p=P) %>% 
      dplyr::select(SNP,A1,A2,freq,b,se,p,N)
    
    fwrite(gsmr_gwas,
           file = paste0(temp_path,"\\",outcome_dir_dat$abr[i],".raw"),
           sep = "\t",col.names = T,row.names = F,quote = F)
  }
  #-----创建暴露与结局的txt文件----------#
  exposure_dir_dat %<>%
    mutate(input_dir=paste0(temp_path,"\\",abr,".raw")) %>% 
    select(abr,input_dir)
  write.table(exposure_dir_dat,
              file = paste0(temp_path,"\\","gsmr_exposure.txt"),
              col.names = F,row.names = F,sep = "\t",quote = F)
  
  outcome_dir_dat %<>%
    mutate(input_dir=paste0(temp_path,"\\",abr,".raw")) %>% 
    select(abr,input_dir)
  write.table(outcome_dir_dat,
              file = paste0(temp_path,"\\","gsmr_outcome.txt"),
              col.names = F,row.names = F,sep = "\t",quote = F)
  #------生成命令------------#
  c1<-paste(gcta_dir,
            "--bfile",gcta_bfile,
            "--gsmr-file ", paste0(temp_path,"\\","gsmr_exposure.txt"),paste0(temp_path,"\\","gsmr_outcome.txt"),
            "--gsmr-direction",gsmr_direction,
            "--out",paste0(out_path,"/",file_name))
  system(c1)
  #-----关闭临时文件-------#
  unlink(temp_path, recursive = TRUE)
  
}
