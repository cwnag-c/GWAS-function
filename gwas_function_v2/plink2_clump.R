ld_clump_qtl<-function(dat,
                       clump_r2=0.01,
                       clump_p=1,
                       bfile,
                       plink_bin,
                       clump_kb=250){
  temp_path <- tempdir()# 创建临时路径
  dat1<-dat %>% select(SNP,pval.exposure) %>% rename(P=pval.exposure)
  write.table(dat1,file = paste0(temp_path,"\\tmp.input.txt"),
              col.names = T,row.names = F,quote = F,sep = "\t")
  
  c<-paste(
    plink_bin,
    "--bfile",bfile,
    "--clump", paste0(temp_path,"\\tmp.input.txt"),
    "--clump-p1",clump_p,
    "--clump-r2",clump_r2,
    "--clump-kb",clump_kb,
    "--out",paste0(temp_path,"\\tmp.output.txt")
  )
  system(c)
  clump_res<-fread(paste0(temp_path,"\\tmp.output.txt.clumps"))
  clump_mis<-read.table(paste0(temp_path,"\\tmp.output.txt.clumps.missing_id"),col.names = F)
  if(!is.null(clump_mis)){
    snp_mis<-as.character(clump_mis$FALSE.)
    dat %<>% filter(!(SNP %in% snp_mis))
    dat1<-dat %>% select(SNP,pval.exposure) %>% rename(P=pval.exposure)
    write.table(dat1,file = paste0(temp_path,"\\tmp.input.txt"),
                col.names = T,row.names = F,quote = F,sep = "\t")
    
    c<-paste(
      plink_bin,
      "--bfile",bfile,
      "--clump", paste0(temp_path,"\\tmp.input.txt"),
      "--clump-p1",clump_p,
      "--clump-r2",clump_r2,
      "--clump-kb",clump_kb,
      "--out",paste0(temp_path,"\\tmp.output.txt")
    )
    system(c)
  }
  
  
}

