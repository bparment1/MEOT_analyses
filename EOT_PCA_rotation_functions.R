#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script analyzes  SST data.
#The goal is to compare EOT and PCA with rotation and without
#This is part of the Time Series analysis project started at IDRISI.
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED: 07/02/2015 
#DATE MODIFIED: 01/09/2016
#
#PROJECT: MEOT/EOT climate variability extraction
#
##################################################################################################
#
###Loading r library and packages

library(raster)                             # loading the raster package
library(gtools)                             # loading ...
library(sp)                                 # spatial objects in R
library(gplots)                             # 
library(rgdal)                              # gdal driver for R
library(RColorBrewer)                       # color scheme, palettes used for plotting
library(gdata)                              # read different format (including .xlsx)
library(plotrix)                            # plot options and functions including plotCI
library(rasterVis)                          # raster visualization
library(colorRamps)                         # contains matlab.like palette
library(lsr)
library(psych)
library(GPArotation)
library(zoo)
library(xts)

#################################################
###### Functions  used in the script  ##########

#[1] "comparison_pca_eot_fun" 
#[2] "cor_series_fun"   
#[3] "cv_test_fun"           
#[4] "pca_to_raster_fun"      
#[5] "run_pca_fun"           

cv_test_fun <- function(x,y) {
  #modified from:
  #http://www.r-bloggers.com/example-8-39-calculating-cramers-v/
  chisq_val <- chisq.test(x, y, correct=FALSE)
  CV = sqrt(chisq_val$statistic /
              (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  print.noquote("Cram?r V / Phi:")
  CV <- as.numeric(CV)
  cv_obj<- list(CV=CV,chisq_val=chisq_val)
  return(cv_obj)
}


pca_to_raster_fun<-function(pc_spdf,ref_raster,NA_val,out_suffix){
  #Input arguments:
  #pc_spdf: spdf with scores components and
  #         must include x,y in the last 2 columns
  #ref_raster: reference raster giving the output extent and NA if present
  #NA_val<- values assigned to NA pixels, if "NA", use default
  
  npc<-ncol(pc_spdf)-2 #number of components 
  pc_scores_lf<-vector("list",npc) 
  for (k in 1:npc){
    pc_scores<-pc_spdf[,k] #vector containing scores for components k
    pc_name<-names(pc_spdf)[k] #name of the current component
    raster_name<-paste("pc_component_",k,"_",out_suffix,".rst",sep="")
    
    if (NA_val=="NA"){
      pc_lag<-rasterize(pc_scores,ref_raster,pc_name,fun=min,overwrite=TRUE,
                        filename=raster_name)
    }
    #Use provided NA value for the background
    if (NA_val!="NA"){
      pc_lag <- rasterize(pc_scores,ref_raster,pc_name,background=NA_val,fun=min,overwrite=TRUE,
                        filename=raster_name)
      NAvalue(pc_lag) <- NA_val
    }
    pc_scores_lf[[k]]<-raster_name
  } 
  return(pc_scores_lf)
}

run_pca_fun <- function(A,mode="T",rotation_opt="none",matrix_val=NULL,npc=1,loadings=TRUE,scores_opt=TRUE){
  
  #Arguments
  #mode: T, s mode or other mode , note that this option is not in use at this moment, use matrix_val
  #rotation_opt: which option to use to rotate pc
  #npc: number of components produced
  #matrix_val: square matrix used in the computation of eigen values, vectors
  #A: data matrix
  
  #######
  
  
  #if(s.null(matrix_val)){
  #  as.matrix(df_val)
  #}
  
  ##Cross product PCA S modes
  
  pca_principal_obj <-principal(r=matrix_val, #
                                nfactors = npc, 
                                residuals = FALSE, 
                                covar=TRUE, #use covar option...
                                rotate = rotation_opt,
                                scores=scores_opt,
                                oblique.scores=TRUE,
                                method="regression")
  
  principal_loadings_df <- as.data.frame(unclass(pca_principal_obj$loadings)) # extract the matrix of ??
  #Add scores...
  principal_scores <- predict(pca_principal_obj,A)
  
  ###
  obj_principal <- list(pca_principal_obj,principal_loadings_df,principal_scores,A)
  names(obj_principal) <- c("pca_principal","loadings","scores","data_matrix")
  return(obj_principal)
  
}

cor_series_fun <- function(ts1,ts2,fig=F,out_suffix){
  
  cor_table <- matrix(data=NA,nrow=length(ts1),ncol=length(ts2))
  #ag_table_lag<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))
  #lag_table_text<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list)) #Formatted table used in the paper
  
  for(i in 1:ncol(ts1)){
    tfs1 <- ts1[,i]
    #retain ccf!!!
    names_tfs1 <- names(ts1)[i]
    
    for(j in 1:ncol(ts2)){
      tfs2 <- ts2[,j]
      names_tfs2 <- names(ts2)[j]
      
      cor_val <- cor(tfs1,tfs2)
      
      if (fig=="TRUE"){
        plot_name<-paste(names_tfs1, "and",names_tfs2,"scatter_plot_cor_analysis",sep="_")
        png(paste(plot_name,"_",out_suffix,".png", sep=""))
        plot(tfs2 ~ tfs1, main= paste(names_tfs1, "and", names_tfs2 ,"cor analysis",sep=" "), 
             ylab=names_tfs2,
             xlab=names_tfs1)
        dev.off()
        
        y_range <- range(c(tfs1,tfs2))
        plot_name<-paste(names_tfs1, "and",names_tfs2,"series_profile_cor_analysis",sep="_")
        plot(tfs1, main= paste(names_tfs1, "and", names_tfs2,"cor analysis",sep=" "), 
             ylab="Series",
             xlab="time steps",
             y_lim=y_range,
             type="l")
        par(new = TRUE)
        plot(tfs2, type = "l", col="blue", axes = FALSE, bty = "n", xlab = "", ylab = "")
        axis(side=4, at = pretty(range(tfs2)))
        mtext("z", side=4, line=3)
        #lines(tfs2, col="blue")
             
        dev.off()
        
      }
      
      ######### NOW FIND THE m
      
      cor_table[i,j] <- format(cor_val,digits=4)
      
      ##Keep ccf lag somewhere
      
    }
  }
  
  cor_table <-as.data.frame(cor_table)
  names(cor_table) <- names(ts2)
  rownames(cor_table)<- names(ts1)
  file_name <- paste("cor_table_",out_suffix,".txt",sep="")
  write.table(cor_table,file=file_name,sep=",")
  
  return(cor_table)
  
}

#This is the main function used in the comparison between EOT and PCA with rotation
comparison_pca_eot_fun<- function(list_no_component,data_matrix,rotation_opt,fig_opt,n_cores,out_suffix){
  ####
  ##  
  #This functions runs pca with or without rotation and compares loadings to a ref_dataset
  #PARAMETERS:
  #data_matrix: data matrix before cross product
  #cross_matrix: matrix used in the PCA
  #rotation_opt: rotation used in pca
  #ref_data: data.frame with loadings or series of the same size to be compared to (e.g. eot)
  #fig_opt: figures to generate
  #list_no_component: components to generate as a list
  #n_cores: number of cores used in the parallel processing, not in used yet
  #out_suffix: output suffix for files generated
  #
  ##### BEGIN SCRIPT ######
  
  out_suffix_str <- paste(rotation_opt,out_suffix,sep="_")
  list_pca <- vector("list",length=length(list_no_component))
  l_loadings <- vector("list",length=length(list_no_component))
  list_cor_df <- vector("list",length=length(list_no_component))
  
  for(i in 1:length(list_no_component)){
    
    no_component <- list_no_component[i]
    out_suffix_str <- paste("pc_components_",no_component,out_suffix_str,sep="")
    list_pca[[i]] <- run_pca_fun(A=data_matrix,mode="T",rotation_opt=rotation_opt,
                                 matrix_val=cross_matrix,
                                 npc=no_component,
                                 loadings=TRUE,
                                 scores_opt=TRUE) 
    l_loadings[[i]] <- list_pca[[i]]$loadings
    list_cor_df[[i]] <- cor_series_fun(ref_data,l_loadings[[i]],fig=fig_opt,out_suffix_str)
    
  }
  
  ### prepare object to return
  comp_pca_eot_obj <- list(list_pca,list_cor_df,list_cor_df)
  names(comp_pca_eot_obj) <- c("list_pca","list_cor","list_cor_df")
  file_name<-paste("cor_analysis_obj_",out_suffix_str,".RData",sep="")
  save(comp_pca_eot_obj,file=file_name)
  
  return(comp_pca_eot_obj)
}

barplot_comparison_eot_pca <- function(comp_eot_pca_rotated_obj,
                                       comp_eot_pca_unrotated_obj,list_no_component){
  #comp_eot_pca_rotated_obj=comp_pca_eot_s_varimax_rotated_obj,
  #comp_eot_pca_unrotated_obj=comp_pca_eot_s_unrotated_obj
  list_cor_df <- vector("list",length=length(list_no_component))
  for(i in 1:length(list_no_component)){
    
    res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
    fig_name <- paste("Figure_barplot_cor_eot_pca_","comparison_with_varimax_and_truncation_",list_no_component[[i]],out_suffix,".png",sep="")
    png(filename= fig_name,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    diag_rotated <- diag(as.matrix(comp_eot_pca_rotated_obj$list_cor[[i]]))
    diag_unrotated <- diag(as.matrix(comp_eot_pca_unrotated_obj$list_cor[[i]]))
    
    comp_trunc1 <- as.vector(cbind(diag_rotated,diag_unrotated))
    cor_df <- as.data.frame(data.frame(cor=comp_trunc1))
    cor_df$type <- c(rep("r",10),rep("u",10))
    names_components <- 1:10
    #barplot(rbind(diag_rotated,diag_unrotated))
    
    barplot(rbind(as.numeric(diag_rotated),as.numeric(diag_unrotated)),main="",
            xlab="component", col=c("darkblue","red"),ylim=c(-0.8,0.8),
            legend = c("rotated","unrotated"), beside=TRUE,
            names.arg=names_components)
    title(paste(main="Correlation between EOT and PCA S mode: trunction ",list_no_component[[i]]))
    dev.off()
    
    list_cor_df[[i]] <- cor_df
  }
  return(list_cor_df)
}





comparison_eot_pca_corr_with_rotation_fun <- function(comp_eot_pca_rotated_obj,comp_eot_pca_unrotated_obj,list_no_component,eot_no){
  
  list_diag_comp_tmp <- vector("list",length=length(list_no_component))
  for(i in 1:length(list_no_component)){
    diag_comp1 <- as.numeric(diag(as.matrix(comp_eot_pca_rotated_obj$list_cor[[i]])))
    diag_comp2 <- as.numeric(diag(as.matrix(comp_eot_pca_unrotated_obj$list_cor[[i]])))
    diag_comp_tmp <- rbind(diag_comp1,diag_comp2)
    list_diag_comp_tmp[[i]] <- diag_comp1
  }
  
  cor_table <- as.data.frame(do.call(rbind,list_diag_comp_tmp))
  names(cor_table)<- paste("EOT",1:eot_no,sep="_")
  #cor_table <- (do.call(rbind,list_diag_comp_tmp))
  
  for(i in 1:eot_no){
    
    #show that the correlation between pca rotated 1 and EOT 1 decreases as the number of pc components
    #increases
    eot_name <- paste("EOT",i,sep="_")
    
    res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
    fig_name1 <- paste("Figure_cor_eot_",i,"_pc_",i,"_","comparison_with_varimax_and_truncation_",out_suffix,".png",sep="")
    png(filename= fig_name1,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    ## Without absolute values
    cor_table$truncation_var <- list_no_component
    plot(((cor_table[[eot_name]]))~truncation_var,ylim=c(-1,1),type="b",
         ylab="correlation",xlab="truncation component",data=cor_table)
    abline(h=abs(as.numeric(diag_comp2[i])),lty="dashed") #correlation for unrotated
    abline(h=0,lty="solid") #correlation for unrotated
    title(paste("Correlation between PC",i, "varimax and",eot_name,"with different truncation",sep=" "))
    dev.off()
    
    res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
    fig_name2 <- paste("Figure_absolute_cor_eot_",i,"_pc_",i,"_","comparison_with_varimax_and_truncation_",out_suffix,".png",sep="")
    png(filename= fig_name2,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    ## now using absolute values
    plot((abs(cor_table[[eot_name]]))~truncation_var,ylim=c(-1,1),type="b",
         ylab="correlation",xlab="truncation component",data=cor_table)
    abline(h=abs(as.numeric(diag_comp2[i])),lty="dashed") #correlation for unrotated
    abline(h=0,lty="solid") #correlation for unrotated
    title(paste("Absolute Correlation between PC",i, "varimax and",eot_name,"with different truncation",sep=" "))
    dev.off()
    
  }
  ##
  return(cor_table)
}



