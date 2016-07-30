#################################    EOT AND PCA ROTATION  #######################################
########################### SPACE-TIME VARIABILITY  #############################################
#This script examines the 1982-2016 dataset downloaded from NOAA:
#http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
#
#The goal is to run EOT, MEOT and PCA on the updated SST dataset using the IDRISI and the R package "remote".
#
#AUTHOR: Benoit Parmentier                                                                       #
#DATE CREATED:07/15/2016 
#DATE MODIFIED: 07/30/2016
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
library(gridExtra)
library(latticeExtra)
library(colorRamps)                         # contains matlab.like palette
library(lsr)
library(psych)
library(GPArotation)
library(zoo)
library(xts)
library(remote)                            # EOT implementation in R/cpp
library(XML)                               # HTML funcitons

#################################################
###### Functions  used in the script  ##########

create_dir_fun <- function(outDir,out_suffix){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

convert_to_numeric <-function(x){
  if(class(x)=="character"){
    x<-as.numeric(x)
  }
  ##
  if(class(x)=="factor"){
    x<-as.numeric(as.character(x))
  }
  return(x)
}


crosscor_lag_analysis_fun<-function(telind,mode_list,d_z,lag_window,fig,out_suffix){
  #This function crosss correlates between two sets of time series given some lag window.
  #The inputs are time series objects of the type zoo.
  #
  #Input Arguments:
  #1)telind: time series 1 as character vector
  #2)modelist: time series 2 as character vector
  #3)d_z: zoo object 
  #4)lag_window:
  #5)fig: TRUE or FALSE, if TRUE generate crossorrrelation figure
  #6)out_suffix
  ##AUTHOR: Benoit Parmentier                                                                       #
  #DATE CREATED:12/15/2013 
  #DATE MODIFIED: 07/29/2016
  ###Comments
  ### Last update on 07/29/2016, this needs improvement!!
  
  lag_table_ext<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))
  lag_table_lag<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list))
  lag_table_text<-matrix(data=NA,nrow=length(telind),ncol=length(mode_list)) #Formatted table used in the paper
  #lag_cross_cor_PCA<-vector("list",length(mode_list))
  lag_m<-seq(-1*lag_window,lag_window,1)
  #retain ccf!!!
  #list_ccf_lag_table
  
  #lag_cross_cor_PCA_m<-array(data=NA,nrow=length(lag_m),ncol=length(mode_list))
  for (i in 1:length(telind)){
    telindex <- telind[i]
    pos1<-match(telindex,names(d_z))
    #retain ccf!!!
    for (j in 1:length(mode_list)){
      mode_n<-mode_list[j]
      pos2<-match(mode_n,names(d_z))
      sum(is.na(d_z[,pos2]))
      
      ccf_obj<-ccf(as.numeric(d_z[,pos1]),as.numeric(d_z[,pos2]), lag=lag_window,na.action=na.pass)  #Note that ccf does not take
      
      lag_m<-seq(-1*lag_window,lag_window,1)
      ccf_obj$lag[,1,1]<-lag_m  #replacing lag values because continuous
      
      if (fig=="TRUE"){
        plot_name<-paste(telindex, "and", mode_n,"lag analysis",sep="_")
        png(paste(plot_name,"_",out_suffix,".png", sep=""))
        plot(ccf_obj, main= paste(telindex, "and", mode_n,"lag analysis",sep=" "), ylab="Cross-correlation",
             xlab="Lag (month)", ylim=c(-1,1))
        dev.off()
      }
      
      ######### NOW FIND THE m
      absext <-max(abs(ccf_obj$acf)) # maximum of the extremum
      pos<-match(absext,ccf_obj$acf) #find the position and lag, if NA it means it was negative
      if (is.na(pos)) {
        pos<-match(absext*-1,ccf_obj$acf)
        absext<-absext*-1   #recover the sign
      } 
      absext_lag<-ccf_obj$lag[pos,1,1] #This is the lag corresponding to the maximum absolute value
      
      lag_table_ext[i,j]<-absext
      lag_table_lag[i,j]<-absext_lag
      #number<-format(absext,digits=3)
      ext<-round(absext,digits=3)
      element<-paste(ext," (",absext_lag,")",sep="")
      lag_table_text[i,j]<-element
      
      ##Keep ccf lag somewhere
    }
  }
  
  lag_table_ext<-as.data.frame(lag_table_ext)
  names(lag_table_ext)<-mode_list
  rownames(lag_table_ext)<-telind
  file_name<-paste("lag_table_extremum_window_", lag_window,"_",out_suffix,".txt",sep="")
  write.table(lag_table_ext,file=file_name,sep=",")
  
  lag_table_lag<-as.data.frame(lag_table_lag)
  names(lag_table_lag)<-mode_list
  rownames(lag_table_lag)<-telind
  file_name<-paste("lag_table_lag_extremum_window_", lag_window,"_",out_suffix,".txt",sep="")
  write.table(lag_table_lag,file=file_name,sep=",")
  
  lag_table_text<-as.data.frame(lag_table_text)
  names(lag_table_text)<-mode_list
  rownames(lag_table_text)<-telind
  file_name<-paste("lag_table_lag_ext_text", lag_window,"_",out_suffix,".txt",sep="")
  write.table(lag_table_text,file=file_name,sep=",")
  
  #create return object
  
  crosscor_obj<-list(lag_table_ext,lag_table_lag,lag_table_text)
  names(crosscor_obj)<-c("extremum","lag_ext","text")
  file_name<-paste("crosscor_obj_lag_analysis_", lag_window,"_",out_suffix,".RData",sep="")
  save(crosscor_obj,file=file_name)
  
  return(crosscor_obj)
}

read_comp_results <- function(data_filename,sheet_name,dseq){
  decomp_dat <- read_ods(data_filename,sheet=sheet_name)
  names(decomp_dat)[1] <- "fnames"
  names(decomp_dat)[2:ncol(decomp_dat)] <- paste0("comp_",1:20)
  test <- subset(decomp_dat,select=paste0("comp_",1:20))
  test<-lapply(test,FUN=convert_to_numeric)
  test<- do.call(cbind,test)
  decomp_dat <- as.data.frame(test)
  decomp_dat_dz <- zoo(decomp_dat,dseq) #create zoo object from data.frame and date sequence object
  time(decomp_dat_dz)  #no time stamp??
  data_obj <- list(decomp_dat_dz,decomp_dat)
  names(data_obj) <- c("dat_dz","dat")
  return(data_obj)
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
      
      #cor_val <- cor(tfs1,tfs2,use="complete.obs")
      cor_val <- cor(tfs1,tfs2,use="na.or.complete") # casewise deletion of NA and assign NA if no complete case for a vector/matrix column
      ### add this as time series later?
      
      if (fig=="TRUE"){
        plot_name<-paste(names_tfs1, "and",names_tfs2,"scatter_plot_cor_analysis",sep="_")
        png(paste(plot_name,"_",out_suffix,".png", sep=""))
        plot(tfs2 ~ tfs1, main= paste(names_tfs1, "and", names_tfs2 ,"cor analysis",sep=" "), 
             ylab=names_tfs2,
             xlab=names_tfs1)
        dev.off()
        
        ## Temporal profiles of both time series
        y_range <- range(c(tfs1,tfs2))
        plot_name<-paste(names_tfs1, "and",names_tfs2,"series_profile_cor_analysis",sep="_")
        
        png(paste(plot_name,"_",out_suffix,".png", sep=""))
        plot(tfs1, main= paste(names_tfs1, "and", names_tfs2,"cor analysis",sep=" "), 
             ylab="Series",
             xlab="time steps",
             y_lim=y_range,
             type="l")
        par(new = TRUE)
        plot(tfs2, type = "l", col="blue", axes = FALSE, bty = "n", xlab = "", ylab = "")
        axis(side=4, at = pretty(range(tfs2)))
        mtext("z", side=4, line=3)
        
        legend("topleft",legend=c(names_tfs1,names_tfs2), cex=0.8, 
               #col=c("blue","darkgreen"),
               col=c("black","blue"),
               lty=c(1,2),lwd=2,
               bty="n")  #lwd=line width
        
        #lines(tfs2, col="blue")
        
        dev.off()
        
      }
      
      ######### NOW FIND THE m
      
      #cor_table[i,j] <- format(cor_val,digits=4)
      cor_table[i,j] <- as.numeric((format(cor_val,digits=4)))
      
      ##Keep ccf lag somewhere
      
    }
  }
  #browser()
  cor_table <-as.data.frame(cor_table)
  names(cor_table) <- names(ts2)
  rownames(cor_table)<- names(ts1)
  file_name <- paste("cor_table_",out_suffix,".txt",sep="")
  write.table(cor_table,file=file_name,sep=",")
  
  return(cor_table)
  
}

generate_barplot_comparison_fun <- function(df1,df2,out_suffix,col_palette=NULL,out_dir=NULL){
  
  #df1: data.frame 1 containing correlation values related to indices/variables
  #df2: data.frame 2 (e.g. pca/eot) containing correlation values related to indices/variables (e.g. telconnections)
  
  #i<-1 #right now works on columns only
  lf_png <- vector("list",length=ncol(df1))
  
  for(i in 1:ncol(df1)){
    
    #names_ref1 <- names(cor_eot_df)[i]
    #names_ref2 <- names(cor_pca_df)[i]
    #names_ind <- rownames(cor_eot_df)
    names_ref1 <- names(df1)[i]
    names_ref2 <- names(df2)[i]
    names_ind <- rownames(df2)  
    data_plot <- (cbind(as.numeric(df1[,i]),as.numeric(df2[,i]))) # make sure it is numeric
    #data_plot <- (cbind(as.numeric(cor_eot_df[,i]),as.numeric(cor_pca_df[,i]))) # make sure it is numeric
    
    #data_plot <- as.numeric((cor_eot_df[,i])) # make sure it is numeric
    
    png_file_name <- paste(names_ref1,"_",names_ref2,"_",out_suffix,".png",sep="")
    
    png(png_file_name)
    
    heights<-as.matrix(t(data_plot))
    barplot(heights,     #data to plot
            main=paste(names_ref1," and ",names_ref2,sep=""),
            #main= names_ref,
            names.arg=names_ind,cex.names=0.8,   #names of the teleconnections indices and size of fonts of axis labes
            beside=TRUE,                         # see two barplots for comparisons...
            xlab="Teleconnection Indices",       # font.lab is 2 to make the font bold
            ylab="Correlation",font.lab=2,
            col=c("blue","red"), ylim=c(-1,1))
    #col=c("blue"), ylim=c(-1,1))
    grid(nx=12,ny=10)      
    
    legend("topright",legend=c(names_ref1,names_ref2), cex=0.9, fill=c("blue","red"),bty="n")
    
    box()
    dev.off()
    
    lf_png[[i]] <- png_file_name
    
  }
  
  ## return list of files created
  return(lf_png)
}


########################## END OF SCRIPT #########################

#data(vdendool) #data of 36 cols and 14 rows! very small

## claculate 2 leading modes
#nh_modes <- eot(x = vdendool, y = NULL, n = 2,
#                standardised = FALSE,
#                verbose = TRUE)
#plot(nh_modes, y = 1, show.bp = TRUE)
#plot(nh_modes, y = 2, show.bp = TRUE)
