########################################       Analyses of space-time dataset       #######################################
########################################### For analyzing MSSA-MEOT and S-T mode PCA results #####################################
#This script analyzes results from synthetic space time datasets in R to test EOT,MEOT,PCA and MSSA.
#AUTHORS:Benoit Parmentier                                             
#DATE: 10/10/2013                                                                                
#PROJECT: Clark Labs, MEOT, time series             
#TO DO:
#Write processing of tsf avl 
###################################################################################################

###Loading R library and packages                                                      

library(gtools)                                        # loading some useful tools 
library(sp)
library(raster)
library(rasterVis)
library(rgdal)
library(vegan) 
library(zoo)

##### Functions used in this script #####

functions_generation_space_time_datasets <- "generation_of_space_time_dataset_functions_10102013v5.R"

script_path<-"/Users/benoitparmentier/Dropbox/Data/MEOT_paper/synthetic_dataset_08262013" #path to script
source(file.path(script_path,functions_generation_space_time_datasets)) #source all functions used in this script.

############# Parameters and arguments ####################

in_dir<- "/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_Working_dir2"
#in_dir<- "/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir1"
#out_dir <- "/Users/benoitparmentier/Dropbox/MEOT"
#in_dir <- "/Users/benoitparmentier/Dropbox/Data/MEOT_paper/synthetic_dataset_09162013"

out_dir <- in_dir
lag_window<-13
out_suffix <- "_meot_10102013"

############ START SCRIPT HERE ###################

list_components_folders <- list.dirs(path=in_dir) #default uses recursive and full.names
list_components_folders <- mixedsort(grep("*.components$",list_components_folders,value=T))

out_suffix_components_folders <- paste("moving_circular_",c(6,12,24,36,48,60,3),sep="")
#create rgf
j<-1
i<-1

in_dir<- "/Users/benoitparmentier/Documents/Data/Benoit/Clark_University/Paper_writings/MSSA_BNP"

list_components_folders <- list.dirs(path=in_dir) #default uses recursive and full.names
list_components_folders <- mixedsort(grep("*.components$",list_components_folders,value=T))
list_components_folders
out_suffix_components_folders <- paste("standing_pattern_",c(12,24,36,48,60,6,3),sep="")

out_dir<- "/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/MEOT_working_dir1"
plot_components_IDRISI_ETM(list_components_folders,out_suffix_components_folders,out_dir,out_suffix,lag_window)


plot_components_IDRISI_ETM <-function(list_components_folders,out_suffix_components_folders,out_dir,out_suffix,lag_window){
  #Plot output components from ETM in IDRISI
  list_raster_name <- vector("list",length=length(list_components_folders))
  for(j in 1:length(list_components_folders)){
    folder_components <-list_components_folders[j]
    out_suffix_s <- out_suffix_components_folders[j]
    #Extract the number of components
    comp_files <- mixedsort(list.files(pattern="*.Comp.*.R.rst$",list_components_folders[1]))
    tmp_str <- strsplit(comp_files,"_")
    comp_str <-unique(unlist(lapply(tmp_str,function(k){grep(pattern="Comp*",k,value=TRUE)})))
    no_comp <- length(comp_str) #number of components in the analysis
    list_comp_s
    list_raster_name <- vector("list",length=length(list_components_folders))
    for(i in 1:length(comp_str)){
      comp_str_processed <-comp_str[i]
      #create_idrisi_rgf(out_suffix_s=paste(comp_str[i],".rst$",sep=""),file_pat=comp_str[i],in_dir=list_meot_folders[j])
      list_r_comp_s <- create_idrisi_rgf(out_suffix_s=paste(comp_str[i],"_R",sep=""),in_dir=folder_components,ending=TRUE)
      r_comp_s <- stack(mixedsort(list_r_comp_s))
      layerNames(r_comp_s)<-paste(comp_str_processed,1:lag_window,sep="_")
      temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
      #levelplot(r_stack,layers=1:48, col.regions=temp.colors)
      p <- levelplot(r_comp_s, col.regions=temp.colors,main=paste(comp_str_processed,"_",out_suffix_s,out_suffix,sep=""))
      
      res_pix<-480
      col_mfrow<-1
      #row_mfrow<-2
      row_mfrow<-1
      
      png_file_name<- paste("Figure_",comp_str_processed,"_",paste(out_suffix_s,out_suffix,sep=""),".png", sep="")
      png(filename=file.path(out_dir,png_file_name),
          width=col_mfrow*res_pix,height=row_mfrow*res_pix)
      par(mfrow=c(row_mfrow,col_mfrow))
      
      print(p)
      dev.off()
    }
  }
}

### Pattern created


out_suffix <- "gen_09162013"

in_dir <- "/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/synthetic_datasets_09162013"
out_dir<- "/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/synthetic_datasets_09162013"

#example for one series:
out_suffix_s <- paste("r_ts1_",out_suffix,sep="")
out_prefix_str <- "test_60c_60r_"

list_r_stack <- create_idrisi_rgf(out_suffix_s=out_suffix_s,file_pat="",
                                  in_dir=file.path(in_dir,out_suffix_s),ending=FALSE)
r_ts_12_standing <- stack(mixedsort(list_r_stack))
temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
#levelplot(r_stack,layers=1:48, col.regions=temp.colors)
p <- levelplot(r_ts_12_standing, layers=12:25,col.regions=temp.colors,main=paste("standing_12","_",out_suffix_s,out_suffix,sep=""))
print(p)

#example for one series:
out_suffix_s <- paste("r_ts6_12_circular_",out_suffix,sep="")
out_prefix_str <- "test_60c_60r_"

list_r_stack <- create_idrisi_rgf(out_suffix_s=out_suffix_s,file_pat="",
                                  in_dir=file.path(in_dir,out_suffix_s),ending=FALSE)
r_ts_12_circular <- stack(mixedsort(list_r_stack))

temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
#levelplot(r_stack,layers=1:48, col.regions=temp.colors)
p <- levelplot(r_ts_12_circular, layers=12:25,col.regions=temp.colors,main=paste("circular_12","_",out_suffix_s,out_suffix,sep=""))
print(p)

## write function for avl
##write function to combine temp and spatial patterns in figure


#### SYNTHETIC DATASET 3x3 and 9x9 pulses

in_dir <- "/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/9x9dataset"
out_suffix<-"t"
#example for one series:
out_suffix_s <- paste("test",sep="")
out_prefix_str <- "test_60c_60r_"

list_r_stack <- create_idrisi_rgf(out_suffix_s=out_suffix_s,file_pat="",
                                  in_dir=file.path(in_dir),ending=FALSE)
test9x9 <- stack(mixedsort(list_r_stack))

temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
p <- levelplot(test9x9, layers=1:25,col.regions=temp.colors,main="test9x9")
print(p)

df9x9 <- as.matrix(as.data.frame(test9x9))
cor_mat <-cor(df9x9)
cor_mat_t <-cor(t(df9x9))

####################

in_dir <- "/Users/benoitparmentier/Dropbox/Data/MEOT_paper/MEOT12272012/3x3dataset"
out_suffix<-"t"
#example for one series:
out_suffix_s <- paste("test",sep="")
out_prefix_str <- "test_60c_60r_"

list_r_stack <- create_idrisi_rgf(out_suffix_s=out_suffix_s,file_pat="",
                                  in_dir=file.path(in_dir),ending=FALSE)
test3x3 <- stack(mixedsort(list_r_stack))

temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
p <- levelplot(test3x3, layers=1:25,col.regions=temp.colors,main="test3x3")
print(p)

df3x3 <- as.matrix(as.data.frame(test3x3))
cor_mat <-cor(df3x3)
cor_mat_t <-cor(t(df3x3))
View(cor_mat)
