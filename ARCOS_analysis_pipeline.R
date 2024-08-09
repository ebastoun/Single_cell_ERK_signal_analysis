## Code for analysis of collective ERK signal with ARCOS package

## set up enviroment
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #setup WD to file path

## Call libraries
library(ARCOS)
library(data.table)
library(ggplot2)
library(sf)
library(dplyr)

## Call matrix from single cell analysis
folder_path = '/Users/Julio/Desktop/20240801_LM_ERK_nuclei/Ctr_4_thld_001'
file_name = 'ARCOS_matrix.csv'
ARCOS_matrix = read.csv(paste(c(folder_path,file_name),collapse = "/"),header = TRUE)
cal_factor = 0.28 # [um/px] Calibration factor
field_of_view = 1608 # [px] Total length of field of view
time_step = 10 # [min] time step between frames

# define arcos object
ARCOS_matrix = arcosTS(data.table(ARCOS_matrix),colPos = c('x','y'), colFrame = 't', colIDobj = 'id', colMeas = 'm')

# Calculate collective events (ERK activation waves)
dcoll = trackColl(ARCOS_matrix[m>0],eps = 70, epsPrev = NULL,  minClSz = as.integer(4), nPrev = as.integer(1)) # calculate events for active objects "m>0" and for 40x magnification 0.28 m/px
# 92.12

# Plot of nuclei and overly of calculated collective events
num_t_steps = unique(dcoll$t)

colna_data_hull = c('t','collid','Area','N_points')
Data_hull = data.frame(matrix(nrow = 0, ncol = length(colna_data_hull)))
colnames(Data_hull) = colna_data_hull

for (i in num_t_steps) { # time point
  num_collid = unique(dcoll[ t == i,4])
  num_collid = as.vector(unlist(num_collid))
  data_plot = ARCOS_matrix[t == i,]
  p = ggplot() + geom_point(data=data_plot,aes(x,y,colour = m),size=0.75) + scale_x_continuous()  + scale_y_continuous() + scale_y_reverse() +
    coord_fixed(ratio=1,xlim = c(0,1608), ylim = rev(c(0,1608)), expand = FALSE) + scale_color_gradient(low="blue", high="red") +
    theme_gray(base_size = 10)
  
  for (j in num_collid ) { # collective events
    # storage: properties of indiviual clusters (collective events)
    dcoll_aux = data.frame(dcoll)
    dcoll_aux = dcoll %>% filter(t == i,collid == j)
    
    # Calculation of concave hull (boundary of cluster) 
    Con_hull = st_concave_hull(st_multipoint(as.matrix(dcoll_aux[,c(5,6)]),dim = "xy"),ratio = 0.5)
    
    if (length(Con_hull[[1]]) == 1 ) {
      next
    }
    
    if (nrow(Con_hull[[1]]) < 4) {
      next
    }
    Con_hull_cen = st_centroid(Con_hull) 
    Area_hull = st_area(Con_hull) # px^2
    
    # Storage cluster characteristics 
    Data_hull = bind_rows(Data_hull,data.frame('t'=i,'collid'=j,'Area'=Area_hull,'N_points'=nrow(dcoll_aux)))
    
    # # Order of concave hull points
    hull_plot = data.frame("x1" = Con_hull[[1]][,1], "y1" = Con_hull[[1]][,2],
                           "x2" = Con_hull[[1]][c(seq(2,nrow(Con_hull[[1]])),1),1],
                           "y2" = Con_hull[[1]][c(seq(2,nrow(Con_hull[[1]])),1),2])
    
    p = p + geom_segment(data = hull_plot,aes(x=x1,y=y1,xend=x2,yend=y2),linewidth = 0.3,color = "magenta")
    
    
  }
  
  # Save image of detected clusters
  ggsave(paste('Frame_',as.character(i),'.png',sep = ""), path = paste(c(folder_path,'Frames',''),collapse = "/"), plot = p,
         width = 10, height = 10, units = "cm", dpi = 300)
  
}

## Save analysis data
write.csv(Data_hull, file = paste(folder_path,"arcos_analysis.csv",sep="",collapse ="/"),row.names = FALSE)

################################################################################

## Load analysis of data 
Analysis_file = paste(c(folder_path,"Frames","arcos_analysis.csv"),sep = "",collapse = "/")
Data_hull = read.csv(Analysis_file)
if (!is.null(Data_hull$X)) {
  Data_hull = Data_hull[,seq(2,5)]
}

## Filter clusters that appear only in one frame
counts_frame = tabulate(Data_hull$collid)
counts_remove = which(counts_frame == 0)
if (length(counts_remove)>0) {
  counts_frame = counts_frame[-counts_remove]
}
cluster_id  = unique(Data_hull$collid)
cluster_remove = which(counts_frame == 1)
if (length(cluster_remove)>0) {
  cluster_keep = cluster_id[-cluster_remove]
} else {
  cluster_keep = cluster_id
}
Data_hull = Data_hull %>% filter(Data_hull$collid %in% cluster_keep)

# Number of events per hour
Total_events = NULL
for (i in unique(Data_hull[,1])) {
  Total_events = rbind(Total_events, nrow(Data_hull %>% filter(t == i)))
}
Total_time = length(unique(Data_hull$t))*(time_step/60) # [h]
Events_per_hour = sum(Total_events)/Total_time/(field_of_view*cal_factor/1000)^2 # [Events/mm^2-h]

# Area per event
Event_total_area = NULL
for (i in unique(Data_hull$collid)) {
  Event_total_area = append(Event_total_area,Data_hull[Data_hull$collid==i,3]) # px^2
}
Event_total_area = Event_total_area*(cal_factor^2) #um2
Median_total_area = median(Event_total_area)

# Results display
print(paste(c('Median number of events:',round(Events_per_hour,digits = 3),"per mm^2-h"), collapse = " "))
print(paste(c('Median event area:',round(Median_total_area,digits = 3),"mm^2"), collapse = " "))
