setwd(dirname(rstudioapi::getSourceEditorContext()$path))#for Rscript inside RStudio


## Settings ------------------------------------------------------------------------
source('config/settings.R')


## Libraries ------------------------------------------------------------------------
library(raster)
library(ggplot2)
library(RColorBrewer)
library(r2stl)
library(reshape2)
library(grid)
library(gridExtra)
library(gtable)
library(rgdal)
library(DescTools)
library(MODIS)
library(cowplot)
library(beepr) 
MODISoptions(MODISserverOrder = c("LPDAAC","LAADS")) #run lpdaacLogin(server = "LPDAAC") first, saves credentials in ~/.netrc, https://www.rdocumentation.org/packages/MODIS/versions/1.1.0/topics/lpdaacLogin



makecalcs <-function(index){
  
  code<-index
  #code<-1 #Milano
  
  #City Name
  city<-cities[code,]$city
  
  #extents (in World Mollweide coords):
  ext<-c(cities[code,]$extent_minX,cities[code,]$extent_maxX,cities[code,]$extent_minY, cities[code,]$extent_maxY) #vector (length=4; order= xmin, xmax, ymin, ymax)

  #start,end points of profile line
  point_A<-c(x = cities[code,]$PointA_X, y = cities[code,]$PointA_Y)
  point_B<-c(x = cities[code,]$PointB_X, y = cities[code,]$PointB_Y)

  OUTPUTDIR<-file.path("output", city) # το directory για την αποθήκευση των διάφορων output
  dir.create(OUTPUTDIR, showWarnings = FALSE)
  
  
  #set CRS definitions
  wgs84<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  moll<-"+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  
  
  ## DMSP/OLS ------------------------------------------------------------------------
  # Crop data
  rast_ols <- crop(raster(OLS_FILE), extent(ext), snap='near',datatype ="INT1U")
  rast_ols[rast_ols>satDN]<-satDN #set all values>saturationDN to saturationDN(63)
  
  
  
  ## DMSP/OLS (radiance calibrated)------------------------------------------------------------------------
  rin <- raster(OLS_FILE_CAL)
  projection(rin) <- CRS(wgs84)
  rout<-"ols_cal.tif"
  rast_ols_cal<-projectRaster(
    rin,
    rast_ols,
    crs       = moll, 
    method    = 'ngb', 
    alignOnly = FALSE, 
    over      = FALSE,
    dataType  = 'INT1U',
    overwrite = TRUE
  )
  plot(rast_ols_cal)
  
  
  
  ## VIIRS ------------------------------------------------------------------------
  rin <- raster(VIIRS_FILE)
  projection(rin) <- CRS(wgs84)
  
  rast_viirs<-projectRaster(
    rin,
    rast_ols,
    crs       = moll, 
    method    = 'ngb', 
    alignOnly = FALSE, 
    over      = FALSE,
    dataType  = 'INT1U',
    overwrite = TRUE
  )
  plot(rast_viirs)
  
  
  
  
  
  ## MODIS/NDVI------------------------------------------------------------------------
  extent.wgs84 <- extent(projectRaster(rast_ols, crs=wgs84))
  
  ## get latest product collection
  cll <- getCollection(product = product, forceCheck = TRUE)
  
  modis.ndvi<-MODIS::runGdal( # run MODIS:::checkTools('GDAL') to check for GDAL library
    job = city,
    product = product,
    extent =  extent.wgs84,
    SDSstring="1", # check with getSds('*.hdf')
    collection = cll,
    begin = ndvistartdate,
    end = ndvienddate,
    outDirPath = "NDVI",
    overwrite= TRUE,
    checkIntegrity = TRUE
  )
  
  
  
  ## Calculate annual NDVI mean ------------------------------------------------------------------------
  s<-stack(unlist(modis.ndvi, recursive = TRUE, use.names = F))
  annualNDVI <- mean(s)
  
  scalefactor<-0.0001 # see.Product page
  annualNDVI <- annualNDVI*scalefactor
  plot(annualNDVI)
  
  
  
  ## Reproject NDVI to  mollweide------------------------------------------------------------------------
  ndvi<-projectRaster(
    annualNDVI,
    rast_ols,
    crs= mollweide,
    method    = 'ngb', 
    alignOnly = FALSE, 
    over      = FALSE
  )
  
  ndvipositive<-ndvi
  ndvipositive[ndvipositive<0]<-NA #valid range for NDVI [0-1],keep only positive values of NDVI, set negative values to NA
  
  # Data Normalization
  norm_rast_ols<-rast_ols/satDN
  norm_rast_ols_cal<- rast_ols_cal/satDN
  norm_rast_viirs<- rast_viirs/satDN
  
  
  
  
  
  ## VANUI ------------------------------------------------------------------------
  
  #' Calculates VANUI
  #'
  #' @param ndvi normalized ndvi raster without negative values 
  #' @param ols  normalized ols raster
  #'
  #' @return VANUI raster
  #' @export
  #'
  #' @examples
  fvanui<-function(ndvi,ols){
    if (length(ndvi[ndvi<0]>0)){
      stop("NDVI contains negative values" )
    }
    
    if (minValue(ols)<0 || maxValue(ols) > 1){
      stop("NTL raster is not [0-1] normalized" )
    }
    
    if (minValue(ndvi)<0 || maxValue(ndvi) > 1){
      stop("NDVI raster is not [0-1] normalized" )
    }
    
    vanui<-(1-ndvi)*(ols)
  }
  
  #rast_vanui <- (1-ndvi)*(rast_ols/satDN)
  rast_vanui <- fvanui(ndvipositive, norm_rast_ols)
  plot(rast_vanui)
  
  
  
  ## Seawinds------------------------------------------------------------------------
  SeaWinds<-raster(SEAWINDS_FILE)
  
  SeaWinds <- projectRaster(
    from = SeaWinds,
    to = rast_ols,
    crs       = moll,
    method    = 'ngb',
    alignOnly = FALSE,
    over      = FALSE
  )
  
  SeaWinds <- SeaWinds/maxValue(SeaWinds)#normalize Seawinds
  
  plot(SeaWinds)
  pdf(file.path(OUTPUTDIR, paste("Seawinds.hist.", city,".pdf",sep = "")))
  hist(SeaWinds, main="SeaWinds/maxValue(SeaWinds)")
  dev.off()
  
  
  
  ## SANUI ------------------------------------------------------------------------
  
  sanui <- (1-ndvipositive)*SeaWinds*(rast_ols/satDN)
  norm_sanui<- sanui
  
  
  ## LST -----------------------------------------------------------------------------------------
  #settings for LST
  product <- "MOD11A2" #Land Surface Temperature and Emissivity 8-Day L3 Global 1km, https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod11a2
  SDSstring <- "100010000000"
  outDirPath <- "LST"
  scale.factor =0.02 
  
  cll <- getCollection(product = product, forceCheck = TRUE)
  
  lst.modis<-MODIS::runGdal(
    job = city,
    product = product,
    extent =  rast_ols,
    SDSstring = SDSstring,
    collection = cll,
    begin = ndvistartdate,
    end = ndvienddate,
    outDirPath = outDirPath,
    overwrite= TRUE,
    checkIntegrity = TRUE,
    wait = 20
  )
  
  #lst.modis$MOD11A2.006<-lst.modis$MOD11A2.006[names(ndvi.modis$MOD13A2.006)] #keep only 16-days, subset with ndvi.modis dates.
  
  lst.rasters<- lapply(lst.modis[[1]], sapply, function(x) raster(x)*scale.factor) #generate rasters (scaled by factor)
  lst.days<- lapply(lst.rasters,function(x) mean(stack(x),na.rm=T))#mean by day
  #lst.days[which(names(lst.days) %in% c("2013-12-19"))] <- NULL #sos: remove some problematic dates
  
  years<-split(x=lst.days, f=format(as.Date(as.character(names(lst.days)), format = "%Y-%m-%d"),"%Y"))#split by year
  months<-lapply(years, function (x) split(x, f=format(as.Date(as.character(names(x)), format = "%Y-%m-%d"),"%m")))
  
  lst.months<-lapply(months, sapply,function(x) mean(stack(unlist(x)),na.rm=T))
  lst.years.mean<-lapply(years, function(x) mean(stack(x),na.rm=T))#calculate annual mean temp.
  
  
  
  
  
  
  
  ## TVANUI -------------------------------------------------------------
  # data normalization 
  lst.norm<-(lst.years.mean$`2013`-minValue(lst.years.mean$`2013`))/(maxValue(lst.years.mean$`2013`)-minValue(lst.years.mean$`2013`)) #norm lst
  
  #' Calculates TVANUI
  #' [Zhang2018] https://www.sciencedirect.com/science/article/pii/S0924271617303611
  #' @param ndvi normalized ndvi raster without negative values 
  #' @param ols  normalized ols raster 
  #' @param lst  normalized lst raster 
  #'
  #' @return tvanui raster
  #' @export
  #'
  #' @examples
  ftvanui<-function(ndvi,ols, lst){
    if (length(ndvi[ndvi<0]>0)){
      stop("NDVI contains negative values..." )
    }
    
    if (minValue(ols)<0 || maxValue(ols) > 1){
      stop("NTL raster is not [0-1] normalized" )
    }
    
    if (minValue(ndvi)<0 || maxValue(ndvi) > 1){
      stop("NDVI raster is not [0-1] normalized" )
    }
    
    if (minValue(lst)<0 || maxValue(lst) > 1){
      stop("NDVI raster is not [0-1] normalized" )
    }
    
    tvanui<-((atan(lst/ndvi))/(pi/2))*ols #return value
  }
  
  rast_tvanui<-ftvanui(ndvipositive,norm_rast_ols,lst.norm)
  
  
  
  
  
  
  
  
  ## Profile lines over OLS/VIRSS ------------------------------------------------------------------------
  #custom function
  getprofile <- function(raster, lines) {
    profile <- raster::extract(raster, lines)
    profile[[1]]
  }
  
  #open pdf device
  pdf(file.path(OUTPUTDIR, paste("profile.", city,".pdf",sep = "")),width=15,height=15,paper='special')
  par(mfrow=c(2,2))
  
  #profile line
  line<-  Line(rbind(point_A,point_B))
  lines = Lines(list(line), ID=1)
  spatlines = SpatialLines(list(lines))
  sp::proj4string(spatlines)<-moll
  df<-SpatialLinesDataFrame(spatlines, data.frame(ID=1:length(spatlines)))
  #writeOGR(df, dsn="." ,layer="spatlines",driver="ESRI Shapefile", overwrite_layer=TRUE) #for debug
  
  
  #map with VIIRS and profile line
  plot(rast_viirs,col=brewer.pal(n = 9, name = "Oranges"), main="Profile line over VIIRS")
  plot(spatlines, add = TRUE)
  text(point_A["x"],point_A["y"],labels=c("A"), adj=1)
  text(point_B["x"],point_B["y"],labels=c("B"), adj=0)
  
  
  #map with OLSand profile line
  plot(rast_ols,col=brewer.pal(n = 9, name = "Oranges"), main="Profile line over OLS (stable lights)")
  plot(spatlines, add = TRUE)
  text(point_A["x"],point_A["y"],labels=c("A"), adj=1)
  text(point_B["x"],point_B["y"],labels=c("B"), adj=0)
  
  
  #map with OLS (radiance calibrated) and profile line
  plot(rast_ols_cal,col=brewer.pal(n = 9, name = "Oranges"), main="Profile line over OLS (radiance calibrated)")
  plot(spatlines, add = TRUE)
  text(point_A["x"],point_A["y"],labels=c("A"), adj=1)
  text(point_B["x"],point_B["y"],labels=c("B"), adj=0)
  



  # #profiles
  # 
  dn_viirs<-getprofile(norm_rast_viirs, spatlines)
  dn_sanui<-getprofile(norm_sanui, spatlines)
  dn_ols<-getprofile(norm_rast_ols, spatlines)
  dn_ols_cal<-getprofile(norm_rast_ols_cal, spatlines)
  dn_lst<-getprofile(lst.norm, spatlines)
  dn_vanui<-getprofile(rast_vanui, spatlines)
  dn_tvanui<-getprofile(rast_tvanui, spatlines)
  dn_ndvi<-1-getprofile(ndvipositive, spatlines)
  dn_seawinds<-getprofile(SeaWinds, spatlines)

  # 
  # #profiles graph
  # 
  # ylimit<-range(c(dn_sanui,dn_ols_cal, dn_ols,dn_viirs, dn_vanui,dn_ndvi, dn_tvanui))#get the range (of all profiles)
  # 
  # plot(dn_viirs,  type = "l" , col="black", ylim=ylimit,xlab="Distance", ylab="DN")
  # 
  # par(new=TRUE)
  # plot(dn_sanui,  type = "l" , col="green", ylim=ylimit, xlab ='',ylab ='')
  # 
  # par(new=TRUE)
  # plot(dn_ols_cal,  type = "l" , col="purple", ylim=ylimit, xlab ='',ylab ='')
  # 
  # par(new=TRUE)
  # plot(dn_vanui,  type = "l" , col="yellow", ylim=ylimit, xlab ='',ylab ='')
  # 
  # par(new=TRUE)
  # plot(dn_tvanui,  type = "l" , col="brown", ylim=ylimit, xlab ='',ylab ='')
  # 
  # par(new=TRUE)
  # plot(dn_ndvi,  type = "l" , col="blue", ylim=ylimit, xlab ='',ylab ='')
  # 
  # par(new=TRUE)
  # plot(dn_ols, type = "l", col="red", ylim=ylimit,  xlab ='',ylab ='',  xaxt = 'n',  yaxt = 'n', main ="Profile Graph from A to B")
  # 
  # abline(h=1, col = "gray60", lty=2)
  # 
  # 
  # legend(x=0,y=ylimit[2],c("OLS (stable lights)","VIIRS", 
  #                          "SANUI","OLSCAL","VANUI","TVANUI","1-NDVI",
  #                          "Saturation"),
  #        lty=c(rep(1,7),2), 
  #        lwd=c(rep(1,5)),
  #        
  #        col=c("red","black","green","purple","yellow","brown","blue",
  #              "gray60") 
  # )
  # 

  dev.off()
  

  
  
  ## Profile lines with ggplot ------------------------------------------------------------------------
  x= seq(1,length(dn_viirs))
  df <- data.frame(x,dn_ols,dn_seawinds,dn_ndvi, dn_tvanui)
  colnames(df)<- c("Distance", "OLS (stable lights)", "SeaWinds", "1-NDVI", "TVANUI") 
  
  dfmelt <- melt(df, id.vars = "Distance")
  cutoff <- data.frame(yintercept=1, cutoff=factor(1))
  p1 <- ggplot(dfmelt, aes(x = Distance, y =value), size=7) +
    #theme_bw() +
    ggtitle("Profile Graph from A to B")+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_line(aes(linetype = variable))+
    scale_linetype_manual(values=c("twodash", "dotted", "dashed","12345678"))+
    theme(legend.title=element_blank())+ 
    xlab("Distance") +
    ylab("DN")+
    theme(legend.key = element_rect(fill = "white"))+
    guides(colour = guide_legend(reverse=T))+
    annotate("text", x=5,y=1.05, label = "OLS saturation")+
    geom_hline(aes(yintercept= 1),linetype="longdash",colour= 'black')
  
  
  #VANUI, SANUI,TVANUI, Mean LST, VIIRS
  
  df2 <- data.frame(x,dn_sanui,dn_vanui,dn_tvanui, dn_lst, dn_viirs)
  colnames(df2)<- c("Distance", "SANUI", "VANUI", "TVANUI", "Mean LST", "VIIRS") 
  dfmelt <- melt(df2, id.vars = "Distance")
  
  p2 <- ggplot(dfmelt, aes(x = Distance, y =value, linetype = variable), size=7) +
    #theme_bw() +
    #ggtitle("Profile Graph from A to B")+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_line()+
    #scale_linetype_manual(values=c("solid", "dotted", "dashed","12345678", "F1"))+
    theme(legend.title=element_blank())+ 
    xlab("Distance") +
    ylab("DN")+
    theme(legend.key = element_rect(fill = "white"))+
    guides(colour = guide_legend(reverse=T))+
    annotate("text", x=5,y=1.08, label = "OLS saturation")+
    geom_hline(aes(yintercept= 1),linetype="longdash",colour= 'black')
  
  vertplot<-plot_grid(p1, p2,  nrow = 2, align = "v")
  save_plot(file.path(OUTPUTDIR, paste("profile.ggplot.", city,".pdf",sep = "")),
            ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            
            dpi = 600,
            
            vertplot,base_aspect_ratio = 1.3)  #ggsave
  
  
  
  
  ## Pearson correlation coefficient------------------------------------------------------------------------
  
  #---- OLS original  --- #
  # ~ VIIRS 
  cor11<-cor(getValues(rast_ols), getValues(rast_viirs),use ="complete") #use = "pairwise.complete.obs"
  # ~ radiance calibrated 
  cor12<-cor(getValues(rast_ols), getValues(rast_ols_cal),use ="complete") #use = "pairwise.complete.obs"
  
  
  
  #---- OLS desat  --- #
  # ~ VIIRS 
  cor21<-cor(getValues(sanui), getValues(rast_viirs) ,use ="complete")
  # ~ radiance calibrated 
  cor22<-cor(getValues(sanui), getValues(rast_ols_cal) ,use ="complete")
  
  
  
  #---- VANUI  --- # 
  #~ VIIRS 
  cor31<-cor(getValues(rast_vanui), getValues(rast_viirs), use ="complete") #use = "pairwise.complete.obs"
  # ~ radiance calibrated 
  cor32<-cor(getValues(rast_vanui), getValues(rast_ols_cal),use ="complete") 
  
  
  
  #---- SEAWINDS --- # 
  #~ VIIRS 
  cor41<-cor(getValues(SeaWinds), getValues(rast_viirs), use ="complete") 
  # ~ radiance calibrated 
  cor42<-cor(getValues(SeaWinds), getValues(rast_ols_cal),use ="complete")
  
  
  #----1-NDVI --- # 
  #~ VIIRS 
  cor51<-cor(1-getValues(ndvi), getValues(rast_viirs), use ="complete") 
  # ~ radiance calibrated 
  cor52<-cor(1-getValues(ndvi), getValues(rast_ols_cal), use ="complete") 
  
  
  #----1- NDVI[0-1] --- # correlation only with positive NDVI values
  #~ VIIRS 
  cor61<-cor(1-ndvi[ndvi>=0&ndvi<=1], rast_viirs[ndvi>=0&ndvi<=1], use ="complete")
  # ~ radiance calibrated 
  cor62<-cor(1-ndvi[ndvi>=0&ndvi<=1], rast_ols_cal[ndvi>=0&ndvi<=1], use ="complete")
  
  
  
  #----TVANUI --- 
  #~ VIIRS 
  cor71<-cor(getValues(rast_tvanui), getValues(rast_viirs), use ="complete") 
  # ~ radiance calibrated 
  cor72<-cor(getValues(rast_tvanui), getValues(rast_ols_cal),use ="complete")
  
  
  
  corVIIRS<-round(c(
    cor11,
    cor21,
    cor31,
    cor41,
    cor51,
    cor61,
    cor71
  ),2)
  
  corOLSCAL<-round(c(
    cor12,
    cor22,
    cor32,
    cor42,
    cor52,
    cor62,
    cor72
  ),2)
  
  df <- data.frame(corVIIRS = corVIIRS,corOLSCAL=corOLSCAL )
  row.names(df) <- c("OLS (stable lights)", 
                     "SANUI", 
                     "VANUI",
                     "SeaWinds",
                     "1-NDVI",
                     "1-NDVI[0-1]",
                     "TVANUI"
                     
  )
  colnames(df) <- c("VIIRS", 
                    "OLS (radiance calibrated)"
  )
  #save to pdf
  pdf(file.path(OUTPUTDIR, paste("corr.table.", city,".pdf",sep = "")), height=11, width=8.5)
  
  
  t1 <- tableGrob(round(df,3))
  
  grid.draw(t1)
  
  grid.text(paste("Pearson Correlation with VIIRS & OLS (radiance callibrated)-",city,sep=""), x=0.4, y=0.9, 
            gp=gpar(fontsize=14))
  
  dev.off()
  write.csv(df,file.path(OUTPUTDIR, paste("corr.table.", city,".csv",sep = "")),row.names=TRUE)
  

    ## Perspective plots------------------------------------------------------------------------
  #custom function
  rast2persp<- function(raster,z_exagg,...){
    r_mat = as.matrix(raster)
    z <- r_mat*z_exagg
    x <- 1:nrow(z)
    y <- 1:ncol(z)
    
    persp(x, y, z,col = "white", scale=F, main= deparse(substitute(raster)), ...)
  }
  
  
  #common persp params
  THETA=90;PHI = 40;EXPAND = 0.5; Z_EXAGG=0.3
  
  
  ## Εξαγωγή σε pdf persp με τον άξονα τομής ------------------------------------------------------------------------
  genPersp<- function(raster,name, z_exagg,zline,...){
    mraster = as.matrix(raster)
    z <- mraster*z_exagg
    x <- 1:nrow(z)
    y <- 1:ncol(z)

    pdf(file.path(OUTPUTDIR, paste("persp." ,name,".pdf",sep = "")),width=8.88,height=8.88, paper='special')
    
    pmat<-persp(x, y, z,theta = 90, phi = 20, expand = 0.5, col = "white", scale=F, main=name)
    
    cell1<-cellFromXY(raster, point_A)
    row1<-rowFromCell(raster, cell1)
    col1<-colFromCell(raster, cell1)
    
    cell2<-cellFromXY(raster, point_B)
    row2<-rowFromCell(raster, cell2)
    col2<-colFromCell(raster, cell2)
    
    
    lines(trans3d(x=c(row1,row2 ), y=c(col1,col2) , z=zline, pmat= pmat), col = "black", lwd=2)
    points(trans3d(x=c(row1,row2 ), y=c(col1,col2) , z=zline, pmat= pmat), col = "black", pch = 16)
    text(trans3d(x=c(row1,row2 ), y=c(col1,col2) , z=zline, pmat= pmat),labels=c("A", "B"),adj=-1, col = "black")
    dev.off()
  }
  
  z_exagg<-5
  
  genPersp(aggregate(norm_rast_ols, fact=2),"OLS (stable lights)",z_exagg, 10)
  genPersp(aggregate(norm_rast_ols_cal, fact=2),"OLS (radiance callibrated)",z_exagg,20)
  genPersp(aggregate(norm_rast_viirs, fact=2),"VIIRS",z_exagg*3,20)
  genPersp(aggregate(SeaWinds, fact=2),"SeaWinds",z_exagg*5,20)
  genPersp(aggregate(lst.norm, fact=2),"LST (normalized)",z_exagg*3,20)
  genPersp(aggregate(rast_tvanui, fact=2),"TVANUI",z_exagg,10)
  genPersp(aggregate(rast_vanui, fact=2),"VANUI",z_exagg*3,20)
  genPersp(aggregate(norm_sanui, fact=2),"SANUI",z_exagg*3,20)
  genPersp(aggregate(1-ndvi, fact=2),"1-NDVI",z_exagg*3,10)
  
  
  ## SoL/SoNDVI ------------------------------------------------------------------------
  #SoL (OLS stable lights)
  sol<-cellStats(rast_ols, stat='sum', na.rm=TRUE)
  sol_perc<-sol/ncell(rast_ols)/satDN # as % of urban core satDN=saturation DN
  sol_perc<-round(sol_perc,digits=2)
  
  # Sum of NDVI[0-1]
  
  #rsondvi<-ndvi
  sondvi<-cellStats(ndvi, stat='sum', na.rm=TRUE)
  sondvi_perc<-sondvi/ncell(ndvi) 
  sondvi_perc<-round(sondvi_perc,digits=2)
  
  df = data.frame("Sol" = c(sol,ncell(rast_ols),sol_perc), "SoNDVI" = c(sondvi,ncell(ndvi),sondvi_perc))
  rownames(df)<-c("Sum","ncells", "Percentage(extra:SoL%/63)" )
  write.csv(df,file.path(OUTPUTDIR, paste("SoL.SoNDVI." ,city,".csv",sep = "")),row.names=TRUE)
  
  
  print (paste0("City: ", city,". Code:", code, ". Completed!"))
  beep(3)
}

 for (i in 10:10){#1:nrow(cities)){#
   makecalcs(i)
}
