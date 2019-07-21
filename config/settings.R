setwd("config")

cities <- read.csv(file="cities.csv", header=TRUE, sep=",",stringsAsFactors = FALSE)

setwd("../")

OLSDIR<-'OLS'
OLS_FILE<-file.path(OLSDIR, 'FIL2013.tif')# DMSP/OLS file, source: dstath
satDN<-63 # saturation DN of DMSP/OLS

OLSCALDIR<-'OLS_CAL'
OLSCAL_URL <- "https://data.ngdc.noaa.gov/instruments/remote-sensing/passive/spectrometers-radiometers/imaging/ols/composites/fg/F16_20100111-20101209_rad_v4.geotiff.tgz"
OLSCALGZ_FILE<-'F16_20100111-20101209_rad_v4.geotiff.tgz'
OLS_FILE_CAL<-'F16_20100111-20101209_rad_v4.avg_vis.tif'

VIIRSDIR<-'VIIRS'
VIIRS_URL<- "https://data.ngdc.noaa.gov/instruments/remote-sensing/passive/spectrometers-radiometers/imaging/viirs/dnb_composites/v10//201401/vcmcfg/SVDNB_npp_20140101-20140131_75N060W_vcmcfg_v10_c201506171538.tgz"
VIIRSGZ_FILE<-'SVDNB_npp_20140101-20140131_75N060W_vcmcfg_v10_c201506171538.tgz'
VIIRS_FILE<-'SVDNB_npp_20140101-20140131_75N060W_vcmcfg_v10_c201506171538.avg_rade9h.tif'


SEAWINDSDIR<-'SEAWINDS'
SEAWINDS_FILE<-'Global_quev_2009_JFM_PR.tif'


#*************** DMSP/OLS Radiance Cal.
#https://www.ngdc.noaa.gov/eog/dmsp/download_radcal.html, file: F16_20100111-20101209_rad_v4
#https://data.ngdc.noaa.gov/instruments/remote-sensing/passive/spectrometers-radiometers/imaging/ols/composites/fg/F16_20100111-20101209_rad_v4.geotiff.tgz

#function to download and unzip files
getUncomprFile <- function(dir, gzfile, URL, file) {
  #check if file already exists else...
  if (!file.exists(file.path(dir ,file))){
    #check if compressed file exists
    if (!file.exists(file.path(dir, gzfile))){
      
      #Download OLSCAL file
      download.file(
        url=URL,
        destfile=file.path(dir ,gzfile),
        method='curl')
    }
    #uncompress file
    untar(file.path(dir, gzfile),files=file, exdir =file.path(dir))
  }
}
getUncomprFile(OLSCALDIR,OLSCALGZ_FILE,OLSCAL_URL,OLS_FILE_CAL)
OLS_FILE_CAL<-file.path(OLSCALDIR,OLS_FILE_CAL)

#*********** VIIRS
#https://www.ngdc.noaa.gov/eog/viirs/download_dnb_composites.html
#01/2014,vcmsl:the radiance values have undergone the stray-light correction procedure : https://data.ngdc.noaa.gov/instruments/remote-sensing/passive/spectrometers-radiometers/imaging/viirs/dnb_composites/v10//201401/vcmcfg/SVDNB_npp_20140101-20140131_75N060W_vcmcfg_v10_c201506171538.tgz
getUncomprFile(VIIRSDIR,VIIRSGZ_FILE,VIIRS_URL,VIIRS_FILE)
VIIRS_FILE<-file.path(VIIRSDIR,VIIRS_FILE)
SEAWINDS_FILE<-file.path(SEAWINDSDIR,SEAWINDS_FILE) # source: dstath

#MODIS
ndvistartdate <-"2013-01-01"
ndvienddate <-"2013-12-31"
product <- "MOD13A3" # check also "MOD13Q1"



