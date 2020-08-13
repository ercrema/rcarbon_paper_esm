nsim = 1000
ncores = 12
run = TRUE 
if(run)
{
print(paste0("Running and Saving Results, using ",nsim," simulations over ",ncores," cores"))
}

if(!run)
{
print(paste0("Plotting Results, using ",nsim," simulations over ",ncores, " cores"))
}


library(rcarbon)
library(maptools)


if (run) #run data
{
  set.seed(1)
  ### Setup ####
print("setup")  
  # IntCal 20 Calibration Curve from Reimer et al 2020
  intcal20 <- read.csv('http://intcal.org/curves/intcal20.14c', encoding="UTF-8",skip=11,header=F)
  colnames(intcal20) <- c("BP","CRA","Error","D14C","Sigma")
  
  # ShCal20 Calibration Curve from Hogg et al 2020
  shcal20 <- read.csv('http://intcal.org/curves/shcal20.14c', encoding="UTF-8",skip=11,header=F)
  colnames(shcal20) <- c("BP","CRA","Error","D14C","Sigma")
  
  # Eastern Meditteranean 14C Dataset from Roberts et al 2018
  data(emedyd) # this is available within the rcarbon package
  emedyd <- subset(emedyd, CRA<=17000&CRA>=8000) 
  southLevant <- subset(emedyd,Region=='1') #
  
  # Sahara 14C Dataset from Manning and Timpson 2014:
  sahara <- read.csv('./data/Manning_and_Timpson2014.csv',header=TRUE)
  sahara <- subset(sahara,Yrs.BP<=17000&Yrs.BP>=8000)
  
  # Brazilian 14C Dataset from Bueno et al 2013:
  brazil <- read.csv('./data/Bueno_etal_2013.csv',header=TRUE)
  brazil <- subset(brazil,OccCRA<=17000&OccCRA>=8000)
  
  # EUROEVOL 14C Datasets from Manning et al 2016:
  data(euroevol) # this is available within the rcarbon package
  data(ewdates) # subset of dates in England and Wales
  data(ewowin) # window of analysis for the English and Welsh subset
  
print("Fig 1 Analysis")  
  ## Figure 1 ####
  grimes <- euroevol[grepl("S2072", euroevol$SiteID),]
  grimesc <- calibrate(x=grimes$C14Age, errors=grimes$C14SD, ids=grimes$C14ID, normalised=FALSE,verbose = FALSE)
  workingstartBP <- 7500
  workingendBP <- 2500
  gspd0 <- spd(x=grimesc, timeRange=c(workingstartBP,workingendBP), datenormalised=TRUE,verbose = FALSE)
  threedates <- grimes[grimes$LabCode %in% c("BM-1022","BM-1033","BM-1066"),]
  threedatesc <- calibrate(x=threedates$C14Age, errors=threedates$C14SD, ids=threedates$C14ID, normalised=FALSE,verbose=FALSE)
  
  ## Thinned version
  n <- 10
  inds <- thinDates(ages=grimes$C14Age,errors=grimes$C14SD,grimes$SiteID,size=n,method="random", seed=99)
  grimesthin <- grimes[inds,]
  grimesthinc <- calibrate(x=grimesthin$C14Age, errors=grimesthin$C14SD, ids=grimesthin$C14ID, normalised=FALSE,verbose = FALSE)
  gspdthin <- spd(x=grimesthinc, timeRange=c(workingstartBP,workingendBP), datenormalised=TRUE,verbose = FALSE)
  
  ## Binned version
  bins50 <- binPrep(sites=grimes$SiteID, ages=grimes$C14Age, h=50)
  bins100 <- binPrep(sites=grimes$SiteID, ages=grimes$C14Age, h=100)
  bins200 <- binPrep(sites=grimes$SiteID, ages=grimes$C14Age, h=200)
  gspd50 <- spd(x=grimesc, bins=bins50, timeRange=c(workingstartBP,workingendBP), datenormalised=TRUE,verbose = FALSE)
  gspd100<- spd(x=grimesc, bins=bins100, timeRange=c(workingstartBP,workingendBP), datenormalised=TRUE,verbose = FALSE)
  gspd200 <- spd(x=grimesc, bins=bins200, timeRange=c(workingstartBP,workingendBP), datenormalised=TRUE,verbose = FALSE)
  
  
  
  
  
  
print("Fig 2 Analysis")  
  ## Figure 2 ####
  
  
  ## Calibrate Dates
  southlevant.x <- calibrate(southLevant$CRA,southLevant$Error,normalised=F,verbose=F,calCurves='intcal20')
  sahara.x <- calibrate(sahara$Yrs.BP,sahara$STD,normalised=F,verbose=F,calCurves='intcal20')
  brazil.x <- calibrate(brazil$CRA,brazil$Error,normalised=F,verbose=F,calCurves='shcal20')
  
  ## Minimal binning to control for sites with many dates.
  southlevant.bins <- binPrep(sites=southLevant$SiteName, ages=southLevant$CRA, h=50)
  sahara.bins <- binPrep(sites=sahara$Sitename, ages=sahara$Yrs.BP, h=50)
  brazil.bins <- binPrep(sites=brazil$SiteName, ages=brazil$CRA, h=50)
  
  ## Summation 
  southlevant.spd.nnorm <- spd(southlevant.x, bins=southlevant.bins, timeRange=c(16000,9000),verbose=F)
  sahara.spd.nnorm <- spd(sahara.x, bins=sahara.bins, timeRange=c(16000,9000),verbose=F)
  brazil.spd.nnorm <- spd(brazil.x, bins=brazil.bins, timeRange=c(16000,9000),verbose=F)
  southlevant.spd.norm <- spd(southlevant.x, bins=southlevant.bins, timeRange=c(16000,9000),verbose=F,datenormalised=T)
  sahara.spd.norm <- spd(sahara.x, bins=sahara.bins, timeRange=c(16000,9000),verbose=F,datenormalised=T)
  brazil.spd.norm <- spd(brazil.x, bins=brazil.bins, timeRange=c(16000,9000),verbose=F,datenormalised=T)
 


print("Fig 3 Analysis")  
  ## Figure 3 ####
  southlevant.x <- calibrate(southLevant$CRA,southLevant$Error,normalised=TRUE,verbose=FALSE)
  southlevant.x2 <- calibrate(southLevant$CRA,southLevant$Error,normalised=FALSE,verbose=FALSE)
  
  # 50 years bin
  sl.bins <- binPrep(sites = southLevant$SiteName, ages = southLevant$CRA, h = 50)
  
  # Fit and simulate theorethical SPDs using the calsample method:
  test.cal.nn = modelTest(southlevant.x2, bins = sl.bins, errors=southLevant$Error,timeRange=c(16000,9000), nsim=nsim, ncores=ncores, method='calsample', normalised = FALSE,runm=50,verbose=FALSE)
  test.cal.n <- modelTest(southlevant.x, bins = sl.bins, errors=southLevant$Error,timeRange=c(16000,9000), nsim=nsim, ncores=ncores, method='calsample', normalised = TRUE,runm=50,verbose=FALSE)
  
  # Fit and simulate theorethical SPDs using the uncalsample method:
  test.uncal.nn <- modelTest(southlevant.x2, bins = sl.bins, errors=southLevant$Error,timeRange=c(16000,9000), nsim=nsim, ncores=ncores, method='uncalsample', normalised = FALSE,runm=50,verbose=FALSE)
  test.uncal.n <- modelTest(southlevant.x, bins = sl.bins, errors=southLevant$Error,timeRange=c(16000,9000), nsim=nsim, ncores=ncores, method='uncalsample', normalised = TRUE,runm=50,verbose=FALSE)
 


print("Fig 4 Analysis")  
  ## Figure 4 ####
  levant.emedyd <- subset(emedyd,Region%in%c(1,2))
  
  alldatesLevant.emedyd <- calibrate(x=levant.emedyd$CRA, errors=levant.emedyd$Error, calCurves='intcal20', normalised=FALSE,verbose=FALSE)
  binsLevant <- binPrep(sites=levant.emedyd$SiteName, ages=levant.emedyd$CRA, h=50)
  alldates.spd <- spd(alldatesLevant.emedyd,bins=binsLevant,timeRange=c(16000,9000),runm=200,spdnormalised=T,verbose=F)
  
  #execute the permutation test
  resultLevant <- permTest(x=alldatesLevant.emedyd, bins=binsLevant, marks=levant.emedyd$Region, timeRange=c(16000,9000), runm=200,nsim=nsim,spdnormalised=T,verbose=FALSE)
  
 

print("Fig 5 Analysis")  
  ## Figure 5 ####
  fig5.x <- calibrate(x=ewdates$C14Age, errors=ewdates$C14SD, normalised=FALSE,verbose=FALSE)
  ## Create centennial timeslices (also with site binning)
  fig5.bins <- binPrep(sites=ewdates$SiteID, ages=ewdates$C14Age, h=50)
  
 
print("Fig 6 Analysis")  
  ## Figure 6 #### 
  library(maptools)
  fig6.x <- calibrate(x=ewdates$C14Age, errors=ewdates$C14SD, normalised=FALSE,verbose=FALSE)
  fig6.bins <- binPrep(sites=ewdates$SiteID, ages=ewdates$C14Age, h=50)
  
  # Create SpatialPoints
  sites <- unique(data.frame(id=ewdates$SiteID,east=ewdates$Eastings,north=ewdates$Northings))
  rownames(sites) <- sites$id
  sites <- sites[,-1]
  sp::coordinates(sites) <- c("east","north")
  sp::proj4string(sites) <- sp::CRS("+init=epsg:27700")
  
  #Compute distance matrix
  d <- sp::spDists(sites,sites,longlat=FALSE)
  #Compute spatial weights 
  w <- spweights(d/1000,h=40)
  # Define temporal blocks
  breaks <- seq(6300,5900,-200) 
  
  # Spatial Permutation Test
  fig6.res <- sptest(calDates=fig6.x,bins=fig6.bins,timeRange=c(6300,5900),locations=sites,permute="locations",nsim=nsim,breaks=breaks,spatialweights=w,verbose=FALSE) 

  ## Save Image ####
  save.image("figres.RData")
}

if (run==FALSE)
{
 load.image("figres.RData") 
}


print("Making figures...")
### FIGURE 1 ####
tiff(file = "figure1.tiff", width = 4, height = 6, units = "in", res = 300)
## Code Figure 1 ####

# plot
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow=TRUE), widths=c(3,3), heights=c(3,2,2))

# Panel a
par(mar=c(4, 4, 1, 0.5)) 
plot(gspd0, calendar="BCAD", fill.p="chocolate1")
plot(gspd0, runm=50, calendar="BCAD", add=TRUE, fill.p=NA, border="chocolate4")
plot(threedatesc[1],add=TRUE,calendar='BCAD', col='black')
plot(threedatesc[2],add=TRUE,calendar='BCAD', col='black')
plot(threedatesc[3],add=TRUE,calendar='BCAD', col='black')

legend("topleft", bty="n", col=c("chocolate1","chocolate4","black"), legend=c("raw SPD", "smoothed (50-yr running mean)", "3 example dates"), lwd=c(5,1,5), cex=0.8)
text(-5400, 0.02, "a", font=2)

# Panel b
par(mar=c(0.5, 3, 0.5, 0.5)) 

plot(gspd0, xaxt="n", yaxt="n", runm=50, fill.p=NA, border="white", spdnormalised=TRUE)
plot(gspdthin, xaxt="n", yaxt="n", runm=50, fill.p="aquamarine", spdnormalised=TRUE, add=TRUE)
plot(gspd0, xaxt="n", yaxt="n", runm=50, fill.p=NA, border="chocolate4", spdnormalised=TRUE, add=TRUE)
legend("topleft", bty="n", col=c("chocolate4","aquamarine"), legend=c("original SPD", paste("Thinned to ", nrow(grimesthin), " random dates", sep="")), lwd=c(1,5), cex=0.8)
text(7200, 0.00015, "b", font=2)

# Panel c
plot(gspd0, xaxt="n", yaxt="n", runm=50, fill.p=NA, border="white", spdnormalised=TRUE)
plot(gspd50, xaxt="n", yaxt="n", runm=50, fill.p="bisque3", spdnormalised=TRUE, add=TRUE)
plot(gspd0, xaxt="n", yaxt="n", runm=50, fill.p=NA, border="chocolate4", spdnormalised=TRUE, add=TRUE)
legend("topleft", bty="n", col=c("chocolate4","bisque3"), legend=c("original SPD", "50 year bins"), lwd=c(1,5), cex=0.8)
text(7200, 0.00015, "c", font=2)

# Panel d
plot(gspd0, xaxt="n", yaxt="n", runm=50, fill.p=NA, border="white", spdnormalised=TRUE)
plot(gspd100, xaxt="n", yaxt="n", runm=50, fill.p="bisque3", spdnormalised=TRUE, add=TRUE)
plot(gspd0, xaxt="n", yaxt="n", runm=50, fill.p=NA, border="chocolate4", spdnormalised=TRUE, add=TRUE)
legend("topleft", bty="n", col=c("chocolate4","bisque3"), legend=c("original SPD", "100 year bins"), lwd=c(1,5), cex=0.8)
text(7200, 0.00015, "d", font=2)

# Panel e
plot(gspd0, xaxt="n", yaxt="n", runm=50, fill.p=NA, border="white", spdnormalised=TRUE)
plot(gspd200, xaxt="n", yaxt="n", runm=50, fill.p="bisque3", spdnormalised=TRUE, add=TRUE)
plot(gspd0, xaxt="n", yaxt="n", runm=50, fill.p=NA, border="chocolate4", spdnormalised=TRUE, add=TRUE)
legend("topleft", bty="n", col=c("chocolate4","bisque3"), legend=c("original SPD", "200 year bins"), lwd=c(1,5), cex=0.8)
text(7200, 0.00015, "e", font=2)

dev.off()


### FIGURE 2 ####
tiff(file = "figure2.tiff", width = 5.5, height = 3.5, units = "in", res = 300,pointsize=8)
## Code Figure 2 ####
fmat <- matrix(c(rep(1,3),rep(2,3),rep(3:5,each=2)),nrow=6,ncol=2)
layout(fmat,widths=c(0.4,0.8),heights=c(1.2,1,0.8,0.8,1,1.2))

# Panel a
x1 <- calibrate(3523,45,normalised = F,verbose=F)
plot(x1,type='auc',cex.lab=1,cex.axis=1)
title("a",line=-2)

# Panel b
x2 <- calibrate(4274,45,normalised = F,verbose=F)
plot(x2,type='auc',cex.lab=1,cex.axis=1)
title("b",line=-2)


# Panel c
par(yaxs="i")
par(xaxs="i")
par(mar=c(0,6,6,1)) 
plot(southlevant.spd.nnorm,fill.p=NA, xlim=c(14000,9000), ylim=c(0,0.25), xaxt="n",yaxt='n',cex.lab=1,cex.axis=0)
rect(x=c(9509.948,10228.445,11241.710,12690)+100,xright=c(9509.948,10228.445,11241.710,12690)-100,ybottom=rep(-100,4),ytop=rep(100,5),border=NA,col='orange')
plot(southlevant.spd.nnorm,fill.p='bisque3',ylim=c(0,0.25), xaxt="n", add=TRUE)
plot(southlevant.spd.norm,type='simple', xaxt="n",add=TRUE)
lines(intcal20[intcal20$BP<=14000 & intcal20$BP>=9000,"BP"],reScale(intcal20[intcal20$BP<=14000 & intcal20$BP>=9000,"CRA"])*par('usr')[4],col="red",lwd=2)
legend("topleft",legend=c("Unnormalised SPD","Normalised SPD", "IntCal20","ShCal20","Areas with spikes"),col=c('bisque3','black',"red","blue","orange"),lwd=c(5,1,1,1,5))
title("c",line=-2)

# Panel d
par(mar=c(0,6,0,1))
plot(sahara.spd.nnorm,fill.p=NA, xlim=c(14000,9000), ylim=c(0,0.35), xaxt="n",yaxt='n',cex.lab=1)
rect(x=c(9509.948,10228.445,11241.710,12690)+100,xright=c(9509.948,10228.445,11241.710,12690)-100,ybottom=rep(-100,4),ytop=rep(100,5),border=NA,col='orange')
plot(sahara.spd.nnorm,fill.p='bisque3',ylim=c(0,0.35), xaxt="n", add=TRUE)
plot(sahara.spd.norm,type = 'simple', xaxt="n",add=TRUE)
lines(intcal20[intcal20$BP<=14000 & intcal20$BP>=9000,"BP"],reScale(intcal20[intcal20$BP<=14000 & intcal20$BP>=9000,"CRA"])*par('usr')[4],col="red",lwd=2)
title("d",line=-2)
title(ylab="Summed Probability",cex.lab=1.3,line=0.5)

#Panel e
par(mar=c(6,6,0,1))
plot(brazil.spd.nnorm,fill.p=NA, xlim=c(14000,9000), ylim=c(0,0.125), xaxt="n",yaxt='n',cex.lab=1,cex.axis=0)
rect(x=c(9509.948,10228.445,11241.710,12690)+100,xright=c(9509.948,10228.445,11241.710,12690)-100,ybottom=rep(-100,4),ytop=rep(100,5),border=NA,col='orange')
plot(brazil.spd.nnorm,fill.p='bisque3',ylim=c(0,0.125), xaxt="n", add=TRUE)
plot(brazil.spd.norm,type = 'simple', xaxt="n", add=TRUE)
axis(side=1, at=seq(14000,9000,-1000), labels=seq(14000,9000,-1000), las=2, cex.axis=1)
lines(shcal20[shcal20$BP<=14000 & shcal20$BP>=9000,"BP"],reScale(shcal20[shcal20$BP<=14000 & shcal20$BP>=9000,"CRA"])*par('usr')[4],col="blue",lwd=2)
title("e",line=-2)

dev.off()









### FIGURE 3 ####
tiff(file = "figure3.tiff", width = 5, height = 5, units = "in", res = 300,pointsize=10)
## Code Figure 3 ####
par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(test.cal.n,lwd.obs = 1.5)
title('a')
plot(test.cal.nn,lwd.obs = 1.5)
title('b')
plot(test.uncal.n,lwd.obs = 1.5)
title('c')
plot(test.uncal.nn,lwd.obs = 1.5)
title('d')
dev.off()













### FIGURE 4 ####
tiff(file = "figure4.tiff", width = 4, height = 8, units = "in", res = 300)
## Code Figure 4 ####
par(mfrow=c(2,1))
plot(resultLevant,focalm='1',lwd=2,main="Southern Levant")
plot(alldates.spd,add=T,type='simple',lty=2,col='grey27')
plot(resultLevant,focalm='2',lwd=2,main="Northern Levant")
plot(alldates.spd,add=T,type='simple',lty=2,col='grey27')
dev.off()






### FIGURE 5 ####
tiff(file = "figure5.tiff", width = 5.5, height = 2, units = "in", res = 300,pointsize=10)
## Code Figure 5 ####
## Plot an example of all four basic outputs for 6000 calBP
stkde1 <- stkde(x=fig5.x, coords=ewdates[,c("Eastings", "Northings")], win=ewowin, sbw=40000, cellres=2000, focalyears=seq(6500, 5000, -100), tbw=100, bins=fig5.bins, backsight=200, outdir=tempdir(),amount=1,verbose=FALSE)
par(mar=c(0.5, 0.5, 2.5, 2))
plot(stkde1, 6000, type="all")
## Add scale bar
arrows(x0 = 478616, x1= 478616+100000, y0=40187, y1=40187, code=3, length=0.01, angle=90)
text(x=478616+50000,y=187,"100km")
dev.off()



### FIGURE 6 ####
tiff(file = "figure6.tiff", width = 5, height = 3.5, units = "in", res = 300,pointsize = 8.5)
## Code Figure 6 ####
base=as(ewowin,"SpatialPolygons")
sp::proj4string(base) <- sp::CRS("+init=epsg:27700")
xlim = bbox(base)[1,]
ylim = bbox(base)[2,]
xlim[1] = xlim[1]-100000
ylim[1] = ylim[1] + 100000
par(mar=c(0,1,1,1),mfrow=c(1,2))


plot(base,col="antiquewhite3",border="antiquewhite3",xlim=xlim,ylim=ylim,yaxs='i')
plot(fig6.res,index=1,option="raw",add=TRUE,legend=TRUE,legSize=0.8,location="topleft",baseSize=0.9)
plot(base,col="antiquewhite3",border="antiquewhite3",xlim=xlim)
plot(fig6.res,index=1,option="test",add=TRUE,legend=TRUE,legSize=0.8,location="topleft",baseSize=0.9)

# Add scale bar
arrows(x0 = 478616, x1= 478616+100000, y0=40187, y1=40187, code=3, length=0.01, angle=90)
text(x=478616+50000,y=187,"100km")

dev.off()




