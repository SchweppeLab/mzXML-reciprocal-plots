# open mzR mzXML reader
library(mzR)
#find desired mzXML file
MSRun=openMSfile(file.choose())
#create variables for the two known spectra to be comapared using scan numbers (using 1000 and 2000 in this case)
RecipPlot=function(mat,spec1,spec2) {
	y=peaks(mat,spec1)
	z=peaks(mat,spec2)
	#set the ranges for both the intensity and m/z axes
	intrange=c(-1,1)
	mzrange=c(0,max(y[,1],z[,1]))
	#plot background axes and plot
	plot(mzrange,intrange,ylab="Normalized intensity",xlab="m/z",type="n",axes=F)
	axis(2,at=c(-1,0,1),labels=c(1,0,1))
	axis(1,at=mzrange,labels=round(mzrange,digits=0))
	par(new=T)
	#plot top of reciprocal plot
	plot(y[,1],y[,2]/max(y[,2]),type="h",axes=F,ylim=range(intrange),xlim=range(mzrange),ylab="",xlab="")
	#add legend showing MS1 or MS2
	legend(1,paste("MS",header(mat,spec1)$msLevel),xjust=0,box.lty=0)
	if(header(mat,spec1)$precursorMZ>0) {
		legend(0.9,paste("precursor:",round(header(mat,spec1)$precursorMZ,digits=2),"m/z"),xjust=0,box.lty=0)
		}
	par(new=T)
	#plot bottom of reciprocal plot
	plot(z[,1],-z[,2]/max(z[,2]),type="h",axes=F,ylim=range(intrange),xlim=range(mzrange),ylab="",xlab="")
	#add legend showing MS1 or MS2
	legend(-0.8,paste("MS",header(mat,spec2)$msLevel),xjust=0,box.lty=0)
	if(header(mat,spec2)$precursorMZ>0) {
		legend(-0.9,paste("precursor:",round(header(mat,spec2)$precursorMZ,digits=2),"m/z"),xjust=0,box.lty=0)
		}
}

RecipPlot(MSRun,1003,1004)

#close the mzXML file to reduce memory load:
close(MSRun)
