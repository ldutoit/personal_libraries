#Useful_function plotting R

library(ggplot2)

#simple plots_of a correlation with a linear model
plot_correlation<- function(x,y,output_file,cex_points=0.5,meth_cor="pearson"){
	pdf(output_file)
	model<-lm(y~x)
	plot(x,y,pch=19,cex=cex_points,sub=paste(meth_cor,"Rsq",summary(model)$r.squared))
	abline(model)
	dev.off()
	summary(model)
}


#Bin and plots a function, bins y by x according to the number of categories specified for y
plot_binning_to_x<- function(x,y,nb_bins,output_file,cex_points=0.5,meth_cor="pearson"){
	x1<-split(x,cut_number(x,n=nb_bins))
	y1<-split(y,cut_number(x,n=nb_bins))
	means_x<-rep(NA,nb_bins)
	means_y<-rep(NA,nb_bins)
	for (i in 1:nb_bins ){
		means_x[i]<-median(unlist(x1[i]),na.rm=T)
		means_y[i]<-median(unlist(y1[i]),na.rm=T)
	}

plot_correlation(means_x,means_y,output_file= output_file ,cex=cex_points,meth_cor=meth_cor)
}



plot_binning_to_y<- function(x,y,nb_bins,output_file,cex_points=0.5,meth_cor="pearson"){
	x1<-split(x,cut_number(y,n=nb_bins))
	y1<-split(y,cut_number(y,n=nb_bins))
	means_x<-rep(NA,nb_bins)
	means_y<-rep(NA,nb_bins)
	for (i in 1:nb_bins ){
		means_x[i]<-median(unlist(x1[i]),na.rm=T)
		means_y[i]<-median(unlist(y1[i]),na.rm=T)
	}

plot_correlation(means_x,means_y,output_file= output_file ,cex=cex_points,meth_cor=meth_cor)

}


plot_quadratic<-function(X,Y){
	#quadratic relation
	model<- lm(Y ~X + I(X^2))
	newX<- seq( min(X), max(X), by=0.00001)

	plot(X,Y)
	lines(xv, predict(model,list(pred=newx)))
}