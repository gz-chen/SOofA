regRange=function(from,to){
	if(from>to){
		return(numeric(0));
	} else {
		return(from:to);
	}
}

fun_plot_test_pch=function(){
	par(mfrow=c(1,1),mar=c(3,3,0.5,0.5),mgp=c(3,1,0));
	plot(x=NULL,xlim=c(0,7),ylim=c(0,8),
		asp=1,main="",ylab="",xlab="",cex.axis=2,cex.lab=2);
	for(i0 in 0:8){
		for(i1 in 0:7){
			points(x=i1,y=i0,pch=i0*8+i1,cex=5);
		}
	}
}
fun_plot_strat_pairone=function(mat_des,vec_stratsize,fileout,display=1){
	n_comp=dim(mat_des)[2];
	list_cex=list(cex_point=5,lwd_point=2.5,lwd_sep=2.5);
	n_seppos=length(vec_stratsize)-1;
	seppos=regRange(1,length(vec_stratsize));
	seppos[length(vec_stratsize)]=n_comp-0.5;
	for(i_sep in 1:n_seppos){
		seppos[length(vec_stratsize)-i_sep]=
			seppos[length(vec_stratsize)-i_sep+1]-vec_stratsize[length(vec_stratsize)-i_sep+1];
	}
	seppos=seppos[1:n_seppos];
	if(display==0){
		postscript(fileout,width=10,height=10,pointsize=16,paper="special",horizontal=FALSE);
	}
	par(mfrow=c(1,1),mar=c(5,5,0.5,0.5),mgp=c(3,1,0));
	plot(x=NULL,xlim=c(0,n_comp-1),ylim=c(0,n_comp-1),
		asp=1,main="",ylab="Component 1",xlab="Component 0",cex.axis=1,cex.lab=1);
	abline(v=seppos,lty=2,lwd=list_cex$lwd_sep,col=rgb(1,0,0));
	abline(h=seppos,lty=2,lwd=list_cex$lwd_sep,col=rgb(1,0,0));
	for(i0 in 0:(n_comp-1)){
		for(i1 in 0:(n_comp-1)){
			if(i0!=i1){
				points(x=i0,y=i1,
					pch=1,cex=list_cex$cex_point,lwd=list_cex$lwd_point,col=rgb(0.5,0.5,0.5));
			}
		}
	}
	points(x=mat_des[,c(1,2)],pch=21,cex=list_cex$cex_point,lwd=list_cex$lwd_point,
		col=rgb(0.5,0.5,0.5),bg=rgb(0,0,1));
	if(display==0){
		dev.off();
	}
}
fun_plot_strat_pairall=function(mat_des,vec_stratsize,fileout,display=1){
	n_comp=dim(mat_des)[2];
	list_cex=list(cex_point=2.5,lwd_point=1.5,lwd_sep=3,cex_text_diag=10);
	for(i in regRange(1,length(list_cex))){
		list_cex[[i]]=list_cex[[i]]/sqrt(n_comp);
	}
	n_seppos=length(vec_stratsize)-1;
	seppos=regRange(1,length(vec_stratsize));
	seppos[length(vec_stratsize)]=n_comp-0.5;
	for(i_sep in 1:n_seppos){
		seppos[length(vec_stratsize)-i_sep]=
			seppos[length(vec_stratsize)-i_sep+1]-vec_stratsize[length(vec_stratsize)-i_sep+1];
	}
	seppos=seppos[1:n_seppos];
	if(display==0){
		postscript(fileout,width=10,height=10,pointsize=16,paper="special",horizontal=FALSE);
	}
	par(mfrow=c(n_comp,n_comp),mgp=c(0,0,0),mar=c(0.25,0.25,0.25,0.25));
	for(i_comp_0 in regRange(1,n_comp)){
		for(i_comp_1 in regRange(1,i_comp_0-1)){
			plot(x=NULL,xlim=c(0,n_comp-1),ylim=c(0,n_comp-1),
				asp=1,main="",ylab="",xlab="",xaxt="n",yaxt="n");
			abline(v=seppos,lty=2,lwd=list_cex$lwd_sep,col=rgb(1,0,0));
			abline(h=seppos,lty=2,lwd=list_cex$lwd_sep,col=rgb(1,0,0));
			for(i0 in regRange(0,n_comp-1)){
				for(i1 in regRange(0,n_comp-1)){
					if(i0!=i1){
						points(x=i0,y=i1,
							pch=1,cex=list_cex$cex_point,lwd=list_cex$lwd_point,
							col=rgb(0.5,0.5,0.5));
					}
				}
			}
			points(x=mat_des[,c(i_comp_1,i_comp_0)],
				pch=21,cex=list_cex$cex_point,lwd=list_cex$lwd_point,
				col=rgb(0.5,0.5,0.5),bg=rgb(0,0,1));
		}
		{
			plot(x=NULL,xlim=c(0,n_comp-1),ylim=c(0,n_comp-1),
				asp=1,main="",ylab="",xlab="",xaxt="n",yaxt="n");
			text(x=(n_comp-1)/2,y=(n_comp-1)/2,labels=sprintf("%d",i_comp_0-1),
				cex=list_cex$cex_text_diag);
		}
		for(i_comp_1 in regRange(i_comp_0+1,n_comp)){
			plot.new();
		}
	}
	if(display==0){
		dev.off();
	}
}
fun_plot_strat_pairall_compare=function(mat_des_0,mat_des_1,vec_stratsize,fileout,display=1){
	n_comp=dim(mat_des_0)[2];
	list_cex=list(cex_point=2.5,lwd_point=1.5,lwd_sep=3,cex_text_diag=10);
	for(i in regRange(1,length(list_cex))){
		list_cex[[i]]=list_cex[[i]]/sqrt(n_comp);
	}
	n_seppos=length(vec_stratsize)-1;
	seppos=regRange(1,length(vec_stratsize));
	seppos[length(vec_stratsize)]=n_comp-0.5;
	for(i_sep in 1:n_seppos){
		seppos[length(vec_stratsize)-i_sep]=
			seppos[length(vec_stratsize)-i_sep+1]-vec_stratsize[length(vec_stratsize)-i_sep+1];
	}
	seppos=seppos[1:n_seppos];
	if(display==0){
		postscript(fileout,width=10,height=10,pointsize=16,paper="special",horizontal=FALSE);
	}
	par(mfrow=c(n_comp,n_comp),mgp=c(0,0,0),mar=c(0.25,0.25,0.25,0.25));
	for(i_comp_0 in regRange(1,n_comp)){
		for(i_comp_1 in regRange(1,i_comp_0-1)){
			plot(x=NULL,xlim=c(0,n_comp-1),ylim=c(0,n_comp-1),
				asp=1,main="",ylab="",xlab="",xaxt="n",yaxt="n");
			abline(v=seppos,lty=2,lwd=list_cex$lwd_sep,col=rgb(1,0,0));
			abline(h=seppos,lty=2,lwd=list_cex$lwd_sep,col=rgb(1,0,0));
			for(i0 in regRange(0,n_comp-1)){
				for(i1 in regRange(0,n_comp-1)){
					if(i0!=i1){
						points(x=i0,y=i1,
							pch=1,cex=list_cex$cex_point,lwd=list_cex$lwd_point,
							col=rgb(0.5,0.5,0.5));
					}
				}
			}
			points(x=mat_des_0[,c(i_comp_1,i_comp_0)],
				pch=21,cex=list_cex$cex_point,lwd=list_cex$lwd_point,
				col=rgb(0.5,0.5,0.5),bg=rgb(0,0,1));
		}
		{
			plot(x=NULL,xlim=c(0,n_comp-1),ylim=c(0,n_comp-1),
				asp=1,main="",ylab="",xlab="",xaxt="n",yaxt="n");
			text(x=(n_comp-1)/2,y=(n_comp-1)/2,labels=sprintf("%d",i_comp_0-1),
				cex=list_cex$cex_text_diag);
		}
		for(i_comp_1 in regRange(i_comp_0+1,n_comp)){
			plot(x=NULL,xlim=c(0,n_comp-1),ylim=c(0,n_comp-1),
				asp=1,main="",ylab="",xlab="",xaxt="n",yaxt="n");
			abline(v=seppos,lty=2,lwd=list_cex$lwd_sep,col=rgb(1,0,0));
			abline(h=seppos,lty=2,lwd=list_cex$lwd_sep,col=rgb(1,0,0));
			for(i0 in regRange(0,n_comp-1)){
				for(i1 in regRange(0,n_comp-1)){
					if(i0!=i1){
						points(x=i0,y=i1,
							pch=1,cex=list_cex$cex_point,lwd=list_cex$lwd_point,
							col=rgb(0.5,0.5,0.5));
					}
				}
			}
			points(x=mat_des_1[,c(i_comp_1,i_comp_0)],
				pch=21,cex=list_cex$cex_point,lwd=list_cex$lwd_point,
				col=rgb(0.5,0.5,0.5),bg=rgb(0,0,1));
		}
	}
	if(display==0){
		dev.off();
	}
}
fun_des_read=function(filename){
	return(as.matrix(
		read.table(file=filename,header=FALSE,comment.char="#",blank.lines.skip=TRUE)));
}
fun_des_perminv=function(mat_des){
	for(i in regRange(1,nrow(mat_des))){
		mat_des[i,]=t(order(as.numeric(mat_des[i,]))-1);
	}
	return(mat_des);
}

fun_plot_eg_motivation=function(fileout,display=1){
	mat_des=fun_des_read("./des_fig/OofAdes_10_5_2_2_1.txt");
	n_comp=dim(mat_des)[2];
	list_cex=list(cex_point=5,lwd_point=2.5,lwd_sep=2.5);
	if(display==0){
		postscript(fileout,width=10,height=10,pointsize=16,paper="special",horizontal=FALSE);
	}
	par(mfrow=c(1,1),mar=c(5,5,0.5,0.5),mgp=c(3,1,0));
	plot(x=NULL,xlim=c(0,n_comp-1),ylim=c(0,n_comp-1),
		asp=1,main="",ylab="Component 1",xlab="Component 0",cex.axis=1,cex.lab=1);
	abline(v=1.5,lty=2,lwd=list_cex$lwd_sep,col=rgb(1,0,0));
	abline(h=1.5,lty=2,lwd=list_cex$lwd_sep,col=rgb(1,0,0));
	abline(v=3.5,lty=3,lwd=1.5,col=rgb(1,0,0));
	abline(h=3.5,lty=3,lwd=1.5,col=rgb(1,0,0));
	for(i0 in 0:(n_comp-1)){
		for(i1 in 0:(n_comp-1)){
			if(i0!=i1){
				points(x=i0,y=i1,
					pch=1,cex=list_cex$cex_point,lwd=list_cex$lwd_point,col=rgb(0.5,0.5,0.5));
			}
		}
	}
	points(x=mat_des[,c(1,2)],pch=21,cex=list_cex$cex_point,lwd=list_cex$lwd_point,
		col=rgb(0.5,0.5,0.5),bg=rgb(0,0,1));
	if(display==0){
		dev.off();
	}
}