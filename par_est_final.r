get.info<-function(str_intv){
	tmp_str_intv=as.character(str_intv)
	my_strsplit=as.numeric(gsub("\\[([0-9]*)-([0-9]*)\\]","\\2",tmp_str_intv))### transfer format,e.g [999-1000] to 1000
	return(my_strsplit)
}

get.range.time<-function(array.search,array.beg,array.end,array.bound=1){
	start.ind<-NULL
	if(all(array.search[array.beg:array.end]<array.bound)){
		start.ind<-array.end
	}
	else{
		start.meet<-which(array.search[array.beg:array.end]>=array.bound)
		if(length(start.meet)==1){
				start.ind<-array.end
		}
		else{
			if(start.meet[2]-start.meet[1]==1){
				start.ind<-start.meet[1]-1+array.beg-1
			}
			else{
				start.ind<-start.meet[2]-1+array.beg-1
			}
		}
	}
	return(start.ind)
}

get.lm<-function(array.y,array.x,range.list){
	array.bg<-range.list[2]
	array.ed<-range.list[3]
	array.fg<-range.list[4]
	if(array.fg==2){
		return(c(range.list[1],array.bg,array.ed,0))
	}
	else{
		y<-array.y[array.bg:array.ed]
		x<-array.x[array.bg:array.ed]
		lm.exp<-round(abs(lm(log(y)~x)$coefficients[2]),2)
		y.pos<-y[which(y<1.0)]
		x.pos<-x[which(y<1.0)]
		if(length(y.pos)==0){
			return(c(lm.exp,array.bg,array.ed,1))
		}
		else{
			lm.det<-round(abs(lm(log((1-y.pos)/y.pos)~x.pos)$coefficients[2]),2)
		}
		return(c(min(lm.exp,lm.det),array.bg,array.ed,1))
	}
}

get.range.tiny<-function(array.search,array.time,array.beg,array.end){
	search.range<-array.end-array.beg+1
	range.stack<-NULL
	ret.tiny<-c()
	if(search.range<2){
		ret.tiny<-rep(0,4)
	}
	else{
		for(j in rev(2:(search.range-1))){
			tmp.stack<-NULL
			for(i in 1:(search.range-j)){
				sub.beg<-array.beg-1+i
				sub.end<-array.beg-1+i+j
				tmp.a<-array.search[sub.end]-array.search[sub.beg]
				tmp.t<-array.time[sub.end]-array.time[sub.beg]
				tmp.s<-tmp.a/tmp.t
				tmp.sc<-c(tmp.s,sub.beg,sub.end)
				tmp.stack<-cbind(tmp.stack,tmp.sc)
			}
			tmp.max<-which.max(abs(tmp.stack[1,]))
			range.stack<-cbind(range.stack,c(tmp.stack[,tmp.max],j))
		}
		range.est<-apply(range.stack,2,get.lm,array.y=array.search,array.x=array.time)
		est.exp.ind<-which.max(range.est[1,])
		ret.tiny<-c(range.est[,est.exp.ind])#,range.est[,est.det.ind])
	}
	name.list<-c('s','bg','ed','tp')
	names(ret.tiny)<-name.list
	return(ret.tiny)		
}

get.start.t<-function(range.info){
	int.tol<-1000
	tree.tol<-4000
	intv.tol<-300
	alp1<-range.info['alp1']
	alp2<-range.info['alp2']
	alp3<-range.info['alp3']
	tau1<-range.info['tau1']
	tau2<-range.info['tau2']
	n1<-range.info['n1']
	n2<-range.info['n2']
	
	coal.size<-(length(range.info)-8)/3
	alp.info<-range.info[paste('w.t',1:coal.size,sep='')]
	tree.info<-range.info[paste('a.t',1:coal.size,sep='')]
	int.info<-range.info[paste('i.t',1:coal.size,sep='')]
	
	start.t<-NULL
	start.ind<-NULL
	start.info<-NULL
	coef.info<-NULL
	if(alp1>alp2&&alp2>alp3){### case1 +/-||+||++, speed slow down
		######################################	
		####_____						     #
		####     |_______		             #
		####			 |_____________      #
		####                                 #
		######################################
		start.ind<-get.range.time(alp.info,n2+1,coal.size,1.0)
		if(int.info[start.ind]>=int.tol){
			start.ind<-start.ind-1
		}
		if(tree.info[start.ind]<tree.tol){
			while(int.info[start.ind]>intv.tol&&start.ind>=n2){
				start.ind=start.ind-1
			}
		}
		start.t<-tree.info[start.ind]
		start.info<-c(start.t,1,ifelse(start.ind==coal.size,0,1))
		if(alp2<1.0){
			coef.info<-get.range.tiny(alp.info,tree.info,n1+1,coal.size)
		}
		else{
			coef.info<-get.range.tiny(alp.info,tree.info,n2+1,coal.size)
		}	
	}
	else if(alp1>alp2&&alp2<alp3){### case2 +/-||++||+/- speed up,convex case
		######################################	
		####	 							 #	
		####_____	 					     #
		####     |			    _________    #		             
		####	 |_____________|     		 #
		####                                 #
		######################################
		#test adjust case2:
		if(alp3<1.0){
			z1<-alp2/(n2-n1)
			z2<-alp3/(coal.size-n2)
			if(z1/(z2*2^(2*n2-n1-coal.size))>1){
				start.ind<-get.range.time(alp.info,n2,coal.size,1.0)
				if(int.info[start.ind]>=int.tol){
					start.ind<-start.ind-1
				}
				if(tree.info[start.ind]<tree.tol){
					while(int.info[start.ind]>intv.tol&&start.ind>=n2){
						start.ind=start.ind-1
					}
				}
				start.t<-tree.info[start.ind]
				start.info<-c(start.t,2,1)
				coef.info<-get.range.tiny(alp.info,tree.info,n2+1,coal.size)		
			}
			else{
				start.ind<-n2+1
				if(int.info[start.ind]>=int.tol){
					start.ind<-start.ind-1
				}
				if(tree.info[start.ind]<tree.tol){
					while(int.info[start.ind]>intv.tol&&start.ind>=n2){
						start.ind=start.ind-1
					}
				}
				start.t<-tree.info[start.ind]
				start.info<-c(start.t,2,1)
				coef.info<-get.range.tiny(alp.info,tree.info,n1+1,n2)
			}
		
		}
		######################################	
		else{#################################
			start.ind<-n2
			if(int.info[start.ind]>=int.tol){
				start.ind<-start.ind-1
			}
			if(tree.info[start.ind]<tree.tol){
				while(int.info[start.ind]>intv.tol&&start.ind>=n2){
					start.ind=start.ind-1
				}
			}
			start.t<-tree.info[start.ind]
			start.info<-c(start.t,2,1)
			coef.info<-get.range.tiny(alp.info,tree.info,n1+1,n2)
		}#####################################
	}
	else if(alp1<alp2&&alp2>alp3){### case3 ++/+||+/-||++/+ speed slow down concave
		######################################	
		####      _____________				 #
		####	 |		 	   |             #
		####_____|			   |             #		             
		####	               |________     #
		####                                 #
		######################################
		if(all(c(alp1,alp2,alp3)<1)){
			start.ind<-get.range.time(alp.info,n2+1,coal.size,1.0)
			if(int.info[start.ind]>=int.tol){
				start.ind<-start.ind-1
			}
			if(tree.info[start.ind]<tree.tol){
				while(int.info[start.ind]>intv.tol&&start.ind>=n2){
					start.ind=start.ind-1
				}
			}
			start.t<-tree.info[start.ind]
			start.info<-c(start.t,3,ifelse(start.ind==coal.size,0,1))
			if(alp1<=alp3){
				coef.info<-get.range.tiny(alp.info,tree.info,1,n1+1)
			}
			else{
				coef.info<-get.range.tiny(alp.info,tree.info,n2+1,coal.size)
			}
		}
		else{	
		if(alp1<1.0&&alp3<1.0){
			z1<-alp1/(n1)
			z2<-alp3/(coal.size-n2)
			if(z1/(z2*2^(n1+n2-coal.size))>1){### accleration in the 3rd block
				start.ind<-get.range.time(alp.info,n2+1,coal.size,1.0)
				if(int.info[start.ind]>=int.tol){
					start.ind<-start.ind-1
				}
				if(tree.info[start.ind]<tree.tol){
					while(int.info[start.ind]>intv.tol&&start.ind>=n2){
						start.ind=start.ind-1
					}
				}
				start.t<-tree.info[start.ind]
				start.info<-c(start.t,3,ifelse(start.ind==coal.size,0,1))
				coef.info<-get.range.tiny(alp.info,tree.info,n2+1,coal.size)
			}
			else{### accleration in the 1st block
				start.ind<-n1
				if(int.info[start.ind]>=int.tol){
					start.ind<-start.ind-1
				}
				if(tree.info[start.ind]<tree.tol){
					while(int.info[start.ind]>intv.tol&&start.ind>=n1){
						start.ind=start.ind-1
					}
				}
				start.t<-tree.info[start.ind]
				start.info<-c(start.t,3,1)
				coef.info<-get.range.tiny(alp.info,tree.info,1,n1)
			}
		}
		else if(alp1<1.0){### in the 1st block
			start.ind<-n1
			if(int.info[start.ind]>=int.tol){
				start.ind<-start.ind-1
			}
			if(tree.info[start.ind]<tree.tol){
				while(int.info[start.ind]>intv.tol&&start.ind>=n1){
					start.ind=start.ind-1
				}
			}
			start.t<-tree.info[start.ind]
			start.info<-c(start.t,3,1)
			coef.info<-get.range.tiny(alp.info,tree.info,1,n1)
		}
		else{### in the 3rd block
			start.ind<-get.range.time(alp.info,n2+1,coal.size,1.0)
			if(int.info[start.ind]>=int.tol){
				start.ind<-start.ind-1
			}
			if(tree.info[start.ind]<tree.tol){
				while(int.info[start.ind]>intv.tol&&start.ind>=n2){
					start.ind=start.ind-1
				}
			}
			start.t<-tree.info[start.ind]
			start.info<-c(start.t,3,ifelse(start.ind==coal.size,0,1))
			coef.info<-get.range.tiny(alp.info,tree.info,n2+1,coal.size)
		}
		}
	}
	else{### case4 +/-||+/-||++/+
		######################################	
		####                   				 #
		####	  		 	    ___________  #
		####      _____________|             #        
		####_____|                           #
		####                                 #
		######################################
		if(alp2<1){
			if(alp3<1){
				start.ind<-get.range.time(alp.info,n2+1,coal.size,1.0)
				if(int.info[start.ind]>=int.tol){
					start.ind<-start.ind-1
				}
				if(tree.info[start.ind]<tree.tol){
					while(int.info[start.ind]>intv.tol&&start.ind>=n2){
						start.ind=start.ind-1
					}
				}
				start.t<-tree.info[start.ind]
				start.info<-c(start.t,3,ifelse(start.ind==coal.size,0,1))
			}
			else{
				start.ind<-n2
				if(int.info[start.ind]>=int.tol){
					start.ind<-start.ind-1
				}
				if(tree.info[start.ind]<tree.tol){
					while(int.info[start.ind]>intv.tol&&start.ind>=n2){
						start.ind=start.ind-1
					}
				}
				start.t<-tree.info[start.ind]
				start.info<-c(start.t,4,1)
			}
		}
		else{
			start.ind<-n1
			if(int.info[start.ind]>=int.tol){
				start.ind<-start.ind-1
			}
			if(tree.info[start.ind]<tree.tol){
				while(int.info[start.ind]>intv.tol&&start.ind>=n1){
					start.ind=start.ind-1
				}
			}
			start.t<-tree.info[start.ind]
			start.info<-c(start.t,4,1)
		}
		coef.info<-get.range.tiny(alp.info,tree.info,1,n1)
	}
	start.info<-start.info[1]
	names(start.info)<-'est.time'
	return(c(start.info,coef.info))
}		

get.filter<-function(res.file,pop,if_an=F){###fileter the file and collect all the information
	#############################################
	#############################################
	cenpos<-1000000
	unit  <-100000
	ret.res<-NULL
	cen_seq<-NULL
	#############################################
	get.power<-function(res.file,pop,if_an){
		res<-read.table(res.file)
		ps.info<-get.info(res[[1]])
		select.var<-switch(pop,YRI=1,1.5)
		cen.bg<-which(ps.info>=cenpos-select.var*unit)[1]
		cen.ed<-which(ps.info>cenpos+select.var*unit)[1]-1
		res<-res[cen.bg:cen.ed,]
		ps.range<-diff(c(cenpos-select.var*unit,ps.info[cen.bg:cen.ed]))
		if(!if_an){
			thres<-switch(pop,CHB=c(0.527,5.9),CEU=c(0.528,6.4),YRI=c(0.5,7.5))
			if(all(res[[2]]>=thres[1])){
				return(0)
			}
			else{
				if(all(res[[7]]-res[[3]]<thres[2])){
					return(0)
				}
				else{
					cen_seq<<-cen.bg:cen.ed
					ret.res<<-cbind(ps.info[cen.bg:cen.ed],res[,c(4:6,8:11)])
					return(1)
				}
			}
		}
		else{
			thres<-switch(pop,CHB=c(0.527,10.8),CEU=c(0.528,11.7),YRI=c(0.5,13.1))
			sgm.index<-which(res[[7]]-res[[3]]>=thres[2])
			sgm.porpt<-sum(ps.range[sgm.index])/(select.var*unit)
			if(any(res[[2]]>thres[1]&res[[2]]<1)&& any(res[[7]]-res[[3]]>=thres[2])&&sgm.porpt>=0.1){
				cen_seq<<-cen.bg:cen.ed
				ret.res<<-cbind(ps.info[cen.bg:cen.ed],res[,c(4:6,8:11)])
				return(1)
			}
			else{
				return(0)
			}
		}
		
	}	
	#############################################
	if(!if_an){
		s.time<-as.numeric(gsub('[a-zA-Z0-9_/]*_(\\d+)_(\\d+\\.\\d+)_.*','\\1',res.file))
		s.coef<-as.numeric(gsub('[a-zA-Z0-9_/]*_(\\d+)_(\\d+\\.\\d+)_.*','\\2',res.file))
	}
	else{
		s.time<-as.numeric(gsub('[a-zA-Z0-9_/.]*_(\\d+\\.\\d+)_(\\d+\\.\\d+)_.*','\\1',res.file))*40000
		s.coef<-as.numeric(gsub('[a-zA-Z0-9_/.]*_(\\d+\\.\\d+)_(\\d+\\.\\d+)_.*','\\2',res.file))
	}
	if(!get.power(res.file,pop,if_an)){
		return(NULL)
	}
	else{
		colnames(ret.res)<-c('pos','alp1','alp2','alp3','tau1','tau2','n1','n2')
		alp.file<-gsub('res','alp',res.file)
		ret.alp<-as.matrix(read.table(alp.file)[cen_seq,-1])
		colnames(ret.alp)<-paste('w.t',1:ncol(ret.alp),sep='')

		tree.file<-gsub('res','tree',res.file)
		ret.tree<-as.matrix(read.table(tree.file)[cen_seq,-1])
		ret.int<-ret.tree
		colnames(ret.int)<-paste('i.t',1:ncol(ret.int),sep='')
		ret.tree<-t(apply(ret.tree,1,cumsum))	
		colnames(ret.tree)<-paste('a.t',1:ncol(ret.tree),sep='')
		range.mat<-cbind(ret.res,ret.alp,ret.tree,ret.int)
		ret.start.t<-t(apply(range.mat,1,get.start.t))
		ret.median<-cbind(s.time,s.coef,median(ret.start.t[,1]),median(ret.start.t[,4]))
		names(ret.median)<-c('time','coef','e.time','e.coef')
		return(ret.median)
	}
}

get.files<-function(indir,pop,if_an=F){
	indir<-gsub('(\\/)?$','\\/',indir)
	files<-list.files(paste(indir,'res',sep=''),full.names=T)
	dt.summary<-do.call(rbind,lapply(files,get.filter,pop=pop,if_an=if_an))
	colnames(dt.summary)<-c('time','coef','e.time','e.coef')
	return(dt.summary)
}

my_args<-commandArgs(T)
indir<-my_args[1]
pop<-my_args[2]
if_an<-as.numeric(my_args[3])
outfiles<-my_args[4]
pop<-get.files(indir,pop,if_an)
write.table(pop,outfiles,row.names=F,quote=F)
q(save='no')	
