ies2rsb<-function(hh_pop1,hh_pop2,popname1=NA,popname2=NA,method="bilateral",allow_different_dimensions=TRUE){

	
	if(!allow_different_dimensions){
  if(!(nrow(hh_pop1)==nrow(hh_pop2))){stop("hh_pop1 and hh_pop2 must have the same dimensions")}
  	ies_1=hh_pop1[,6] ; ies_2=hh_pop2[,6]
  if(sum(hh_pop1[,2]==hh_pop2[,2])<nrow(hh_pop1)){stop("SNP position in hh_pop1 and hh_pop2 must be the same")}
	}else{
		consenesus_rows = intersect(hh_pop1[,2],hh_pop2[,2])
		updated_pop = hh_pop1[match(consenesus_rows,hh_pop1[,2]),]
		ies_1 = hh_pop1[match(consenesus_rows,hh_pop1[,2]),6]
		ies_2 = hh_pop2[match(consenesus_rows,hh_pop2[,2]),6]
	}
	print(ies_1[1000])
	print(ies_2[1000])
  tmp_rsbnc=log(ies_1/ies_2) ;
  tmp_rsbnc[is.infinite(tmp_rsbnc)|is.nan(tmp_rsbnc)] = NA
  tmp_med=median(tmp_rsbnc,na.rm=T) ; tmp_sd=sd(tmp_rsbnc,na.rm=T)
  rsbcor=(tmp_rsbnc-tmp_med)/tmp_sd
  tmp_pval=rsbcor*0
  if(method=="bilateral"){tmp_pval=-1*log10(1-2*abs(pnorm(rsbcor)-0.5))}
  if(method=="unilateral"){tmp_pval=-1*log10(1-pnorm(rsbcor))}
  tmp_pval2=tmp_pval ; tmp_pval2[tmp_pval2=="Inf"]=NA  
  tmp_pval[tmp_pval=="Inf"]=max(tmp_pval2,na.rm=TRUE) + 1 
  res.rsb=cbind(updated_pop[,1:2],rsbcor,tmp_pval)
  colnames(res.rsb)[3]=paste("rSB (",popname1," vs ",popname2,")",sep="")
  colnames(res.rsb)[4]=paste("Pvalue (",method,")",sep="")

  return(list(res.rsb=res.rsb))
}
