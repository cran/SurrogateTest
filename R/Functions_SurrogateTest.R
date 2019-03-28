early.delta.test = function(Axzero, Adeltazero, Aszero, Bxzero, Bdeltazero, Bszero, Bxone, Bdeltaone, Bsone, t, landmark, perturb = T, extrapolate = T, transform = F) {
	n0.A = length(Axzero)
	n1.B = length(Bxone)
	n0.B = length(Bxzero)
	n.B = n1.B+n0.B
	if(n0.A < 200 | n1.B < 200 | n0.B < 200) {print("Warning: Samples sizes are small; this procedure may not produce stable results.")}
	 if(transform){	 	    	
    	mean.o= mean(c(Aszero[Axzero>landmark], Bszero[Bxzero>landmark], Bsone[Bxone>landmark]), na.rm = T)
  		sd.o = sd(c(Aszero[Axzero>landmark], Bszero[Bxzero>landmark], Bsone[Bxone>landmark]), na.rm = T)
    	Aszero.new = pnorm((Aszero[Axzero>landmark] - mean.o)/sd.o)
    	Bszero.new = pnorm((Bszero[Bxzero>landmark] - mean.o)/sd.o)
    	Bsone.new = pnorm((Bsone[Bxone>landmark] - mean.o)/sd.o)
    	sA0.new.all = rep(NA, length(Aszero))
		sB0.new.all = rep(NA, length(Bszero))
		sB1.new.all = rep(NA, length(Bsone))		
		sA0.new.all[Axzero > landmark] = Aszero.new
		sB0.new.all[Bxzero > landmark] = Bszero.new
		sB1.new.all[Bxone > landmark] = Bsone.new
    	Aszero = sA0.new.all
    	Bszero = sB0.new.all
    	Bsone = sB1.new.all
	} 
	A0.range = range(Aszero, na.rm = T)
	B1.range = range(Bsone, na.rm = T); B0.range = range(Bszero, na.rm = T)
	if((B1.range[1]< A0.range[1]) | (B0.range[1]< A0.range[1]) |  (B1.range[2]> A0.range[2]) |  (B0.range[2]> A0.range[2]) ) { print("Warning: Surrogate value ranges do not appear to overlap; consider transforming the surrogate values.")}
	delta.eb = delta.eb.single(Axzero = Axzero, Adeltazero = Adeltazero, Aszero = Aszero, Bxzero = Bxzero, Bdeltazero = Bdeltazero, Bszero = Bszero, Bxone = Bxone, Bdeltaone = Bdeltaone, Bsone = Bsone, t = t, landmark=landmark, extrapolate = extrapolate)
	var.closed.eb = var.delta.eb(Axzero = Axzero, Adeltazero = Adeltazero, Aszero = Aszero, Bxzero = Bxzero, Bdeltazero = Bdeltazero, Bszero = Bszero, Bxone = Bxone, Bdeltaone = Bdeltaone, Bsone = Bsone, t = t, landmark=landmark, extrapolate = extrapolate)
	se.closed.deltaeb = sqrt(var.closed.eb)/sqrt(n.B)
	Z.closed.eb = sqrt(n.B)*delta.eb/sqrt(var.closed.eb)
	p.closed = (1-pnorm(Z.closed.eb))*2
	conf.closed.norm = c(delta.eb - 1.96*sqrt(var.closed.eb)/sqrt(n.B), delta.eb + 1.96*sqrt(var.closed.eb)/sqrt(n.B))
	if(perturb){
		weightA.mat = matrix(rexp(500*(n0.A), rate=1), ncol = 500)
		weightB.mat = matrix(rexp(500*(n1.B+n0.B), rate=1), ncol = 500)
		weight.together = rbind(weightA.mat, weightB.mat)
		delta.eb.p.vec = apply(weight.together, 2, delta.eb.single, Axzero = Axzero, Adeltazero = Adeltazero, Aszero = Aszero, Bxzero = Bxzero, Bdeltazero = Bdeltazero, Bszero = Bszero, Bxone = Bxone, Bdeltaone = Bdeltaone, Bsone = Bsone, t = t, landmark=landmark, extrapolate = extrapolate, weightA=NULL, weightB = NULL)
		var.perturb.eb = n.B*var(delta.eb.p.vec)
		se.perturb.deltaeb = sd(delta.eb.p.vec)
		conf.quantile.delta.eb = c(quantile(delta.eb.p.vec, 0.025 ), quantile(delta.eb.p.vec, 0.975 ))
		Z.perturb.eb = sqrt(n.B)*delta.eb/sqrt(var.perturb.eb)
		p.perturb = (1-pnorm(Z.perturb.eb))*2
		conf.perturb.norm = c(delta.eb - 1.96*sqrt(var.perturb.eb)/sqrt(n.B), delta.eb + 1.96*sqrt(var.perturb.eb)/sqrt(n.B))
	}
	if(!perturb){
		return(list("delta.eb" = delta.eb, "se.closed" = se.closed.deltaeb,"Z.closed" = Z.closed.eb, "p.value.closed" = p.closed, "conf.closed.norm" = conf.closed.norm))}
	if(perturb){
		return(list("delta.eb" = delta.eb, "se.closed" = se.closed.deltaeb,"Z.closed" = Z.closed.eb, "p.value.closed" = p.closed, "conf.closed.norm" = conf.closed.norm, "se.perturb" = se.perturb.deltaeb,"Z.perturb" = Z.perturb.eb, "p.value.perturb" = p.perturb, "conf.perturb.norm" = conf.perturb.norm, "delta.eb.CI" = conf.quantile.delta.eb))}
}

delta.eb.single = function(Axzero, Adeltazero, Aszero, Bxzero, Bdeltazero, Bszero, Bxone, Bdeltaone, Bsone, t, landmark, weightA = NULL, weightB = NULL, weight.both = NULL,  extrapolate) {
	if(!is.null(weight.both)) {
		weightA = weight.both[c(1:length(Axzero))]
		weightB = weight.both[-c(1:length(Axzero))]	
	}
	if(is.null(weightA)) {weightA = rep(1,length(Axzero))}
	if(is.null(weightB)) {weightB = rep(1,length(Bxone)+length(Bxzero))}
	weightB.group1 = weightB[1:length(Bxone)]
	weightB.group0 = weightB[(1+length(Bxone)):(length(Bxone)+length(Bxzero))]
	mu.s = pred.smooth.surv.new(Axzero.f=Axzero[Axzero>landmark], Adeltazero.f=Adeltazero[Axzero>landmark], Aszero.f=Aszero[Axzero>landmark], Bsnew.f=Bsone[Bxone>landmark], Bsnew2.f = Bszero[Bxzero>landmark], myt=t, weight.pred = weightA[Axzero>landmark], extrapolate = extrapolate)
	censorB1.landmark = censor.weight(Bxone, Bdeltaone, landmark, weight = weightB.group1)
	censorB0.landmark = censor.weight(Bxzero, Bdeltazero, landmark, weight = weightB.group0)
	G1 = sum(weightB.group1[Bxone>landmark]*mu.s$Phat.ss.1)/sum(weightB.group1*censorB1.landmark) 
	G0 = sum(weightB.group0[Bxzero>landmark]*mu.s$Phat.ss.2)/sum(weightB.group0*censorB0.landmark) 
	delta.eb = G1-G0	
	return(delta.eb)
}

pred.smooth.surv.new = function(Axzero.f, Adeltazero.f, Aszero.f, Bsnew.f, Bsnew2.f=NULL, myt, bw = NULL, weight.pred, extrapolate)
   #Takes group 0 from study A and creates the estimator, then returns the predictions/approximations for Bsnew; or you can send both groups from Study B at the same time using Bsnew2 in which case this returns a list with a vector for each group
  { 
  	
    if(is.null(bw))
      {
        bwini = bw.nrd(Aszero.f)
        n.s = length(Aszero.f)
        bw <- bwini/(n.s^0.11)
      }      
    kerni.ss = Kern.FUN(zz=Aszero.f,zi=Bsnew.f,bw)           
    tmpind = (Axzero.f<=myt)&(Adeltazero.f==1); tj = Axzero.f[tmpind]; 
    kerni.1 = t(weight.pred*t(kerni.ss))  #weights for study A, group 0
    pihamyt0.tj.ss = helper.si(tj, "<=", Axzero.f, Vi=t(kerni.1)) ## n.tj x n.ss matrix ##   
    dLamhat.tj.ss = t((kerni.1[,tmpind]))/pihamyt0.tj.ss; 
    ret = apply(dLamhat.tj.ss,2,sum)
    Phat.ss  =exp(-ret)
    if(sum(is.na(Phat.ss))>0 & extrapolate){
    	c.mat = cbind(Bsnew.f, Phat.ss)
    	for(o in 1:length(Phat.ss)) {
    		if(is.na(Phat.ss[o])){
    			distance = abs(Bsnew.f - Bsnew.f[o])
    			c.temp = cbind(c.mat, distance)
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where predication is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			Phat.ss[o] = new.est
    	}
    }
	}
	if(is.null(Bsnew2.f)) { return(Phat.ss)}
	if(!is.null(Bsnew2.f)){
		Phat.ss.1 = Phat.ss
		kerni.ss = Kern.FUN(zz=Aszero.f,zi=Bsnew2.f,bw)           
    	tmpind = (Axzero.f<=myt)&(Adeltazero.f==1); tj = Axzero.f[tmpind]; 
    	kerni.1 = t(weight.pred*t(kerni.ss))  #weights for study A, group 0
    	pihamyt0.tj.ss = helper.si(tj, "<=", Axzero.f, Vi=t(kerni.1)) ## n.tj x n.ss matrix ##   
    	dLamhat.tj.ss = t((kerni.1[,tmpind]))/pihamyt0.tj.ss; 
    	ret = apply(dLamhat.tj.ss,2,sum)
    	Phat.ss  =exp(-ret)
    	if(sum(is.na(Phat.ss))>0 & extrapolate){
    		c.mat = cbind(Bsnew2.f, Phat.ss)
    		for(o in 1:length(Phat.ss)) {
    			if(is.na(Phat.ss[o])){
    				distance = abs(Bsnew2.f - Bsnew2.f[o])
    				c.temp = cbind(c.mat, distance)
    				c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where predication is not na
    				new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    				Phat.ss[o] = new.est
    			}	
    		}
			}
		return(list("Phat.ss.1" = Phat.ss.1, "Phat.ss.2" = Phat.ss)	)
	}   
    }


censor.weight = function(data.x, data.delta, t, weight=NULL) {
	if(is.null(weight)) {weight = rep(1,length(data.x))}
	S.KM = survfit(Surv(data.x,1-data.delta)~1, weights = weight)
	S.t.KM = approx(S.KM$time,S.KM$surv,t)$y
	return(S.t.KM)
}

cumsum2 <- function(mydat)     #cumsum by row, col remains the same
  {
    if(is.null(dim(mydat))) return(cumsum(mydat))
    else{
      out <- matrix(cumsum(mydat), nrow=nrow(mydat))
      out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
      return(out)
    }
  }

helper.si <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
  {  
    if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
    if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
    pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
    if(is.null(Vi)){return(pos)}else{
      Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
      out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
      out[pos!=0,] <- Vi[pos,]
      if(is.null(dim(Vi))) out <- c(out)
      return(out) ## n.y x p
    }
  } 



VTM=function(vc, dm){
	#takes vc and makes it the repeated row of a matrix, repeats it dm times
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }
    
Kern.FUN = function(zz,zi,bw,kern0="gauss") ## returns an (n x nz) matrix ##
  { 
    out = (VTM(zz,length(zi))- zi)/bw
    switch(kern0,
            "epan"= 0.75*(1-out^2)*(abs(out)<=1)/bw,
            "gauss"= dnorm(out)/bw
           )
  }

recover.B= function(Axzero, Adeltazero, Aszero, Axone, Adeltaone, Asone, Bxzero, Bdeltazero, Bszero, Bxone, Bdeltaone, Bsone, t, landmark, extrapolate = T, transform = F){
	n0.A = length(Axzero)
	n1.A = length(Axone)
	n1.B = length(Bxone)
	n0.B = length(Bxzero)
	if(n0.A < 200 | n1.B < 200 | n0.B < 200 | n1.A < 200) {print("Warning: Samples sizes are small; this procedure may not produce stable results.")}
	if(transform){	 	
    	mean.o= mean(c(Aszero[Axzero>landmark], Asone[Axone>landmark], Bszero[Bxzero>landmark], Bsone[Bxone>landmark]), na.rm = T)
  		sd.o = sd(c(Aszero[Axzero>landmark], Asone[Axone>landmark], Bszero[Bxzero>landmark], Bsone[Bxone>landmark]), na.rm = T)
    	Aszero.new = pnorm((Aszero[Axzero>landmark] - mean.o)/sd.o)
    	Asone.new = pnorm((Asone[Axone>landmark] - mean.o)/sd.o)
    	Bszero.new = pnorm((Bszero[Bxzero>landmark] - mean.o)/sd.o)
    	Bsone.new = pnorm((Bsone[Bxone>landmark] - mean.o)/sd.o)
    	sA0.new.all = rep(NA, length(Aszero))
		sA1.new.all = rep(NA, length(Asone))
		sB0.new.all = rep(NA, length(Bszero))
		sB1.new.all = rep(NA, length(Bsone))		
		sA0.new.all[Axzero > landmark] = Aszero.new
		sA1.new.all[Axone > landmark] = Asone.new
		sB0.new.all[Bxzero > landmark] = Bszero.new
		sB1.new.all[Bxone > landmark] = Bsone.new
    	Aszero = sA0.new.all
    	Asone = sA1.new.all
    	Bszero = sB0.new.all
    	Bsone = sB1.new.all
	}
	A1.range = range(Asone, na.rm = T); A0.range = range(Aszero, na.rm = T)
	B1.range = range(Bsone, na.rm = T); B0.range = range(Bszero, na.rm = T)
	if((A1.range[1]< A0.range[1]) | (A1.range[2]> A0.range[2]) | (B1.range[1]< A0.range[1]) | (B0.range[1]< A0.range[1]) |  (B1.range[2]> A0.range[2]) |  (B0.range[2]> A0.range[2]) ) { print("Warning: Surrogate value ranges do not appear to overlap; consider transforming the surrogate values.")}
	n.B = n1.B+n0.B
	delta.eb = delta.eb.single(Axzero = Axzero, Adeltazero = Adeltazero, Aszero = Aszero, Bxzero = Bxzero, Bdeltazero = Bdeltazero, Bszero = Bszero, Bxone = Bxone, Bdeltaone = Bdeltaone, Bsone = Bsone, t = t, landmark=landmark, extrapolate = extrapolate)
	delta.A = delta.estimate(xone = Axone, xzero = Axzero, deltaone = Adeltaone, deltazero = Adeltazero, t = t)
	delta.ea = delta.ea.single(Axzero = Axzero, Adeltazero = Adeltazero, Aszero = Aszero, Bxzero = Axzero, Bdeltazero = Adeltazero, Bszero = Aszero, Bxone = Axone, Bdeltaone = Adeltaone, Bsone = Asone, t = t, landmark=landmark, extrapolate = extrapolate)
	Ra = delta.ea/delta.A
	recovered.deltaB = delta.eb/Ra
	weightA.mat = matrix(rexp(500*(n1.A+n0.A), rate=1), ncol = 500)
	weightB.mat = matrix(rexp(500*(n1.B+n0.B), rate=1), ncol = 500)
	delta.A.p = apply(weightA.mat, 2, delta.estimate, xone = Axone, xzero = Axzero, deltaone = Adeltaone, deltazero = Adeltazero, t = t)
	weight.together = rbind(weightA.mat[c((n1.A+1):(n1.A+n0.A)),], weightB.mat)
	delta.eb.p = apply(weight.together, 2, delta.eb.single, Axzero = Axzero, Adeltazero = Adeltazero, Aszero = Aszero, Bxzero = Bxzero, Bdeltazero = Bdeltazero, Bszero = Bszero, Bxone = Bxone, Bdeltaone = Bdeltaone, Bsone = Bsone, t = t, landmark=landmark, extrapolate = extrapolate, weightA=NULL, weightB = NULL)
	weight.together = rbind(weightA.mat[c((n1.A+1):(n1.A+n0.A)),], weightA.mat)
	delta.ab.p = apply(weight.together, 2, delta.ea.single, Axzero = Axzero, Adeltazero = Adeltazero, Aszero = Aszero, Bxzero = Axzero, Bdeltazero = Adeltazero, Bszero = Aszero, Bxone = Axone, Bdeltaone = Adeltaone, Bsone = Asone,t = t, landmark=landmark, extrapolate = extrapolate, weightA=NULL, weightB = NULL)
	Ra.p = delta.ab.p/delta.A.p
	recovered.deltab.p = delta.eb.p/(delta.ab.p/delta.A.p)
	sd.recovered.deltaB  = sd(recovered.deltab.p)
	conf.quantile.recovered.deltaB = c(quantile(recovered.deltab.p, 0.025 ), quantile(recovered.deltab.p, 0.975 ))
	return(list("recovered.deltaB" = recovered.deltaB, "se.recovered.deltaB" = sd.recovered.deltaB,"conf.quantile.recovered.deltaB" = conf.quantile.recovered.deltaB))
}

delta.estimate= function(xone,xzero, deltaone, deltazero, t, weight = NULL, KM = FALSE) {
	if(is.null(weight)) {weight = rep(1,length(xone)+length(xzero))}
	censor1.t = censor.weight(xone, deltaone, t, weight = weight[1:length(xone)])
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight[(1+length(xone)):(length(xone)+length(xzero))])
	delta = sum(1*(xone > t)*weight[1:length(xone)])/sum(weight[1:length(xone)]*censor1.t) - sum(1*(xzero > t)*weight[(1+length(xone)):(length(xone)+length(xzero))])/sum(weight[(1+length(xone)):(length(xone)+length(xzero))]* censor0.t)
	S.KM.1 = survfit(Surv(xone,deltaone)~1, weights = weight[1:length(xone)])
	S.t.KM.1 = approx(S.KM.1$time,S.KM.1$surv,t)$y
	S.KM.0 = survfit(Surv(xzero,deltazero)~1, weights = weight[(1+length(xone)):(length(xone)+length(xzero))])
	S.t.KM.0 = approx(S.KM.0$time,S.KM.0$surv,t)$y
	if(!KM) {return(delta)}
	if(KM) {return(S.t.KM.1 - S.t.KM.0)}

}

delta.ea.single = function(Axzero, Adeltazero, Aszero, Bxzero, Bdeltazero, Bszero, Bxone, Bdeltaone, Bsone,  t, landmark, weightA = NULL, weightB = NULL, weight.both = NULL, extrapolate) {
	if(!is.null(weight.both)) {
		weightA = weight.both[c(1:length(Axzero))]
		weightB = weight.both[-c(1:length(Axzero))]	
	}
	if(is.null(weightA)) {weightA = rep(1,length(Axzero))}
	if(is.null(weightB)) {weightB = rep(1,length(Bxone)+length(Bxzero))}
	weightB.group1 = weightB[1:length(Bxone)]
	weightB.group0 = weightB[(1+length(Bxone)):(length(Bxone)+length(Bxzero))]
	mu.s = pred.smooth.surv.new(Axzero.f=Axzero[Axzero>landmark], Adeltazero.f=Adeltazero[Axzero>landmark], Aszero.f=Aszero[Axzero>landmark], Bsnew.f=Bsone[Bxone>landmark], Bsnew2.f = Bszero[Bxzero>landmark], myt=t, weight.pred = weightA[Axzero>landmark], extrapolate = extrapolate)
	censorB1.landmark = censor.weight(Bxone, Bdeltaone, landmark, weight = weightB.group1)
	censorB0.landmark = censor.weight(Bxzero, Bdeltazero, landmark, weight = weightB.group0)
	G1 = sum(weightB.group1[Bxone>landmark]*mu.s$Phat.ss.1)/sum(weightB.group1*censorB1.landmark) 
	censorB0.t = censor.weight(Bxzero, Bdeltazero, t, weight = weightB.group0)
	G0= sum(1*(Bxzero > t)*weightB.group0)/sum(weightB.group0* censorB0.t) 
	delta.eb = G1-G0	
	return(delta.eb)
}

design.study= function(Axzero, Adeltazero, Aszero, Axone=NULL, Adeltaone=NULL, Asone=NULL, delta.ea = NULL, psi = NULL, R.A.given = NULL, t, landmark, extrapolate = T, adjustment = F,n=NULL, power=NULL,pi.1=0.5,pi.0=0.5, cens.rate, transform = F){
	#checking if adjustment is needed
	S.KM = survfit(Surv(Axzero,Adeltazero)~1)
	S.t.KM = approx(S.KM$time,S.KM$surv,t)$y
	if(S.t.KM>=0.90) {adjustment = T}
	if(!is.null(Axone) & !is.null(Adeltaone)){
		S.KM = survfit(Surv(Axone,Adeltaone)~1)
		S.t.KM = approx(S.KM$time,S.KM$surv,t)$y
		if(S.t.KM>=0.90) {adjustment = T}
	}
	n0.A = length(Axzero)
	if(!is.null(Axone)) {n1.A = length(Axone)}
	if(is.null(Axone)) {n1.A = 1000}
	if(n0.A < 200 | n1.A < 200) {print("Warning: Samples sizes are small; this procedure may not produce stable results.")}
	if(transform){	 	
    	mean.o= mean(c(Aszero[Axzero>landmark], Asone[Axone>landmark]), na.rm = T)
  		sd.o = sd(c(Aszero[Axzero>landmark], Asone[Axone>landmark]), na.rm = T)
    	Aszero.new = pnorm((Aszero[Axzero>landmark] - mean.o)/sd.o)
    	Asone.new = pnorm((Asone[Axone>landmark] - mean.o)/sd.o)
    	s0.new.all = rep(NA, length(Aszero))
		s1.new.all = rep(NA, length(Asone))
		s0.new.all[Axzero > landmark] = Aszero.new 
		s1.new.all[Axone > landmark] = Asone.new
    	Aszero = s0.new.all
    	Asone = s1.new.all
	} 
	if(!is.null(Asone)) {
		A1.range = range(Asone, na.rm = T); A0.range = range(Aszero, na.rm = T)
		if((A1.range[1]< A0.range[1]) | (A1.range[2]> A0.range[2])) { print("Warning: Surrogate value ranges do not appear to overlap; consider transforming the surrogate values.")}
	}
	param.need = FALSE
	if(!is.null(Axone) & !is.null(Adeltaone) & !is.null(Asone)){ param.need = TRUE}
	if(!is.null(delta.ea) & !is.null(psi)) { param.need = TRUE}
	if(!is.null(R.A.given) & !is.null(psi)) { param.need = TRUE}
	if(!param.need) {print("Error: Need to provide either (1) data from treatment arm in Study A or (2) hypothesized values for delta.ea (or R.A) and psi or (3) data from treatment arm in Study A and hypothesized psi (if different from observed treatment effect at t in Study A).")}
	if(!is.null(delta.ea) & !is.null(psi)) {
		if(delta.ea > psi) {print("Warning: delta.ea must not be larger than psi.")}
	}
	if(!is.null(R.A.given)) {
		if(R.A.given < 0 | R.A.given > 1) {print("Warning: R must be between 0 and 1.")}
	}
	if(is.null(n) & is.null(power)) {
		param.need = FALSE; print("Error: either n or power should be provided.")
	}
	if(param.need){
		if(is.null(R.A.given) & is.null(delta.ea) & is.null(psi)){
			delta.ea = delta.ea.single(Axzero = Axzero, Adeltazero = Adeltazero, Aszero = Aszero, Bxzero = Axzero, Bdeltazero = Adeltazero, Bszero = Aszero, Bxone = Axone, Bdeltaone = Adeltaone, Bsone = Asone, t = t, landmark=landmark, extrapolate = extrapolate)
			psi = delta.estimate(xone = Axone, xzero = Axzero, deltaone = Adeltaone, deltazero = Adeltazero, t = t)
			R.A = delta.ea/psi
			if( delta.ea> psi) {print("Warning: delta.ea (the early treatment effect at t0) should not be larger than psi (the treatment effect at t) in Study A.")}
			}  else if(is.null(R.A.given) & is.null(delta.ea) & !is.null(psi)){
			delta.ea = delta.ea.single(Axzero = Axzero, Adeltazero = Adeltazero, Aszero = Aszero, Bxzero = Axzero, Bdeltazero = Adeltazero, Bszero = Aszero, Bxone = Axone, Bdeltaone = Adeltaone, Bsone = Asone, t = t, landmark=landmark, extrapolate = extrapolate)
			delta.A = delta.estimate(xone = Axone, xzero = Axzero, deltaone = Adeltaone, deltazero = Adeltazero, t = t)
			if( delta.ea> psi) {print("Warning: delta.ea (the early treatment effect at t0) should not be larger than psi (the treatment effect at t) in Study A.")}
			R.A = delta.ea/delta.A
		}  else if(!is.null(delta.ea) & !is.null(psi) & is.null(Axone) & is.null(Adeltaone) & is.null(Asone)){
			R.A = delta.ea/psi
		}
		if(!is.null(R.A.given)) {R.A = R.A.given}
		mu.s = pred.smooth.surv.new(Axzero.f=Axzero[Axzero>landmark], Adeltazero.f=Adeltazero[Axzero>landmark], Aszero.f=Aszero[Axzero>landmark], Bsnew.f=Aszero[Axzero>landmark], myt=t, weight.pred = rep(1, sum(Axzero>landmark)), extrapolate = extrapolate)
		survB0.t = censor.weight(Axzero, 1-Adeltazero, t)
		censorA0.landmark = censor.weight(Axzero, Adeltazero, landmark)
		mu1hatzero=sum(mu.s)/censorA0.landmark/length(Axzero)

		censorB0.landmark = 1-pexp(landmark, rate = cens.rate)
		term2 = mu1hatzero^2*(1+int.hazc.plan(landmark = landmark, delta = Adeltazero, x = Axzero, cens.rate= cens.rate, B=10000))
		
		term1 = (1/length(Axzero))*sum(mu.s^2/(censorA0.landmark*censorB0.landmark))
		var.term1 = (1/(pi.0*pi.1))*(term1-term2)
		if(!adjustment){
			var.term=var.term1
		}
		if(adjustment){
			p0=censor.weight(Axzero, 1-Adeltazero, landmark)
			p1=p0+ R.A*psi/mu1hatzero
			f = (p1*(1-p1)*pi.0)/(p0*(1-p0)*pi.1)
			var.term =var.term1*(1+f)/(1+pi.0/pi.1)
		}
		if(is.null(n)){
			n = (sqrt(var.term)*(1.96-qnorm(1-power))/(R.A*psi))^2
			return(list("n" = n))
		}
		if(is.null(power)){
			power = 1-pnorm(1.96 - (sqrt(n)*R.A*psi)/(sqrt(var.term)))
			return(list("power" = power))
		}
	}
}


int.hazc = function(landmark, delta, x)     
  {

    ######### sort the failure time ######
	n=length(x)
	id=order(x)
	x=x[id]
	delta=delta[id]
	#####################################

	########## estimate the hazard function ##
	hazardc=rep(0, n)
	wx=rep(0, n)
	for(hh in c(1:n)){
	hazardc[hh]=(1-delta[hh])/sum(x>=x[hh])
	wx[hh]=mean(x>=x[hh])
	}
	return(sum(hazardc/wx*(x<=landmark)))
    
  }

int.hazc.plan = function(landmark, delta, x, cens.rate, B=10000)

{
	tt=seq(0, landmark, length=B)
	survt=rep(0, B)
	S.KM = survfit(Surv(x,delta)~1)
	for(hh in c(1:B)){
		survt[hh]=min(c(1,S.KM$surv)[c(0, S.KM$time)<=tt[hh]])}
	survc=exp(-tt*cens.rate)

	return(landmark*mean(cens.rate/survc/survt))
	}
	

var.delta.eb = function(Axzero, Adeltazero, Aszero, Bxone, Bdeltaone, Bsone, Bxzero, Bdeltazero, Bszero, t, landmark=landmark, extrapolate) {
	weightA = rep(1,length(Axzero))
	weightB = rep(1,length(Bxone)+length(Bxzero))
	weightB.group1 = weightB[1:length(Bxone)]
	weightB.group0 = weightB[(1+length(Bxone)):(length(Bxone)+length(Bxzero))]
	mu.s = pred.smooth.surv.new(Axzero.f=Axzero[Axzero>landmark], Adeltazero.f=Adeltazero[Axzero>landmark], Aszero.f=Aszero[Axzero>landmark], Bsnew.f=Bsone[Bxone>landmark], Bsnew2.f = Bszero[Bxzero>landmark], myt=t, weight.pred = weightA[Axzero>landmark], extrapolate = extrapolate)
	censorB1.landmark = censor.weight(Bxone, Bdeltaone, landmark, weight = weightB.group1)
	censorB0.landmark = censor.weight(Bxzero, Bdeltazero, landmark, weight = weightB.group0)
	n1 = length(Bxone)
	n0 = length(Bxzero)
	n.total = n1+n0
	pi.0 = n0/n.total
	pi.1 = n1/n.total
	mu0b2 = (1/n0)*sum((mu.s$Phat.ss.2^2)/censorB0.landmark) 
	mu0b1 = (1/n0)*sum(mu.s$Phat.ss.2/censorB0.landmark) 
	term0 = (1/pi.0)*(mu0b2/censorB0.landmark - mu0b1^2*(1+int.hazc(landmark = landmark, delta = Bdeltazero, x = Bxzero)))
	mu1b2 = (1/n1)*sum((mu.s$Phat.ss.1^2)/censorB1.landmark) 
	mu1b1 = (1/n1)*sum(mu.s$Phat.ss.1/censorB1.landmark) 
	term1 = (1/pi.1)*(mu1b2/censorB1.landmark - mu1b1^2*(1+int.hazc(landmark = landmark, delta = Bdeltaone, x = Bxone)))
	var.closed = term0+term1	
	return(var.closed)
}
