# Post-fit uncertainty only; no fitting or FPCA recomputation.
band_assert<-function(x,msg)if(!isTRUE(x))stop(msg,call.=FALSE)
band_nonnegative<-function(x,label) {
 band_assert(all(is.finite(x)),paste('Nonfinite variance:',label))
 band_assert(min(x)>=-1e-9*max(1,max(abs(x))),paste('Negative variance:',label))
 pmax(x,0)
}
band_quadratic<-function(C,V)band_nonnegative(rowSums((C%*%V)*C),'quadratic covariance')
band_many_cov<-function(C,covs) {
 k<-ncol(C);Q<-do.call(cbind,lapply(seq_len(k),function(b)C*C[,b]))
 band_nonnegative(Q%*%do.call(cbind,lapply(covs,as.vector)),'coefficient variances')
}
band_loading_moments<-function(f,block,s=1) {
 if(block=='shared'){a<-f$mu_q_a_original;second<-sweep(f$term_a,1,f$response_scale^2,'*');mu<-f$mu_q_normal_a;vv<-f$Sigma_q_normal_a;pp<-f$mu_q_gamma_a}
 else {a<-f$mu_q_b_specific_original[[s]];second<-sweep(f$term_b_specific[[s]],1,f$response_scale^2,'*');mu<-f$mu_q_normal_b[[s]];vv<-f$Sigma_q_normal_b[[s]];pp<-f$mu_q_gamma_b[[s]]}
 band_assert(max(abs(second-sweep(pp*(vv+mu^2),1,f$response_scale^2,'*')))<1e-8,'Loading second-moment binding')
 band_assert(max(abs(a-sweep(pp*mu,1,f$response_scale,'*')))<1e-8,'Loading mean binding')
 list(mean=a,second=second,variance=band_nonnegative(second-a^2,'effective loading'),
      slab_mean=sweep(mu,1,f$response_scale,'*'),slab_variance=sweep(vv,1,f$response_scale^2,'*'),ppi=pp)
}
band_cache<-function(pack,t) {
 f<-pack$estimate;C<-stage7m_C(t,pack$basis)
 list(t=t,C=C,
      mean=lapply(f$mu_q_nu_mu_original,function(x)C%*%do.call(cbind,x)),
      mean_var=lapply(f$Sigma_q_nu_mu_original,function(x)band_many_cov(C,x)),
      beta=lapply(seq_len(f$d),function(r)C%*%do.call(cbind,lapply(f$mu_q_nu_beta_original,`[[`,r))),
      beta_var=lapply(seq_len(f$d),function(r)band_many_cov(C,lapply(f$Sigma_q_nu_beta_original,`[[`,r))),
      shared_u=lapply(f$mu_q_nu_phi,function(x)C%*%x),
      shared_v=lapply(f$Sigma_q_nu_phi,function(x)band_many_cov(C,x)),
      specific_u=lapply(f$mu_q_nu_psi,function(x)lapply(x,function(y)C%*%y)),
      specific_v=lapply(f$Sigma_q_nu_psi,function(x)lapply(x,function(y)band_many_cov(C,y))),
      shared_load=band_loading_moments(f,'shared'),
      specific_load=lapply(seq_len(f$S),function(s)band_loading_moments(f,'specific',s)))
}
band_process<-function(U,V,m,Sigma) {
 m<-as.numeric(m);mean<-as.vector(U%*%m)
 score_var<-band_quadratic(U,Sigma)
 time_var<-as.vector(V%*%(diag(Sigma)+m^2))
 list(mean=mean,variance=band_nonnegative(score_var+time_var,'factor process'),
      score_var=score_var,time_var=time_var)
}
band_subject<-function(pack,cache,s,i) {
 f<-pack$estimate;nr<-length(cache$t);p<-f$p;zero<-matrix(0,nr,p)
 shared<-specific<-shared_var<-specific_var<-zero
 decomp<-list(mean=cache$mean_var[[s]],beta=zero,shared_score=zero,shared_time=zero,
              shared_loading=zero,specific_score=zero,specific_time=zero,specific_loading=zero)
 beta<-zero
 for(r in seq_len(f$d)){z<-pack$observation$Z[[s]][i,r];beta<-beta+z*cache$beta[[r]];decomp$beta<-decomp$beta+z^2*cache$beta_var[[r]]}
 for(l in seq_along(cache$shared_u)) {
  a<-band_process(cache$shared_u[[l]],cache$shared_v[[l]],f$mu_q_zeta[[s]][[l]][i,],f$Sigma_q_zeta[[s]][[l]][[i]])
  lm<-cache$shared_load
  shared<-shared+outer(a$mean,lm$mean[,l])
  decomp$shared_score<-decomp$shared_score+outer(a$score_var,lm$second[,l])
  decomp$shared_time<-decomp$shared_time+outer(a$time_var,lm$second[,l])
  decomp$shared_loading<-decomp$shared_loading+outer(a$mean^2,lm$variance[,l])
 }
 for(l in seq_along(cache$specific_u[[s]])) {
  a<-band_process(cache$specific_u[[s]][[l]],cache$specific_v[[s]][[l]],f$mu_q_xi[[s]][[l]][i,],f$Sigma_q_xi[[s]][[l]][[i]])
  lm<-cache$specific_load[[s]]
  specific<-specific+outer(a$mean,lm$mean[,l])
  decomp$specific_score<-decomp$specific_score+outer(a$score_var,lm$second[,l])
  decomp$specific_time<-decomp$specific_time+outer(a$time_var,lm$second[,l])
  decomp$specific_loading<-decomp$specific_loading+outer(a$mean^2,lm$variance[,l])
 }
 shared_var<-decomp$shared_score+decomp$shared_time+decomp$shared_loading
 specific_var<-decomp$specific_score+decomp$specific_time+decomp$specific_loading
 mean<-list(signal=cache$mean[[s]]+beta+shared+specific,factor=shared+specific,
            shared_contribution=shared,specific_contribution=specific,mean=cache$mean[[s]],beta_contribution=beta)
 variance<-list(signal=cache$mean_var[[s]]+decomp$beta+shared_var+specific_var,
                 factor=shared_var+specific_var,shared_contribution=shared_var,
                 specific_contribution=specific_var,mean=cache$mean_var[[s]],beta_contribution=decomp$beta)
 variance<-lapply(variance,band_nonnegative,label='summed process')
 noise<-matrix(f$sigsq_eps_hat[s,],nr,p,byrow=TRUE)
 band_assert(all(is.finite(noise))&&all(noise>0),'Invalid stored E_q(noise variance)')
 variance$response<-variance$signal+noise;mean$response<-mean$signal;decomp$noise<-noise
 list(mean=mean,variance=variance,decomposition=decomp)
}
band_spike_quantile<-function(prob,ppi,mean,sd) {
 band_assert(all(ppi>=0&ppi<=1)&all(sd>=0),'Invalid spike/slab distribution')
 out<-rep(0,length(ppi));nz<-ppi>0&sd>0
 lo<-ppi*pnorm(0,mean,sd);hi<-1-ppi+lo
 left<-nz&prob<lo;right<-nz&prob>hi
 out[left]<-qnorm(prob/ppi[left],mean[left],sd[left])
 out[right]<-qnorm((prob-(1-ppi[right]))/ppi[right],mean[right],sd[right])
 # A degenerate slab is a two-point distribution.
 at<-ppi>0&sd==0
 out[at&mean<0&prob<ppi]<-mean[at&mean<0&prob<ppi]
 out[at&mean>0&prob>1-ppi]<-mean[at&mean>0&prob>1-ppi]
 out
}
band_summary<-function(estimate,variance,target,w,z,cols=seq_len(ncol(estimate))) {
 estimate<-as.matrix(estimate)[,cols,drop=FALSE];target<-as.matrix(target)[,cols,drop=FALSE]
 variance<-as.matrix(variance)[,cols,drop=FALSE]
 band_assert(identical(dim(estimate),dim(target))&&identical(dim(estimate),dim(variance)),'Band metric dimensions')
 if(!length(cols)||!length(w)||sum(w)<=0)return(c(coverage=NA,width=NA,rmse=NA,all_grid_covered=NA,standardized_squared_error=NA,points=0))
 err<-estimate-target;sd<-sqrt(variance);cover<-abs(err)<=z*sd
 ww<-w/sum(w);coverage<-mean(colSums(ww*cover));width<-mean(colSums(ww*(2*z*sd)))
 nz<-variance>1e-20;std<-if(all(nz))mean(colSums(ww*(err^2/variance)))else NA_real_
 c(coverage=coverage,width=width,rmse=sqrt(mean(colSums(ww*err^2))),
   all_grid_covered=mean(colSums(!cover)==0),standardized_squared_error=std,points=length(estimate))
}
