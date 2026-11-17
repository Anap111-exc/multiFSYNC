# Reconstruction only. No fitting, truth-dependent scaling, or cloud access.
stage7m_weights<-function(t){dt<-diff(t);c(dt[1]/2,(head(dt,-1)+tail(dt,-1))/2,tail(dt,1)/2)}
stage7m_basis<-function(time_obs,K) {
 flat<-unlist(time_obs,recursive=FALSE)
 knots<-unname(quantile(sort(unique(unlist(flat))),seq(0,1,length.out=K)[-c(1,K)]))
 allKnots<-c(rep(0,4),knots,rep(1,4));nk<-length(knots);nn<-3*(nk+8)
 xtilde<-(rep(allKnots,each=3)[-c(1,nn-1,nn)]+rep(allKnots,each=3)[-c(1,2,nn)])/2
 wts<-rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,nk+7)
 Bdd<-splines::spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),outer.ok=TRUE)$design
 Omega<-crossprod(Bdd*wts,Bdd);sv<-svd(Omega);ix<-seq_len(nk+2)
 LZ<-t(t(sv$u[,ix,drop=FALSE])/sqrt(sv$d[ix]))
 list(schema='stage7m_train_only_osullivan_basis_v1',K=K,internal_knots=knots,boundary=c(0,1),degree=3L,all_knots=allKnots,penalty=Omega,transform=LZ,training_time_sha256=digest::digest(time_obs,algo='sha256'),training_grid=seq(0,1,length.out=51))
}
stage7m_C<-function(t,basis) {
 B<-splines::spline.des(basis$all_knots,as.numeric(t),outer.ok=TRUE)$design
 cbind(1,as.numeric(t),B%*%basis$transform)
}
stage7m_reported_coefficients<-function(fit) {
 project<-function(x) {
  coef<-qr.solve(fit$C_g,x);err<-max(abs(fit$C_g%*%coef-x))
  assert(err<1e-9,'FPCA output not in fitted basis span');coef
 }
 list(shared=lapply(fit$list_Phi_hat,project),specific=lapply(fit$list_Phi_hat_spec,function(x)lapply(x,project)))
}
stage7m_estimate<-function(fit,basis,s,i,t,representation='raw',reported=NULL) {
 C<-stage7m_C(t,basis);p<-fit$p;d<-fit$d
 mu<-C%*%do.call(cbind,fit$mu_q_nu_mu_original[[s]])
 beta<-lapply(seq_len(d),function(r)C%*%do.call(cbind,lapply(fit$mu_q_nu_beta_original,`[[`,r)))
 if(representation=='raw') {
  shared<-lapply(seq_len(fit$L_f),function(l)as.vector(C%*%fit$mu_q_nu_phi[[l]]%*%fit$mu_q_zeta[[s]][[l]][i,]))
  specific<-lapply(seq_along(fit$mu_q_nu_psi[[s]]),function(l)as.vector(C%*%fit$mu_q_nu_psi[[s]][[l]]%*%fit$mu_q_xi[[s]][[l]][i,]))
  A<-fit$mu_q_a_original;B<-fit$mu_q_b_specific_original[[s]]
 } else {
  assert(representation=='reported','Unknown estimate representation')
  if(is.null(reported))reported<-stage7m_reported_coefficients(fit)
  idx<-sum(head(fit$n_s,s-1))+i
  shared<-lapply(seq_along(reported$shared),function(l)as.vector(C%*%reported$shared[[l]]%*%fit$list_Zeta_hat[[l]][idx,]))
  specific<-lapply(seq_along(reported$specific[[s]]),function(l)as.vector(C%*%reported$specific[[s]][[l]]%*%fit$list_Zeta_hat_spec[[s]][[l]][i,]))
  A<-fit$mu_q_a_hat;B<-fit$mu_q_b_specific_hat[[s]]
 }
 shared_contribution<-matrix(0,length(t),p);specific_contribution<-shared_contribution
 for(l in seq_along(shared))shared_contribution<-shared_contribution+outer(shared[[l]],A[,l])
 for(l in seq_along(specific))specific_contribution<-specific_contribution+outer(specific[[l]],B[,l])
 list(mean=mu,beta=beta,shared=shared,specific=specific,shared_contribution=shared_contribution,specific_contribution=specific_contribution,factor=shared_contribution+specific_contribution)
}
stage7m_prediction<-function(fit,obs,basis,test) {
 reported<-stage7m_reported_coefficients(fit)
 predict_rep<-function(rep)lapply(seq_len(fit$S),function(s)lapply(seq_len(fit$n_s[s]),function(i){
  a<-stage7m_estimate(fit,basis,s,i,test$time[[s]][[i]],rep,reported)
  ans<-a$mean+a$factor
  for(r in seq_len(fit$d))ans<-ans+obs$Z[[s]][i,r]*a$beta[[r]]
  ans
 }))
 list(schema='stage7m_blind_predictions_v1',test_time=test$time,raw=predict_rep('raw'),reported=predict_rep('reported'),truth_read=FALSE)
}
stage7m_truth_basis<-function(truth) {
 make_block<-function(meta) {
  system<-Splinets::splinet(meta$knots,degree=meta$degree,periodic=meta$periodic,norm=TRUE)
  co<-meta$stiefel_coefficients
  list(os=system$os,coefficients=co)
 }
 m<-truth$mean_common_meta;os<-Splinets::splinet(m$knots,degree=m$degree,periodic=m$periodic,norm=TRUE)$os
 list(shared=lapply(truth$phi_basis_meta,make_block),specific=lapply(truth$psi_basis_meta,function(x)lapply(x,make_block)),mean_os=os)
}
stage7m_truth_functions<-function(truth,t,truth_basis=NULL) {
 if(is.null(truth_basis))truth_basis<-stage7m_truth_basis(truth)
 values<-function(block)Splinets::evspline(block$os,x=t)[,-1,drop=FALSE]%*%block$coefficients
 shared<-lapply(truth_basis$shared,values);specific<-lapply(truth_basis$specific,function(x)lapply(x,values))
 means<-Splinets::evspline(truth_basis$mean_os,x=t)[,-1,drop=FALSE]
 mean<-lapply(truth$mean_basis_meta,function(x)means%*%x$basis_coefficients)
 beta<-lapply(seq_len(truth$d_use),function(r){
  if(identical(truth$beta_truth_method,'zero_effect_d2'))return(matrix(0,length(t),truth$p))
  B<-splines::bs(t,degree=3,knots=c(.25,.5,.75),Boundary.knots=c(0,1),intercept=TRUE)
  raw<-B%*%truth$beta_bspline_coefficients[,,r]
  ans<-sweep(sweep(raw,2,truth$beta_bspline_raw_mean[,r],'-'),2,truth$beta_bspline_scale[,r],'*')
  sweep(ans,2,truth$beta_bspline_target_mean[,r],'+')
 })
 list(shared=shared,specific=specific,mean=mean,beta=beta)
}
stage7m_true_parts<-function(truth,obs,s,i,t,truth_basis=NULL,functions=NULL) {
 a<-functions%||%stage7m_truth_functions(truth,t,truth_basis)
 shared<-lapply(seq_along(a$shared),function(l)as.vector(a$shared[[l]]%*%truth$zeta_fpca_true[[s]][[l]][i,]))
 specific<-lapply(seq_along(a$specific[[s]]),function(l)as.vector(a$specific[[s]][[l]]%*%truth$xi_fpca_true[[s]][[l]][i,]))
 # zeta_fpca_true / xi_fpca_true are the Scheme-1 variance-bearing true FPCA scores.
 gc<-matrix(0,length(t),truth$p);sc<-gc
 for(l in seq_along(shared))gc<-gc+outer(shared[[l]],truth$a_true[,l])
 for(l in seq_along(specific))sc<-sc+outer(specific[[l]],truth$b_true[[s]][,l])
 bc<-matrix(0,length(t),truth$p)
 for(r in seq_len(truth$d_use))bc<-bc+obs$Z[[s]][i,r]*a$beta[[r]]
 list(mean=a$mean[[s]],beta=bc,shared=shared,specific=specific,shared_contribution=gc,specific_contribution=sc,factor=gc+sc,signal=a$mean[[s]]+bc+gc+sc,functions=a)
}
