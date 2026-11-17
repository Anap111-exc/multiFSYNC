# Projection only: no basis reconstruction, fit update, initialization or RNG.
# Whole retained objects (including Linux knots/transforms/C) are preserved.
project_unseen_fit_inputs<-function(observation,metadata){
  require_fields<-function(x,names,label){
    stopifnot(is.list(x),!anyDuplicated(base::names(x)))
    missing<-setdiff(names,base::names(x));if(length(missing))stop(label,' lacks ',paste(missing,collapse=','),call.=FALSE)
  }
  obs_fields<-c('bundle_class','experiment_id','scenario_id','data_id','Y','Z','time_obs','dimensions',
    'generated_d','truth_available_to_fit','test_response_available_to_fit')
  dim_fields<-c('S','n_s','p','d','K')
  fixed_fields<-c('schema','basis','grid','weights','C','fit_metadata','aggregation','truth_read')
  fitmeta_fields<-c('S','n_s','p','d','K','K_total','C_g','time_g','list_cp_C','list_hyper')
  require_fields(observation,obs_fields,'observation');require_fields(observation$dimensions,dim_fields,'dimensions')
  require_fields(metadata,fixed_fields,'metadata');require_fields(metadata$fit_metadata,fitmeta_fields,'fit metadata')
  stopifnot(is.null(observation$Z),observation$dimensions$d==0L,isFALSE(observation$truth_available_to_fit),isFALSE(metadata$truth_read))
  obs<-observation[obs_fields];obs$dimensions<-observation$dimensions[dim_fields]
  fixed<-metadata[fixed_fields];fixed$fit_metadata<-metadata$fit_metadata[fitmeta_fields]
  stopifnot(identical(obs$Y,observation$Y),identical(obs$Z,observation$Z),identical(obs$time_obs,observation$time_obs),
    identical(fixed$basis,metadata$basis),identical(fixed$C,metadata$C),
    identical(fixed$fit_metadata$C_g,metadata$fit_metadata$C_g),identical(fixed$fit_metadata$time_g,metadata$fit_metadata$time_g),
    !any(c('L_f','L_s','L_s_by_study','M','M_f','M_s')%in%names(obs$dimensions)),
    !any(c('L_f','L_s','L_s_by_study','M','M_f','M_s')%in%names(fixed$fit_metadata)))
  list(observation=obs,metadata=fixed,audit=list(
    observation_fields=obs_fields,dimension_fields=dim_fields,metadata_fields=fixed_fields,fit_metadata_fields=fitmeta_fields,
    removed_observation_fields=setdiff(names(observation),obs_fields),
    removed_dimension_fields=setdiff(names(observation$dimensions),dim_fields),
    removed_metadata_fields=setdiff(names(metadata),fixed_fields),
    removed_fit_metadata_fields=setdiff(names(metadata$fit_metadata),fitmeta_fields),
    observation_values_unchanged=TRUE,basis_knots_transform_C_Cg_time_unchanged=TRUE))
}
