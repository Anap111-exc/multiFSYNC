# Isolated NC experiment I/O; no model updates are defined here.
options(warn=1,stringsAsFactors=FALSE)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
           NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',BLIS_NUM_THREADS='1')
assert<-function(ok,message) if(!isTRUE(ok))stop(message,call.=FALSE)
`%||%`<-function(x,y)if(is.null(x))y else x
sha<-function(path)digest::digest(file=path,algo='sha256',serialize=FALSE)
check_hash<-function(path,expected)assert(file.exists(path)&&identical(tolower(sha(path)),tolower(expected)),paste('BINDING_CONFLICT',path))
read_csv<-function(path)read.csv(path,check.names=FALSE)
new_rds<-function(x,path){
 assert(!file.exists(path)&&!file.exists(paste0(path,'.part')),paste('Refusing overwrite',path))
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE,mode='0700')
 tmp<-paste0(path,'.part');saveRDS(x,tmp,version=3L)
 check<-readRDS(tmp);assert(identical(x,check),paste('Serialization roundtrip failed',path));rm(check)
 assert(file.rename(tmp,path),paste('Atomic rename failed',path))
}
new_csv<-function(x,path){
 assert(!file.exists(path)&&!file.exists(paste0(path,'.part')),paste('Refusing overwrite',path))
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 tmp<-paste0(path,'.part');write.csv(x,tmp,row.names=FALSE,na='')
 assert(nrow(read_csv(tmp))==nrow(x),'CSV row-count mismatch');assert(file.rename(tmp,path),'CSV atomic rename failed')
}
new_lines<-function(x,path){
 assert(!file.exists(path)&&!file.exists(paste0(path,'.part')),paste('Refusing overwrite',path))
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 tmp<-paste0(path,'.part');writeLines(x,tmp,useBytes=TRUE);assert(file.rename(tmp,path),'Text atomic rename failed')
}
save_new<-new_rds;csv_new<-new_csv;lines_new<-new_lines
extract_functions<-function(path,target,envir=.GlobalEnv){
 found<-character()
 for(expr in parse(path,keep.source=FALSE)){
  if(!is.call(expr)||!identical(expr[[1]],as.name('<-'))||!is.symbol(expr[[2]])||!is.call(expr[[3]])||!identical(expr[[3]][[1]],as.name('function')))next
  nm<-as.character(expr[[2]]);if(nm%in%target){eval(expr,envir);found<-c(found,nm)}
 }
 assert(setequal(found,target)&&!anyDuplicated(found),'Function extraction failed')
}
rbind_fill<-function(rows){
 rows<-Filter(function(x)is.data.frame(x)&&nrow(x)>0L,rows)
 if(!length(rows))return(data.frame())
 cols<-unique(unlist(lapply(rows,names)))
 do.call(rbind,lapply(rows,function(x){for(nm in setdiff(cols,names(x)))x[[nm]]<-NA;x[,cols,drop=FALSE]}))
}
all_finite<-function(x){if(is.numeric(x))return(all(is.finite(x)));if(is.list(x))return(all(vapply(x,all_finite,logical(1))));TRUE}
