# Isolated protocol I/O. No fitting, response or truth values loaded on source.
iv_assert <- function(ok,msg) if(!isTRUE(ok))stop(msg,call.=FALSE)
iv_sha <- function(path)digest::digest(file=path,algo='sha256',serialize=FALSE)
iv_hash <- function(x)digest::digest(x,algo='sha256',serialize=TRUE)
iv_check <- function(path,expected)iv_assert(file.exists(path)&&identical(iv_sha(path),tolower(expected)),paste('SHA mismatch:',path))
iv_json <- function(x,path) {
  iv_assert(!file.exists(path)&&!file.exists(paste0(path,'.part')),'Refuse overwrite')
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  jsonlite::write_json(x,paste0(path,'.part'),auto_unbox=TRUE,pretty=TRUE,null='null',na='null',digits=17)
  iv_assert(file.rename(paste0(path,'.part'),path),'Atomic JSON rename failed')
}
iv_csv <- function(x,path) {
  iv_assert(!file.exists(path),'Refuse CSV overwrite');dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  write.csv(x,path,row.names=FALSE,na='')
}
iv_rds <- function(x,path) {
  iv_assert(!file.exists(path)&&!file.exists(paste0(path,'.part')),'Refuse RDS overwrite')
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE,mode='0700')
  saveRDS(x,paste0(path,'.part'),version=3,compress='gzip')
  iv_assert(identical(x,readRDS(paste0(path,'.part'))),'RDS roundtrip failed')
  iv_assert(file.rename(paste0(path,'.part'),path),'Atomic RDS rename failed')
}
iv_verify_package <- function(ex) {
  p<-file.path(ex,'CODE_SHA256SUMS.txt');lines<-readLines(p,warn=FALSE)
  iv_assert(length(lines)>0L&&all(grepl('^[0-9a-f]{64}  .+',lines)),'Invalid source manifest')
  rel<-substring(lines,67L)
  iv_assert(!anyDuplicated(rel)&&!any(grepl('(^/|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$))',rel)),'Unsafe package manifest')
  for(i in seq_along(lines))iv_check(file.path(ex,rel[i]),substr(lines[i],1,64))
  iv_sha(p)
}
iv_data_binding <- function(out,data_id,role) {
  b<-read.csv(file.path(out,'DATA_BINDINGS.csv'),stringsAsFactors=FALSE,check.names=FALSE)
  b<-b[b$data_id==data_id&b$role==role,,drop=FALSE]
  iv_assert(nrow(b)==1L,'Missing/duplicate data binding')
  iv_assert(!grepl('(^/|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$))',b$relative_path),'Unsafe data relative path')
  b$path<-file.path(out,b$relative_path);b
}
iv_read_binding <- function(binding) {iv_check(binding$path,binding$sha256);readRDS(binding$path)}
iv_threads <- function() {
  for(k in c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'))
    iv_assert(identical(Sys.getenv(k),'1'),paste('Numerical threads must be 1:',k))
}
iv_rng <- function()list(kind=RNGkind(),present=exists('.Random.seed',.GlobalEnv,inherits=FALSE),
  value=if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))get('.Random.seed',.GlobalEnv)else NULL)
iv_functions <- function(path,wanted,env) {
  found<-character()
  for(e in parse(path,keep.source=FALSE))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&
    is.symbol(e[[2]])&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function'))&&as.character(e[[2]])%in%wanted) {
    eval(e,env);found<-c(found,as.character(e[[2]]))
  }
  iv_assert(setequal(found,wanted)&&!anyDuplicated(found),'Frozen function extraction anchors differ')
}
iv_splinets <- function(ex) {
  dep<-jsonlite::read_json(file.path(ex,'new_subjects/frozen/DEPENDENCY_BINDING.json'),simplifyVector=FALSE)
  iv_assert(Sys.info()[['sysname']]=='Linux'&&as.character(getRversion())=='4.5.3','Registered Linux R 4.5.3 required')
  .libPaths(unique(c(dirname(dep$package_root),dep$R_LIBS_USER,.libPaths())))
  iv_assert(as.character(utils::packageVersion('Splinets'))=='1.5.1','Frozen Splinets version differs')
  iv_assert(identical(normalizePath(getNamespaceInfo(asNamespace('Splinets'),'path'),mustWork=TRUE),
    normalizePath(dep$package_root,mustWork=TRUE)),'Frozen Splinets path differs')
  for(entry in dep$package_files)iv_check(file.path(dep$package_root,entry$relative_path),entry$sha256)
  invisible(dep)
}
