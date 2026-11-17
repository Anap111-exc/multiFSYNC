# Evaluation-only fix: exact numeric PPI ties; no decimal label conversion.
stage7m_auc<-function(ppi,support) {
 np<-sum(support);nn<-sum(!support)
 auc<-if(np&&nn)(sum(rank(ppi,ties.method='average')[support])-np*(np+1)/2)/(np*nn)else NA_real_
 ord<-order(-ppi)
 group_id<-match(ppi[ord],unique(ppi[ord]))
 groups<-split(seq_along(ppi),factor(group_id,levels=seq_len(max(group_id))))
 tp<-fp<-0;ap<-0
 for(ix in groups){delta<-sum(support[ord[ix]]);tp<-tp+delta;fp<-fp+length(ix)-delta;if(np)ap<-ap+(delta/np)*tp/(tp+fp)}
 selected<-ppi>=.5;tp<-sum(selected&support);fp<-sum(selected&!support);fn<-sum(!selected&support)
 precision<-if(tp+fp)tp/(tp+fp)else 0;recall<-if(np)tp/np else NA_real_
 data.frame(roc_auc=auc,pr_average_precision=if(np)ap else NA_real_,precision_ppi05=precision,recall_ppi05=recall,f1_ppi05=if(is.finite(recall)&&precision+recall>0)2*precision*recall/(precision+recall)else 0,support_size=np,selected_size=sum(selected),true_positive=tp,false_positive=fp,false_negative=fn)
}

