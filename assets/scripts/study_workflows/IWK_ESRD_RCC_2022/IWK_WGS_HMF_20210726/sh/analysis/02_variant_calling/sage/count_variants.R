library(data.table)

load("/Users/toddjohnson/eclipse-workspace/RemoteSystemsTempFiles/SLOGIN.HGC.JP/yshare2/ZETTAI_path_WA_slash_home_KARA/home/tjohnson/workspace/runs/IWK_WGS_HMF_20210726/result_summaries/sage_variant_snpEff_summaries.Rdata")

germline.genomewide.variant.type.counts[,Total:=gsub(pattern=',', replacement='', Total, fixed=TRUE)]
germline.genomewide.variant.type.counts[,Total:=as.numeric(Total)]
germline.variant.cts <- germline.genomewide.variant.type.counts[,list(total=sum(Total)), by=list(subject.id)]

mean(germline.variant.cts$total)
[1] 4496808
> sd(germline.variant.cts$total)
[1] 100430.7