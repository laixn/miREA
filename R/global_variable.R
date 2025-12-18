default.method_list = c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score", "Edge_ORA", "Edge_Score", "Edge_2Ddist", "Edge_Topology", "Edge_Network")

default.scoreFun = c("none", "rank", "sigmoid", "normalize")

default.pAdjMethod = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none")

default.pathway_list = c("Reactome", "hallmark", "TP")

default.TP_cancer = c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "KIRC",
                      "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD",
                      "READ", "STAD", "THCA", "UCEC")
default.pvalueType = c("pos", "std", "neg")
default.GGI_source = c("Omnipath", "Reactome")

combine_parallel = function(res) {
  obj = data.frame()
  for(i in 1:length(res)) {
    obj = rbind(obj, as.data.frame(res[[i]]))
  }
  return(obj)
}
