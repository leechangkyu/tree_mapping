
# install scrattch.hicat@dev_zy from github
library(devtools)
devtools::install_github("AllenInstitute/scrattch.hicat", ref="dev_zy")
library(scrattch.hicat)


#
# mapping_cal_atleastone_de.genes(cl, norm.dat)
# = calculate DEG so that there are at least one DEG's in all pairwise comparisions
# 
# - input
# cl       : assigned cluster id
# norm.dat : count matrix,  log2(1 + CPM)
#
mapping_cal_atleastone_de.genes <- function(cl, norm.dat) {

   print("calculating de.genes with set of parameters")
   q1thr = c(0.4, 0.3, 0.0)
   pthr  = c(0.01, 0.01, 1)
   dethr = c(120, 100, 0)
   qdiffthr = c(0.7, 0.2, 0.1)
   flag.0DE = TRUE
   t=0
   while (flag.0DE && t<3) {
      t=t+1
      print(paste0("pthr=", pthr[t], ", q1thr=", q1thr[t], ", de.thr=", dethr[t], ", q.diff.thr=", qdiffthr[t]))
      de.genes = display_cl(cl=cl, norm.dat=norm.dat,  max.cl.size=100, mc.cores=10, 
                            de.param = de_param(q1.th=q1thr[t], de.score.th=dethr[t], 
                                                padj.th=pthr[t], q.diff.th=qdiffthr[t]))$de.genes
      flag.0DE=FALSE
      print("!!! Checking whether any pairwise cluster comparison gives 0 DEG's")
      for (i in 1:length(de.genes)) {
         if (de.genes[[i]]$num < 1) {
            print(paste(i, names(de.genes)[i], de.genes[[i]]$num))
            flag.0DE=TRUE
         }
      }
      if (flag.0DE) print("!!! there are 0 DEG in pariwise comparison. finding DEG with relaxed de.param criteria :")
   }
   return(de.genes)
}


#
# mapping_build_reference (cl, norm.dat, dend, de.genes, n.markers)
# = build the referene
#
# - input
# cl         : Seurat-assigned cluster (numeric preferred, if character, without "_")
# norm.dat   : log2(1+CPM)
# dend       : dendrogram of assigned cluster
# de.genes   : DEG in pairwise cluster comparison
# na.markers : use upto n.markers genes in pairwise DEG
#
mapping_build_reference <- function(cl, norm.dat, dend, de.genes=NULL, n.makers=30) {

   # if de.genes is not provided, 
   # then calculate DEG so that there is at least one DEG in all pairwise comparisions
   if (is.null(de.genes)) mapping_cal_atleastone_de.genes(cl, norm.dat) 
       
   reference = build_reference(cl=cl, norm.dat=norm.dat, dend=dend, de.genes=de.genes, n.markers=n.markers)

   return(reference)
}


#
# mapping_tree_mapping(reference, query.dat.norm, query.dat.cells=colnames(query.dat.norm))
# = do the tree mapping based on reference for query.dat.norm 
# 
# - input
# reference      : reference built by scrattch.hicat::build_reference()
# query.dat.norm : count matrix : log2(1 + CPM)
#
# - output 
# mapping.df
# $cl    : mapped cluster
# $score : membership score       
# other...
#
mapping_tree_mapping <- function(reference, query.dat.norm, query.dat.cells=colnames(query.dat.norm)) {

   # dendrogram & median gene expression for each cluster
   dend=reference$dend
   cl.med=reference$cl.dat

   # tree mapping : 
   # calculate the membership (probability) of each data point belonging to each cluster 
   # (terminal node of dendrogram)
   memb = map_dend_membership(dend, cl.dat=cl.med,
                              map.dat=query.dat.norm, map.cells=query.dat.cells,
                              mc.cores=10, bs.num=100, p=0.7, low.th=0.15)
   # 
   mapping.df = summarize_cl(dend, memb, query.dat.norm)
   return(mapping.df)
}

