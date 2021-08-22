#####################################################################
#
# build_patchseq_reference.R
#
#####################################################################

# Load libraries
library(dendextend)
library(dplyr)
library(feather)
library(scrattch.hicat)
library(scrattch.io)
library(matrixStats)
library(patchseqtools)  # devtools::install_github("AllenInstitute/patchseqtools")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/Utils.R")
options(stringsAsFactors = FALSE)
options(error=recover)


#refFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/changkyul/Mapping_PatchSeq/Ref_test/mouse_patchseq_VISp_20210412_collapsed40_cpm/"
#completeDendFile = "dend.RData"
#refDendFile = "dend_ref.RData"

# set input arguments
print(paste0("refFolder = ", refFolder))
print(paste0("completeDendFile = ", completeDendFile))
print (paste0("refDendFile = ", refDendFile))

if (TRUE ||
    !file.exists(file.path(refFolder, "reference.rda")) ||
    !file.exists(file.path(refFolder, "membership_information_reference.rda")) ||
    !file.exists(file.path(refFolder, "keep_ref.KL.mapping.rda")) ) {

   if (species=="mouse" && taxonomy_str=="ForeBrain") {
      ##working directory
      d = workFolder #"/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/foreBrain/"
      ##Shiny
      anno = read_tome_anno("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/tomes/facs/mouse_ForeBrain_20201204/smrt.tome")
      load(file.path(workFolder,"dend.rda"))
      load(file.path(workFolder,"cl.final.rda"))
      file.copy(file.path(workFolder, "cl.final.rda"), file.path(refFolder, "cl.final.rda"))
   
      ###Decide which clusters should be included in the reference
      select.roi = c(unique(grep("TH", anno$roi_label,value=T)), unique(grep("STR[v|d]", anno$roi_label,value=TRUE)), unique(grep("PAL_GP",anno$roi_label,value=TRUE)))
      tb=with(anno, table(cluster_label, roi_label))
      tb.frac= tb/rowSums(tb)
      select.counts = rowSums(tb[,select.roi])
      select.ratio=rowSums(tb.frac[,select.roi])
      select.cl_label = intersect(names(select.ratio)[select.ratio > 0.1 & select.counts>=4], cl.df[cl.df$class_label %in% c("Ex","Inh","Hybrid"),"cluster_label"])
      select.cl = droplevels(cl[cl %in% (cl.df %>% filter(cluster_label %in% select.cl_label) %>% pull(cl))])
      select.sample_name = names(select.cl)
      #select2.cl = select.cl
   
      ###Prune dendrogram with only selected clusters
      select.dend = prune_dend(dend, setdiff(labels(dend), select.cl))
      select.dend = reorder_dend(select.dend, with(cl.df[labels(select.dend),],setNames(cluster_id, cl)))
      labeled.dend = select.dend
   
      dend_ref = select.dend
      cluster_id.dend_ref = labels(dend_ref)
      cluster_label.dend_ref = cl.df[match(cluster_id.dend_ref, cl.df$cl), "cluster_label"]
      names(cluster_label.dend_ref) = cluster_id.dend_ref
       
      labels(dend_ref) = cluster_label.dend_ref
      dend = dend_ref

      load(file.path(d, "norm.dat.rda"))
      norm.dat = norm.dat[, match(select.sample_name, colnames(norm.dat))]
      thisdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_ForeBrain_20201203"
      save(norm.dat, file=file.path(thisdir, "norm.sel.rda"))
      print("save norm.sel.rda for mapping")

   } else {
      ##
      # Read/Load in anno, data, dendrogram for reference
      ##
      anno_ref = feather(file.path(refFolder,"anno.feather"))
      dat_ref  = feather(file.path(refFolder,"data.feather"))
      dat_ref  = as.data.frame(dat_ref)
      if (file.exists(file.path(refFolder, "cl.df.rda"))) { 
         load(file.path(refFolder, "cl.df.rda")) 
      } else {
         cl.df = anno_ref[match(labels(dend), anno_ref$cluster_label), ]
      }
      if (length(rois)> 0) {
         sel_idx = which(anno_ref$Region %in% rois)
         anno_ref = anno_ref[sel_idx,]
         dat_ref  = dat_ref[sel_idx,]
      }
      if ("cl" %in% colnames(anno_ref)) cl = as.factor(anno_ref$cl)
      if ("cluster_id" %in% colnames(anno_ref)) cl = as.factor(anno_ref$cluster_id)
      if ("sample_id" %in% colnames(anno_ref)) names(cl) = anno_ref$sample_id
      if ("exp_component_name" %in% colnames(anno_ref)) names(cl) = anno_ref$exp_component_name

      a = load(file.path(refFolder,refDendFile)) 
      eval(parse(text=paste("dend <-",a)))
      
      # maker sure sample_id is used
      colnames(dat_ref) <- gsub("sample_name","sample_id",colnames(dat_ref))
      
      print("Define normalized data matrix")
      counts   = as.matrix(dat_ref[,colnames(dat_ref)!="sample_id"])
      
      counts   = t(counts)  
      colnames(counts)  <- dat_ref$sample_id
      norm.dat = logCPM(counts) 
      rm(counts)
      colnames(norm.dat) <- dat_ref$sample_id
      rm(dat_ref)
       
      ###Decide which clusters should be included in the reference
      #if ("roi_label" %in% colnames(anno)) {
      #   tb=with(anno, table(cluster_label, roi_label))
      #   tb.frac= tb/rowSums(tb)
      #   select.counts = rowSums(tb)
      #   select.ratio=rowSums(tb.frac)
      #   select.cl_label = names(select.ratio)[select.ratio > 0.1 & select.counts>=4]#, cl.df[cl.df$class_label %in% c("Ex","Inh","Hybrid"),"cluster_label"])
      #   select.cl = droplevels(cl[cl %in% (cl.df %>% filter(cluster_label %in% select.cl_label) %>% pull(cl))])
      #   select.sample_name = names(select.cl)
      #} else {
         select.counts = table(anno$cluster_label)
         select.cl_label = names(select.counts)[select.counts>=4]#, cl.df[cl.df$class_label %in% c("Ex","Inh","Hybrid"),"cluster_label"])
         select.cl = droplevels(cl[cl %in% (cl.df %>% filter(cluster_label %in% select.cl_label) %>% pull(cl))])
         select.sample_name = names(select.cl)
      #}
      print("###Prune dendrogram with only selected clusters")
      select.dend = prune_dend(dend, setdiff(labels(dend), select.cl_label))
      select.dend = reorder_dend(select.dend, with(cl.df[match(labels(select.dend),cl.df$cluster_label),],setNames(cluster_id, cl)))
      labeled.dend = select.dend
   
      dend_ref = select.dend
      cluster_label.dend_ref = labels(dend_ref)
      cluster_id.dend_ref = cl.df[match(cluster_label.dend_ref, cl.df$cluster_label), "cluster_id"]
      cluster_cl.dend_ref = cl.df[match(cluster_label.dend_ref, cl.df$cluster_label), "cl"]
      names(cluster_label.dend_ref) = cluster_id.dend_ref
      names(cluster_cl.dend_ref) = cluster_id.dend_ref
       
      dend = dend_ref
      labels(dend) = cluster_id.dend_ref

   }
   
   
   
   # Subsample prior to feeding into display_cl to avoid function crashing
   if (TRUE || !file.exists(file.path(refFolder, "norm.dat2.rda"))) {
      max.cl.size <- 100
      kpSamp     <- sample_cells(select.cl,max.cl.size)
      norm.dat2  <- norm.dat[,kpSamp]
      select2.cl <- select.cl[kpSamp]
      levels(select2.cl) <- labels(dend)
      save(norm.dat2, select2.cl, max.cl.size, file=file.path(refFolder, "norm.dat2.rda"))
   } else {
      load(file.path(refFolder, "norm.dat2.rda"))
   }
      
   if (FALSE && file.exists(file.path(taxoFolder_FACS, "de.genes.rda"))) {
      file.copy(file.path(taxoFolder_FACS, "de.genes.rda"), file.path(refFolder, "de.genes.rda"))
      load(file.path(taxoFolder_FACS, "de.genes.rda"))
   } else {
      print("calculating de.genes")
      q1thr = c(0.4, 0.3, 0.0)
      pthr  = c(0.01, 0.01, 1)
      dethr = c(120, 100, 0)
      qdiffthr = c(0.7, 0.2, 0.1)
      flag.0DE = TRUE
      t=0
      #select.cl = droplevels(select.cl)
      while (flag.0DE && t<3) {
         t=t+1
         print(paste0("!!! p=", pthr[t], ", q1=", q1thr[t], ", de.thr=", dethr[t], ", q.diff.thr=", qdiffthr[t]))
         de.genes = display_cl(cl=select2.cl, norm.dat=norm.dat2,  max.cl.size=max.cl.size, mc.cores=10, 
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
         if (flag.0DE) {
            print("!!! there are 0 DEG in pariwise comparison. finding DEG with relaxed de.param criteria :")
         }
      }
      save(de.genes, file=file.path(refFolder, "de.genes.rda"))
   }
   
   
   print("Define gene scores")
   if (FALSE && file.exists(file.path(taxoFolder_FACS, "gene.score.rda"))) {
      file.copy(file.path(taxoFolder_FACS, "gene.score.rda"), file.path(refFolder, "gene.score.rda"))
      load(file.path(taxoFolder_FACS, "gene.score.rda"))
   } else {
      if (!exists("de.genes")) load(file.path(refFolder, "de.genes.rda"))
      gene.score = get_gene_score(de.genes)
      save(gene.score, file=file.path(refFolder, "gene.score.rda"))
   }   
   
   print("Define select.markerss")
   if (TRUE || !file.exists(file.path(refFolder, "select.markers.rda"))) {
      labels(dend_ref) = cluster_id.dend_ref
      #select.markers = select_markers( norm.dat, cl=select.cl, n.markers = 20, de.genes = de.genes)
      select.markers = select_markers( norm.dat2, cl=select2.cl, n.markers = 20, de.genes = de.genes)
      save(select.markers, file=file.path(refFolder, "select.markers.rda"))
   } else {
      load(file.path(refFolder, "select.markers.rda"))
   }
   
   print("Now, Ready to Build the reference")
   
   ##u# Build the Reference
   ##
   #print("!!!! In Buliding Reference, if you get this error message ")
   #print("------------------------------------------------------------------------------------------")
   #print(" Error in select_dend_markers(dend, norm.dat = norm.dat, cl = cl, de.genes = de.genes,  : ")
   #print(" no loop for break/next, jumping to top level ")
   #print("------------------------------------------------------------------------------------------")
   #print("!!!! Please relax thresholds in de.param() to detect DEG in ")
   #print("!!!! pairwise cluster comparison and rerun")
   
   if (FALSE && file.exists(file.path(refFolder, "reference.rda"))) {
      load(file.path(refFolder, "reference.rda")) 
   
   } else {
      #reference = build_reference(cl=select2.cl, norm.dat=norm.dat2, dend=dend,de.genes=de.genes, 
      #reference = build_reference(cl=select.cl, norm.dat=norm.dat, dend=dend,de.genes=de.genes, 
      # make sure select.cl, labels(dends) are consistent.
      if (any(is.na(match(unique(select2.cl),labels(dend))))) labels(dend) = names(labels(dend)) 
      reference = build_reference(cl=select2.cl, norm.dat=norm.dat2, dend=dend,de.genes=de.genes, 
                               up.gene.score=gene.score$up.gene.score, cl.label=NULL,
                               down.gene.score=gene.score$down.gene.score, n.markers=20)
   
      print("before")
      mylabel0 = labels(reference$dend)
      labels(reference$dend) <- setNames(colnames(reference$cl.dat),colnames(reference$cl.dat))
      #save(reference, file="reference.rda")
      print("after")
      mylabel1 = labels(reference$dend)
      cbind(mylabel0, mylabel1)
         
      print("This is needed if the starting dendrogram is from the nomenclature GitHub ")
      revert_dend_label <- function(dend, value, attribute="label") {
         if(attr(dend, attribute)=="")
             attr(dend, attribute) <- value[attr(dend,"original_label")]
         if (length(dend)>1) for(i in 1:length(dend))
             dend[[i]]=revert_dend_label(dend[[i]], value=value, attribute)
         return(dend)
      }
      if(sum(!is.na(get_nodes_attr(reference$dend, "original_label"))>0)){
         reference$dend <- revert_dend_label(reference$dend,get_nodes_attr(reference$dend, "original_label"),"label")
      }
      
      #
      # update dend_ref and cluster_label.dend_ref in keep_ref.KL.mapping.rda
      #
      
      dend_ref  = reference$dend
      cl.dat = reference$cl.dat
   
      cluster_id.dend_ref = labels(dend_ref)
      cluster_label.dend_ref = as.character(cl.df[match(cluster_id.dend_ref, cl.df$cl), "cluster_label"])
      labels(dend_ref) = cluster_label.dend_ref
      labels(reference$dend) = cluster_label.dend_ref
   
      reference$cl.dat = cl.dat[, cluster_id.dend_ref]
      colnames(reference$cl.dat) = cluster_label.dend_ref
   
      save(reference, file=file.path(refFolder, "reference.rda"))
   }
   
   if (!file.exists(file.path(refFolder, "membership_information_reference.rda"))) {
      ##
      # Membership table of reference
      ##
      print("Build membership table of reference vs. reference for use with patch-seq mapping")
      memb.ref  = map_dend_membership( dend=reference$dend, 
                                       cl.dat=reference$cl.dat, 
                                       map.dat=norm.dat2,
                                       map.cells=names(select2.cl),
                                       mc.cores=10, bs.num=100, p=0.75, low.th=0.15 )
     #                                  mc.cores=1, bs.num=3, p=0.75, low.th=0.15 )
      save( memb.ref, file=file.path(refFolder, "membership_information_reference.rda"))
      print("now summarize it ")
      map.df.ref = summarize_cl(reference$dend, memb.ref, norm.dat2)
       
      my.cl = cluster_label.dend_ref[select.cl]
      names(my.cl) = names(select.cl)
      myselect <- cluster_label.dend_ref
      memb.ref = memb.ref[, labels(reference$dend)]
      Tree_Map_Prob = compute_map_prob( memb = memb.ref,
                                  ref.cl= my.cl,
                                  select.cells =colnames(norm.dat2),
                                  select.cl = myselect)
   
   
      save( memb.ref, map.df.ref, Tree_Map_Prob, 
            file=file.path(refFolder, "membership_information_reference.rda"))
   }
   
   save( cl.df, select2.cl, select.sample_name, dend_ref, select.markers,
         cluster_label.dend_ref, cluster_id.dend_ref, 
         file=file.path(refFolder, "keep_ref.KL.mapping.rda") )
   
} else {
   print("reference are already built!")
}

browser()




