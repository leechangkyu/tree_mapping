##
# setup_patchseq_reference.R
##

# Load libraries
library(dendextend)
library(dplyr)
library(feather)
library(scrattch.io)
library(scrattch.hicat)
library(patchseqtools)  # devtools::install_github("AllenInstitute/patchseqtools")
library(stringr)
library(matrixStats)
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/Utils.R")
options(error=recover)
options(stringsAsFactors = FALSE)

# set input arguments
print(paste0("file.type = ", file.type))  
print(paste0("taxoFolder_FACS = ", taxoFolder_FACS))
print(paste0("refFolder = ", refFolder))
print(paste0("qcCellsPerCluster = ", qcCellsPerCluster))
print(paste0("completeDendFile = ", completeDendFile))
print(paste0("refDendFile = ", refDendFile))

if (file.type=="tome") {
   print(paste0("tomeFN = ", tome.FN))
}

## 
# Data clean-up and saving
# -- Make sure label names match between dendrogram and final clusters
# -- Exclude all non-neuronal types from the tree
# -- Subset and correct annotation and data files 
# -----(anno.feather, desc.feather, dend.RData, data.feather, tsne.feather, tsne_desc.feather)
# -- Transfer updated data files to the patch-seq folder
##
print("Read in the data and metadata")

if (file.type=="tome") {
   anno.tome <- read_tome_anno(tome.FN)
   anno <- as.data.frame(anno.tome)
   colnames(anno) <- gsub("sample_name","sample_id",colnames(anno))  # Make sure it's called sample_id
   load(file.path(taxoFolder_FACS, "norm.dat.rda"))

   this.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_ForeBrain_20201203"
   #anno <- feather(file.path(this.dir, "anno.feather"))
   dat  <- feather(file.path(this.dir, "data.feather"))
} else {
   anno <- feather(file.path(taxoFolder_FACS,"anno.feather"))
   dat  <- feather(file.path(taxoFolder_FACS, "data.feather"))
   anno <- as.data.frame(anno)
   colnames(anno) <- gsub("sample_name","sample_id",colnames(anno))  # Make sure it's called sample_id
}

cn   <- c("sample_id","cluster_color","cluster_id","cluster_label")
anno <- anno[,c(cn,setdiff(colnames(anno),cn))]

print("Reorder anno to match dat")
if (file.type=="tome") {
   samples = colnames(norm.dat)
} else {
   if("sample_id" %in% colnames(dat)){
     samples <- as.character(dat$sample_id)
   } else {
     samples <- as.character(dat$sample_name)
   }
}

common.samples = intersect(samples, anno$sample_id)

anno <- anno[match(common.samples,anno$sample_id),]  # Sometimes called sample_name
dat <- dat[match(common.samples,samples),]  # Sometimes called sample_name

##
# Subset FACS Data According to Each Reference
##
print("Subset dendrogram to only include neurons")
print("========================================================")
print("    HERE, CODE SHOULD BE CHANGED FOR EACH DENDROGRAM    ")
print("========================================================")
if (species=="mouse" && taxonomy_str=="ForeBrain") { #THIS IS FOR ForeBrain
   load(file.path(workFolder, "cl.final.rda"))
   #select.roi = unique(anno$roi_label)
   select.roi = c( unique(grep("TH", anno$roi_label,value=T)), 
                   unique(grep("STR[v|d]", anno$roi_label,value=TRUE)), 
                   unique(grep("PAL_GP",anno$roi_label,value=TRUE)))

   tb=with(anno, table(cluster_label, roi_label))
   tb.frac= tb/rowSums(tb)
   select.counts = rowSums(tb[,select.roi])
   select.ratio=rowSums(tb.frac[,select.roi])
   select.cl_label = intersect(names(select.ratio)[select.ratio > 0.1 & select.counts>=4], cl.df[cl.df$class_label %in% c("Ex","Inh","Hybrid"),"cluster_label"])
   select.cl = droplevels(cl[cl %in% (cl.df %>% filter(cluster_label %in% select.cl_label) %>% pull(cl))])
   select.sample_name = names(select.cl)

   kpSub = anno$sample_id %in% select.sample_name
  
   load(file.path(workFolder, "dend.rda"))

} else {
   dend = readRDS(file.path(taxoFolder_FACS,"dend.RData")) # Usually called "dend.RData"
   common.cluster_label = intersect(labels(dend), anno$cluster_label) 
   cl.df = anno[match(common.cluster_label, anno$cluster_label),]
   cl.df$cl = cl.df$cluster_id
   save(cl.df, file=file.path(refFolder, "cl.df.rda"))

   cl = as.factor(anno$cluster_id)
   names(cl) = anno$sample_id

   select.roi = rois
   if (length(select.roi)>0) {
      tb=with(anno, table(cluster_label, Region))
      tb.frac= tb/rowSums(tb)
      select.counts = rowSums(tb[,select.roi])
      select.ratio=rowSums(tb.frac[,select.roi])
      select.cl_label = names(select.ratio)[select.ratio > 0.1 & select.counts>=4]
   } else {
      tb=table(anno$cluster_label)
      select.cl_label = names(tb)[tb>=4]
   }
   select.cl = droplevels(cl[cl %in% (cl.df %>% filter(cluster_label %in% select.cl_label) %>% pull(cl))])
   select.sample_name = names(select.cl)

   kpSub = anno$sample_id %in% select.sample_name
   set.seed(1)
   if (str_detect(reference_str, "neuron")) {
      kpNeuron= !(anno$class_label %in% c("Non-Neuronal","Non-neuronal"))
      kpSub = kpSub & kpNeuron & subsampleCells(clusters = anno$cluster_label, subSamp = cellsPerCluster)
   } else {
      kpSub    = kpSub & subsampleCells(clusters = anno$cluster_label, subSamp = cellsPerCluster)
   }
   # incude Non-neuronal
   
   select.cl = anno$cluster_label[kpSub]
   names(select.cl) = anno$sample_id[kpSub]
   print(table(select.cl))
   print("check kpSub")
   print(sum(kpSub))

}   

print(" # Prune dendrogram with only selected clusters ")
select.dend = prune_dend(dend, setdiff(labels(dend), select.cl))
all.label = labels(dend)
prune.label = labels(select.dend)
lrank = order(match(prune.label, all.label))
select.dend = reorder_dend(select.dend, lrank)
print("(== Ignore NA NA here)")
dend_ref = select.dend
##
# confirm the label
##
if (!any(is.na(match(labels(dend_ref), cl.df$cl)))) {
   labels(dend_ref) = cl.df[match(labels(dend_ref), cl.df$cl), "cluster_label"]
}
print("Check that label names match and change if needed (ideally unnecessary)")
#if ("cl" %in% colnames(cl.df)) labels(dend_ref) = cl.df[labels(dend_ref), "cluster_label"]
cls <- unique(anno$cluster_label[order(anno$cluster_id)])[1:length(labels(dend_ref))]
new <- setdiff(labels(dend_ref),anno$cluster_label)
old <- setdiff(cls,labels(dend_ref))
if (length(new)>0)
  for (i in 1:length(new))
    anno$cluster_label[anno$cluster_label==old[i]] <- new[i]
rm(cls,old,new)

print("Subset the annotation and data files to only include the tree samples")

kpSamp   <- kpSub&is.element(anno$cluster_label,labels(dend_ref))

anno_ref <- anno[kpSamp, ]
dat_ref  <- dat[kpSamp, ]
cl_ref = cl[kpSamp]
select.sample_name = anno$sample_id[kpSamp]

select.cl = as.factor(anno$cluster_id[kpSamp])
names(select.cl) = anno$sample_id[kpSamp]
cl.label = labels(dend_ref)

print("Write reference feather")
print(table(select.cl))
##
# WRITE taxonomy  anno_ref, dat_ref, tSNEinto feather to reference folder
##
print("Write anno, data, tSNE into feather")
if(!file.exists(refFolder))  dir.create(refFolder)  # Create folder if needed
write_feather(anno_ref, file.path(refFolder,"anno.feather"))
anno_desc <- create_desc(as.data.frame(anno_ref), use_label_columns = TRUE)
write_feather(anno_desc, file.path(refFolder,"desc.feather"))
write_feather(dat_ref, file.path(refFolder,"data.feather"))

print("Subset and copy tsne.  (Revisit: this may not be necessary)")
if (file.exists(file.path(taxoFolder_FACS,"tsne.feather"))){
   tsne     <- feather(file.path(taxoFolder_FACS,"tsne.feather"))
   tsne_ref <- tsne[match(select.sample_name, tsne$sample_id),]
   write_feather(tsne_ref, file.path(refFolder,"tsne.feather"))
} else {
   print(paste("tsne.feather is missing from ", taxoFolder_FACS))
}
if (file.exists(file.path(taxoFolder_FACS,"tsne_desc.feather"))){
   file.copy(file.path(taxoFolder_FACS,"tsne_desc.feather"),file.path(refFolder,"tsne_desc.feather"))
} else {
   print(paste("tsne_desc.feather is missing from ", taxoFolder_FACS))
}

print("Write the subsetted dendrogram")
save(dend, file=file.path(refFolder, completeDendFile))
save(dend_ref, file=file.path(refFolder, refDendFile))


## 
# Save marker genes for patchseqQC 
#
# -- NOTE: THIS SECTION WILL LIKELY NEED TO BE UPDATED FOR EACH REFERENCE
# Required inputs:
# ----- Subclass calls for each cell
# ----- Broad class class calls for each cell
# ----- Distinction of neuron (e.g., mappable type) vs. non-neuron (e.g., contamination type)
#
##
#dend.label = cl.label[labels(dend_ref)]

print("Define conversion table (with cluster_label and subclass_label variables)")
conversions <- as.data.frame(anno_ref[match(labels(dend_ref),anno_ref$cluster_label),])

qcCellsPerCluster=100
print("Subsample and filter metadata and data")
goodSamp <- (anno_ref$subclass_label!="Low Quality")  & !is.na(anno_ref$class_label)  
kpSamp2  <- subsampleCells(anno_ref$subclass_label, qcCellsPerCluster )
kpSamp2  <- kpSamp2&goodSamp
annoQC   <- anno_ref[kpSamp2,]
annoQC$subclass_label = make.names(annoQC$subclass_label)

datQC    <- as.matrix(dat_ref[kpSamp2,names(dat_ref)!="sample_id"])
rownames(datQC) <- annoQC$sample_id
datQC    <- t(datQC)

print("Define class and subclass labels")
# --- For class we wrap neurons by class but retain non-neurons by subclass
classBr   <- annoQC$subclass_label
NN.labels = c("Non-Neuronal", "Non-neuronal", "glia", "Glia", "NN")
classBr[!(annoQC$class_label %in% NN.labels)] = annoQC$class_label[!(annoQC$class_label %in% NN.labels)]
classBr   <- factor(classBr)
subclassF <- factor(annoQC$subclass_label)
if (species=="mouse" && taxonomy_str=="ForeBrain") { #THIS IS FOR ForeBrain
   subclassF <- annoQC$subclass_label
   subclassF[grep("^HY", subclassF)] = "HY"
   subclassF[grep("^TH", subclassF)] = "TH"
   subclassF[grep("sAMY", subclassF)] = "sAMY"
   subclassF[grep("^STR", subclassF)] = "STR"
   subclassF[grep("^Sst", subclassF)] = "Sst"
   subclassF[grep("^OLF", subclassF)] = "OLF"
   subclassF = factor(subclassF) 
}


print("Define and output marker genes for each broad class and contamination class  (use 50 for now) ") 
# -- These are selected using some reasonable approach that could probably improved, if needed.    
save(datQC, subclassF,classBr, annoQC, NN.labels, file="Debug.rda")
markers    <- defineClassMarkers(datQC,subclassF,classBr,numMarkers = 50)
allMarkers <- unique(unlist(markers))
rownames(datQC) <- make.names(rownames(datQC))
countsQC   <- datQC[allMarkers,]
cpmQC      <- cpm(datQC)[allMarkers,]
normQC     <- log2(cpmQC + 1)
browser()
if (species=="mouse" && taxonomy_str=="ForeBrain") {
   allMarkers = intersect(allMarkers, rownames(norm.dat))
   normQC = norm.dat[allMarkers, kpSamp2]
   cpmQC  = 2^normQC - 1
}
save(markers, annoQC, countsQC, cpmQC, normQC, classBr, subclassF, allMarkers, file=file.path(refFolder, "QC_markers.RData"))
print("=========== OK, Set Up is done : Ready for build_patchseq_reference.R")
