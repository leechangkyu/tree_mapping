#############################################
#
# run_patchseq_reference.R
#
#
# input args
# args[1] = refFolder_FACS /allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_ForeBrain1250_07222020
# args[2] = refFolder /allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_ForeBrain1250_07222020/
# args[3] = cellsPerCluster : Number of cells per cluster to keep in reference (set to Inf if you don't want to subsample)
# args[4] = completeDendFile 
# args[5] = file.type "tome" or "feather" 
# if (args[5] =="tome")
# args[6] = tome.FN (if file.type=="tome")
# args[7] = workFolder /allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_ForeBrain1250_07222020/
#
# CTX-HPF
# /allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/CTX-HIP-20200915 
#############################################

##
# assign input args to variables 
##
args              = commandArgs(TRUE)

species = "mouse"
taxonomy_str = "ForeBrain"
#taxonomy_str = "CTX-HPF"
#taxonomy_str = "V1-ALM"

#Mouse V1
if (species == "mouse" && taxonomy_str == "V1-ALM") {
   args[1]="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520"
   args[2]="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/reference_taxonomies/mouse/V1-ALM"
   args[3]=500
   args[4]="dend.RData"
   args[5]= "feather"
   args[6]= ""
   args[7]= ""
   reference_str = "V1-ALM" 
   rois = c()
}

#Human MTG
if (species == "human" && taxonomy_str == "MTG") {
   args[1]="/allen/programs/celltypes/workgroups/rnaseqanalysis/Nik/Analyses_for_great_apes_paper/Shiny_obj/human/"
   args[2]="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/reference_taxonomies/human/MTG_great_ape_neuron_test0721"
   #args[2]="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/reference_taxonomies/human/MTG_great_ape_all_unified"
   args[3]=500
   args[4]="dend.RData"
   args[5]= "feather"
   args[6]= ""
   args[7]= ""
   #reference_str = "MTG_great_ape_all" 
   reference_str = "MTG_great_ape_neuron" 
   rois = c()
}

#Muse CTX-HIP
if ( species == "mouse" && taxonomy_str == "CTX-HPF" ) {
   #args[1]="/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/Cortex_HIP"
   args[1]="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_CTX_HIP_20200915"
   args[2]="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/reference_taxonomies/mouse/CTX-HPF/"
   args[3]= 500
   args[4]="dend.RData"
   args[5]= "feather"
   args[6]= ""
   args[7]= ""
   reference_str = "CTX-HPF" 
   #rois = c("VISp","TCx","FCx","MOp","TEa","HIPCA1")
   #rois = c("VISp","MOp","TEa-PERI-ECT","HIP")
   rois = c()

}

#Mouse FB
if (species == "mouse" && taxonomy_str == "ForeBrain") {
   #args[1]= "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_ForeBrain_20201203"
   args[1]= "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/foreBrain/"
   #args[2]="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/reference_taxonomies/mouse/ForeBrain/"
   args[2]= "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_forebrain_20201204"
   args[3]= 500
   args[4]= "dend.rda"
   args[5]= "tome"
   args[6]= "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/tomes/facs/mouse_ForeBrain_20201204/smrt.tome"
   args[7]= "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/foreBrain/"
   flag_setup_skip = FALSE
   rois = c()
}

source("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/Utils.R")
set.seed(100)

taxoFolder_FACS   = args[1] 
refFolder         = args[2] 
cellsPerCluster   = as.numeric(args[3])
completeDendFile  = args[4] 
file.type         = args[5]
tome.FN        = args[6] 
workFolder     = args[7] 

refDendFile       = "dend.reference.RData"
qcCellsPerCluster = 100 


options(error=recover)
print(" === Prepare for patchseq reference  ")
source("setup_patchseq_reference.R")
browser()
print(" === Build patchseq reference        ")
source("build_patchseq_reference.R")

print(" REFERRENCE is BUILT in ")
print(refFolder)
print(dir(refFolder))
