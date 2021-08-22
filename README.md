# tree_mapping
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
