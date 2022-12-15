
CFmatrix_file = system.file("extdata", "ground_truth_tree.CFMatrix", package = "converTree")
t_all = cf2treedata(CFmatrix_file)
nt2=t_all
nttre2 <- t_all %>% ape::as.phylo()

## fix label info
index_label <- function(node) {

  nt2[nt2$node==node,]$label

}

nttre2[["tip.label"]] <- Map(index_label,nttre2[["tip.label"]]) %>% unlist()

treeio::write.tree(nttre2)

