#' Stack mutations to node and make branch length to number of mutation
#'
#' stack_mutationTree is a naive function for stack mutation node in branch to
#' node in a mutation tree of class "phylo".
#'
#' @param tre a mutation tree
#' @importFrom treeio as_tibble
#' @importFrom dplyr bind_rows pull distinct mutate rename group_modify group_by arrange n full_join select
#' @importFrom ape as.phylo
#'
#' @return
#' a tree of class "phylo"
#' @export
#'
#' @examples
#' library(ggtree)
#' nwk_text = "((((((((((8)7,(1)12)9)10)2)4)3,5,(((17)11)16)14)15)6)13)18;"
#' tre <- treeio::read.newick(text = nwk_text)
#' ggtree(tre) +
#'     geom_nodepoint(size = 1,color = "red",shape = 20) +
#'     geom_nodelab(color = "red",nudge_y = 0.1)+
#'     geom_nodelab(aes(label = node))+
#'     geom_tippoint(size = 1,color = "red",shape = 20)+
#'     geom_tiplab(color = "red",nudge_y = 0.1)+
#'     geom_tiplab(aes(label = node))
#' stre <- stack_mutationTree(tre)
#' ggtree(stre)+
#'   geom_nodepoint(size = 1,color = "red",shape = 20) +
#'   geom_nodelab(color = "red",nudge_y = 0.1,nudge_x = -0.4)+
#'   geom_nodelab(aes(label = node))+
#'   geom_tippoint(size = 1,color = "red",shape = 20)+
#'   geom_tiplab(color = "red",nudge_y = 0.1,nudge_x = -0.4)+
#'   geom_tiplab(aes(label = node))+
#'   theme_tree2()+
#'   scale_x_continuous(breaks = c(0:10), limits = c(0,12))
stack_mutationTree <- function(tre) {

  # tree in data frame format
  tre_dat <-  tre %>%
    as_tibble()

  # nodes
  tre_dat2 <- tre_dat %>%
    group_by(parent) %>%
    mutate(cn = n())

  node_index <- factor(tre_dat2$node, labels = tre_dat2$label)
  tre_dat2$parent_label <- node_index[tre_dat2$parent]

  # split node
  tre_dat2$group <- NA

  tre_dat2$group[tre_dat2$node == tre_dat2$parent] = "root"
  tre_dat2$group[!tre_dat2$node %in% tre_dat2$parent] = "leaf"
  root_df = tre_dat2[tre_dat2$group == "root",] %>% na.omit()
  leaf_df = tre_dat2[tre_dat2$group == "leaf",] %>% na.omit()

  real_nodes = tre_dat2$parent[tre_dat2$cn >= 2 & is.na(tre_dat2$group)] %>% unique()
  real_nodes = real_nodes[real_nodes != root_df$node]
  tre_dat2$group[tre_dat2$node %in% real_nodes] = "nodes"

  tre_dat2$group[is.na(tre_dat2$group)] = "branch"

  # node recount
  in_df = tre_dat2[tre_dat2$group == "nodes",]
  br_df = tre_dat2[tre_dat2$group == "branch",]

  root_df$parentN <- NULL
  root_df$parentN <- root_df$parent
  root_df$nodeN <- root_df$node
  root_node_ = root_df$parent[1]

  ## real parent node should be recount from root
  in_df$parentN <- NULL
  in_df$nodeN <- factor(x=in_df$node,labels = order(unique(in_df$node))) %>% as.numeric()
  in_df$nodeN <- in_df$nodeN+root_node_
  ## real child node should be recount from root

  ## for parent in branch
  find_real_fp <- function(fp){

    while (!(fp %in% root_df$parent | fp %in% in_df$node)) {
      fp <- tre_dat2[tre_dat2$node == fp,]$parent
    }
    return(fp)

  }

  in_index <- function(x) {
    in_df[in_df$node == x,]$nodeN
  }



  in_df$parentN <- apply(in_df, 1, function(i){
    tmp = as.numeric(i["parent"])
    t <- find_real_fp(tmp)

    if (t == root_node_) {
      root_node_
    } else {
      in_index(t)
    }

  }) %>% unlist()



  leaf_df$parentN <- NULL
  ## parent is not real node or root
  index1 <- leaf_df$parent %in% br_df$parent
  tmp1 <- leaf_df[index1,]

  ## parent in real node
  tmp1$parentN <- sapply(tmp1$parent,in_index)
  tmp1$nodeN <- tmp1$node

  ## parent in branch
  index2 <- !index1
  tmp2 <- leaf_df[index2,]

  ## for parent in branch
  find_real_leafsP <- function(leafsP){
    while (!(leafsP %in% root_df$parent | leafsP %in% in_df$node)) {
      leafsP <- tre_dat2[tre_dat2$node == leafsP,]$parent
    }
    return(leafsP)
  }

  ## find real node
  leaf_index <- function(x) {

    tmp_df <- data.frame(
      nodes = c(in_df$parent, in_df$node),
      N = c(in_df$parentN, in_df$nodeN)
    ) %>% unique()
    tmp_df[tmp_df$nodes == x,]$N

  }

  tmp2$parentN <- apply(tmp2, 1, function(i){
    tmp = as.numeric(i["parent"])
    t <- find_real_leafsP(tmp)

    if (t == root_node_) {
      root_node_
    } else {
      leaf_index(t)
    }

  }) %>% unlist()

  tmp2$nodeN <- tmp2$node

  leaf_df <- dplyr::bind_rows(tmp1,tmp2)

  ## br_df
  ## find real parent
  find_real_p <- function(p){

    while (!(p %in% root_df$parent | p %in% in_df$node)) {
      p <- tre_dat2[tre_dat2$node == p,]$parent
    }
    return(p)

  }
  ## find
  br_df$parentN <- NULL
  br_df$parentN <- apply(br_df, 1, function(x){

    tmp = as.numeric(x["parent"])
    t <- find_real_p(tmp)

    ## real parent new number
    if (t == root_node_) {
      root_node_
    } else {
      leaf_index(t)
    }

  }) %>% unlist()


  ## find real node
  find_real_l <- function(l){

    while (!(l %in% leaf_df$node | l %in% in_df$node)) {
      l <- tre_dat2[tre_dat2$parent == l,]$node
    }
    return(l)

  }

  leaf_index2 <- function(x) {

    tmp_df <- data.frame(
      nodes = c(in_df$parent, in_df$node, leaf_df$parent,leaf_df$node),
      N = c(in_df$parentN, in_df$nodeN, leaf_df$parentN,leaf_df$nodeN)
    ) %>% unique()
    tmp_df[tmp_df$nodes == x,]$N

  }

  ## find
  br_df$nodeN <- apply(br_df, 1, function(x){
    tmp = as.numeric(x["node"])
    t <- find_real_l(tmp)

    ## real parent new number
    if (t == root_node_) {
      root_node_
    } else {
      leaf_index2(t)
    }

  }) %>% unlist()

  ## recombined
  new_tree <- dplyr::bind_rows(root_df,leaf_df,in_df,br_df)

  ## to tree ----
  ## add label with stack mutations

  stackMutationLabel <- function(dt) {

    dt$labelN = NULL
    dt2 <- dplyr::arrange(dt,parent)
    dt2$labelN = paste0(dt2$label,collapse = "|")
    return(dt2)

  }

  t22 <- new_tree %>%
    dplyr::group_by(parentN,nodeN) %>%
    dplyr::group_modify(~ {
      .x %>%
        stackMutationLabel
    })

  ## add branch length
  nt2 <- t22[,c("parentN","nodeN","labelN")] %>% dplyr::rename(parent = parentN, node = nodeN, label = labelN) %>%
    dplyr::group_by(parent, node) %>%
    dplyr::mutate(branch.length = dplyr::n())%>%
    dplyr::distinct()

  ## root should be NA
  nt2[nt2$parent == nt2$node,]["branch.length"] = NA
  nt2[nt2$parent == nt2$node,]["label"] = NA

  ## convert to phylo format
  nttre2 <- nt2 %>% ape::as.phylo(length = "branch.length")
  ## fix label info
  index_label <- function(node) {

    nt2[nt2$node==node,]$label

  }

  nttre2[["node.label"]] <- Map(index_label,nttre2[["node.label"]]) %>% unlist()
  nttre2[["tip.label"]] <- Map(index_label,nttre2[["tip.label"]]) %>% unlist()

  return(nttre2)

}
