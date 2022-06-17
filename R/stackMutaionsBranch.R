#' Stack mutations to node and make branch length to number of mutation
#'
#' stack_mutationTree is a naive function for stack mutation node in branch to
#' node in a mutation tree of class "phylo".
#'
#' @param tre a mutation tree
#' @importFrom treeio as_tibble
#' @importFrom dplyr bind_rows pull distinct mutate rename group_modify group_by arrange n
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

  ## tree to data ----
  tre_dat <- tre %>% treeio::as_tibble()

  ## subset tree data by node and tip ----
  root_df <- tre_dat[tre_dat$parent == tre_dat$node,]
  leaf_df <- tre_dat[!tre_dat$node %in% tre_dat$parent,]
  # node internal
  dt_tmp <- as.data.frame(table(tre_dat$parent))
  in_node <- dt_tmp$Var1[dt_tmp$Freq >= 2 & dt_tmp$Var1 != root_df$parent]
  in_df <- tre_dat[tre_dat$parent %in% in_node & !(tre_dat$node %in% leaf_df$node),]
  # node only in branch
  b_df <- tre_dat[!tre_dat$node %in% root_df$node & !tre_dat$parent %in% in_df$parent & !tre_dat$node %in% leaf_df$node,]

  # function for find new number ----
  indf_new <- function(in_df) {

    ## for parent
    in_df$parentN <- ""
    in_df$parentN <- factor(x=in_df$parent,labels = order(unique(in_df$parent))) %>% as.numeric()
    in_df$parentN <- in_df$parentN+root_df$parent[1]

    ## node
    ## for node with leaf
    in_df_o <- in_df[in_df$node %in% leaf_df$node,]
    ## for node with leaf branch
    in_df_m <- in_df[!in_df$node %in% leaf_df$node,]
    index_branch_leaf <- apply(in_df_m, 1, function(x){

      all(as.numeric(x["node"]) > in_df_m$parent)

    })
    in_df_i <- in_df_m[!index_branch_leaf,]
    in_df_b <- in_df_m[index_branch_leaf,]
    # internal node
    in_df_i$nodeN <- in_df_i$parentN+1
    # link to leaf inderectly
    # find node to leaf
    find_real_leafs <- function(leafsP){

      while (!leafsP %in% leaf_df$node) {
        leafsP <- tre_dat[tre_dat$parent == leafsP,]$node
      }
      return(leafsP)
    }
    ## find
    in_df_b$nodeN <- apply(in_df_b, 1, function(x){
      tmp = x["node"]
      t <- find_real_leafs(tmp)
      return(t)
    })

    dplyr::bind_rows(in_df_o,in_df_i,in_df_b)

  }

  leaf_new <- function(leaf_df) {

    leaf_df$parentN <- ""
    ## parent is not real node or root
    index1 <- leaf_df$parent %in% in_df$parent | leaf_df$parent %in% root_df$parent
    tmp1 <- leaf_df[index1,]
    ## parent in branch
    index2 <- !(leaf_df$parent %in% in_df$parent | leaf_df$parent %in% root_df$parent)
    tmp2 <- leaf_df[index2,]

    ## for parent not in branch
    tmp1$parentN <- unique(in_df[in_df$parent == tmp1$parent,]$parentN)
    ## for parent in branch
    find_real_leafsP <- function(leafsP){
      while (!(leafsP %in% root_df$parent | leafsP %in% in_df$parent)) {
        leafsP <- tre_dat[tre_dat$node == leafsP,]$parent
      }
      return(leafsP)
    }
    ## find
    tmp2$parentN <- apply(tmp2, 1, function(x){
      tmp = x["parent"]
      t <- find_real_leafsP(tmp)
      ## real parent new number
      tt <- unique(in_df[in_df$parent == t,]$parentN)
      return(tt)
    })

    ndat <- dplyr::bind_rows(tmp1,tmp2)
    ndat$nodeN <- ndat$node
    return(ndat)
  }

  b_new <- function(b_df) {

    ## find real parent
    find_real_p <- function(p){

      while (!(p %in% root_df$parent | p %in% in_df$parent)) {
        p <- tre_dat[tre_dat$node == p,]$parent
      }
      return(p)

    }
    ## find
    b_df$parentN <- ""
    b_df$parentN <- apply(b_df, 1, function(x){

      tmp = as.numeric(x["parent"])
      t <- find_real_p(tmp)
      t <- as.numeric(t)

      ## real parent new number
      tt <- unique(in_df[in_df$parent == t,]$parentN)

      if(identical(numeric(0),tt)){
        tt <- unique(root_df[root_df$parent == t,]$parentN)
      }

      return(tt)
    })


    ## find real node
    find_real_l <- function(l){

      while (!(l %in% leaf_df$node | l %in% in_df$parent)) {
        l <- tre_dat[tre_dat$parent == l,]$node
      }
      return(l)

    }
    ## find
    b_df$nodeN <- apply(b_df, 1, function(x){
      tmp = as.numeric(x["node"])
      t <- find_real_l(tmp)
      t <- as.numeric(t)
      ## real parent new number
      tt <- unique(in_df[in_df$parent == t,]$parentN)
      if(identical(numeric(0),tt)){
        tt <- unique(leaf_df[leaf_df$node == t,]$nodeN)
      }
      return(tt)
    })

    return(b_df)

  }


  ## end of functions -----

  # find new node and parent ----
  root_df$parentN <- ""
  root_df$parentN <- root_df$parent
  root_df$nodeN <- root_df$node

  in_df <- indf_new(in_df)
  leaf_df <- leaf_new(leaf_df)
  b_df <- b_new(b_df)

  new_tree <- dplyr::bind_rows(root_df,leaf_df,in_df,b_df)

  ## to tree ----
  ## add label with stack mutations
  t2 <- merge(tre_dat,new_tree)

  stackMutationLabel <- function(dt) {

    dt$labelN = ""
    dt2 <- dplyr::arrange(dt,parent)
    dt2$labelN = paste0(dt2$label,collapse = "|")
    return(dt2)

  }

  t22 <- t2 %>%
    dplyr::group_by(nodeN) %>%
    dplyr::group_modify(~ {
      .x %>%
        stackMutationLabel
    })

  ## add branch length
  nt2 <- t22[,c(5,1,6)] %>% dplyr::rename(parent = parentN, node = nodeN, label = labelN) %>%
    dplyr::group_by(parent, node) %>%
    dplyr::mutate(branch.length = dplyr::n())%>%
    dplyr::distinct()

  ## root should be NA
  nt2[nt2$parent == nt2$node,]["branch.length"] = NA
  nt2[nt2$parent == nt2$node,]["label"] = NA

  ## convert to phylo format
  nttre2 <- nt2 %>% ape::as.phylo(length = "branch.length")
  ## fix label info
  nttre2[["node.label"]] <- dplyr::pull(nt2[nt2$node %in% nttre2[["node.label"]],],label)
  nttre2[["tip.label"]] <- dplyr::pull(nt2[nt2$node %in% nttre2[["tip.label"]],],label)
  return(nttre2)

}
