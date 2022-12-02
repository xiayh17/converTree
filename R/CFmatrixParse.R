#' Convert CF matrix to tree data
#'
#' CF matrix is conflict-free matrix, tree data is a data table contains
#' `from`, `to` and `label` columns at least.
#'
#' @param CFmatrix_file a file path to CFMatrix file.
#'
#' @importFrom readr read_delim
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr rename group_by slice
#' @importFrom plyr mapvalues
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' CFmatrix_file = system.file("extdata", "ground_truth_tree.CFMatrix", package = "converTree")
#' cf2treedata(CFmatrix_file)
cf2treedata <- function(CFmatrix_file) {

  cf_mat <- readr::read_delim(CFmatrix_file,
                delim = "\t", escape_double = FALSE,
                trim_ws = TRUE)

  # duplicate process -------------------------------------------------------
  ## duplicate mutations and cells detecting
  dp_index_mut <- duplicated(t(cf_mat))
  dp_index_cell <- duplicated(cf_mat2[,-1])
  cf_mat2 <- cf_mat[,!dp_index_mut]
  cf_mat3 <- cf_mat2[!dp_index_cell,]

  ## index duplicate mutations and cells
  seed_mut = colnames(cf_mat3)[-1]
  seed_cell = unlist(cf_mat3[,1],use.names = F)

  ## duplicate mutations group
  same_columns_l = list()
  loosed_mut = colnames(cf_mat[,dp_index_mut])
  for (i in seq_along(seed_mut)) { # for every mutation remains,find loosed same mutation

    mut_i = seed_mut[i]
    # print(i)
    same_columns <- vector()
    for (j in loosed_mut) { ## from loosed mutation only
      # print(j)

      if (!(mut_i %in% unlist(same_columns_l)) & !(j %in% unlist(same_columns_l))) {

        if (all(cf_mat[,mut_i] == cf_mat[,j])) {

          same_columns <- c(same_columns,mut_i, j)

        }}

    }

    same_columns_l[[mut_i]] = unique(same_columns)

  } ## duplicate mutations group

  ## duplicate cells group
  same_cells_l = list()
  loose_cell = cf_mat[dp_index_cell,]$cellIDxmutID
  for (i in seq_along(seed_cell)) {
    cell_i = seed_cell[i]
    # print(i)
    same_cells <- vector()
    for (j in loose_cell) {
      # print(j)
      if (!(cell_i %in% unlist(same_cells_l)) & !(j %in% unlist(same_cells_l))) {
        if (all(cf_mat[cf_mat$cellIDxmutID==cell_i,-1] == cf_mat[cf_mat$cellIDxmutID==j,-1])) {
          same_cells <- c(same_cells,cell_i, j)
        }
      }
    }

    same_cells_l[[cell_i]] = unique(same_cells)

  } ## duplicate cells group
  # duplicate process -------------------------------------------------------

  # rebuild links from matrix -----------------------------------------------
  ## build links
  cell2mut <- cf_mat3 %>%
    pivot_longer(cols = 2:ncol(.),
                 names_to = "mutID") %>%
    rename(cellID=cellIDxmutID)

  cell2mut2 <- cell2mut[cell2mut$value != 0,]

  ## sort mutations for each cell according to mutations sum
  cell2mut2 <- cell2mut2 %>%
    group_by(mutID) %>%
    mutate(count=n()) %>%
    group_by(cellID) %>%
    mutate(rank = order(count,decreasing = TRUE))

  ## join cells mutations path
  ## a fake root
  ## if a top mutation dose exist in other cell and still top (True Top)
  ## it's parent node should to be fake root
  ## Different true top means subclone for these cells
  top_dat <- cell2mut2[cell2mut2$rank == 1,]
  ## show all true top
  unique(top_dat$mutID)

  ## leaves node should be buttom mutations
  b_dat <- cell2mut2 %>%
    group_by(cellID) %>%
    slice(which.max(rank))

  ## internal node
  in_dat <- setdiff(cell2mut2,top_dat)
  in_dat <- setdiff(in_dat,b_dat)

  ## some cells maybe only one mutation
  ## canbe find by
  # intersect(top_dat,b_dat)

  ## build tree data
  ## add fake root
  ## a tree that leaves numbers equal to cells
  ncells=unique(cell2mut2$cellID) %>% length()

  ## fake root
  t_root<-data.frame(
    parent=ncells+1,
    node=ncells+1,
    label=NA_character_
  )
  root_node = t_root$parent
  m_root<-merge(t_root$parent,
                unique(top_dat$mutID))
  m_root$label =unique(top_dat$mutID)
  colnames(m_root) <- c("parent","node","label")

  t_root <- rbind(t_root,m_root)

  ## leaves
  ## contains node number for each cell (leaf)
  t_leaves <- data.frame(
    parent=b_dat$mutID,
    node=1:length(b_dat$cellID),
    label=b_dat$cellID
  )

  link_mut <- function(x) {
    min = min(x)
    max = max(x)

    if (length(x)>2) {
      md = x[!x %in% c(min,max)]
      vec=c(min, rep(sort(md),each=2), max)
    } else {
      vec=c(min,max)
    }
    mat=matrix(vec,nrow=2,ncol=length(vec)/2)
    return(mat)
  }

  ## parse node number for each mutation (node)
  for (i in seq_along(t_leaves$label)) {
    cell = t_leaves$label[i]

    df = cell2mut2[cell2mut2$cellID %in% cell,]

    if (i==1) {
      tmp_df=data.frame(
        mutIn = df$mutID,
        node = df$rank+root_node
      )
      link_o=link_mut(tmp_df$node)

      link_df=data.frame(
        parent=link_o[1,],
        node=link_o[2,],
        label=link_o[2,]
      )
    } else {
      mut_df = cell2mut2[cell2mut2$mutID %in% tmp_df$mutIn,]

      mutNow = df$mutID
      mutNew = setdiff(mutNow,tmp_df$mutIn)
      newNode = seq_along(mutNew) + max(tmp_df$node)
      print(mutNew)
      print(newNode)
      if (length(mutNew)>0) {
        tmp_df=data.frame(
          mutIn = c(tmp_df$mutIn,mutNew),
          node = c(tmp_df$node,newNode)
        )
      }

      if (length(mutNew)==1) {
        pairRank=df[df$mutID == mutNew,]$rank-1
        if (pairRank != 0) {
          pairMut=df[df$rank == pairRank,]$mutID
          pairNode=tmp_df[tmp_df$mutIn==pairMut,]$node
          pairNodes=c(pairNode,newNode)
          link_o=link_mut(pairNodes)
          link_df=data.frame(
            parent=c(link_df$parent,link_o[1,]),
            node=c(link_df$node,link_o[2,]),
            label=c(link_df$node,link_o[2,])
          )
        }
      }

      if (length(mutNew)>=2) {
        checkRank=df[df$mutID == mutNew[1],]$rank-1
        if (checkRank != 0) {
          pairMut=df[df$rank == checkRank,]$mutID
          pairNode=tmp_df[tmp_df$mutIn==pairMut,]$node
          pairNodes=c(pairNode,newNode)
          link_o=link_mut(pairNodes)
          link_df=data.frame(
            parent=c(link_df$parent,link_o[1,]),
            node=c(link_df$node,link_o[2,]),
            label=c(link_df$node,link_o[2,])
          )
        } else {
          link_o=link_mut(newNode)
          link_df=data.frame(
            parent=c(link_df$parent,link_o[1,]),
            node=c(link_df$node,link_o[2,]),
            label=c(link_df$node,link_o[2,])
          )
        }

      }

    }


  }

  link_df<-unique(link_df)

  link_df$label <- plyr::mapvalues(link_df$label,
                                   from = tmp_df$node,
                                   to = tmp_df$mutIn,
                                   warn_missing = FALSE)

  t_all <- rbind(t_root,t_leaves,link_df)

  if (identical(cell2mut2$mutID[!cell2mut2$mutID %in% t_all$label[grep("M",t_all$label)]] %>% unique(),character(0))) {
    message("Tree data bone seems ok")
  }
  # rebuild links from matrix -----------------------------------------------

  # insert duplicate back ---------------------------------------------------
  # for each label, is there any duplications?
  same_cells_l
  same_columns_l
  # mutation
  # mutation is easy as modify label of tree dat
  for (i in seq_along(same_columns_l)) {
    ele=same_columns_l[[i]]
    if (length(ele) >0 ) {
      mut=names(same_columns_l)[i]
      new_label=paste0(ele,collapse = "|")
      t_all$label[t_all$label == mut] = new_label
    }
  }

  # cell
  # cell is a little complex for add leaves to tree
  # maybe like mutation is a good idea
  for (i in seq_along(same_cells_l)) {
    ele=same_cells_l[[i]]
    if (length(ele) >0 ) {
      cell=names(same_cells_l)[i]
      new_label=paste0(ele,collapse = "|")
      t_all$label[t_all$label == cell] = new_label
    }
  }
  # insert duplicate back ---------------------------------------------------

  return(t_all)

}

# cf2gvfile <- function(CFmatrix_file) {
#
# }

