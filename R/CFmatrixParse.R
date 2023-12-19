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

  ## is CF
  if (!isCFM(CFmatrix_file)) {
    stop("Not a conflict free matrix!")
  }

  ## import data
  cf_mat <- readr::read_delim(CFmatrix_file,
                              delim = "\t", escape_double = FALSE,
                              trim_ws = TRUE, show_col_types = FALSE)

  ## count mutations
  sum_mut=colSums(cf_mat[,-1])
  sum_cell=rowSums(cf_mat[,-1])
  mut_names=colnames(cf_mat)[-1]

  ## if there is a column all = 0
  ## drop it
  zero_muts_index=sum_mut==0
  zero_mut=mut_names[zero_muts_index]
  if (length(zero_mut)>0) {
    warning(c(paste0(zero_mut,collapse = ",")," Not be contained in any cell. Removed!"))
  }

  ## if there is a row all = 0
  ## drop it
  zero_cells_index=sum_cell==0
  if (sum(zero_cells_index)>0) {
    warning(c(length(zero_cells_index)," Cells Not contain any mutation. Removed!"))
  }

  ## remove
  cf_mat=cf_mat[!zero_cells_index,c(TRUE,!zero_muts_index)]

  ## if there is a column all = 1
  ## set as root
  all_mut_index=sum_mut==nrow(cf_mat)
  all_mut=mut_names[all_mut_index]
  if (length(all_mut)>0) {
    warning(c(paste0(all_mut,collapse = ",")," found in all cells."))
  }

  # duplicate process -------------------------------------------------------
  ## duplicate mutations and cells detecting
  dp_index_mut <- duplicated(t(cf_mat))
  dp_index_cell <- duplicated(cf_mat[,-1])
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
    mutate(rank = rank(-count))

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
  # in_dat <- dplyr::setdiff(cell2mut2,top_dat)
  # in_dat <- dplyr::setdiff(in_dat,b_dat)

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
  # if (length(all_mut)>0) {
  #   t_root$label=paste0(all_mut,collapse = "|")
  # }
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

  ## from leaves and root
  ## number every mutation is possible
  ## in each cell, rank convert to node number according to root and other cell
  link_mut <- function(x) {

    min = min(x)
    max = max(x)

    if (length(x)>2) {
      md = x[!x %in% c(min,max)]
      vec=c(min, rep(sort(md),each=2), max)
    } else if (length(x)==2) {
      vec=c(min,max)
    } else if (length(x)==1) {
      stop("vec must be of a vector of length greater than or equal to 2!")
    }
    # vec_link=c(min-1,min,vec)
    vec_link=vec
    mat=matrix(vec_link,nrow=2,ncol=length(vec_link)/2)
    return(mat)
  }

  ## parse node number for each mutation (node)
  cells = t_leaves$label
  for (i in seq_along(cells)) {
    # for (i in 1:8) {
    cell = cells[i]
    df = cell2mut2[cell2mut2$cellID %in% cell,]
    # print(cell)
    ## in the first cell, only according to root
    ## rank add root node number directly
    if (i==1) {
      # print(df)
      tmp_df=data.frame(
        mutIn = df$mutID,
        node = df$rank+root_node
      )
      # print(tmp_df)
      ## a link of node is also can be infered from node number
      newNode = tmp_df$node
      if (length(newNode) == 1) {
        pairsNodes = c(root_node,newNode)
        link_o=link_mut(pairsNodes)
      } else {
        link_o=link_mut(newNode)
      }
      ## first row be parent column
      ## second row be node column
      link_df=data.frame(
        parent=link_o[1,],
        node=link_o[2,],
        label=link_o[2,]
      )
      # print(link_df)
      ## above tmp_df storage the mutations and node number
      ## link_df storage the link of nodes
    } else {
      # print(df)
      ## according to other cells node number storage in tmp_df
      ## check current cell mutations is show before
      mutNow = df$mutID
      ## for these not include before
      mutNew = setdiff(mutNow,tmp_df$mutIn)
      # print(mutNew)
      # if there newNode exists
      if (length(mutNew)>0) {

        ## add to the node number continuly
        ## should rank as the same
        before_end = max(tmp_df$node)
        # newNode_s = before_end+1
        # newNode_e = before_end+before_end
        Ranks = df[df$mutID %in% mutNew,]$rank
        newNode = Ranks-(min(Ranks)-1) + before_end

        # storage New node number with olders
        tmp_df=data.frame(
          mutIn = c(tmp_df$mutIn,mutNew),
          node = c(tmp_df$node,newNode)
        )

        # for the rank 1 mutation
        # it should be link to root, which already in t_root, skip
        # newNode rank
        if (length(newNode)==1) {

          newRank=df[df$mutID == mutNew,]$rank

          if (newRank!=1) {

            ## find rank-1 node as pair
            pairMut=df[df$rank == (newRank-1),]$mutID
            pairNode=tmp_df[tmp_df$mutIn==pairMut,]$node
            pairNodes=c(pairNode,newNode)
            # print(pairNodes)
            link_o=link_mut(pairNodes)
            link_df=data.frame(
              parent=c(link_df$parent,link_o[1,]),
              node=c(link_df$node,link_o[2,]),
              label=c(link_df$node,link_o[2,])
            )

          }

        } else if (length(newNode)>1) {

          for (k in seq_along(mutNew)) {

            newRank=df[df$mutID == mutNew[k],]$rank

            if (newRank!=1) {

              ## find rank-1 node as pair
              # infer link of new nodes
              pairMut=df[df$rank == (newRank-1),]$mutID
              pairNode=tmp_df[tmp_df$mutIn==pairMut,]$node
              pairNodes=c(pairNode,newNode[k])
              # print(pairNodes)
              link_o=link_mut(pairNodes)
              # storage New link with olders
              link_df=data.frame(
                parent=c(link_df$parent,link_o[1,]),
                node=c(link_df$node,link_o[2,]),
                label=c(link_df$node,link_o[2,])
              )

            }

          }

        }

      }



    }

  }

  link_df<-unique(link_df)

  # link_df[link_df$parent==link_df$label,]
  # unlist(link_df) |> unique() |> sort()

  ## map node number to tree data
  t_root$node <- plyr::mapvalues(t_root$node,
                                 from = tmp_df$mutIn,
                                 to = tmp_df$node,
                                 warn_missing = FALSE) %>% as.numeric()

  link_df$label <- plyr::mapvalues(link_df$label,
                                   from = tmp_df$node,
                                   to = tmp_df$mutIn,
                                   warn_missing = FALSE)

  t_leaves$parent <- plyr::mapvalues(t_leaves$parent,
                                     from = tmp_df$mutIn,
                                     to = tmp_df$node,
                                     warn_missing = FALSE) %>% as.numeric()

  t_all <- unique(rbind(t_root,t_leaves,link_df))

  if (identical(cell2mut2$mutID[!cell2mut2$mutID %in% t_all$label[grep("M",t_all$label)]] %>% unique(),character(0))) {
    message("Tree data bone seems ok")
  }
  # rebuild links from matrix -----------------------------------------------

  # insert duplicate back ---------------------------------------------------
  # for each label, is there any duplications?
  # same_cells_l
  # same_columns_l
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
  # know which cell will insert others
  for (i in seq_along(same_cells_l)) {
    ele=same_cells_l[[i]]
    if (length(ele) >0 ) {
      cell=names(same_cells_l)[i]
      keynode=as.numeric(na.omit(t_all$node[t_all$label == cell]))
      newcell=setdiff(ele,cell)
      newnode=keynode+(1:length(newcell))
      newparent=na.omit(t_all$parent[t_all$label == cell])+length(newcell)
      new_df = data.frame(
        parent=rep(newparent,length(newcell)),
        node=newnode,
        label=newcell
      )
      t_all$parent=t_all$parent+length(newcell)
      t_all$node[t_all$node>keynode]=t_all$node[t_all$node>keynode]+length(newcell)
      t_all=rbind(t_all,new_df)
    }
  }
  # insert duplicate back ---------------------------------------------------

  return(t_all)

}

#' Convert tree data to gv file
#'
#' [gv](https://graphviz.org/about/) file is a graph format in dot language.
#'
#' @param treedata tree data from cf2treedata
#' @param gvfile a file name to save result
#' @param highlight which mutation or cell name to highlight
#' @param highcolor highlight color
#' @param lightNodeKey mutation or cell name to target a node
#' @param lightNodeColor color of targeted node
#' @param is_cell `auto_name` or `auto_position` or a vector list of bool
#' @param cell_prefix if `auto_name` used in is_cell, a common prefix for cell name
#' @param cell_style more in [here](https://graphviz.org/docs/attrs/style/)
#' @param cell_shape more in [here](https://graphviz.org/doc/info/shapes.html)
#' @param cell_fill more in [here](https://graphviz.org/docs/attrs/fillcolor/)
#' @param mut_style more in [here](https://graphviz.org/docs/attrs/style/)
#' @param mut_shape more in [here](https://graphviz.org/doc/info/shapes.html)
#' @param mut_fill more in [here](https://graphviz.org/docs/attrs/fillcolor/)
#'
#' @importFrom igraph graph_from_data_frame vertex.attributes write_graph vertex.attributes<-
#'
#' @return a file
#' @export
#'
#' @examples
#' CFmatrix_file = system.file("extdata", "ground_truth_tree.CFMatrix", package = "converTree")
#' t_all = cf2treedata(CFmatrix_file)
#' treedata2gv(t_all,gvfile="tree.gv",highlight = c("M3","M42","M1"))
treedata2gv<-function(treedata,gvfile="tree.gv",
                      highlight=NULL,
                      highcolor="red",
                      lightNodeKey=NULL,
                      lightNodeColor="red",
                      is_cell="auto_name",cell_prefix="sc",
                      cell_style="filled",cell_shape="box",cell_fill="white",
                      mut_style="filled",mut_shape="ellipse",mut_fill="grey82") {
  edges=treedata[!treedata$parent==treedata$node,1:2]
  colnames(edges)=c("from","to")
  nodes=treedata$node
  treedata$label[is.na(treedata$label)]="-"
  node_labels=formatNodelabel(treedata$label,highlight=highlight,highcolor=highcolor)
  ## distinguish cells and mutations
  if (is_cell=="auto_name") {
    ## auto_name
    is_cell = startsWith(treedata$label,"sc")
  } else if (is_cell=="auto_position") {
    ## auto_position
    is_cell = ifelse(treedata$node < root_,TRUE,FALSE)
  } else if (is_cell=="manually") {
    ## manually
    is_cell = is_cell
  }

  root_ = treedata$node[treedata$parent==treedata$node]
  g <- igraph::graph_from_data_frame(edges, directed=TRUE, vertices=nodes)
  vertex.attributes(g)$name=NULL
  vertex.attributes(g)$label=node_labels
  vertex.attributes(g)$style=ifelse(is_cell,cell_style,mut_style)
  vertex.attributes(g)$shape=ifelse(is_cell,cell_shape,mut_shape)
  vertex.attributes(g)$fillcolor=ifelse(is_cell,cell_fill,mut_fill)
  vertex.attributes(g)$fillcolor=lightNode(vertex.attributes(g)$label,
                                           vertex.attributes(g)$fillcolor,
                                           lightNodeKey,lightNodeColor)
  igraph::write_graph(g, gvfile, "dot")
  # read the text file
  text <- readLines(gvfile)

  # search for the line containing "<NA>" in the label
  line_index <- which(grepl("label=.*", text))

  # remove the quotation marks from the line
  text[line_index] <- gsub("label=\"", "label=", text[line_index])
  text[line_index] <- gsub("\"$", "", text[line_index])
  text <- gsub('\\\\\"', '\"',text)

  # save the modified text to a new file
  writeLines(text, gvfile)
}

#' Conflict matrix or Not
#'
#' A logical return, TRUE for conflict free.
#'
#' @param CFmatrix_file a file path to CFMatrix file.
#'
#' @return a bool value
#' @export
#'
#' @examples
#' CFmatrix_file = system.file("extdata", "ground_truth_tree.CFMatrix", package = "converTree")
#' isCFM(CFmatrix_file)
isCFM <- function(CFmatrix_file) {

  cf_mat <- readr::read_delim(CFmatrix_file,
                              delim = "\t", escape_double = FALSE,
                              trim_ws = TRUE, show_col_types = FALSE)

  num_mat<-as.matrix(cf_mat[,-1])

  res_bool <- isConflictFree(num_mat)

  return(res_bool)

}

# R version code too slow
# rewrite as function `isConflictFree` in is.conflict_free.cpp
# is.conflict <- function(CFmatrix_file) {
#
#   cf_mat <- readr::read_delim(CFmatrix_file,
#                               delim = "\t", escape_double = FALSE,
#                               trim_ws = TRUE,show_col_types = FALSE)
#
#   conflict_free = TRUE
#
#   for (p in 2:ncol(cf_mat)) {
#
#     for (q in (p+1):ncol(cf_mat)) {
#
#       if (q <= ncol(cf_mat)) {
#         oneone = FALSE
#         zeroone = FALSE
#         onezero = FALSE
#         for (r in 1:nrow(cf_mat)) {
#           print(c(r,p,q))
#           print(c(unlist(cf_mat[r, p]),unlist(cf_mat[r, q])))
#           if (cf_mat[r, p] == 1 & cf_mat[r, q] == 1) oneone = TRUE
#           if (cf_mat[r, p] == 0 & cf_mat[r, q] == 1) zeroone = TRUE
#           if (cf_mat[r, p] == 1 & cf_mat[r, q] == 0) onezero = TRUE
#           print(c(oneone,zeroone,onezero))
#         }
#         if (all(oneone,zeroone,onezero)) {
#           conflict_free = FALSE
#           print(conflict_free)
#           return(conflict_free)
#         }
#       }
#
#     }
#   }
#
# }

## light
lightNode <- function(labels,fillcolor,lightNodeKey=NULL,lightNodeColor="red"){

  if (!is.null(lightNodeKey)) {

    for (i in seq_along(lightNodeKey)) {
      key=lightNodeKey[i]

      serch_regrex=paste0("(?<=^|[^A-Za-z0-9])",key,"(?=$|[^A-Za-z0-9])")

      loc_index <- grep(serch_regrex,labels,perl = TRUE)

      if (length(lightNodeColor)==1) {
        fillcolor[loc_index] = lightNodeColor
      } else if (length(lightNodeColor)==length(lightNodeKey)) {
        fillcolor[loc_index] = lightNodeColor[i]
      } else {
        warning("lightNodeColor is dont match to lightNodeKey length!")
        fillcolor[loc_index] = lightNodeColor[1]
      }

    }

  }

  return(fillcolor)

}

## a help function
formatNodelabel <- function(labels,highlight=NULL,highcolor="red"){

  if (!is.null(highlight)) {

    for (i in seq_along(highlight)) {
      key=highlight[i]

      serch_regrex=paste0("(?<=^|[^A-Za-z0-9])",key,"(?=$|[^A-Za-z0-9])")

      if (length(highcolor)==1) {
        replace_regrex=paste0('<FONT COLOR="',highcolor,'">',key,'</FONT> ')
      } else if (length(highcolor)==length(highlight)) {
        replace_regrex=paste0('<FONT COLOR="',highcolor[i],'">',key,'</FONT> ')
      } else {
        warning("highcolor is dont match to highlight length!")
        replace_regrex=paste0('<FONT COLOR="',highcolor[1],'">',key,'</FONT> ')
      }

      labels=gsub(serch_regrex,replace_regrex,labels,perl = TRUE)

    }

  }
  labels=paste0("<",labels,">")
  labels=gsub("\\|",'<br/>',labels)
  return(labels)
}
