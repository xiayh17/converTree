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

  t_all <- rbind(t_root,t_leaves,link_df)

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
        parent=newparent,
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
treedata2gv<-function(treedata,gvfile="tree.gv",highlight=NULL,highcolor="red",
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

## a help function
formatNodelabel <- function(labels,highlight=NULL,highcolor="red"){
  if (!is.null(highlight)) {

    for (i in seq_along(highlight)) {
      key=highlight[i]
      serch_regrex=paste0("(?<=^|[^A-Za-z0-9])",key,"(?=$|[^A-Za-z0-9])")
      replace_regrex=paste0('<FONT COLOR="',highcolor,'">',key,'</FONT> ')
      labels=gsub(serch_regrex,replace_regrex,labels,perl = TRUE)
    }

  }
  labels=paste0("<",labels,">")
  labels=gsub("\\|",'<br/>',labels)
  return(labels)
}
