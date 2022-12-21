#' Convert treedata to CFMatrix
#'
#' Convert treedata to CFMatrix, storage result to CFmatrix file and return a datafram
#'
#' @param treedata a treedata format dataframe
#' @param filename file name to save CFMatrix, default is `test.CFMatrix`
#'
#' @importFrom utils write.table
#' @importFrom stringr str_sort str_split
#'
#' @return a dataframe and a file
#' @export
#'
#' @examples
#' CFmatrix_file = system.file("extdata", "ground_truth_tree.CFMatrix", package = "converTree")
#' tre_dat <- cf2treedata(CFmatrix_file)
#' cf_mat<-treedata2cf(tre_dat)
#' filename="test.CFMatrix"
#' write.table(cf_mat,file = filename,quote = FALSE,sep = "\t",row.names = FALSE)
treedata2cf <- function(treedata,filename="test.CFMatrix") {

  tre_dat = as.data.frame(treedata)
  ## to CFmatrix
  root_node <- as.numeric(na.omit(tre_dat[tre_dat$parent == tre_dat$node,"node"]))
  leaves_dat <- tre_dat[tre_dat$node < root_node,]
  mut_dat <- tre_dat[tre_dat$node > root_node,]

  cells <- leaves_dat$label
  ncell <- length(cells)
  mutations <- mut_dat$label
  main_mut <- gsub("\\|.*","",mut_dat$label)
  nmmut = length(main_mut)

  ## matrix from main mutations
  ## init matrix
  init_cf_mat <- matrix(data = 0,nrow = ncell,ncol = nmmut)
  colnames(init_cf_mat) <- main_mut
  rownames(init_cf_mat) <- cells
  tre_dat2 <- tre_dat
  tre_dat2$label <- gsub("\\|.*","",tre_dat2$label)

  ## find mutatins for each cell
  # tre_dat[tre_dat$label,]
  findParent <- function(node) {
    p<-tre_dat2[tre_dat2$node == node,"parent"]
    p <- as.numeric(p)
    return(p)
  }

  for (i in seq_along(cells)) {

    cell=cells[i]
    # print(cell)
    ktmp=c()
    k0=as.numeric(na.omit(leaves_dat[leaves_dat$label==cell,"node"]))
    # print(k0)
    k1=findParent(k0)
    # print(k1)
    cell_tmp=c()

    if (!k1 %in% ktmp) {
      ktmp=c(ktmp,k0)
      mut_tmp=c()
      k=k0
      while (k!=root_node) {
        k=findParent(k)
        # print(k)
        if (k!=root_node) {
          new_tmp=as.character(na.omit(tre_dat2[tre_dat2$node==k,"label"]))
          # print(new_tmp)
          mut_tmp=c(mut_tmp,new_tmp)
        }
      }
      # print(mut_tmp)
      init_cf_mat[cell,mut_tmp] = 1
      cell_tmp=c(cell_tmp,cell)
    } else {
      ## for cells have the same parent
      ## skip for they have the same value
      k_dup = tre_dat2[tre_dat2$label==cell,"node"]
      cell_dups = tre_dat2[tre_dat2$node==k_dup,"label"]
      cell_dup = cell_dups[cell_dups %in% cell_tmp]
      init_cf_mat[cell,] = init_cf_mat[cell_dup,]
    }


  }

  ## return duplicate mutations
  # main_mut
  all_mut <- mut_dat$label

  for (i in seq_along(main_mut)) {

    mut = main_mut[i]
    pat <- paste0("^",mut,"\\|")
    # print(pat)
    muts <- all_mut[grep(pat,all_mut,perl = TRUE)]
    # print(muts)
    res <- stringr::str_split(muts,"\\|")
    if (length(res)==0) {
      res <- mut
    } else {
      res <- res[[1]]
    }
    a_mut = setdiff(res,mut)
    new_m <- matrix(rep(0,length(a_mut)*ncell), nrow = ncell, ncol = length(a_mut))
    colnames(new_m) <- a_mut
    init_cf_mat <- cbind(init_cf_mat,new_m)
    init_cf_mat[,a_mut] = init_cf_mat[,mut]


  }

  ## sort cell and mut
  init_cf_mat <- init_cf_mat[stringr::str_sort(rownames(init_cf_mat),numeric = TRUE),stringr::str_sort(colnames(init_cf_mat),numeric = TRUE)]

  ## rowname to first col
  cell_column <- matrix(data = rownames(init_cf_mat),nrow = ncell,ncol = 1,dimnames = list(NULL,"cellIDxmutID"))
  cf_mat_g <- cbind(cell_column,init_cf_mat)
  cf_mat_g <- as.data.frame(cf_mat_g)
  utils::write.table(cf_mat_g,file = filename,quote = FALSE,sep = "\t",row.names = FALSE)
  return(cf_mat_g)

}

#' Convert nwk to CFMatrix
#'
#' @param nwk the string of newick format tree.
#' @param filename file name to save CFMatrix, default is `test.CFMatrix`
#'
#' @return a dataframe and a file
#' @export
#'
#' @examples
#' nwk_text="((((((((((8)7,(1)12)9)10)2)4)3,5,(((17)11)16)14)15)6)13)18;"
#' nwk2cf(nwk=nwk_text)
nwk2cf <- function(nwk,filename="test.CFMatrix") {

  tre <- treeio::read.newick(text = nwk)
  stre <- stack_mutationTree(tre)
  tre_dat <- stre %>% treeio::as_tibble()
  tre_dat <- tre_dat[,-3]
  cf <- treedata2cf(tre_dat,filename)
  return(cf)

}
