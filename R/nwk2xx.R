#' Convert the newick format tree to an ancestor matrix
#'
#' @param nwk the string of newick format tree.
#'
#' @importFrom dplyr select pull
#' @importFrom treeio as_tibble read.tree
#'
#' @return
#' the string of newick format tree.
#' @export
#'
#' @examples
#' p = c(11, 2, 3, 14, 14, 16, 8, 6, 9, 1, 15, 8, 10, 14, 5, 13, 17)
#' nwk <- parentVector2nwk(parents=p,n = 16)
#' nwk2anc(nwk)
nwk2anc <- function(nwk) {

  nwk = read.tree(text = nwk)

  ## nwk to data
  tree_data <- nwk %>% as_tibble() %>%
    select(parent,node,label)

  ## filter root
  noroot_data <- tree_data[!tree_data$parent == tree_data$node,]

  ## replace by labels
  test_fact <- with(noroot_data,factor(label,
                                       levels = unique(label),
                                       labels = unique(label)))

  ## init  matrix
  dims <- pull(noroot_data,node) %>% length()
  fm = matrix(data = 0, nrow = dims,ncol = dims)
  diag(fm) = 1

  ## retrieve tree to matrix
  for (i in 1:nrow(noroot_data)) {

    # print(i)
    x = noroot_data[i,1] %>% as.numeric() %>% test_fact[.]
    y = noroot_data[i,2] %>% as.numeric() %>% test_fact[.]
    fm[x,y] <- 1

  }

  ## return matrix
  return(fm)

}

#' Convert the newick format tree to a parent vector format
#'
#' @param nwk the string of newick format tree.
#' @param n A single integer. The number of mutations.
#'
#' @return
#' a integer vector; parent vector format tree.
#' @export
#'
#' @examples
#' p = c(11, 2, 3, 14, 14, 16, 8, 6, 9, 1, 15, 8, 10, 14, 5, 13, 17)
#' nwk = parentVector2nwk(parents=p,n = 16)
#' nwk2parentVector(nwk,n=16)
nwk2parentVector <- function(nwk,n) {

  ##to matrix
  tm <- nwk2anc(nwk)

  ##matrix to pv
  pv <- anc2parentVector(anc=tm,n=n)

  return(pv)

}
