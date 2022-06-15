#' Converts a tree given as lists of children to the Newick tree format
#'
#' @param list A list. the list of children of tree.
#' @param n A single integer. The number of mutations.
#'
#' @return
#' the string of newick format tree.
#' @export
#'
#' @examples
#' p = c(11, 2, 3, 14, 14, 16, 8, 6, 9, 1, 15, 8, 10, 14, 5, 13, 17)
#' cl = getChildListFromParentVector(p,16)
#' childList2nwk(list = cl,n = 16)
childList2nwk <- function(list,n) {

  tmp <- getNewickCode(list = list,root=n)
  tmp <- paste0(tmp,";")
  return(tmp)

}

#' Convert a parent vector format to the newick format tree
#'
#' @param parents A integer vector. Parent vector format tree.
#' @param n A single integer. The number of mutations.
#'
#' @return
#' the string of newick format tree.
#' @export
#'
#' @examples
#' p = c(11, 2, 3, 14, 14, 16, 8, 6, 9, 1, 15, 8, 10, 14, 5, 13, 17)
#' parentVector2nwk(parents=p,n = 16)
parentVector2nwk <- function(parents,n) {
  childList <- getChildListFromParentVector(parents=parents,n=n)
  nwk <- childList2nwk(list = childList, n = n)
  return(nwk)
}

#' Converts a tree given as parent vector format to an ancestor matrix
#'
#' @param parents A integer vector. Parent vector format tree.
#' @param n A single integer. The number of mutations.
#'
#' @return
#' a logical matrix;0 and 1
#' @export
#'
#' @examples
#' p = c(11, 2, 3, 14, 14, 16, 8, 6, 9, 1, 15, 8, 10, 14, 5, 13, 17)
#' parentVector2anc(parents=p,n = 16)
parentVector2anc <- function(parents,n) {

  ancBool <- parentVector2ancMatrix(parents=parents,n=n)
  ancBool[isTRUE(ancBool)] = 1
  return(ancBool)

}

#' Converts an ancestor matrix to a tree in parent vector format
#'
#' @param anc A matrix. An ancestor matrix.
#' @param n A single integer. The number of mutations.
#'
#' @return
#' a integer vector; parent vector format tree.
#' @export
#'
#' @examples
#' p = c(11, 2, 3, 14, 14, 16, 8, 6, 9, 1, 15, 8, 10, 14, 5, 13, 17)
#' anc <- parentVector2anc(parents=p,n = 16)
#' anc2parentVector(anc=anc,n=16)
anc2parentVector <- function(anc,n) {

  pv <- ancMatrixToParVector(anc=anc,n=n)
  return(pv)

}

#' Converts a tree given as parent vector format to an ancestor matrix
#'
#' @param anc A matrix. An ancestor matrix.
#' @param n A single integer. The number of mutations.
#'
#' @return
#' a string of newick format tree.
#' @export
#'
#' @examples
#' p = c(11, 2, 3, 14, 14, 16, 8, 6, 9, 1, 15, 8, 10, 14, 5, 13, 17)
#' anc <- parentVector2anc(parents=p,n = 16)
#' anc2nwk(anc=anc,n=16)
anc2nwk <- function(anc,n){
  pv <- ancMatrixToParVector(anc=anc,n=n)
  nwk <- parentVector2nwk(parents = pv,n=n)
  return(nwk)
}
