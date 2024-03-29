% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/improveCore.R
\name{childList2nwk}
\alias{childList2nwk}
\title{Converts a tree given as lists of children to the Newick tree format}
\usage{
childList2nwk(list, n)
}
\arguments{
\item{list}{A list. the list of children of tree.}

\item{n}{A single integer. The number of mutations.}
}
\value{
the string of newick format tree.
}
\description{
Converts a tree given as lists of children to the Newick tree format
}
\examples{
p = c(11, 2, 3, 14, 14, 16, 8, 6, 9, 1, 15, 8, 10, 14, 5, 13, 17)
cl = getChildListFromParentVector(p,16)
childList2nwk(list = cl,n = 16)
}
