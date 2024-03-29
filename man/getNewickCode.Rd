% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getNewickCode}
\alias{getNewickCode}
\title{Converts a tree given as lists of children to the Newick tree format}
\usage{
getNewickCode(list, root)
}
\arguments{
\item{list}{A list. the list of children of tree.}

\item{root}{A single integer. The number of mutations.}
}
\value{
the string of newick format tree.
}
\description{
Converts a tree given as lists of children to the Newick tree format
}
