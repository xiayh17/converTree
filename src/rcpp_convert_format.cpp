#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;

//
// These function using Rcpp that convert tree format
// modified from https://github.com/cbg-ethz/infSCITE..
//

//' Convert a parent vector format to the list of children
//'
//' @param parents A integer vector. Parent vector format tree.
//' @param n A single integer. The number of mutations.
//' @export
//' @return
//' the list of children of tree.
// [[Rcpp::export]]
std::vector<std::vector<int> > getChildListFromParentVector(IntegerVector parents, int n){

  std::vector<std::vector<int> > childList(n+1);
  for(int i=0; i<n; i++){
    //cout << "child lists at " << parents[i] << " out of 0 to " << n << "\n";
    childList.at(parents[i]).push_back(i);
  }
  return childList;
}

//' Converts a tree given as lists of children to the Newick tree format
//'
//' @param list A list. the list of children of tree.
//' @param root A single integer. The number of mutations.
//' @return
//' the string of newick format tree.
// [[Rcpp::export]]
std::string getNewickCode(std::vector<std::vector<int> > list, int root){
  std::stringstream newick;
  std::vector<int> rootChilds = list.at(root);
  if(!rootChilds.empty()){
    newick << "(";
    bool first = true;
    for(int i=0; i<rootChilds.size(); i++){
      if(!first){
        newick << ",";
      }
      first = false;
      newick << getNewickCode(list, rootChilds.at(i));
    }
    newick << ")";
  }
  newick << root+1;
  return newick.str();
}

// [[Rcpp::export]]
LogicalMatrix init_boolMatrix(int n, int m, bool value){

  LogicalMatrix matrix = LogicalMatrix(n, m);     // allocate

  for (int i=0; i<n; ++i)             // initialize
  {
    for (int j=0; j<m; ++j)
    {
      matrix(i,j) = value;
    }
  }
  return matrix;
}

//' Converts a tree given as parent vector format to an ancestor matrix
//'
//' @param parents A integer vector. Parent vector format tree.
//' @param n A single integer. The number of mutations.
//' @return
//' a logical matrix;TRUE and FALSE
// [[Rcpp::export]]
LogicalMatrix parentVector2ancMatrix(IntegerVector parents, int n){
  LogicalMatrix ancMatrix = init_boolMatrix(n, n, false);
  int root = n;
  for(int i=0; i<n; i++){
    int anc = i;
    int its =0;
    while(anc < root){                              // if the ancestor is the root node, it is not represented in the adjacency matrix
      if(parents[anc]<n){
        ancMatrix(parents[anc],i) = true;
      }

      anc = parents[anc];
      its++;
    }
  }
  for(int i=0; i<n; i++){
    ancMatrix(i,i) = true;
  }
  return ancMatrix;
}

//' Converts a an ancestor matrix to tree in parent vector format
//'
//' @param anc A matrix. An ancestor matrix.
//' @param n A single integer. The number of mutations.
//' @return
//' a integer vector; parent vector format tree.
// [[Rcpp::export]]
IntegerVector ancMatrixToParVector(LogicalMatrix anc, int n){
  IntegerVector parVec = IntegerVector(n);
  for(int i=0; i<n; i++){
    parVec(i) = n;
  }
  for(int i=0; i<n; i++){
    for(int k=0; k<n; k++){
      if(k!=i && anc(k,i)==true){  // k is true ancestor of i
        bool cand = true;
        for(int l=0; l<n; l++){
          if(l!=i && l!=k && anc(l,i)==true && anc(k,l)==true){   // k is ancestor of l, and l is ancestor of i
            cand = false;                                        // k is no longer candidate for being parent of i
            break;
          }
        }
        if(cand==true){           // no ancestor of i is descendant of k -> k is parent of i
          parVec(i) = k;
        }
      }

    }
  }
  return parVec;
}
