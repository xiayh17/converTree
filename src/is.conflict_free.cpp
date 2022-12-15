#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
bool isConflictFree(NumericMatrix cf_mat) {
  int nrow = cf_mat.nrow(), ncol = cf_mat.ncol();

  for (int p = 0; p < ncol; ++p) {
    for (int q = p + 1; q < ncol; ++q) {
      bool oneone = false, zeroone = false, onezero = false;
      for (int r = 0; r < nrow; ++r) {
        if (cf_mat(r, p) == 1 && cf_mat(r, q) == 1) oneone = true;
        if (cf_mat(r, p) == 0 && cf_mat(r, q) == 1) zeroone = true;
        if (cf_mat(r, p) == 1 && cf_mat(r, q) == 0) onezero = true;
        if (oneone && zeroone && onezero) return false;
      }
    }
  }

  return true;
}
