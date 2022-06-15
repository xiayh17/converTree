p = c(11, 2, 3, 14, 14, 16, 8, 6, 9, 1, 15, 8, 10, 14, 5, 13, 17)
getChildListFromParentVector(p,16)
cl = getChildListFromParentVector(p,16)
getNewickCode(cl,16)
tmp = getNewickCode(cl,16)
anc = parentVector2ancMatrix(p,16)
anc[isTRUE(anc)] = 1
ancMatrixToParVector(anc,17)
