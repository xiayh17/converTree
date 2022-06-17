## add license
usethis::use_mit_license("YongheXia")

## unload dll file
.onUnload <- function (libpath) {
  library.dynam.unload("converTree", libpath)
}



## update document
devtools::document(".")

## load all
devtools::load_all(".")

## build
devtools::build(".")

## add dependecies
# pak::pkg_install("attachment")
attachment::att_amend_desc()

## update readme
devtools::build_readme()

## more document
# usethis::use_pkgdown()

## push to github
# usethis::use_github()
