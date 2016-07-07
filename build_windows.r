library("devtools")
library("git2r")

document(".")  

# Check that working directory is clean
git_repo <- repository(".")
if (any(sapply(status(git_repo), length) != 0)) {
  stop("Working directory is not clean")
}

# Use git to make a more useful version number:
# "Release"."number of commit since release"-"SHA of current commit"
list_commits <- commits(git_repo)
last_release <- rev(tags(git_repo))[1]
release_name <- names(last_release)
n_commits <- which(sapply(list_commits , identical, y = last_release[[1]])) - 1
short_sha <- sub("^(.{7})(.*)$", "\\1", list_commits[[1]]@sha)
build_version <- paste0('"',release_name,".", n_commits,"-", short_sha,'"')

build_year <- format(Sys.Date(), "%Y")

# Update CITATION file
citation  <- readLines("./inst/CITATION")
citation  <- gsub(pattern = "build_version", replace = build_version, citation)
citation  <- gsub(pattern = "build_year", replace = build_year, citation)
writeLines(citation, "./inst/CITATION")

try({
# Perform build
devtools::build(binary = TRUE, args = c('--preclean'))
build()
})

# Revert changes made to CITATION file
citation  <- readLines("./inst/CITATION")
citation  <- gsub(pattern = build_version, replace = "build_version", citation)
citation  <- gsub(pattern = build_year, replace = "build_year", citation)
writeLines(citation, "./inst/CITATION")
