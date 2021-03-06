# The exstra_db class
# Includes information on STRs, whether they are disease-causing or of a general nature


# check if the object is of this class
#' @import data.table
#' @import stringr
#' @import xlsx
#' @import testit
#' 
#' @export
is.exstra_db <- function(x) inherits(x, "exstra_db")

# Create a new object of this class (not for the user)
#
exstra_db_new_ <- function(strd, input_type = NULL) {
  # Transforms a data.frame or data.table into a exstra_db object
  if (!is.data.frame(strd)) stop("strd must be data.frame")
  strd <- data.table(strd)
  strd <- strd[!(is.na(chrom) | is.na(chromStart) | is.na(chromEnd))]
  strd$input_order <- seq(1, dim(strd)[1], 1)
  # TODO this should always just be locus, hack away!
  if(!is.null(strd$disease.symbol)) {
    #setkey(strd, "disease.symbol")
    setnames(strd, "disease.symbol", "locus")
  } 
  setkey(strd, "locus")
  structure(list(db = strd, input_type = input_type), class = c("exstra_db"))
}


#' @export
print.exstra_db <- function(x, ...) {
  cat(class(x)[1], " object with ", dim(x$db)[1], " loci ($db) of type ",  x$input_type, "\n",
    sep = "")
}

#' @export
exstra_db_text <- function(file) {
  if (!is.character(file)) stop("file must be character")
  stop("Text exstra_db reading not yet implemented")
}

#' @export
loci.exstra_db <- function(exstra_db) {
  # Give the loci names
  loci <- exstra_db$db[order(input_order), locus]
  assert("Could not identify the str loci", !is.null(loci))
  loci
}

# TODO: mae these into generics

#' @export
loci_text_info.exstra_db <- function(x, locus) {
  # gives text info for the locus, usually used in plot titles
  # TODO: modify this:
  assert("The class of x must be exstra_db", is.exstra_db(x))
  locus.in <- locus
  if(x$input_type == "named") {
    x.info <- x$db[locus.in == locus]
    #TODO: this is wrong
    assert(paste("The locus", locus, "was not found"), dim(x.info)[1] >= 1)
    assert(paste("There were multiple entries for locus", locus), dim(x.info)[1] <= 1)
    rs.len <- with(x.info, nchar(as.character(Repeat.sequence)))
    normal.copyNum <- with(x.info, ifelse(is.null(read_detect_size), floor(copyNum), floor(read_detect_size / rs.len)))
    normal.size.bp <- with(x.info, ifelse(is.null(read_detect_size), floor(copyNum * rs.len), read_detect_size))
    return(with(x.info,  
      paste0(locus, " (", 
        Location.of.repeat.within.gene, " ", Repeat.sequence, ") norm: ", normal.copyNum, 
        " (", normal.size.bp, "bp) , exp: ", rn.unst.low, " (", 
        floor(rn.unst.low * rs.len), "bp)"))
    )
  } else if (x$input_type == "ucsc") {
    x.info <- x$db[locus.in == locus] 
    #TODO: this is wrong
    assert(paste("The locus", locus, "was not found"), dim(x.info)[1] >= 1)
    assert(paste("There were multiple entries for locus", locus), dim(x.info)[1] <= 1)
    rs.len <- with(x.info, nchar(as.character(Repeat.sequence)))
    #normal.copyNum <- with(x.info, ifelse(is.null(read_detect_size), floor(copyNum), floor(read_detect_size / rs.len)))
    #normal.size.bp <- with(x.info, ifelse(is.null(read_detect_size), floor(copyNum * rs.len), read_detect_size))
    return(with(x.info,  
      paste0(locus))
    )
  } else {
    stop("Unrecognised input_type in exstra_db. Got ", x$input_type)
  }
}



#' @export
loci_normal <- function(x, locus) {
  loci_normal_exp (x, locus)[1]
}

#' @export
loci_min_exp <- function(x, locus) {
  # Give the minimum expanded STR in bp
  loci_normal_exp (x, locus)[2]
}

#' @export
`[.exstra_db` <- function(x, fil) {
  assert("locus not the key of x$db", key(x$db)[1] == "locus")
  x$db <- x$db[eval(substitute(fil))]
  x
}

# copy data.table inside
#' @export
copy.exstra_db <- function(x) {
  x$db <- copy(x$db)
  x
}


# I think the following was code that was left over from another time
#Y <- exstra_db_read("/Users/tankard/Documents/Research/repeats/disease_repeats/repeat_disorders.xlsx")
#class(Y)


# TODO method for seqnames(exstra_db)

#TODO length(exstra_db) =  number of loci

#TODO dim(exstra_db) = dim(exstra_db$db)