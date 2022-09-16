fromList <- function (input) {
  # source: https://github.com/hms-dbmi/UpSetR/issues/85#issuecomment-327900647
  # Same as original UpSetR::fromList()...
  # function creates an upset data frame where each column represents a sample and each row represents an observation.
  # Rows are values 0/1 for presence/absence of the observation within the sample.
  # the behavior of this function differs from UpSetR::fromList() in that this one preserves the names of the observations as row names in the final DF.
  # this is helpful for ComplexUpset plots so that the intersection bar charts can be colored by metadata associated with the observation.
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
      x <- as.vector(match(elements, x))
      }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
 }

read_gather <- function(path){
  # read in output files from sourmash gather.
  # specifies each col type to facilitate binding.
  # coltypes can be mis-read when data is missing. 
  library(readr)
  gather <- read_csv(path, col_types = "ddddddddcccddddcccdc")
}
