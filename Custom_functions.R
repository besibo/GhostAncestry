# Custom_functions.R

# ---- Needed to create the Observation file ----

# This function counts the number of distinct bases found in a (vector of) string(s) of bases.
# By default, "Ns" are removed before counting and all bases are separated by the same empty string
n_distinct_bases <- function(string, remove_N = TRUE, sep = "") {
  
  require(tidyverse)
  if (remove_N) string <- str_remove_all(string, "N")
  tmp_str <- str_split(string = string,
                       pattern = sep)
  
  sapply(tmp_str, n_distinct)
}



# ---- Needed to create the Weight file ----

# This function returns a 3-column table (chr, start and end) and a variable number of rows:
# - if an uncallable region falls inside a predefined window, then only one row is returned.
# - if an uncallable region extends over several windows, then it is split in as many rows
# as necessary: one row for each window.
# The function needs :
# chr: name of the scaffold
# start: position of the start of an uncallable region
# end: position of the end of an uncallable region
# window_start: lower bound of the predefined window in which "start" falls
# window_end: lower bound of the predefined window in which "end" falls
add_lines <- function(chr, start, end, window_start, window_end, ..., window_width = 1000) {
  
  if (window_start != window_end) {
    seq_1 <- c(start, seq(window_start, window_end, by = window_width)[-1])
    seq_2 <- c(seq(window_start, window_end, by = window_width)[-1], end)
    bed_tbl <- tibble(chr = chr, start = seq_1, end = seq_2)
  } else {
    bed_tbl <- tibble(chr = chr, start = start, end = end)
  }
  if (start == end) {
    bed_tbl <- tibble(chr = chr, start = start-1, end = end)
  }
  return(bed_tbl)
}



# ---- Needed to create the Weight file ----

# This function replaces the value in the last row
# and column "col" of a tibble "tbl" by "val and returns 
# The updated tibble
replace_last <- function(tbl, val, col = 3) {
  last <- nrow(tbl)
  tbl[last,col] <- tbl[last,2] / val
  tbl
}

