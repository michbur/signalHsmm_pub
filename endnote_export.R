library(dplyr)

bib_lines <- readLines("lokalizom.txt")

starts <- c(grep("@", bib_lines), length(bib_lines) + 1)

bib_names <- sapply(1L:(length(starts) - 1), function(i) {
  pub_lines <- bib_lines[starts[i]:((starts[i + 1]) - 1)]
  author_line <- pub_lines[grep("author = ", pub_lines)]
  if(length(author_line) > 0) {
    author_name <- tolower(strsplit(strsplit(author_line, "= {", fixed = TRUE)[[1]][2], split = ", ")[[1]][1])
    splitted_author <- strsplit(author_name, "")[[1]]
    #removes any non-letters
    final_author <- splitted_author[splitted_author %in% letters] %>% paste0(collapse = "")
  } else {
    final_author <- c()
  }
  
  year_line <- pub_lines[grep("year = ", pub_lines)]
  
  if(length(year_line) > 0) {
    final_year <- strsplit(strsplit(year_line, "= {", fixed = TRUE)[[1]][2], split = "}", fixed = TRUE)[[1]]
  } else {
    final_year <- c()
  }
  
  title_line <- pub_lines[grep("title = ", pub_lines)]
  if(length(title_line) > 0) {
    title_name <- tolower(strsplit(strsplit(title_line, "= {", fixed = TRUE)[[1]][2], split = " ")[[1]][1])
    splitted_title <- strsplit(title_name, "")[[1]]
    final_title <- splitted_title[splitted_title %in% letters] %>% paste0(collapse = "")
  } else {
    final_title <- c()
  }
  
  paste0(final_year, final_author, final_title)
})

bib_names <- paste0(bib_names, ",")
bib_starts <- grep("@", bib_lines)
for(i in 1L:length(bib_starts)) {
  bib_lines[bib_starts[i]] <- paste0(bib_lines[bib_starts[i]], bib_names[i])
}
writeLines(bib_lines, "lokalizom.bib")
