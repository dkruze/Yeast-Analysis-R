library(tidyverse)

# The following function takes a data.table object and a string and uses them to export a table into a CSV file
# It takes a string and sets it equal to a global fileName variable, then uses an interpolator to make a filename string out of the given string
# Then, that filename string is used to write a new table using the data from the input table
# Expected input: data.table "tableVar", string "y"
# Expected output: external CSV table "[filename].csv"
fileName <- ""
export <- function(tableVar, y) {
  fileName <- y
  tableName <- stringr::str_interp("${fileName}.csv")
  write.table(tableVar, file = tableName, append = FALSE, sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
}