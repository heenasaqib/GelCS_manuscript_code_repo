# Date Nov 14, 2023
# The purpose of this script is to convert diff_sign_pathways_with_score.csv file into a datatable
# (Yes...Super non-generalized...)

library(DT)

options("DT.TOJSON_ARGS" = list(na = "string"))
options(stringsAsFactors = FALSE)

# Setup Directory
# setwd("~/Desktop/School Work/Stanford/Yang Lab/2023 Spring rotation_Heena") #change this line

# Read CSV
df <- read.csv("diff_sign_pathways_with_score.csv", header = TRUE, stringsAsFactors = FALSE)

datatable(df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in M1-like Mø',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 20,
                         lengthMenu = c("10", "25", "50", "100"),
                         columnDefs = list(list(targets = 7,  # Adjust column index as needed
                         render = JS(
                            "function(data, type, row, meta) {
                                if (type === 'sort' || type === 'type') {
                                    return (data === 'Inf' ? 1e6 : data);
                                    // Replace Infinity with a large number for sorting
                                }
                                return data; // For display, return the data as is
                            }"
                         )))),
          editable = TRUE) %>%
  formatRound(columns=c(3,7), digits=4) %>% #formats column style - rounding
  formatSignif(columns=c(2,6), digits=3)    #formats column style - scientific notation
