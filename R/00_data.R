# 00_data.R -- load and prepare the Venice data
# Run from the repository root. Creates: venice, years, V2, N, tt

library(ismev)
library(dplyr)

if (!exists("year_range")) year_range <- c(1931, 1981)  # years analysed in the paper; NULL = all rows

if (file.exists("data/venice.rda")) load("data/venice.rda") else data(venice)
venice <- as.data.frame(venice)

# find the year column (any capitalisation); fall back to the first column
yr_col <- names(venice)[tolower(names(venice)) == "year"][1]
if (is.na(yr_col)) yr_col <- names(venice)[1]

if (!is.null(year_range)) {
  venice <- venice %>% filter(.data[[yr_col]] %in% year_range[1]:year_range[2])
}

years <- venice[[yr_col]]                                  # year labels
V2    <- venice %>% select(-all_of(yr_col)) %>% select(1:10)  # 10 largest values per year
N     <- nrow(V2)                                          # number of years (51 for 1931-1981)
tt    <- seq_len(N)                                        # time index (trend covariate)

# The data are rounded to the nearest cm, so adjacent order statistics are often tied
# (87 ties in 1931-1981). This is why the GoF analysis perturbs the data (see 02_gof_bootstrap.R).
n_ties <- sum(V2[, 1:9] == V2[2:10], na.rm = TRUE)
