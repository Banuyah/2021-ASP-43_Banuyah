# Obtain GPL6244 platform annotation. Prefer the committed local copy so the
# workflow can run without repeatedly downloading the platform table.
suppressPackageStartupMessages(library(GEOquery))

if (!file.exists("data/GPL6244_annotation.csv")) {
  gpl <- GEOquery::getGEO("GPL6244")
  annotation <- GEOquery::Table(gpl)
  dir.create("data", showWarnings = FALSE)
  write.csv(annotation, "data/GPL6244_annotation.csv", row.names = FALSE)
} else {
  annotation <- read.csv(
    "data/GPL6244_annotation.csv",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

cat("GPL6244 annotation rows:", nrow(annotation), "\n")
cat("GPL6244 annotation columns:\n")
print(names(annotation))
