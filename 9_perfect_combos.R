# Read file
txt <- readLines("plot/panel4_CAT1.txt")
txt_all <- paste(txt, collapse = "\n")

# Split into combo blocks
blocks <- unlist(strsplit(txt_all, "={20,}"))

extract_combo_info <- function(block) {
  
  # Combo number
  combo <- regmatches(block, regexpr("Combo:\\s*\\d+", block))
  combo <- gsub("Combo:\\s*", "", combo)
  
  # Genes (single line or multiline safe)
  genes <- regmatches(block,
                      regexpr("Genes:\\s*[^\n]+", block))
  genes <- gsub("Genes:\\s*", "", genes)
  
  # Count Accuracy = 1 occurrences
  acc <- regmatches(block, gregexpr("Accuracy\\s*:\\s*1(\\s|$)", block))[[1]]
  
  # Keep only if both training & validation are perfect
  if (length(acc) >= 2 && combo != "") {
    return(data.frame(
      Combo = combo,
      Genes = genes,
      stringsAsFactors = FALSE
    ))
  } else {
    return(NULL)
  }
}

# Apply to blocks
result <- do.call(rbind, lapply(blocks, extract_combo_info))

write.table(result,
            file = "plot/perfect_combos_CAT1.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
