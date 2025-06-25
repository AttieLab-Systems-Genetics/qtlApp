# Test pattern matching for plasma metabolites
dataset_names <- c(
  "HC_HF Plasma plasma_metabolite, additive",
  "HC_HF Plasma plasma_metabolite, interactive (Sex)",
  "HC_HF Liver Genes, additive",
  "HC_HF Systemic Clinical Traits, additive"
)

pattern <- "HC_HF.*Plasma.*plasma_metabolite"

for (name in dataset_names) {
  result <- grepl(pattern, name, ignore.case = TRUE)
  cat(sprintf("Dataset: '%s' -> Matches: %s\n", name, result))
}
