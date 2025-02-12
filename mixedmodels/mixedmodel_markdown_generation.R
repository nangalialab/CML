
library(bookdown)

dir.create(path = "../figures/v1/", recursive = T)

# Telomere
# Tel code as html
rmarkdown::render(
  "Telomere_mixed_models.Rmd", 
  output_format = "html_document",
  output_file = paste0("Telomere_mixed_models.html"),
  params = list(code = TRUE, boot = TRUE)
)


# SNV
# SNV code as html
rmarkdown::render(
  "SNV_mixed_models.Rmd", 
  output_format = "html_document",
  output_file = paste0("SNV_mixed_models.html"),
  params = list(code = TRUE, boot = TRUE)
)




