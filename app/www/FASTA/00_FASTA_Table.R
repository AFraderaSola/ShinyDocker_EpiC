species <- sort(unique(c(names(DNA_5mCG_quantified_proteins_df), names(DNA_5mCHG_quantified_proteins_df),
                  names(RNA_Nuclear_quantified_proteins_df), names(RNA_Nuclear_quantified_proteins_df),
                  names(Histone_H3K4me3_quantified_proteins_df), names(Histone_H3K9me2_quantified_proteins_df))))

DNA_5mCG <- gsub(pattern = FALSE, 
             replacement = "",
             x = gsub(pattern = TRUE, 
                      replacement = "X", 
                      x = species %in% names(DNA_5mCG_quantified_proteins_df)))

DNA_5mCHG <- gsub(pattern = FALSE, 
              replacement = "",
              x = gsub(pattern = TRUE, 
                       replacement = "X", 
                       x = species %in% names(DNA_5mCHG_quantified_proteins_df)))

RNA_nuclear <- gsub(pattern = FALSE, 
                    replacement = "",
                    x = gsub(pattern = TRUE, 
                             replacement = "X", 
                             x = species %in% names(RNA_Nuclear_quantified_proteins_df)))

RNA_total <- gsub(pattern = FALSE, 
                    replacement = "",
                    x = gsub(pattern = TRUE, 
                             replacement = "X", 
                             x = species %in% names(RNA_Total_quantified_proteins_df)))

Histone_H3K4me3 <- gsub(pattern = FALSE, 
                  replacement = "",
                  x = gsub(pattern = TRUE, 
                           replacement = "X", 
                           x = species %in% names(Histone_H3K4me3_quantified_proteins_df)))

Histone_H3K9me2 <- gsub(pattern = FALSE, 
                        replacement = "",
                        x = gsub(pattern = TRUE, 
                                 replacement = "X", 
                                 x = species %in% names(Histone_H3K9me2_quantified_proteins_df)))

FASTA <- list.files(pattern = ".*fasta")

FASTA_df <- tibble("Species" = species,
                          "DNA - 5mCG" = DNA_5mCG,
                          "DNA - 5mCHG" = DNA_5mCHG,
                          "RNA - m6A Nuclear" = RNA_nuclear,
                          "RNA - m6A Total" = RNA_total,
                          "Histone - H3K4me3" = Histone_H3K4me3,
                          "Histone - H3K9me2" = Histone_H3K9me2,
                          "FASTA file" = FASTA)

save(FASTA_df, file = "00_FASTA.RData")

