#' \code{compare_codon_usage} Screen out the codons with usage preferences in the two sets of sequences.
#' @param cds_list_one: A DNA sequence to analyze.
#' @param cds_list_two: Another DNA sequence to analyze.
#' @param gcid: A parameter in get_codon_table() to get the required codon table for mapping codons to amino acids.
#' @returns A data frame containing the computed metrics.
#' @examples Comparing two sets of genes in yeast
#' result <- compare_codon_usage(yeast_cds[1:50], yeast_cds[51:100], gcid = "1" )
#' print(result)
library(Biostrings)
library(dplyr)
compare_codon_usage <- function(cds_list_one, cds_list_two, gcid = NULL) {
  
  yeast_cds_qc_one <- check_cds(cds_list_one)
  yeast_cds_qc_two <- check_cds(cds_list_two) 

  # Count the total number of each codon
  yeast_cf_one <- count_codons(yeast_cds_qc_one)
  yeast_cf_two <- count_codons(yeast_cds_qc_two)
  codon_totals_1 <- colSums(yeast_cf_one)
  codon_totals_2 <- colSums(yeast_cf_two)

  # Combined data
  data <- cbind(codon_totals_1, codon_totals_2)
  data_df <- as.data.frame(data)
  data_df$codon <- rownames(data_df)
  rownames(data_df) <- NULL

  # Mapping codons to amino acids
  codon_table <- get_codon_table(gcid = gcid)
  data_df <- merge(data_df, codon_table)

  # Filter out stop codons and single codons
  stop_codons <- c("TAA", "TAG", "TGA")
  single_codon_amino_acids <- c("TGG", "ATG")
  data_tr <- data_df[!(data_df$codon %in% c(stop_codons, single_codon_amino_acids)),]
  
  # Fisher test
  by_amino_acid <- split(data_tr, data_tr$amino_acid)
  results <- data.frame(AminoAcid = character(),
                        Codon = character(),
                        Codon_aa_code = character(),
                        Codon_subfam = character(),
                        P_Value = numeric(),
                        stringsAsFactors = FALSE)
  
  for (amino_acid in names(by_amino_acid)) {
    codon_group <- by_amino_acid[[amino_acid]]

    for (i in 1:nrow(codon_group)) {
      row <- codon_group[i, ]
      codon_total_1 <- as.numeric(row[['codon_totals_1']])
      codon_total_2 <- as.numeric(row[['codon_totals_2']])
      other_codons_totals_1 <- sum(codon_group[['codon_totals_1']]) - codon_total_1
      other_codons_totals_2 <- sum(codon_group[['codon_totals_2']]) - codon_total_2

      # Form the 2x2 contingency table
      contingency_table <- matrix(c(codon_total_1, other_codons_totals_1,
                                    codon_total_2, other_codons_totals_2), nrow = 2, byrow = TRUE)

      # Check for non-positive entries explicitly
      if (any(contingency_table <= 0, na.rm = TRUE)) {
        cat("Non-positive entry detected in contingency table for amino acid", amino_acid, "and codon", row['codon'], "\n")
        print(contingency_table)
        next # Skip this codon
      }
      fisher_result <- fisher.test(contingency_table)
      
      # Append the results
      results <- rbind(results, data.frame(AminoAcid = amino_acid,
                                           Codon = as.character(row[['codon']]),
                                           Codon_aa_code = as.character(row[['aa_code']]),
                                           Codon_subfam = as.character(row[['subfam']]),
                                           P_Value = fisher_result$p.value,
                                           stringsAsFactors = FALSE))
      }
    }
  
  # Apply BH correction
  results$P_Adjusted_BH <- p.adjust(results$P_Value, method = "BH")

  return(results)
}

