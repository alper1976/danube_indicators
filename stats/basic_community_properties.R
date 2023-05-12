#' @title Determine basic community properties
#'
#' @description This function determines different community properties from phyloseq object
#' @param phylo_seq_object A phyloseq object with an OTU table and taxonomy table.
#' @details
#' Community properties
#' @return A cumulative sum plot.
#' @return A prevalence plot.
#' @return A prevalence histogram.

#' @return A prevalence total counts plot.


#' @keywords phyloseq
#' @export
#' @examples
#' basic_community_properties()

basic_community_properties <- function(phylo_seq_object){


	tdt = data.table(tax_table(phylo_seq_object),
	                 TotalCounts = taxa_sums(phylo_seq_object),
	                 OTU = taxa_names(phylo_seq_object))

	taxcumsum = tdt[, .N, by = TotalCounts]
	setkey(taxcumsum, TotalCounts)
	taxcumsum[, CumSum := cumsum(N)]

	# Prevalence
	mdt = fast_melt(phylo_seq_object)
	prevdt = mdt[, list(Prevalence = sum(count > 0),
	                    TotalCounts = sum(count)),
	             by = TaxaID]
	prevcumsum = prevdt[, .N, by = Prevalence]
	setkey(prevcumsum, Prevalence)
	prevcumsum[, CumSum := cumsum(N)]

	dim(prevdt)

	ordered_prevdt = prevdt[order(prevdt$Prevalence, decreasing = TRUE),]

	data_taxa = as.data.frame(tax_table(phylo_seq_object))
	data_taxa$TaxaID = rownames(data_taxa)

	ordered_prevdt_taxa = dplyr::full_join(ordered_prevdt, data_taxa, by = "TaxaID")
	ordered_prevdt_taxa[1:20,]

	abundorder_prevdt_taxa = ordered_prevdt_taxa[order(ordered_prevdt_taxa$TotalCounts, decreasing = TRUE),]
	abundorder_prevdt_taxa[1:20,]

	## get most diverse groups

	most_diverse <- summary(tax_table(phylo_seq_object))

	## get number of sequences not assigned to genera and family
	genus_na_reads <- abundorder_prevdt_taxa %>% group_by(Genus) %>% summarize(n()) %>% filter(is.na(Genus))
	genus_na_reads$'n()'/sum(abundorder_prevdt_taxa$TotalCounts, na.rm = TRUE)
	family_na_reads <- abundorder_prevdt_taxa %>% group_by(Family) %>% summarize(n()) %>% filter(is.na(Family))
	family_na_reads$'n()'/sum(abundorder_prevdt_taxa$TotalCounts, na.rm = TRUE)
	return(list(most_diverse, genus_na_reads, family_na_reads))
}
