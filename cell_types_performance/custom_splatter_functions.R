splatPopSimConditionalEffects <- function (key, means.pop, conditions, cde.groupSpecific, sample_heterogeneity) 
{
	condition.list <- unique(conditions)
	if (is.null(cde.groupSpecific)) {
		for (c in condition.list) {
			c.samples <- names(conditions[conditions == c])
			c.de <- key[[paste0("ConditionDE.", c)]]
			if (is.list(means.pop)) {
				for (i in names(means.pop)) {
					if (is.null(sample_heterogeneity)) {
						means.pop[[i]][, c.samples] <- mapply("*", means.pop[[i]][, 
																				  c.samples], c.de)
					} else {
						for (sample in c.samples) {
							non1Genes <- which(c.de != 1)
							nExcludedGenes <- round(length(non1Genes) * sample_heterogeneity)
							excludedGenes <- sample(non1Genes, nExcludedGenes)
							c.de.sample <- c.de
							c.de.sample[excludedGenes] <- 1.0
							means.pop[[i]][,sample] <- means.pop[[i]][, sample] * c.de.sample
						}
					}
					means.pop[[i]][means.pop[[i]] < 0] <- 1e-05
				}
			}
			else {
				means.pop[, c.samples] <- mapply("*", means.pop[, 
																c.samples], c.de)
				means.pop[means.pop <= 0] <- 1e-05
			}
		}
	} else{
		for (c in condition.list) {
			c.samples <- names(conditions[conditions == c])
			c.de <- key[[paste0("ConditionDE.", c)]]
			if (is.list(means.pop)) {
				for (i in names(means.pop)) {
					c.de.group = c.de - 1
					c.de.group = c.de.group * cde.groupSpecific[[i]]
					c.de.group = c.de.group + 1
					if (is.null(sample_heterogeneity)) {
						means.pop[[i]][, c.samples] <- mapply("*", means.pop[[i]][, 
																				  c.samples], c.de.group)
					} else {
						for (sample in c.samples) {
							non1Genes <- which(c.de.group != 1)
							nExcludedGenes <- round(length(non1Genes) * sample_heterogeneity)
							excludedGenes <- sample(non1Genes, nExcludedGenes)
							c.de.group.sample <- c.de.group
							c.de.group.sample[excludedGenes] <- 1.0
							means.pop[[i]][,sample] <- means.pop[[i]][, sample] * c.de.group.sample
						}
					}

					means.pop[[i]][means.pop[[i]] < 0] <- 1e-05
				}
			}
			else {
				means.pop[, c.samples] <- mapply("*", means.pop[, 
																c.samples], c.de)
				means.pop[means.pop <= 0] <- 1e-05
			}
		}
	}

	return(means.pop)
}

splatPopSimulateMeans <- function (vcf = mockVCF(), params = newSplatPopParams(nGenes = 1000), 
								   verbose = TRUE, key = NULL, gff = NULL, eqtl = NULL, means = NULL,
								   cde.groupSpecific = NULL, sample_heterogeneity = NULL,
								   ...) 
{
	set.seed(getParam(params, "seed"))
	nGroups <- getParam(params, "nGroups")
	quant.norm <- getParam(params, "pop.quant.norm")
	vcf <- splatter:::splatPopParseVCF(vcf, params)
	group.names <- paste0("Group", seq_len(nGroups))
	samples <- colnames(VariantAnnotation::geno(vcf)$GT)
	conditions <- splatter:::splatPopDesignConditions(params, samples)
	if (is.null(key)) {
		if (is.null(eqtl) || is.null(means)) {
			if (verbose) {
				message("Simulating data for genes in GFF...")
			}
			key <- splatter:::splatPopParseGenes(params, gff)
		}
		else {
			if (verbose) {
				message("Using base gene means from data provided...")
			}
			key <- splatter:::splatPopParseEmpirical(vcf = vcf, gff = gff, 
										  eqtl = eqtl, means = means, params = params)
			params <- setParams(params, nGenes = nrow(key))
		}
	}
	else {
		if (verbose) {
			message("Simulating data for genes in key...")
		}
	}
	if (!all(c("meanSampled", "cvSampled") %in% names(key))) {
		key <- splatter:::splatPopAssignMeans(params, key)
	}
	if (!all(c("eQTL.group", "eSNP.ID", "eQTL.EffectSize") %in% 
			 names(key))) {
		key <- splatter:::splatPopeQTLEffects(params, key, vcf)
	}
	if (length(group.names) > 1) {
		key <- splatter:::splatPopGroupEffects(params, key, group.names)
	}
	if (!all(c("eQTL.condition", "ConditionDE.Condition1") %in% 
			 names(key))) {
		key <- splatter:::splatPopConditionEffects(params, key, conditions)
	}
	if (verbose) {
		message("Simulating gene means for population...")
	}
	means.pop <- splatter:::splatPopSimMeans(vcf, key, means)
	if (quant.norm && ncol(means.pop) > 4) {
		means.pop <- splatter:::splatPopQuantNorm(params, means.pop)
		key <- splatter:::splatPopQuantNormKey(key, means.pop)
	}
	eMeansPop <- splatter:::splatPopSimEffects("global", key, conditions, 
									vcf, means.pop)
	if (length(group.names) > 1) {
		eMeansPopq.groups <- list()
		for (id in group.names) {
			eMeansPop.g <- splatter:::splatPopSimEffects(id, key, conditions, 
											  vcf, eMeansPop)
			eMeansPop.g[eMeansPop.g <= 0] <- 1e-05
			eMeansPopq.groups[[id]] <- eMeansPop.g
		}
		eMeansPop <- eMeansPopq.groups
	}
	sim.means <- splatPopSimConditionalEffects(key, eMeansPop, 
											   conditions, cde.groupSpecific,
											   sample_heterogeneity)
	return(list(means = sim.means, key = key, conditions = conditions))
}