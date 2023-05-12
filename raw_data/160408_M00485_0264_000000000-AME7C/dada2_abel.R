#!/usr/bin/env Rscript

# sequencing run
# 160408_M00485_0264_000000000-AME7C

##load packages
if (!require("optparse")) {
   install.packages("optparse", dependencies = TRUE)
   library(optparse)
   }

if (!require("dada2")) {
   source("https://bioconductor.org/biocLite.R"); biocLite("dada2")
   library(dada2); packageVersion("dada2")
   }

if (!require("phyloseq")) {
   source("https://bioconductor.org/biocLite.R"); biocLite("phyloseq", dependencies = TRUE)
   library(phyloseq); packageVersion("phyloseq")
   }

if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2); packageVersion("ggplot2")
   }

if (!require("vegan")) {
   install.packages("vegan", dependencies = TRUE)
   library(vegan)
   }

## set parameters

option_list = list(

    make_option(c("-tf", "--truncate_length_fwd"),
        type="integer",
        default=200,
        help="Truncate length of forward reads [default %default]",
        metavar="number"),

    make_option(c("-tr", "--truncate_length_rev"),
        type="integer",
        default=160,
        help="Truncate length of reverse reads [default %default]",
        metavar="number"),

    make_option(c("-tq", "--truncate_quality"),
        type="integer",
        default=10,
        help="Truncates based on quality score [default %default]",
        metavar="character"),

    make_option(c("-p", "--path"),
        type="character",
        default="/work/users/alexaei/output/160408_M00485_0264_000000000-AME7C/AdaptersRemoved",
        help="path to input",
        metavar="character"),

    make_option(c("-o", "--output"),
        type="character",
        default="/work/users/alexaei/output/160408_M00485_0264_000000000-AME7C",
        help="path to output",
        metavar="character"),

    make_option(c("-s", "--path_to_silva"),
        type="character",
        default="/work/users/alexaei/databases/silva/silva_nr_v132_train_set.fa.gz",                      ### needs to be fixed if necessary
        help="path to database",
        metavar="character"),

    make_option(c("-n", "--max_n"),
        type="integer",
        default=0,
        help="maximum of ambiguous bases (Ns) [default %default]",
        metavar="number")


    )

opt = parse_args(OptionParser(option_list=option_list))


    #####################################################################
    # clean and match
    #####################################################################

    ## Removes all sequences containing any ambiguous bases (N).
    ## Tries to match up all pairs of forward and reverse reads by ID,
    ## Removes all sequences that have been orphaned by filtering.
    ## Creates new directory for clean and matched files

    path_to_input = opt$path
    path_to_output = opt$output

    path_to_fwd = list.files(path_to_input, pattern = "R1_001.fastq.gz", all.files = TRUE, full.names = TRUE, recursive = FALSE)
    path_to_rev = list.files(path_to_input, pattern = "R2_001.fastq.gz", all.files = TRUE, full.names = TRUE, recursive = FALSE)

    fastq_fwd <- basename(path_to_fwd)
    fastq_rev <- basename(path_to_rev)
    filtered_path <- file.path(path_to_output, "filtered")
    filtered_fwd <- file.path(filtered_path, fastq_fwd)
    filtered_rev <- file.path(filtered_path, fastq_rev)

    ##Plot quality profiles


  plot_fwd_qual <- plotQualityProfile(file.path(path_to_fwd)[1:10])
  pdf(file.path(path_to_output,"plot_fwd_qual.pdf"))
  plot_fwd_qual
  dev.off()

  plot_rev_qual <- plotQualityProfile(file.path(path_to_rev)[1:10])
  pdf(file.path(path_to_output, "plot_rev_qual.pdf"))
  plot_rev_qual
  dev.off()


  if(length(path_to_fwd) != length(path_to_fwd))
        # DevNote: Returns error to stderr
        write("The number of input files looks dodgy", stderr())

   # write(file.path(path_to_fwd), stdout())
   # write(file.path(filtered_fwd), stdout())
   # write(file.path(path_to_rev), stdout())
   # write(file.path(filtered_rev), stdout())

    # Filtering: these parameters need to be set by user
    out <- filterAndTrim(fwd = file.path(path_to_fwd), filt = file.path(filtered_fwd),
                         rev = file.path(path_to_rev), filt.rev = file.path(filtered_rev),
                         truncLen = c(opt$truncate_length_fwd, opt$truncate_length_rev),
                         maxEE = 2,
                         truncQ = opt$truncate_quality,
                         maxN = opt$max_n,
                         rm.phix=TRUE,
                         compress=TRUE,
                         verbose=TRUE,
                         multithread=FALSE)

    if (length(path_to_fwd)*2 != length(list.files(filtered_path))) {
        write(paste0("Filtering did not result in any sequences from some files.\n",
                     "The number of filtered sequence files is ",
                     length(list.files(filtered_path)),
                     " \n",
                     "while the number of raw files is ",
                     length(path_to_fwd) + length(path_to_rev),
                     " \n"),
                     stdout())
    }

  # Infer Sequence variance

    filtered_fwd = list.files(filtered_path, pattern = "R1_001.fastq.gz", all.files = TRUE, full.names = TRUE, recursive = TRUE)
    filtered_rev = list.files(filtered_path, pattern = "R2_001.fastq.gz", all.files = TRUE, full.names = TRUE, recursive = TRUE)

    sample.names <- basename(filtered_fwd)  # Assumes filename = samplename_XXX.fastq.gz
    sample.names.rev <- basename(filtered_rev) # Assumes filename = samplename_XXX.fastq.gz
    names(filtered_fwd) <- sample.names
    names(filtered_rev) <- sample.names
    set.seed(100)

    # Learn forward error rates
    error_fwd <- learnErrors(filtered_fwd, nreads=1e7, randomize = TRUE, multithread=FALSE)
    # Learn reverse error rates
    error_rev <- learnErrors(filtered_rev, nreads=1e7, randomize = TRUE, multithread=FALSE)

    # Sample inference and merger of paired-end reads
    mergers <- vector("list", length(sample.names))
    names(mergers) <- sample.names
    dada_fwd <- vector("list", length(sample.names))
    dada_rev <- vector("list", length(sample.names))

    for(sam in sample.names) {
      cat("Processing:", sam, "\n")
        derepF <- derepFastq(filtered_fwd[[sam]])
        ddF <- dada(derepF, err = error_fwd, multithread=FALSE)
        dada_fwd[[sam]] <- ddF
        derepR <- derepFastq(filtered_rev[[sam]])
        ddR <- dada(derepR, err = error_rev, multithread=FALSE)
        dada_rev[[sam]] <- ddR
        merger <- mergePairs(ddF, derepF, ddR, derepR)
        mergers[[sam]] <- merger
    }
    rm(derepF); rm(derepR); rm(merger); rm(ddF); rm(ddR)


    # Construct sequence table and remove chimeras
    seqtab <- makeSequenceTable(mergers)
    write(paste0("Merging resulted in ", ncol(seqtab), " sequences."), stdout())


    if(ncol(seqtab) == 0){
        # DevNote: Returns error to stderr
        write(paste0("Merging did not result in any sequences.\n", "You may consider a different truncation for forward and reverse reads. For now the program will continue running the concatenate option instead."),
        stdout())
        rm(seqtab)
        for(sam in sample.names) {
            cat("Processing:", sam, "\n")
            derepF <- derepFastq(filtered_fwd[[sam]])
            ddF <- dada(derepF, err = error_fwd, multithread=FALSE)
            dada_fwd[[sam]] <- ddF
            derepR <- derepFastq(filtered_rev[[sam]])
            ddR <- dada(derepR, err = error_rev, multithread=FALSE)
            dada_rev[[sam]] <- ddR
            concatenate <- mergePairs(ddF, derepF, ddR, derepR, justConcatenate=TRUE)
            mergers[[sam]] <- concatenate
        }
        seqtab <- makeSequenceTable(mergers)
        write(paste0("Concatenating the sequences resulted in ", ncol(seqtab), " sequences."))
        rm(derepF); rm(derepR); rm(concatenate); rm(ddF); rm(ddR)
        seqtab.final <- seqtab
    } else {
        seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)  ### should be used for plotting taxonomy
        saveRDS(seqtab.nochim, file.path(path_to_output, "seqtab_nochim.rds"))
        seqtab.final <- seqtab.nochim
    }

    saveRDS(seqtab.final, file.path(path_to_output, "seqtab.rds"))

    # Assign taxonomy
    write(paste0("Starting taxonomy analysis."), stdout())

    tax <- assignTaxonomy(seqtab.final, opt$path_to_silva, multithread=FALSE)  ### should be used for plotting taxonomy

    # Write to file
    saveRDS(tax, file.path(path_to_output, "tax_final.rds"))
    # tax = readRDS(file.path(path_to_output, "tax_final.rds"))

    write(paste0("Finalized taxonomy analysis."), stdout())

    #####################################################################
    # track sequences through analysis
    #####################################################################
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out,
                   sapply(mergers, getN),
                   rowSums(seqtab.final))
    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    colnames(track) <- c("input", "filtered", "merged", "final")
    rownames(track) <- sample.names

    # write to file
    write.csv(track, file.path(path_to_output, "track_sequence_stats.csv"))
    write(paste0("Write tracking file to output."), stdout())

    write(paste0("Analysis run to the end."), stdout())

