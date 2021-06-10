rm(list = ls())
#!/usr/bin/envRscript
library("optparse")
cat("\t\t\t****************************************\n")
cat("\t\t\t*                 GESim                *\n")
cat("\t\t\t*      Gene Expression Simulator       *\n")
cat("\t\t\t*             V0.1 - 2021              *\n")
cat("\t\t\t*       author: Ziad Al Bkhetan        *\n")
cat("\t\t\t*      ziad.albkhetan@gmail.com        *\n")
cat("\t\t\t*   https://github.com/ziadbkh/GESim   *\n")
cat("\t\t\t****************************************\n")

option_list <- list(
  make_option(
    c("-i", "--vcf"),
    type = "character",
    default = NULL,
    help = "Haplotype/Genotype file path (.vcf).
    \t\tIf simulations based on haplotypes are reqiured, the alleles must be phased and '|' separated.",
    metavar = "character"
  ),
  
  make_option(
    c("-s", "--snp"),
    type = "character",
    default = NULL,
    help = "SNP file path.\n\t\tIt contain one column with the SNP index (the order of the SNP) in the VCF file.
    \t\tSee the sample files for the format",
    metavar = "character"
  ),
  
  make_option(
    c("--pair_a"),
    type = "character",
    default = NULL,
    help = "SNP pairs file path to be used for  simulations based on additive impact of a SNP pair.
    \t\tIt contains two columns (tab-separated) containing the SNP index (the order of the SNP) in the VCF file.
    \t\tSee the sample files for the format",
    metavar = "character"
  ),
  make_option(
    c("--pair_i"),
    type = "character",
    default = NULL,
    help = "SNP pairs file path to be used for SNP interaction-based simulations.
    \t\tIt contains two columns (tab-separated) containing the SNP index (the order of the SNP) in the VCF file.
    \t\tSee the sample files for the format",
    metavar = "character"
  ),
  make_option(
    c("--hap"),
    type = "character",
    default = NULL,
    help = "Haplotype file path to be used for haplotype-based simulations.
    \t\tIt contains three columns (tab-separated) as follows: the SNP determninng the begining of the haplotype, the SNP determninng the end of the haplotype, then the haplotype stretch used for encoding.
    \t\tSee the sample files for the format.",
    metavar = "character"
  ),
  make_option(
    c("--h2"),
    type = "double",
    default = 0.05,
    help = "heritability value between 0 and 1.
    \t\tIt referes to the proportion of the expression variation caused by the genetic architecture.",
    metavar = "numeric"
  ),
  make_option(
    c("--random"),
    type = "integer",
    default = 0,
    help = "Number of simulations with no causal genetic architecture."
  ),
  make_option(
    c("-o", "--out"),
    type = "character",
    default = "",
    help = "output file path with the prefix of the files names.",
    metavar = "character"
  )
)


opt_parser <- OptionParser(option_list = option_list)

opt <- parse_args(opt_parser)

start_time <- Sys.time()

cat(paste("Analysis started at:", start_time), quote = F)


if (is.null(opt$vcf) | is.null(opt$out)) {
  cat("Inpout and output files are reqiured! Please use --help command for more details.\n",
      quote = FALSE)
  stop("Exiting..")
}

if (is.null(opt$hap) &
    is.null(opt$snp)  & is.null(opt$pair) & opt$random == 0) {
  cat("At least one parameter of --snp, --pair_i, --pair_a, or --hap is reqiured. Please use --help command for more details.\n",
      quote = FALSE)
  stop("Exiting..")
}

if (file.access(opt$vcf, mode = 4) != 0) {
  stop("Variant file can not be accessed !", call. = FALSE)
}

if (!is.null(opt$snp))
{
  if (file.access(opt$vcf, mode = 4) != 0) {
    stop("Snp file can not be accessed !", call. = FALSE)
  }
}

if (!is.null(opt$pair))
{
  if (file.access(opt$pair, mode = 4) != 0) {
    stop("SNP pairs file can not be accessed !", call. = FALSE)
  }
}

if (!is.null(opt$hap))
{
  if (file.access(opt$hap, mode = 4) != 0) {
    stop("haplotype file can not be accessed !", call. = FALSE)
  }
}


simulate_no_causal <- function(individial_count) {
  return (rnorm(individial_count))
}

simulate_gene_expression_one_snp <-
  function(genotypes, heritability = 0.05) {
    y_g <- genotypes
    var_y_g <- var(y_g)
    y_e <-
      sqrt((var_y_g * (1 - heritability) / heritability)) * rnorm(length(y_g))
    varYe <- var(y_e)
    Y <- y_g + y_e
    sm_gene_expression <- (Y - mean(Y)) / sd(Y)
    return (sm_gene_expression)
  }

simulate_gene_expression_two_snps <-
  function(genotypes1,
           genotypes2,
           type = "additive",
           heritability = 0.05) {
    y_g1 <- genotypes1
    y_g2 <- genotypes2
    
    r2 <- cor(y_g1, y_g2)
    if (r2 * r2 > 0.8)
    {
      cat(
        "Warning: the correlation between the input SNPs is greater than 0.8, we recommend to use SNPs with less correlation for simulations.\n"
      )
    }
    
    if (type == "additive")
    {
      y_g <- y_g1 + y_g2
    }
    else if (type == "interaction")
    {
      y_g <- y_g1 * y_g2
    }
    
    if (length(unique(y_g)) == 1)
    {
      return (c("skipped"))
    }
    r2 <- cor(y_g1, y_g)
    if (r2 * r2 > 0.8)
    {
      cat(
        paste0(
          "Warning: the correlation between the first SNP and the combined pair (",
          type,
          ") is greater than 0.8, we recommend using other SNPs.\n"
        )
      )
    }
    
    r2 <- cor(y_g, y_g2)
    if (r2 * r2 > 0.8)
    {
      cat(
        paste0(
          "Warning: the correlation between the second SNP and the combined pair (",
          type,
          ") is greater than 0.8, we recommend using other SNPs.\n"
        )
      )
    }
    
    var_y_g <- var(y_g)
    y_e <-
      sqrt((var_y_g * (1 - heritability) / heritability)) * rnorm(length(y_g))
    varYe <- var(y_e)
    Y <- y_g1 + y_e
    sm_gene_expression <- (Y - mean(Y)) / sd(Y)
    return (sm_gene_expression)
  }


get_genotype <- function(alleles) {
  x <- strsplit(alleles, ":")[[1]][[1]]
  if (grepl("[.]", x))
    return (0)
  
  if (grepl("/", x))
    return (as.numeric(strsplit(x, "/")[[1]][[1]]) + as.numeric(strsplit(x, "/")[[1]][[2]]))
  else if (grepl("[|]", x))
    return (as.numeric(strsplit(x, "[|]")[[1]][[1]]) + as.numeric(strsplit(x, "[|]")[[1]][[2]]))
  else
    stop (paste0("Something wrong with allele separator ", alleles))
}

get_haplotype <- function(alleles) {
  x <- strsplit(alleles, ":")[[1]][[1]]
  x <- gsub("[.]", "0", x)
  if (grepl("[|]", x))
    return (c(as.numeric(strsplit(x, "[|]")[[1]][[1]]) , as.numeric(strsplit(x, "[|]")[[1]][[2]])))
  else
    stop (paste0("Something wrong with allele separator ", alleles))
}


get_block_haplotypes <- function(block) {
  hap1 <- c()
  hap2 <- c()
  individual <- 1
  for (individual in 1:length(block))
  {
    hap1 <-
      c(hap1, paste0(sapply (block[, individual], get_haplotype)[1, ], collapse = ""))
    hap2 <-
      c(hap2, paste0(sapply (block[, individual], get_haplotype)[2, ], collapse = ""))
  }
  return (list("hap1" = hap1, "hap2" = hap2))
}


header_row <- 1
variant_data <- file(opt$vcf, "r")
while (TRUE) {
  line <- readLines(variant_data, n = 1)
  #print(line)
  if (startsWith(line, "#CHROM")) {
    break
  }
  header_row <- header_row + 1
}
close(variant_data)

variant_data <-
  read.table(
    opt$vcf,
    skip = header_row - 1,
    comment.char = "",
    header = T,
    stringsAsFactors = F
  )

colnames(variant_data)[[1]] <- "CHROM"

if (!is.null(opt$snp))
{
  snp_data <- read.table(opt$snp)
  snp_out <- data.frame()
  for (snp_id in snp_data[, 1])
  {
    genotypes <-
      sapply(variant_data[snp_id, 10:length(variant_data)], get_genotype)
    if (length(unique(genotypes)) == 1)
    {
      cat(
        paste0(
          "Simulation is skipped - There is no variance in the genotype of the snp (",
          snp_id ,
          ").\n"
        )
      )
    }
    else
    {
      sm <- simulate_gene_expression_one_snp(genotypes, opt$h2)
      snp_out <- rbind(snp_out, c(snp_id, sm))
    }
  }
  colnames(snp_out) <-
    c("snp_id", colnames(variant_data)[10:length(variant_data)])
  write.table(
    snp_out,
    paste0(opt$out, "_snps_simulations.txt"),
    quote = F,
    row.names = F,
    sep = "\t"
  )
  cat(paste0(
    "SNP-based simulations are saved at:",
    opt$out,
    "_snps_simulations.txt\n"
  ))
  rm(list = c("snp_out", "snp_data"))
}

if (!is.null(opt$pair_a))
{
  pair_data <- read.table(opt$pair_a)
  pair_out <- data.frame()
  i <- 1
  for (i in 1:nrow(pair_data))
  {
    genotypes1 <-
      sapply(variant_data[pair_data[i, 1], 10:length(variant_data)], get_genotype)
    genotypes2 <-
      sapply(variant_data[pair_data[i, 2], 10:length(variant_data)], get_genotype)
    if (length(unique(genotypes1)) == 1 |
        length(unique(genotypes2)) == 1)
    {
      cat(
        paste0(
          "Simulation is skipped - There is no variance in the genotype of the snp pair (",
          pair_data[i, 1],
          ", ",
          pair_data[i, 2] ,
          ").\n"
        )
      )
    }
    else
    {
      sm <-
        simulate_gene_expression_two_snps(genotypes1,
                                          genotypes2,
                                          type = "additive",
                                          heritability = opt$h2)
      if (sm[[1]] == "skipped")
      {
        cat(
          paste0(
            "Simulation is skipped - There is no variance in the combined impact (additive) of the snp pair (",
            pair_data[i, 1],
            ", ",
            pair_data[i, 2] ,
            ").\n"
          )
        )
        
      } else
        pair_out <-
          rbind(pair_out, c(pair_data[i, 1], pair_data[i, 2], sm))
    }
    
  }
  colnames(pair_out) <-
    c("snp_id1", "snp_id2", colnames(variant_data)[10:length(variant_data)])
  write.table(
    pair_out,
    paste0(opt$out, "_pair_additive_simulations.txt"),
    quote = F,
    row.names = F,
    sep = "\t"
  )
  cat(
    paste0(
      "Simulations based on additive impact of a SNP pair are saved at: ",
      opt$out,
      "_pair_additive_simulations.txt\n"
    )
  )
  
  rm(list = c("pair_data", "pair_out"))
}

if (!is.null(opt$pair_i))
{
  pair_data <- read.table(opt$pair_i)
  pair_out <- data.frame()
  i <- 1
  for (i in 1:nrow(pair_data))
  {
    genotypes1 <-
      sapply(variant_data[pair_data[i, 1], 10:length(variant_data)], get_genotype)
    genotypes2 <-
      sapply(variant_data[pair_data[i, 2], 10:length(variant_data)], get_genotype)
    if (length(unique(genotypes1)) == 1 |
        length(unique(genotypes2)) == 1)
    {
      cat(
        paste0(
          "Simulation is skipped - There is no variance in the genotype of the snp pair (",
          pair_data[i, 1],
          ", ",
          pair_data[i, 2] ,
          ").\n"
        )
      )
    }
    else
    {
      sm <-
        simulate_gene_expression_two_snps(genotypes1,
                                          genotypes2,
                                          type = "interaction",
                                          heritability = opt$h2)
      if (sm[[1]] == "skipped")
      {
        cat(
          paste0(
            "Simulation is skipped - There is no variance in the combined impact (interaction) of the snp pair (",
            pair_data[i, 1],
            ", ",
            pair_data[i, 2] ,
            ").\n"
          )
        )
        
      } else
        pair_out <-
          rbind(pair_out, c(pair_data[i, 1], pair_data[i, 2], sm))
    }
  }
  colnames(pair_out) <-
    c("snp_id1", "snp_id2", colnames(variant_data)[10:length(variant_data)])
  write.table(
    pair_out,
    paste0(opt$out, "_pair_interaction_simulations.txt"),
    quote = F,
    row.names = F,
    sep = "\t"
  )
  cat(
    paste0(
      "Simulations based on interaction impact of a SNP pair are saved at: ",
      opt$out,
      "_pair_interaction_simulations.txt\n"
    )
  )
  rm(list = c("pair_data", "pair_out"))
}

if (opt$random > 0)
{
  no_causal_out <- data.frame()
  for (i in 1:opt$random)
  {
    sm <- simulate_no_causal(length(variant_data) - 9)
    no_causal_out <- rbind(no_causal_out, c(i, sm))
  }
  colnames(no_causal_out) <-
    c("sim_id", colnames(variant_data)[10:length(variant_data)])
  write.table(
    no_causal_out,
    paste0(opt$out, "_no_causal_simulations.txt"),
    quote = F,
    row.names = F,
    sep = "\t"
  )
  cat(paste0(
    "Simulations without causal are saved at:",
    opt$out,
    "_no_causal.txt\n"
  ))
  rm(list = c("no_causal_out"))
}

if (!is.null(opt$hap))
{
  hap_data <-
    read.table(opt$hap, colClasses = c("numeric", "numeric", "character"))
  hap_out <- data.frame()
  i <- 1
  for (i in 1:nrow(hap_data))
  {
    block <-
      variant_data[hap_data[i, 1]:hap_data[i, 2], 10:length(variant_data)]
    res <- get_block_haplotypes(block)
    b_hap1 <- res[["hap1"]]
    b_hap2 <- res[["hap2"]]
    
    block_encoding <-
      ifelse(b_hap1 == hap_data[i, 3], 1, 0) + ifelse(b_hap2 == hap_data[i, 3], 1, 0)
    
    if (length(unique(block_encoding)) == 1)
    {
      cat(
        paste0(
          "Simulation is skipped - There is no variance in the haplotype block (",
          i ,
          ").\n"
        )
      )
    }
    else
    {
      sm <-
        simulate_gene_expression_one_snp(block_encoding, heritability = opt$h2)
      hap_out <- rbind(hap_out, c(i, sm))
    }
  }
  colnames(hap_out) <-
    c("block_id", colnames(variant_data)[10:length(variant_data)])
  write.table(
    hap_out,
    paste0(opt$out, "_haplotype_simulations.txt"),
    quote = F,
    row.names = F,
    sep = "\t"
  )
  cat(
    paste0(
      "Simulations based on haplotypes are saved at: ",
      opt$out,
      "_haplotype_simulations.txt\n"
    )
  )
  rm(list = c("hap_out", "hap_data"))
}

end_time <- Sys.time()

cat(paste("Analysis finished at:", end_time, "\n"))

duration <- difftime(end_time, start_time, units = "mins")

cat(paste("Analysis duration is:", duration , "miutes\n"))
rm(list = ls())
