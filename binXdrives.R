#!/usr/bin/env R

# Description ------------------------------------------------------------------
# The script models the impact of a Bipartite Expression Drive (BED) on a 
# population. The BED may be intended to spread 1) toxic proteins that sterilize
# BED-carrying females or 2) lethal seminal proteins that BED-carrying males 
# transfer to females, killing them right after mating.


# Libraries --------------------------------------------------------------------
suppressPackageStartupMessages({
  library("optparse")
  library("stringr")
  library("zeallot")
  library("rlist")
  library("dplyr")
  library("bettermc")
  library("logger")
  library("ggplot2")
})


# Define logging format --------------------------------------------------------
logger <- layout_glue_generator(format = paste(
  "[{crayon::italic(format(time, \"%Y-%m-%d %H:%M:%S\"))}]", 
  "{crayon::bold(colorize_by_log_level(level, levelr))}",
  "{grayscale_by_log_level(msg, levelr)}"))
log_layout(logger)


# Auxiliary Functions ----------------------------------------------------------
validate_argument_required <- function(options, argument) {
  if (is.null(options[[argument]])) {
    log_error("--{argument} must be provided. Run the script with --help for ",
              "more details.")
    quit()
  }
}

validate_argument_positive <- function(options, argument) {
  if (options[[argument]] < 0) {
    log_error("--{argument} must be a positive value: {options[[argument]]}. ", 
              "Run the script with --help for more details.")
    quit()
  }
}

validate_argument_in_list <- function(options, argument, items) {
  if (!(options[[argument]] %in% items)) {
    log_error("--{argument} doesn't match an expected value: ", 
              "{options[[argument]]} [{paste0(items, collapse = \", \")}]. ", 
              "Run the script with --help for more details.")
    quit()
  }
}

validate_argument_in_range <- function(options, argument, left, right) {
  if (!between(options[[argument]], left, right)) {
    log_error("--{argument} is out of range: {options[[argument]]} ",
              "[{left}, {right}]. Run the script with --help for more details.")
    quit()
  }
}

validate_argument_integer_in_range <- function(options, argument, left, right) {
  if (options[[argument]] %% 1 != 0) {
    log_error("--{argument} must be an integer value: {options[[argument]]}. ",
              "Run the script with --help for more details.")
    quit()
  }
  if (!between(options[[argument]], left, right)) {
    log_error("--{argument} is out of range: {options[[argument]]} ",
              "[{left}, {right}]. Run the script with --help for more details.")
    quit()
  }
}

get_gamete_haplotypes <- function(value, vector, allele1, allele2, homing_efficiency, resistance_formation) {
  #
  # Params
  #   value: character
  #     genotype of the gametes provider
  #   vector: array(double)
  #     initial haplotypes probabilities
  #   allele1: character
  #     maternal allele
  #   allele2: character
  #     paternal allele
  #   homing_efficiency: double [0, 1]
  #     drive conversion probability
  #   resistance_formation: double [0, 1]
  #     probability of resistance allele formation
  #
  # Output
  #   Returns a probability vector which is added to the initial haplotype 
  #   probabilities.
  
  a <- vector[1]
  b <- vector[2]
  c <- vector[3]
  
  # Create the M, m, N, n versions
  l1u <- toupper(allele1)
  l1l <- tolower(allele1)
  l2u <- toupper(allele2)
  l2l <- tolower(allele2)
  
  # Prepare a list with all the combinations to use in the if clause
  options <- list(
    "Mm"=paste0(l1u, l1l),
    "mM"=paste0(l1l, l1u),
    "MM"=paste0(l1u, l1u),
    "mm"=paste0(l1l, l1l),
    "Mn"=paste0(l1u, l2l),
    "nM"=paste0(l2l, l1u),
    "mn"=paste0(l1l, l2l),
    "nm"=paste0(l2l, l1l),
    "nn"=paste0(l2l, l2l)
  )
  
  if (value %in% c(options$Mm, options$mM)) {
    out <- c(
      a + (1 + homing_efficiency) / 2,
      b + resistance_formation / 2,
      c + 1 - (1 + homing_efficiency) / 2 - resistance_formation / 2
    )
  } else if (value %in% c(options$MM)) {
    out <- c(a + 1, b, c)
  } else if (value %in% c(options$mm)) {
    out <- c(a, b, c + 1)
  } else if (value %in% c(options$Mn, options$nM)) {
    out <- c(a + 0.5, b + 0.5, c)
  } else if (value %in% c(options$mn, options$nm)) {
    out <- c(a, b + 0.5, c + 0.5)
  } else if (value %in% c(options$nn)) {
    out <- c(a, b + 1, c)
  }
  
  return(out)
}

initialize_simulation <- function(initial_population_size = initial_population_size) {
  #
  # Params
  #   initial_population_size: integer
  #     initial number of virgin adults (males plus females)
  #       it is a constant (see --carrying capacity)
  #
  # Notes
  #   Frequencies/numbers refers to virgin adults right after emerging from pupae
  #   Each generation starts with the egg stage
  #   First raw of the output table represent the last generation before release
  #   We refers to this generation as generation zero, so it is the SECOND raw the
  #   one that reflects the first generation after release
  #
  # Output
  #   Returns a data frame to record each generation results
  #     each row represents a generation (total number of generation is 'gen + 1')
  #     initial_pop_size: initial number of virgin adults (males and females)
  #     mothers: number of inseminated females that will produce next generation eggs
  #     terminated_females: number of females killed by terminators
  #     freq_drive_A: drive A allelic frequency
  #     freq_drive_B: drive B allelic frequency
  #     freq_A_B_: frequency of AB carrying genotypes (A_B_)
  #     freq_A_resistance: A_resistance allelic frequency
  #     freq_B_resistance: B_resistance allelic frequency

  output <-
    data.frame(
      generation         = 0,
      drive_release      = TRUE,
      virgin_adults      = initial_population_size,
      virgin_females     = round(initial_population_size / 2),
      unmated_females    = 0,
      terminated_females = 0,
      mothers            = round(initial_population_size / 2),
      freq_A_B_          = 0,
      freq_drive_A       = 0,
      freq_drive_B       = 0,
      freq_A_resistance  = 0,
      freq_B_resistance  = 0,
      genetic_load       = 0
    )
  if (bed_design == "no") { # freq_drive_A correction for single drives
  	output$freq_drive_A[1] <- 1 
  }
  
  return(output)

}

initialize_genotypes <- function(initial_population_size,
                                 wt_genotype,
                                 A_genotype,
                                 release_A,
                                 B_genotype,
                                 release_B
                        ) {
  #
  # Params
  #  initial_population_size: integer
  #     initial number of virgin adults (males plus females)
  #       it is a constant (see --carrying capacity)
  #  wt_genotype: string
  #     homozygous wt genotype
  #       wt_genotype is a fixed parameter ("aabb")
  #  A_genotype: string
  #     genotype of transgenic males carrying the A drive
  #       A_genotype is a constant ("AAbb")
  #  release_A: integer
  #     number of A-carrying males to be released in the first generation
  #       it is a constant
  #  B_genotype: string
  #     genotype of transgenic males carrying the B drive
  #       B_genotype is a constant ("aaBB")
  #  release_B: integer
  #     number of B-carrying males to be released in the first generation
  #       it is a constant
  #
  # Output
  #     Returns the genotypes and numbers of virgin adults in the generation zero
  #     (in the generation zero number of females and males are equal
  #     until transgenic males are added)
  
  out <- vector("list", length = 5)
  # number of virgin females
  out[[1]] <- round(initial_population_size / 2)
  # virgin females genotypes
  out[[2]] <- rep(wt_genotype, out[[1]])
  # virgin males genotypes
  out[[3]] <- c(out[[2]], rep(A_genotype, release_A), rep(B_genotype, release_B))
  # after-release number of virgin males
  out[[4]] <- length(out[[3]])
   # after-release number of adults
  out[[5]] <- out[[1]] + out[[4]]

  return(out)
}

make_rerelease_decision <- function(generation,
                                    males,
                                    estimated_freqA_B_,
                                    estimated_freqA,
                                    estimated_freqB,
                                    simulation_table,
                                    re_release_limit,
                                    re_release_lower_limit,
                                    critical_difference,
                                    A_genotype,
                                    release_A,
                                    B_genotype,
                                    release_B
                           ) {
  #
  # Params
  #  estimated_freqA_B_: double (0, 1)
  #     frequency of A_B_ genotypes of a random sample of adults
  #  estimated_freqA: double (0, 1)
  #     A allelic frequency of a random sample of adults
  #  estimated_freqB: double (0, 1)
  #     B allelic frequency of a random sample of adults
  #  simulation_table: dataframe
  #     table that records simulation results per generation
  #       see initialize_simulation function
  #  re_release_limit: double (0, 1)
  #     upper limit of estimated_freqA_B_ for re-release to occur
  #       it is a constant
  #  re_release_lower_limit: double (0, 1)
  #     lower limit of estimated_freqA_B_ for re-release to occur 
  #       it is a constant
  #  critical_difference: double (0, 1)
  #     lower limit of the difference between estimated allelic frequencies
  #     for re-release to occur
  #       it is a constant
  #  A_genotype: string
  #     genotype of transgenic males carrying the A drive
  #       A_genotype is a constant ("AAbb")
  #  release_A: integer
  #     number of A-carrying males to be released in the first generation
  #       it is a constant
  #  B_genotype: string
  #     genotype of transgenic males carrying the B drive
  #       B_genotype is a constant ("aaBB")
  #  release_B: integer
  #     number of B-carrying males to be released in the first generation
  #       it is a constant
  #
  # Output
  #     Returns the decision of re-release transgenic males carrying the
  #     less frequent drive and, if re-release occurs, it also updates
  #     after-re-release genotypes and numbers of adult males and females
  
  simulation_table$drive_release[generation] <- FALSE
  if (estimated_freqA_B_ > re_release_lower_limit && estimated_freqA_B_ < re_release_limit) {
    if (estimated_freqA > (1 + critical_difference) * estimated_freqB) {
      males <- c(males, rep(B_genotype, release_B))
      simulation_table$drive_release[generation] <- TRUE
    } else if ((1 + critical_difference) * estimated_freqA < estimated_freqB) {
      males <- c(males, rep(A_genotype, release_A))
      simulation_table$drive_release[generation] <- TRUE
    }
  }
  out <- vector("list", length = 2)
  out[[1]] <- simulation_table
  out[[2]] <- males
  
  return(out)
}

get_allee_effect <- function(population_size, critical_population_size, sensitivity) {
  #
  # Params
  #   population_size: integer
  #     number of virgin adults (males plus females)
  #   critical_population_size: double [>= 0]
  #     number of virgin adults (males plus females) below which Allee effect occurs
  #
  # Note
  #   This function models a mate-finding Allee effect that makes mating oportunities
  #   less likely
  #   The chances that virgin females achieve mating (and become mothers), as well
  #   as the number of mates per female, is modeled as a sigmoid function of the
  #   population size, yielding a 50% reduction at the critical population size.
  #  
  # Output
  #   Returns a proportion number [0, 1] which can be interpreted as both
  #   1) the probability that a virgin female find any mate and
  #   2) the post-'mate-finding Allee effect' maintained fraction of lambda for
  #      the Poisson expected number of mates per mother.
  
  out <- 1 - 1 / (1 + exp(sensitivity * (population_size - critical_population_size)))
  return(out)

}

get_mate_finding_outcome <- function(virgin_females_number,
                                     allee_effect,
                                     females,
                                     simulation_table,
                                     generation
                                     ) {
  #
  # Params
  #  virgin_females_number: integer
  #     number of virgin females
  #  allee_effect: double (0, 1)
  #     see get_allee_effect function
  #  females
  #     concatenated females' genotypes
  #  simulation_table: dataframe
  #     table that records simulation results per generation
  #       see initialize_simulation function
  #
  # Note
  #  This function is only called when population size is sufficiently
  #  low and an Alee-effect occurs
  #
  # Output
  #   Returns a list of three objects after removing females that wont mate.
  #   [[1]] updates females object (female genotypes)
  #   [[2]] updates simulation object (recorded variables)

  out <- vector("list", length = 2)
  mate_finding <- sample(c(0, 1),
                         virgin_females_number,
                         replace = T,
                         prob = c(1 - allee_effect, allee_effect)
                  )
  out[[1]] <- females[mate_finding == 1]
  simulation_table$unmated_females[generation] <- virgin_females_number - length(out[[1]])
  out[[2]] <- simulation_table
  
  return(out)
}

get_mates_probability_density <- function(lambda, lower_bound, upper_bound) {
  #
  # Params
  #   lambda: double [> 0]
  #     Poisson parameter for the expected number of mates per female
  #   lower_bound: integer
  #     minimum admitted number of mates per female
  #   upper_bound: integer [> lower_bound]
  #     maximum admitted number of mates per female
  #  
  # Output
  #   Returns a probability vector where:
  #   first element is the probability of having 'lower_bound' mates per female,
  #   element n is the probability of having 'lower_bound' plus n minus 1 mates per female, and
  #   last element is the probability of having 'upper_bound' mates per female.

  out <- vapply(lower_bound:upper_bound, function(i) {
        (lambda)^i * exp(-lambda) / factorial(i)
      }, double(1))
  out <- out / sum(out)
  
  return(out)
}

get_male_mating_success <- function(unintended_reproductive_cost_A,
                                    unintended_reproductive_cost_B,
                                    bed_design,
                                    males
                           ) {
  #
  # Params
  #  unintended_reproductive_cost_A: double (0, 1)
  #     unintended dominant mating success cost for drive A
  #       if it is 0.1, A-carrying males will exhibit a 10% reduction
  #  unintended_reproductive_cost_B: double (0, 1)
  #     unintended dominant mating success cost for drive B
  #       if it is 0.1, B-carrying males will exhibit a 10% reduction
  #  bed_design: string
  #     "yes" (for BED designs) or "no" (for SDs)
  #  males: string
  #     a 1-dim vector containing the genotypes of the adult virgin males
  #
  # Note
  #  Unintended reproductive cost of each drive reduces male mating success
  #  of drive-carrying males. This cost is dominant.
  #
  # Output
  #   Returns a 1-dim vector with the mating success probability for all males
  
  males_number <- length(males)
  out <- rep(1, males_number)
  for (i in 1:males_number) { 
        if (grepl("A", males[i]) && bed_design == "yes") {
          out[i] <- out[i] - unintended_reproductive_cost_A
        }
        if (grepl("B", males[i])) {
          out[i] <- out[i] - unintended_reproductive_cost_B
        }
      }

  return(out)

}

get_potential_mates <- function(polyandry,
                                virgin_females_number,
                                males,
                                male_mating_success
                       ) {
  #
  # Params
  #  polyandry: integer [> 0]
  #     maximum number of mates per mother (see get_mates_probability_density)
  #  virgin_females_number: integer [> 0]
  #     number of adult virgin females that will mate
  #  males: string
  #     a 1-dim vector containing the genotypes of the adult virgin males
  #  male_mating_success: double (0, 1)
  #     a 1-dim vector with the mating success probability for all males
  #
  # Output
  #   Returns a dataframe containing the genotypes of the males sampled
  #   as potential mates. Each row represents a mother and each column her
  #   potential mates

  out <- data.frame(row.names = 1:virgin_females_number)
  for (i in 1:polyandry) {
    mates <- data.frame(male = sample(males,
                                      virgin_females_number,
                                      replace = T,
                                      prob = male_mating_success
                               )
             )
    out <- cbind(out, mates)
  }

  return(out)
}

get_summary_input <- function(results_list, last_generation = 36, type = "undefined_design") {
  simulations_number <- length(results_list)
  generations <- last_generation + 1
  for (i in 1:simulations_number) {
    simulation_i <- results_list[[i]]
    simulation_i$virgin_adults <-
      simulation_i$virgin_adults / simulation_i$virgin_adults[1]
    number_of_raws_i <- length(simulation_i[,1])
    difference <- generations - number_of_raws_i
    if (difference > 0) {
      new_raws <- data.frame(
                    generation         = 0,
                    drive_release      = FALSE,
                    virgin_adults      = 0,
                    virgin_females     = 0,
                    unmated_females    = 0,
                    terminated_females = 0,
                    mothers            = 0,
                    freq_A_B_          = simulation_i$freq_A_B_[number_of_raws_i],
                    freq_drive_A       = simulation_i$freq_drive_A[number_of_raws_i],
                    freq_drive_B       = simulation_i$freq_drive_B[number_of_raws_i],
                    freq_A_resistance  = simulation_i$freq_A_resistance[number_of_raws_i],
                    freq_B_resistance  = simulation_i$freq_B_resistance[number_of_raws_i],
                    genetic_load       = 0
                  )
      new_raws <- do.call("rbind", replicate(difference, new_raws, simplify = FALSE))
      new_raws$generation <- number_of_raws_i:last_generation 
      simulation_i <- rbind(simulation_i, new_raws)				
    } else if (difference < 0) {
      simulation_i <- simulation_i[1:generations,]
    }
    simulation_i$simulation <- i
    results_list[[i]] <- simulation_i
  }
  
  simulations_summary <- results_list[[1]]
  if (simulations_number > 1) {
    for (i in 2:simulations_number) {
      simulations_summary <- rbind(simulations_summary, results_list[[i]])	
    }
  }
  
  return(simulations_summary)
}

print_crash <- function(simulation_table, generation) {
  print(paste0(
          "population crash! ",
          "females: ",
          simulation_table$virgin_females[generation],
          ", mothers: ",
          simulation_table$mothers[generation]
    )
  )
}

save_crash_generation <- function(simulation_table,
                                  generation,
                                  stage
                         ) {
  #
  # Params
  #  simulation_table: dataframe
  #     table that records simulation results per generation
  #       see initialize_simulation function
  #  generation: integer [> 0]
  #     generation number: generation 1 starts with eggs after first release
  #  stage: string
  #     "eggs", "larvae" or "adults"
  #
  # Note
  #   This function is called when population crashes
  #
  # Output
  #   Updates the simulation table adding the last generation values

  out <- simulation_table
  out[generation + 1, names(out) == "generation"] <- generation
  out$drive_release[generation + 1]      <- FALSE
  out$virgin_adults[generation + 1]      <- 0
  out$virgin_females[generation + 1]     <- 0
  out$unmated_females[generation + 1]    <- 0
  out$terminated_females[generation + 1] <- 0
  out$mothers[generation + 1]            <- 0
  out$freq_A_B_[generation + 1]          <- out$freq_A_B_[generation]
  out$freq_drive_A[generation + 1]       <- out$freq_drive_A[generation]
  out$freq_drive_B[generation + 1]       <- out$freq_drive_B[generation]
  out$freq_A_resistance[generation + 1]  <- out$freq_A_resistance[generation]
  out$freq_B_resistance[generation + 1]  <- out$freq_B_resistance[generation]
  if (stage != "adults") {
    out$genetic_load[generation + 1] <- out$genetic_load[generation]
  }
  
  return(out)
}

save_simulation_results <- function(simulation_table,
                                    generation,
                                    total_adults_number,
                                    virgin_females_number,
                                    adults
                           ) {
  #
  # Params
  #  simulation_table: dataframe
  #     table that records simulation results per generation
  #       see initialize_simulation function
  #  generation: integer [> 0]
  #     generation number: generation 1 starts with eggs after first release
  #  total_adults_number: integer [> 0]
  #     number of emerged adults (males plus females)
  #  virgin_females_number: integer [> 0]
  #     number of emerged females
  #  adults: string
  #     a 1-dim vector containing the genotypes of the emerged adults
  #  stage: string
  #     "eggs", "larvae" or "adults"
  #
  # Note
  #   This function is called at last generation
  #
  # Output
  #   Updates the simulation table adding the last generation values
  
  out <- simulation_table
  out[generation + 1, names(out) == "generation"] <- generation
  out$virgin_adults[generation + 1] <- total_adults_number
  out$virgin_females[generation + 1] <- virgin_females_number
  out$freq_drive_A[generation + 1] <-
    as.numeric(format(
      sum(str_count(adults, "A")) / (2 * total_adults_number),
      scientific = F,
      digits = 9,
      nsmall = 9
    )          )
  out$freq_A_resistance[generation + 1] <-
    as.numeric(format(
      sum(str_count(adults, "x")) / (2 * total_adults_number),
      scientific = F,
      digits = 9,
      nsmall = 9
    )          )
  out$freq_drive_B[generation + 1] <-
    as.numeric(format(
      sum(str_count(adults, "B")) / (2 * total_adults_number),
      scientific = F,
      digits = 9,
      nsmall = 9
    )          )
  out$freq_B_resistance[generation + 1] <-
    as.numeric(format(
      sum(str_count(adults, "y")) / (2 * total_adults_number),
      scientific = F,
      digits = 9,
      nsmall = 9
    )          )
  out$freq_A_B_[generation + 1] <-
    as.numeric(format(
      sum(grepl("A", adults) & grepl("B", adults)) / total_adults_number,
      scientific = F,
      digits = 9,
      nsmall = 9
    )          )

    return(out)
}



# Arguments Definition ---------------------------------------------------------
option_list <- list(
  
  # Simulation parameters
  make_option(
    c("--simulations", "-s"),
    type = "integer",
    metavar = "INTEGER",
    default = 100,
    help = "number of simulations/replicates/populations [default: %default]"
  ),
  make_option(
    c("--generations", "-g"),
    type = "integer",
    metavar = "INTEGER",
    default = 30,
    help = "number of generations to be simulated [default: %default]"
  ),
  
  # Population parameters
  make_option(
    c("--carrying_capacity", "-k"),
    type = "integer",
    metavar = "INTEGER",
    default = NULL,
    help = paste("[REQUIRED] carrying capacity of adult population, which equals ",
                 "the initial adult population size")
  ),
  make_option(
    c("--fecundity", "-f"),
    type = "double",
    metavar = "DOUBLE",
    default = 48,
    help = "female fecundity [default: %default]"
  ),
  make_option(
    c("--larval_survival"),
    type = "double",
    metavar = "DOUBLE",
    default = 0.5,
    help = "pre-density larval to adult survival rate [default: %default]"
  ),
  make_option(
    c("--critical_population_size"),
    type = "double",
    metavar = "DOUBLE",
    default = 250,
    help = paste("mate-finding Allee effect occurs if virgin adult population ",
                 "size falls below this critical population size ",
                 "[default: %default]")
  ),
  make_option(
    c("--polyandry", "-p"),
    type = "integer",
    metavar = "INTEGER",
    default = 3,
    help = "maximum number of mates per mother [default: %default]"
  ),  
  
  # Intervention parameters
  make_option(
    c("--release_A", "-a"),
    type = "integer",
    metavar = "INTEGER",
    default = NULL,
    help = "[REQUIRED] number of released A-carrying males (AAbb)"
  ),
  make_option(
    c("--release_B", "-b"),
    type = "integer",
    metavar = "INTEGER",
    default = NULL,
    help = paste("[REQUIRED] number of released B-carrying males, ",
                 "homozygous aaBB for BEDs, heterozygous AABb for single drives")
  ),
  make_option(
    c("--release_critical_difference"),
    type = "double",
    metavar = "DOUBLE",
    default = 0.3,
    help = paste("next generation re-release of modified males of the less ",
                 "frequent drive (D1) is only allowed if ",
                 "(1 + release_critical_difference) * freq(D1) < freq(D2) ",
                 "[default: %default]")
  ),
  make_option(
    c("--re_release_limit"),
    type = "double",
    metavar = "DOUBLE",
    default = 0.85,
    help = paste("prevent next generation re-release of modified males of the ",
                 "less frequent drive when freq(A_B_) is higher than this ",
                 "limit [default: %default]")
  ),
  
  # System parameters
  make_option(
    c("--bed_design"),
    type = "character",
    metavar = "STRING",
    default = "yes",
    help = paste("modeling a Bipartite Expression Drive (BED)? ",
                 "yes for BED, no for single suppression drives",
                 "[default: %default]")
  ),
  make_option(
    c("--conversion_efficiency", "-c"),
    type = "double",
    metavar = "DOUBLE",
    default = 0.9,
    help = paste("drive convertion efficiency, i.e. homing probability ",
                 "[default: %default]")
  ),
  make_option(
    c("--resistance_formation"),
    type = "double",
    metavar = "DOUBLE",
    default = 0.01,
    help = paste("formation probability of resistance alleles, i.e., allelles ",
                 "that prevent drive homing [default: %default]")
  ),
  make_option(
    c("--resistance_functionality"),
    type = "double",
    metavar = "DOUBLE",
    default = 0,
    help = paste("In BED designs (where target genes are essential for development) ",
                 "it is the probability that viable eggs homozygous for resistant ",
                 "alleles of any drive develop into a larva. In SD designs, where ",
                 "the target gene is essential for female fertility, it is the ",
                 "fraction of Fecundity preserved in mothers homozygous for ",
                 "resistant alleles. ",
                 "Use 1 for r1 resistance alleles that preserve target function ",
                 "[default: %default]")
  ),
  make_option(
    c("--dominance"),
    type = "double",
    metavar = "DOUBLE",
    default = 0.6,
    help = paste0("dominance of intended effects: [default: %default]\n ",
                  "  1: complete dominance (AABB effect = A_B_ effect > aabb effect = 0)\n ",
                  "  0 < x < 1: incomplete dominance (AaBb effect = x * AABB effect ",
                  "  and AABb effect = AaBB effect = [x + (1 - x) / 2] * AABB effect)\n ",
                  "  0: no dominance, i.e., intended effects are recessive\n ",
                  "  (no effects for non-AABB genotypes)")
  ),  
  make_option(
    c("--terminator_efficiency", "-e"),
    type = "double",
    metavar = "DOUBLE",
    default = 0.95,
    help = paste("terminator efficiency, i.e. killing probability during ",
                 "mating [default: %default]")
  ),
  make_option(
    c("--intended_fecundity_cost"),
    type = "double",
    metavar = "DOUBLE",
    default = 0,
    help = "intended fertility/fecundity cost [default: %default]"
  ),
  make_option(
    c("--unintended_reproductive_cost_A"),
    type = "double",
    metavar = "DOUBLE",
    default = 0.025,
    help = paste("unintended dominant fertility/fecundity/mating cost for drive A ",
                 "[default: %default]")
  ),
  make_option(
    c("--unintended_reproductive_cost_B"),
    type = "double",
    metavar = "DOUBLE",
    default = 0.025,
    help = paste("unintended dominant fertility/fecundity/mating cost for drive B ",
                 "[default: %default]")
  ),  
  make_option(
    c("--intended_viability_cost"),
    type = "double",
    metavar = "DOUBLE",
    default = 0,
    help = "intended egg-viability cost [default: %default]"
  ),
  make_option(
    c("--unintended_viability_cost_A"),
    type = "double",
    metavar = "DOUBLE",
    default = 0.025,
    help = paste("unintended dominant egg-viability cost of drive A ",
                 "[default: %default]")
  ),
  make_option(
    c("--unintended_viability_cost_B"),
    type = "double",
    metavar = "DOUBLE",
    default = 0.025,
    help = paste("unintended dominant egg-viability cost of drive B ",
                 "[default: %default]")
  ),

  # Output parameters
  make_option(
    c("--output", "-o"),
    type = "character",
    metavar = "STRING",
    default = "results",
    help = "output directory name [default: %default]"
  ),
  
  # Processing parameters
  make_option(
    c("--threads", "-t"),
    type = "integer",
    metavar = "INTEGER",
    default = 1,
    help = "number of threads to use [default: %default]"
  ),
  make_option(
    c("--seed"),
    type = "integer",
    metavar = "INTEGER",
    default = 2,
    help = "seed value for reproducibility [default: %default]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))


# Arguments Validation ---------------------------------------------------------
# Validation of required arguments
validate_argument_required(opt, "carrying_capacity")
validate_argument_required(opt, "release_A")
validate_argument_required(opt, "release_B")

# Validation of [0,1] arguments
validate_argument_in_range(opt, "larval_survival", 0, 1)
validate_argument_in_range(opt, "re_release_limit", 0, 1)
validate_argument_in_range(opt, "conversion_efficiency", 0, 1)
validate_argument_in_range(opt, "resistance_formation", 0, 1)
validate_argument_in_range(opt, "resistance_functionality", 0, 1)
validate_argument_in_range(opt, "dominance", 0, 1)
validate_argument_in_range(opt, "terminator_efficiency", 0, 1)
validate_argument_in_range(opt, "intended_fecundity_cost", 0, 1)
validate_argument_in_range(opt, "unintended_reproductive_cost_A", 0, 1)
validate_argument_in_range(opt, "unintended_reproductive_cost_B", 0, 1)
validate_argument_in_range(opt, "intended_viability_cost", 0, 1)
validate_argument_in_range(opt, "unintended_viability_cost_A", 0, 1)
validate_argument_in_range(opt, "unintended_viability_cost_B", 0, 1)

# Validate of integer ranged arguments
validate_argument_integer_in_range(opt, "threads", 1, parallel::detectCores())

# Validation of positive arguments
validate_argument_positive(opt, "simulations")
validate_argument_positive(opt, "generations")
validate_argument_positive(opt, "carrying_capacity")
validate_argument_positive(opt, "fecundity")
validate_argument_positive(opt, "critical_population_size")
validate_argument_positive(opt, "polyandry")
validate_argument_positive(opt, "release_A")
validate_argument_positive(opt, "release_B")
validate_argument_positive(opt, "release_critical_difference")
validate_argument_positive(opt, "threads")
validate_argument_positive(opt, "seed")

# Validation against lists
validate_argument_in_list(opt, "bed_design", c("yes", "no"))

# Special validations
if (opt$conversion_efficiency + opt$resistance_formation > 1) {
  log_error("--resistance_formation can not be higher than ",
            "--conversion_efficiency's complement: ", 
            "resistance_formation < 1 - conversion_efficiency. ",
            "Run the script with --help for more details.")
  quit()
}
# Output directory
if (file.exists(opt$output)) {
  log_error("--output points to an already existing folder: {opt$output}. ",
            "Remove the folder first and try again.")
  quit()
}
dir.create(opt$output)


# Constants Definition ---------------------------------------------------------
num_simulations              <- opt$simulations                  # s
num_generations              <- opt$generations                  # g
initial_population_size      <- opt$carrying_capacity            # k
fecundity                    <- opt$fecundity                    # f
larval_survival              <- opt$larval_survival
critical_population_size     <- opt$critical_population_size
polyandry                    <- opt$polyandry                    # p
release_A                    <- opt$release_A                    # a
release_B                    <- opt$release_B                    # b
critical_difference          <- opt$release_critical_difference
re_release_limit             <- opt$re_release_limit
bed_design                   <- opt$bed_design
homing_efficiency            <- opt$conversion_efficiency        # c
resistance_formation         <- opt$resistance_formation
resistance_functionality     <- opt$resistance_functionality
dominance                    <- opt$dominance
terminator                   <- opt$terminator_efficiency        # e
intended_fecundity_cost      <- opt$intended_fecundity_cost
unintended_reproductive_cost_A  <- opt$unintended_reproductive_cost_A
unintended_reproductive_cost_B  <- opt$unintended_reproductive_cost_B
intended_viability_cost      <- opt$intended_viability_cost
unintended_viability_cost_A  <- opt$unintended_viability_cost_A
unintended_viability_cost_B  <- opt$unintended_viability_cost_B

output                       <- opt$output                       # o
threads                      <- opt$threads                      # t
seed                         <- opt$seed


# Set the seed for reproducibility
set.seed(seed)


# Fixed Parameters
wt_genotype <- "aabb"         # Wild-type genotype
A_genotype  <- "AAbb"         # Genotype of released drive A-carrying males
B_genotype  <- "aaBB"         # Genotype of released drive B-carrying males
allee_sensitivity <- 0.01     # Density dependence of the Allee effect
re_release_lower_limit <- 0   # Minimum A_B_ frequency for re-release activation
two_alleles <- dominance      # Proportion of intended system effect in AaBb genotype
three_alleles <-
  dominance + (1 - dominance) / 2  # Proportion of intended system effect
                                   # in AABb or AaBB genotype is intermediate
                                   # between dominance and 1.
if (resistance_functionality < 1) {
  resistance_type <- "r2"
} else {
  resistance_type <- "r1"
}
costly_resistance_genotypes <- "xx|yy"  # r2 non-functional genotypes for A or B target
                                        # genes

# Some parameters need to be adjusted for SD designs
if (bed_design == "no") {
  wt_genotype <- "AAbb"
  B_genotype  <- "AABb"
  re_release_lower_limit <- 1
}
  
maximum_fecundity <- 2 * fecundity  # Maximum admitted fecundity (upper bound of
                                    # Poisson number of produced eggs p/mother)
# Larval population size that reduces survival rate to 50%
population_density_limit <- initial_population_size / (larval_survival - 2 / fecundity)
# Low density (or intrinsic) population growth rate
intrinsic_increase <- fecundity * larval_survival / 2
# Probability density of the number of mates per female
#   [1] - probability of 1 mate
#   [2] - probability of 2 mates
#   [n] - probability of 'polyandry' mates
mates_probability_density <- 
  get_mates_probability_density(
    lambda = 2 * polyandry,
    lower_bound = 1,
    upper_bound = polyandry
  )
# Larvae sample size for estimating genetic frequencies each generation
sample_size <- 100


# Main Program -----------------------------------------------------------------

run_simulation <- function(
  round,
  num_generations, 
  population_density_limit, 
  fecundity, 
  larval_survival, 
  critical_population_size, 
  polyandry, 
  release_A, 
  release_B, 
  critical_difference, 
  re_release_limit, 
  homing_efficiency, 
  resistance_formation,
  resistance_functionality,
  dominance, 
  terminator, 
  intended_fecundity_cost, 
  unintended_reproductive_cost_A, 
  unintended_reproductive_cost_B, 
  intended_viability_cost, 
  unintended_viability_cost_A, 
  unintended_viability_cost_B, 
  output, 
  wt_genotype, 
  A_genotype, 
  B_genotype, 
  allee_sensitivity, 
  re_release_lower_limit, 
  two_alleles, 
  three_alleles,
  resistance_type,
  maximum_fecundity, 
  initial_population_size, 
  intrinsic_increase, 
  mates_probability_density,
  bed_design,
  costly_resistance_genotypes,
  sample_size
) {
  
  simulation <- initialize_simulation(initial_population_size)
  print(paste0("simulation ", round, "..."))
  
  for (gen in 1:(num_generations + 1)) { # gen=1 corresponds to 2nd raw on summary table
    # Initializing virgin adult female and male genotypes/numbers
    if (gen == 1) {
      initial_genotypes     <- initialize_genotypes(initial_population_size,
                                                    wt_genotype,
                                                    A_genotype,
                                                    release_A,
                                                    B_genotype,
                                                    release_B)
      virgin_females_number <- initial_genotypes[[1]]
      females               <- initial_genotypes[[2]]
      males                 <- initial_genotypes[[3]]
      virgin_males_number   <- initial_genotypes[[4]]
      total_adults_number   <- initial_genotypes[[5]]
      release_on            <- TRUE
    # updating virgin adult male genotypes/numbers in case of re-release
    } else {
      rerelease  <- make_rerelease_decision(generation = gen,
                                            males = males,
                                            estimated_freqA_B_ = estimated_freqA_B_,
                                            estimated_freqA = estimated_freqA,
                                            estimated_freqB = estimated_freqA,
                                            simulation_table = simulation,
                                            re_release_limit = re_release_limit,
                                            re_release_lower_limit = re_release_lower_limit,
                                            critical_difference = critical_difference,
                                            A_genotype = A_genotype,
                                            release_A = release_A,
                                            B_genotype = B_genotype,
                                            release_B = release_B)
      simulation <- rerelease[[1]]
      males      <- rerelease[[2]]
      rm(rerelease)
    }
    
    # mate-finding Allee effect occurs if population size is low
    if (total_adults_number < 2 * critical_population_size) { # TODO: in a next version replace condition with 'virgin_males_number < critical_population_size'
      allee_effect <- get_allee_effect(total_adults_number, critical_population_size, allee_sensitivity)
      # some females wont find mate and wont pass to the next stage (mothers)
      mate_finding_outcome <- get_mate_finding_outcome(
                                virgin_females_number = virgin_females_number,
                                allee_effect = allee_effect,
                                females = females,
                                simulation_table = simulation,
                                generation = gen)
      females <- mate_finding_outcome[[1]]
      simulation <- mate_finding_outcome[[2]]
      rm(mate_finding_outcome)
      virgin_females_number <- length(females)
      
      # population crashes if no female finds mate
      if (virgin_females_number == 0) {
        simulation$terminated_females[gen] <- 0
        simulation$mothers[gen] <- 0
        rm(females, males)
        print_crash(simulation_table = simulation, generation = gen)
        break
      }
      # mating females will have a reduced number of mates
      adjusted_mates_probability_density <- get_mates_probability_density(
        lambda = 2 * polyandry * allee_effect,
        lower_bound = 1,
        upper_bound = polyandry
      )
    } else {
      adjusted_mates_probability_density <- mates_probability_density
      simulation$unmated_females[gen] <- 0
    }
    
    # defining male mating success
    male_mating_success <- get_male_mating_success(unintended_reproductive_cost_A,
                                                   unintended_reproductive_cost_B,
                                                   bed_design,
                                                   males)
    # polyandry represents the maximum number of mates per mother
    # determining all potential mates (columns) per female (rows)
    matings <- get_potential_mates(polyandry,
                                   virgin_females_number,
                                   males,
                                   male_mating_success)
    rm(males, male_mating_success)
    
    # females may be killed by terminator males during mating
    # surv saves whether a female survives to matings
    surv <- rep(1, virgin_females_number)
    # pandry saves the number of mates per female
    pandry <- sample(1:polyandry, 
                     virgin_females_number, 
                     replace = T, 
                     prob = mates_probability_density)
    # filling surv (TODO: esto se podria poner dentro de un 'if (bed_design == "yes")' porque todas las hembras sobreviven al mating cuando bed_design == "no")
    for (i in 1:virgin_females_number) {
      for (j in 1:pandry[i]) {
        # TODO: implementar funcion get_intended_effect
        # surv[i] <- get_intended_effect (matings[i, j], terminator, two_alleles, three_alleles)
        # if (surv[i] == 0) { break }
        if (dominance == 1) {
          # under a dominant system
          if (grepl("A", matings[i, j]) && grepl("B", matings[i, j])) {
            surv[i] <- sample(c(0, 1), 1, prob = c(terminator, 1 - terminator))
            if (surv[i] == 0) {
              break
            }
          }

        } else if (dominance == 0) {
          # under a recessive system
          if (matings[i, j] == "AABB") {
            surv[i] <- sample(c(0, 1), 1, prob = c(terminator, 1 - terminator))
            if (surv[i] == 0) {
              break
            }
          }

        } else {
          # under a semidominant system
          if (matings[i, j] == "AABB") { 
            surv[i] <- sample(c(0, 1), 1, prob = c(terminator, 1 - terminator))
            if (surv[i] == 0) {
              break
            }
          } else if ( (grepl("Aa|aA|Ax|xA", matings[i, j]) && grepl("BB", matings[i, j])) ||
                      (grepl("AA", matings[i, j]) && grepl("Bb|bB|By|yB", matings[i, j])) ) {
            surv[i] <- sample(c(0, 1), 1, prob = c(terminator * three_alleles, 1 - terminator * three_alleles))
            if (surv[i] == 0) {
              break
            }
          } else if (grepl("Aa|aA|Ax|xA", matings[i, j]) && grepl("Bb|bB|By|yB", matings[i, j])) {
            surv[i] <- sample(c(0, 1), 1, prob = c(terminator * two_alleles, 1 - terminator * two_alleles))
            if (surv[i] == 0) {
              break
            }
          } 
        }
      }
    }
    
    # Beginning of mothers stage -----------------------------------------------
    surviving_females_number <- sum(surv)
    simulation$terminated_females[gen] <- virgin_females_number - surviving_females_number
    simulation$mothers[gen] <- surviving_females_number
    
    # population crashes if no female survives to mating
    if (surviving_females_number == 0) {
    	 # genetic load is considered 1 if all females are killed by terminators
    	 # simulation[gen + 1, names(simulation) == "genetic_load"] <- 1 EN REVISION
      rm(females, matings, pandry)
      print_crash(simulation_table = simulation, generation = gen)
      break
    }
    
    # calculation of virgin-to-mother survival chances 
    # this will the be used for calculating genetic_load
    mother_survival <- surviving_females_number / virgin_females_number
    
    # only mothers will move on to produce eggs
    females <- females[surv == 1]
    matings <- matings[surv == 1, , drop = FALSE]
    pandry  <- pandry[surv == 1]
    
    # homing drives activity during meiosis (conversion and resistance allele formation)
    # and fecundity costs  
    eggs <- character(0)
    remaining_fecundity <- numeric(surviving_females_number)
    for (i in 1:surviving_females_number) {
      # female A gametes
      female_A_gametes <- get_gamete_haplotypes(
        substr(females[i], 1, 2), c(0, 0, 0), 
        "a", "x", homing_efficiency, resistance_formation
      )
      # female B gametes
      female_B_gametes <- get_gamete_haplotypes(
        substr(females[i], 3, 4), c(0, 0, 0), 
        "b", "y", homing_efficiency, resistance_formation
      )
      
      # male gametes
      male_A_gametes <- c(0, 0, 0)
      male_B_gametes <- c(0, 0, 0)
      for (j in 1:pandry[i]) {
        # male A gametes
        male_A_gametes <- get_gamete_haplotypes(
          substr(matings[i, j], 1, 2), 
          male_A_gametes, "a", "x", homing_efficiency, resistance_formation
        )
        # male B gametes
        male_B_gametes <- get_gamete_haplotypes(
          substr(matings[i, j], 3, 4), 
          male_B_gametes, "b", "y", homing_efficiency, resistance_formation
        )
      }
      male_A_gametes <- male_A_gametes / pandry[i]
      male_B_gametes <- male_B_gametes / pandry[i]
      
      # unintended dominant fertility/fecundity costs of drives on females
      post_costs_fecundity <- fecundity
      if (grepl("A", females[i])) {
        post_costs_fecundity <- post_costs_fecundity * (1 - unintended_reproductive_cost_A)
      }
      if (grepl("B", females[i])) {
        post_costs_fecundity <- post_costs_fecundity * (1 - unintended_reproductive_cost_B)
      }
      
      # intended fertility/fecundity costs of drives on females
      # two_alleles: proportion of remaining cost for AaBb genotype of semidominant systems
      # three_alleles: proportion of remaining cost for AABb or AaBB genotype of semidominant systems
      if (dominance == 1) {
        # under a dominant system
        if ( (grepl("A", females[i]) && grepl("B", females[i])) ||
             (bed_design == "no" && grepl("y", females[i])) ) {
          post_costs_fecundity <- post_costs_fecundity * (1 - intended_fecundity_cost)
        }
        
      } else if (dominance == 0) {
        # under a recessive system
        if (females[i] == "AABB") {
          post_costs_fecundity <- post_costs_fecundity * (1 - intended_fecundity_cost)
        } else if (bed_design == "no" && grepl("yy|By|yB", females[i])) {
          post_costs_fecundity <- post_costs_fecundity * (1 - intended_fecundity_cost * (1 - resistance_functionality))
        }
      
      } else {
        # under a semidominant system
        if (females[i] == "AABB") {
          post_costs_fecundity <- post_costs_fecundity * (1 - intended_fecundity_cost)
        } else if (bed_design == "no") {
          if (grepl("yy", females[i])) {
          	post_costs_fecundity <- post_costs_fecundity * (1 - intended_fecundity_cost * (1 - resistance_functionality))
          } else if (grepl("By|yB", females[i])) {
          	post_costs_fecundity <- post_costs_fecundity * (1 + intended_fecundity_cost * (resistance_functionality * (1 - two_alleles) - 1)) # fecundity of 'yB' genotype is between 'BB' and 'Bb' genotypes depending linearly on the 'y' allele resistance_cost (the complement of its functionality)
          } else if (grepl("Bb|bB", females[i])) {
          	post_costs_fecundity <- post_costs_fecundity * (1 - intended_fecundity_cost * two_alleles)
          } else if (grepl("by|yb", females[i])) {
          	post_costs_fecundity <- post_costs_fecundity * (1 - intended_fecundity_cost * two_alleles * (1 - resistance_functionality)) # fecundity of 'yb' genotype is between 'Bb' and 'bb' genotypes depending linearly on the 'y' allele resistance_cost (the complement of its functionality)  
          }
        } else if ( (grepl("Aa|aA|Ax|xA", females[i]) && grepl("BB", females[i])) ||
                    (grepl("AA", females[i]) && grepl("Bb|bB|By|yB", females[i])) ) {
          post_costs_fecundity <- post_costs_fecundity * (1 - intended_fecundity_cost * three_alleles)
        } else if (grepl("Aa|aA|Ax|xA", females[i]) && grepl("Bb|bB|By|yB", females[i])) {
          post_costs_fecundity <- post_costs_fecundity * (1 - intended_fecundity_cost * two_alleles)
        }

      }
      
      # defining post-costs fecundity and eggs genotypes
      # maximum_fecundity is a constant
      produced_eggs <- min(rpois(1, post_costs_fecundity), maximum_fecundity)
      fA <- sample(c("A", "x", "a"), produced_eggs, replace = T, prob = female_A_gametes)
      mA <- sample(c("A", "x", "a"), produced_eggs, replace = T, prob = male_A_gametes)
      fB <- sample(c("B", "y", "b"), produced_eggs, replace = T, prob = female_B_gametes)
      mB <- sample(c("B", "y", "b"), produced_eggs, replace = T, prob = male_B_gametes)
      eggs <- c(eggs, paste0(fA, mA, fB, mB))
      
      # calculation of "remaining fecundity" to later calculate genetic load
      # "remaning fecundity" is the proportion of fecundity that remains after
      # the drives effects on females fecundity
      remaining_fecundity[i] <- post_costs_fecundity / fecundity
    }
    rm(females, matings, surv, pandry, fA, mA, fB, mB)
    
    # Remaining fecundity (for later calculation of genetic load)
    remaining_fecundity <- mean(remaining_fecundity)

    # Beginning of eggs stage --------------------------------------------------
    eggs_number <- length(eggs)
    
    # population crashes if no fertile eggs are produced
    if (eggs_number == 0) {
      if (gen < num_generations + 1) {
        # last generation results are added to the simulation table
        simulation <- save_crash_generation(simulation_table = simulation,
                                            generation = gen,
                                            stage = "eggs")
      }
      print_crash(simulation_table = simulation, generation = gen)
      break
    }

    # viability costs of drives on females
    surv <- rep(1, eggs_number)
    for (i in 1:eggs_number) {
      
      # TODO: METER TODO ESTO EN UNA FUNCION PROBABLEMENTE LA MISMA QUE PARA FECUNDITY
      
      # viability costs imposed by r2 alleles in BEDs (in SDs resistance cost affect female fecundity)
      if (resistance_type == "r2" && bed_design == "yes" && grepl(costly_resistance_genotypes, eggs[i])) {
        surv[i] <- resistance_functionality
      }
      if (surv[i] != 0) {
        # unintended dominant viability costs
        if (grepl("A", eggs[i])) {
          surv[i] <- surv[i] * (1 - unintended_viability_cost_A)
        }
        if (grepl("B", eggs[i])) {
          surv[i] <- surv[i] * (1 - unintended_viability_cost_B)
        }
        # intended viability costs
        if (dominance == 1) {
          # under a dominant system
          if ( (grepl("A", eggs[i]) && grepl("B", eggs[i])) ||
               (bed_design == "no" && grepl("y", eggs[i])) ) {
            surv[i] <- surv[i] * (1 - intended_viability_cost)
          }

        } else if (dominance == 0) {
          # under a recessive system
          if (eggs[i] == "AABB") {
            surv[i] <- surv[i] * (1 - intended_viability_cost)
          } else if (bed_design == "no" && grepl("yy", eggs[i])) {
            surv[i] <- surv[i] * (1 - intended_viability_cost * (1 - resistance_functionality))
          }

        } else {
          # under a semidominant system
          if (eggs[i] == "AABB") {
            surv[i] <- surv[i] * (1 - intended_viability_cost)
          } else if (bed_design == "no") {
            if (grepl("yy", eggs[i])) {
          	  surv[i] <- surv[i] * (1 - intended_viability_cost * (1 - resistance_functionality))
            } else if (grepl("By|yB", eggs[i])) {
          	  surv[i] <- surv[i] * (1 + intended_viability_cost * (resistance_functionality * (1 - two_alleles) - 1)) # survival of 'yB' genotype is between 'BB' and 'Bb' genotypes depending linearly on the 'y' allele resistance_cost (the complement of its functionality)
            } else if (grepl("Bb|bB", eggs[i])) {
              surv[i] <- surv[i] * (1 - intended_viability_cost * two_alleles)
            } else if (grepl("by|yb", eggs[i])) {
          	  surv[i] <- surv[i] * (1 - intended_viability_cost * two_alleles * (1 - resistance_functionality)) # survival of 'yb' genotype is between 'Bb' and 'bb' genotypes depending linearly on the 'y' allele resistance_cost (the complement of its functionality)
          	}
          } else if ( (grepl("Aa|aA|Ax|xA", eggs[i]) && grepl("BB", eggs[i])) ||
                      (grepl("AA", eggs[i]) && grepl("Bb|bB|By|yB", eggs[i])) ) {
            surv[i] <- surv[i] * (1 - intended_viability_cost * three_alleles)
          } else if (grepl("Aa|aA|Ax|xA", eggs[i]) && grepl("Bb|bB|By|yB", eggs[i])) {
            surv[i] <- surv[i] * (1 - intended_viability_cost * two_alleles)
          }

        }
      }
      surv[i] <- sample(c(0, 1), 1, prob = c(1 - surv[i], surv[i]))
    }
    
    # Remaining viability for later calculation of genetic load
    remaining_viability <- mean(surv)
    
    # Beginning of larvae stage ------------------------------------------------
    larvae_number <- sum(surv)
    
    # population crashes if no eggs hatches into larva
    if (larvae_number == 0) {
      if (gen < num_generations + 1) {
        # last generation results are added to the simulation table
        simulation <- save_crash_generation(simulation_table = simulation,
                                            generation = gen,
                                            stage = "larvae")
      }
      rm(eggs, surv)
      print_crash(simulation_table = simulation, generation = gen)
      break
    }
    
    # genetic load calculation
    simulation[gen + 1, names(simulation) == "genetic_load"] <- 1 - mother_survival * remaining_fecundity * remaining_viability
    
    larvae <- eggs[surv == 1]
    rm(eggs, surv)
    
    # estimating freq_drive_A, freq_drive_B, and freq_A_B_
    sampled_larvae <- sample(larvae, sample_size, replace = TRUE)
    estimated_freqA    <- sum(str_count(sampled_larvae, "A")) / (2 * sample_size)
    estimated_freqB    <- sum(str_count(sampled_larvae, "B")) / (2 * sample_size)
    estimated_freqA_B_ <- sum(grepl("A", sampled_larvae) & grepl("B", sampled_larvae)) / sample_size
    
    # density-dependent mortality
    # population_density_limit: population size for which survival probability is reduced by 50%
    # larval survival: pre-density (density independent) survival probability
    surv <- rep(1, larvae_number)
    survival_prob <- larval_survival * population_density_limit / (population_density_limit + larvae_number)
    surv <- sample(c(0, 1), larvae_number, replace = T, prob = c(1 - survival_prob, survival_prob))
    
    # Beginning of adults stage ------------------------------------------------
    total_adults_number <- sum(surv)
    
    # population crashes if no larva survives to adult
    if (total_adults_number == 0) {
      if (gen < num_generations + 1) { # TODO: pasar a funcion 'save_crash_generation'
        # last generation results are added to the simulation table
        simulation <- save_crash_generation(simulation_table = simulation,
                                            generation = gen,
                                            stage = "adults")
      }
      
      rm(larvae, surv)
      print_crash(simulation_table = simulation, generation = gen)
      break
    }
    
    adults <- larvae[surv == 1]
    rm(larvae)
    
    # adults sex determination
    sex <- sample(c("female", "male"), total_adults_number, replace = T)
    females <- adults[sex == "female"]
    virgin_females_number <- length(females)
    males <- adults[sex == "male"]
    virgin_males_number <- length(males)
    
    # population crashes if only males or females reach adulthood
    if (virgin_females_number == 0 || virgin_males_number == 0) {
      if (gen < num_generations + 1) {
        # last generation results are added to the simulation table
        simulation <- save_crash_generation(simulation_table = simulation,
                                            generation = gen,
                                            stage = "adults")
      }
      
      print_crash(simulation_table = simulation, generation = gen)
      break
    }
    
    # saving last generation results into table
    if (gen < num_generations + 1) {
      simulation <- save_simulation_results(simulation_table = simulation,
                                            generation = gen,
                                            total_adults_number = total_adults_number,
                                            virgin_females_number = virgin_females_number,
                                            adults = adults)
    }
    
    # printing simulation progress (TODO: pasar a funcion 'get_simulation_progress')
    previous_generation_progress <- floor((gen - 1) * 10 / (num_generations + 1))
    current_generation_progress <- floor(gen * 10 / (num_generations + 1))
    current_percent <- 10 * current_generation_progress
    if (current_generation_progress > previous_generation_progress) {
      print(
        paste0(
          "progress: ",
          current_percent,
          "%, ",
          "females: ",
          simulation$virgin_females[gen],
          ", mothers: ",
          simulation$mothers[gen]
        )
      )
    }
  }
  
  # removing unuseful objects
  if (exists("females")) {
    rm(females)
  }
  if (exists("males")) {
    rm(males)
  }
  if (exists("adults")) {
    rm(adults)
  }
  
  # completing postcrash-generations of the result table (simulation) (TODO: pasar a funcion 'complete_postcrash_generations')
  last_generation <- length(simulation[,1]) - 1
  lacking_generations <- num_generations - last_generation
  if (lacking_generations > 0) {
    new_raws <- data.frame(
                  generation         = 0,
                  drive_release      = FALSE,
                  virgin_adults      = 0,
                  virgin_females     = 0,
                  unmated_females    = 0,
                  terminated_females = 0,
                  mothers            = 0,
                  freq_A_B_          = simulation$freq_A_B_[last_generation + 1],
                  freq_drive_A       = simulation$freq_drive_A[last_generation + 1],
                  freq_drive_B       = simulation$freq_drive_B[last_generation + 1],
                  freq_A_resistance  = simulation$freq_A_resistance[last_generation + 1],
                  freq_B_resistance  = simulation$freq_B_resistance[last_generation + 1],
                  genetic_load       = simulation$genetic_load[last_generation + 1]
    )
    new_raws <- do.call("rbind", replicate(lacking_generations, new_raws, simplify = FALSE))
    new_raws$generation <- (last_generation + 1):num_generations 
    simulation <- rbind(simulation, new_raws)  	
  }
  
  # writing per-generation results (TODO: pasar a funcion 'write_simulation')
  write.table(
    simulation,
    file = paste0(output, "/simulation_", round, ".txt"),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )
  
  return(simulation)
}

# TODO: meter este print en las nuevas funciones y agregar la impresion de todos los parametros.

# List to store each simulation dataframe
# simulations <- vector("list", length=num_simulations)

# Printing population features.
print(paste0(
  "Intrinsic rate of increase (Rm) = ", intrinsic_increase)
)

# Simulations
simulations <- bettermc::mclapply(
  1:num_simulations, 
  run_simulation, 
  num_generations, 
  population_density_limit, 
  fecundity, 
  larval_survival, 
  critical_population_size, 
  polyandry, 
  release_A, 
  release_B, 
  critical_difference, 
  re_release_limit, 
  homing_efficiency, 
  resistance_formation,
  resistance_functionality,
  dominance, 
  terminator, 
  intended_fecundity_cost, 
  unintended_reproductive_cost_A, 
  unintended_reproductive_cost_B, 
  intended_viability_cost, 
  unintended_viability_cost_A, 
  unintended_viability_cost_B, 
  output, 
  wt_genotype, 
  A_genotype, 
  B_genotype, 
  allee_sensitivity, 
  re_release_lower_limit, 
  two_alleles, 
  three_alleles,
  resistance_type,
  maximum_fecundity, 
  initial_population_size, 
  intrinsic_increase, 
  mates_probability_density,
  bed_design,
  costly_resistance_genotypes,
  sample_size,
  mc.cores = threads, 
  mc.preschedule = FALSE, 
  mc.progress = FALSE
)

# writing results
list.save(simulations, paste0(output, "/all_simulations.rds"))

# summary output
summary_input  <- get_summary_input(simulations, last_generation = num_generations)

# suppression results
# creating an empty table for suppression results (TODO: pasar a funcion 'initialize_suppression_results')
suppression_results <- data.frame(simulation       = 1:num_simulations,
                                  elimination      = logical(num_simulations),
                                  crash_generation = rep(NA, num_simulations)
                       )

# filling summary table (TODO: pasar a funcion 'get_suppression_results')
for (i in 1:num_simulations) {
	already_crashed <- summary_input$mothers[summary_input$simulation == i] == 0
	if (sum(already_crashed) > 0) {
		suppression_results$elimination[i] <- TRUE
		suppression_results$crash_generation[i] <- min(which(already_crashed)) - 1
	}
}

# writing summary (TODO: pasar a funcion 'write_summary')
write.table(
    suppression_results,
    file = paste0(output, "/suppression.txt"),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
)

# printing final results (TODO: pasar este print a funcion e incluirlo en otro lado)
print(paste0("Elimination rate: ",
             mean(suppression_results$elimination) * 100,
             "%.  ",
             "Mean time for colapse: ",
             mean(suppression_results$crash_generation, na.rm = TRUE),
             " generations."
      )
)

# summary
# creating an empty table for summary (TODO: pasar a funcion 'initialize_summary')
total_generations <- num_generations + 1 # TODO: pasar a Fixed Parameters
variables <- c("virgin_adults", "freq_drive_A", "freq_drive_B", "freq_A_B_", "freq_A_resistance", "freq_B_resistance", "genetic_load", "drive_release")
variables_number <- length(variables)
summary <- data.frame(generation = rep(0:num_generations, variables_number),
                      variable   = c(rep(variables[1], total_generations),
                                     rep(variables[2], total_generations),
                                     rep(variables[3], total_generations),
                                     rep(variables[4], total_generations),
                                     rep(variables[5], total_generations),
                                     rep(variables[6], total_generations),
                                     rep(variables[7], total_generations),
                                     rep(variables[8], total_generations)
                                    ),
                      min        = numeric(variables_number * total_generations),
                      mean       = numeric(variables_number * total_generations),
                      max        = numeric(variables_number * total_generations)
           )

# filling summary table (TODO: pasar a funcion 'get_summary')
for (i in 0:num_generations) {
  # i indexing
  summary_positions_i <- summary$generation == i
  input_raws_i        <- summary_input$generation == i
	
  for (j in 1:variables_number) {
  	# j indexing
	summary_positions_j <- summary$variable == variables[j]
    input_column_j      <- names(summary_input) == variables[j]
	
	# filling table
	summary$min[summary_positions_i & summary_positions_j] <-
      min(summary_input[input_raws_i, input_column_j])
    summary$mean[summary_positions_i & summary_positions_j] <-
      mean(summary_input[input_raws_i, input_column_j])
    summary$max[summary_positions_i & summary_positions_j] <-
      max(summary_input[input_raws_i, input_column_j])
  }
}

# defining the order of variable levels for plot (esto podria ir dentro de 'get_summary' o de 'get_plot' (ver mas abajo))
summary$variable <- factor(summary$variable, levels = unique(summary$variable))

# writing summary (TODO: pasar a funcion 'write_summary')
write.table(
    summary,
    file = paste0(output, "/summary.txt"),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
)

# output plot
# TODO: pasar a Fixed Parameters
plot_colors <- c("#7460CB", "#D18150", "#EAE646", "#ECD1E0", "#6DE0E7", "#94D998")
# plot_colors_2 <- c("#6841A2", "#E8E048", "#62BADF", "#E0A8D4", "#865F5D", "#828DA5")
# plot_colors_3 <- c("#6841A2", "#E8E048", "#62BADF", "#E0A8D4", "#865F5D", "#F44371")
plot_labels <- c("virgin adults", "drive A", "drive B", "A_B_ genotype", "A-resistance alleles", "B-resistance alleles")

# plot (TODO: pasar a funcion 'get_plot')
plot_data <- summary[summary$variable == "virgin_adults" |
                     summary$variable == "freq_drive_A" |
                     summary$variable == "freq_drive_B" |
                     summary$variable == "freq_A_B_" |
                     summary$variable == "freq_A_resistance" |
                     summary$variable == "freq_B_resistance",]

# tiff(paste0(output, "/plot.tiff"),
#      units = "mm",
#      width = 219,
#      height = 95,
#      res = 300
# )

pdf(paste0(output, "/plot.pdf"), width = 219, height = 95)

ggplot(data=plot_data, aes(x = generation, y = mean, fill = variable)) +
  geom_line(aes(colour = variable))+
  geom_ribbon(aes(ymin = min, ymax = max, fill = variable), alpha = 0.3) +
  labs(x = "Generations", y = "Frequency") +
  scale_fill_manual(labels = plot_labels, values = plot_colors) +
  scale_color_manual(labels = plot_labels, values = plot_colors) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title     = element_blank(),
        legend.position  = "right",
        legend.direction = "vertical",
        panel.grid.major = element_line(colour = "grey90", size = 0.1),
        panel.grid.minor = element_line(colour = "grey90", size = 0.1)
  )

dev.off()
