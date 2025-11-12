version 1.0

# Forked from commit faa824e0322e2ab455ef1cb88bc82c47d753338c of https://github.com/broadinstitute/palantir-workflows

import "Structs.wdl"

task ScoreVcf {
  input {
    File    vcf
    String  basename
    File    weights
    String? extra_args
    File?   sites
    String? chromosome_encoding
    String  docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    Int     addldisk     = 20
    Int     mem_size     = 8
    Int     preemptible  = 1
  }

  Int    base_memory     = if mem_size < 8 then 8 else mem_size
  Int    plink_memory    = base_memory * 1000
  Int    runtime_memory  = base_memory + 2
  Int    file_size       = 3 * ceil(size(vcf, "GB"))
  Int    final_disk_size = file_size + addldisk
  String devdir          = 'DEV'
  String inputsdir       = devdir + '/INPUTS'
  String outputsdir      = devdir + '/OUTPUTS'
  Array[String] inputs   = if defined(sites)
                           then [inputsdir + '/vcf',
                                 inputsdir + '/weights',
                                 inputsdir + '/sites']
                           else [inputsdir + '/vcf',
                                 inputsdir + '/weights']

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset
    set -o xtrace

    mkdir --parents '~{inputsdir}' '~{outputsdir}'
    cp '~{vcf}'     "~{inputsdir}/vcf"
    cp '~{weights}' "~{inputsdir}/weights"

    if '~{if defined(sites) then "true" else "false"}'; then
        cp '~{sites}' "~{inputsdir}/sites"
    fi

    COLUMNS='maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums'
    /plink2                                      \
        --allow-extra-chr ~{extra_args}          \
        ~{"--extract " + sites}                  \
        --memory ~{plink_memory}                 \
        --new-id-max-allele-len 1000 missing     \
        --out ~{basename}                        \
        ~{"--output-chr " + chromosome_encoding} \
        --set-all-var-ids '@:#:$r:$a'            \
        --score '~{weights}'                     \
          header                                 \
          ignore-dup-ids                         \
          list-variants                          \
          no-mean-imputation                     \
          cols="${COLUMNS}"                      \
        --vcf '~{vcf}'
  >>>

  output {
    File        score        = "~{basename}.sscore"
    File        log          = "~{basename}.log"
    File        sites_scored = "~{basename}.sscore.vars"
    Array[File] INPUTS       = inputs
  }

  runtime {
    docker: "~{docker_image}"
    disks: "local-disk ~{final_disk_size} HDD"
    memory: "~{runtime_memory} GB"
    preemptible: preemptible
  }
}

task CheckWeightsCoverSitesUsedInTraining {
  input {
    File      sites_used_in_training
    WeightSet weight_set
    String    docker_image = "python:3.9.10"
    Int       addldisk     = 5
    Int       mem_size     = 4
    Int       preemptible  = 1
  }

  Int file_size = ceil(size(sites_used_in_training, "GB"))
  Int final_disk_size = addldisk + file_size

  command <<<
    python3 << "EOF"
    import csv
    import sys

    with open("~{sites_used_in_training}") as f_sites_used_in_training:
      sites_used_in_training = {s.strip() for s in f_sites_used_in_training}

    sites_in_weight_set = set()
    with open("~{weight_set.linear_weights}") as f_linear_weights:
      linear_weights_reader = csv.reader(f_linear_weights, delimiter='\t')
      next(linear_weights_reader)
      for line in linear_weights_reader:
        sites_in_weight_set.add(line[0])

    if ~{if (defined(weight_set.interaction_weights)) then "True" else "False"}:
      with open("~{weight_set.interaction_weights}") as f_interaction_weights:
        for line in csv.DictReader(f_interaction_weights, delimiter='\t'):
          sites_in_weight_set.add(line['id_1'])
          sites_in_weight_set.add(line['id_2'])

    sites_missing_from_weight_set = sites_used_in_training - sites_in_weight_set

    if len(sites_missing_from_weight_set) > 0:
      sys.exit(f"Error: {len(sites_missing_from_weight_set)} sites used in model training are missing from weights files.")
    EOF
  >>>

  runtime {
    docker : "~{docker_image}"
    disks: "local-disk ~{final_disk_size} SSD"
    memory : "~{mem_size} GB"
    preemptible: preemptible
  }
}

task TrainAncestryModel {
  input {
    File   population_pcs
    File   population_scores
    String output_basename
    String docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
    Int    disk_size    = 100
    Int    mem_size     = 2
    Int    preemptible  = 1
  }

  command <<<
    Rscript -<< "EOF"
      library(dplyr)
      library(readr)
      library(tibble)
      population_pcs = read_tsv("~{population_pcs}")
      population_scores = read_tsv("~{population_scores}")

      population_data = inner_join(population_pcs, population_scores, by=c("IID" = "#IID"))

      # generate the linear model from the population data using the first 4 PCs
      population_model = glm(SCORE1_SUM ~ PC1 + PC2 + PC3 + PC4, data = population_data, family = "gaussian")

      population_data <- population_data %>% mutate(residual_score2 = resid(population_model)^2)

      # generate the linear model for the variance of the score using the first 4 PCs
      population_var_model <- glm(residual_score2 ~ PC1 + PC2 + PC3 + PC4, data = population_data, family = Gamma(link = "log"))

      # use linear model to fit full likelihood model

      # linear transformation to predict variance
      f_sigma2 <- function(t, theta) {
      PC1 = t %>% pull(PC1)
      PC2 = t %>% pull(PC2)
      PC3 = t %>% pull(PC3)
      PC4 = t %>% pull(PC4)
      PC5 = t %>% pull(PC5)
      sigma2 <- exp(theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4)
      }

      # linear transformation to predict mean
      f_mu <- function(t, theta) {
      PC1 = t %>% pull(PC1)
      PC2 = t %>% pull(PC2)
      PC3 = t %>% pull(PC3)
      PC4 = t %>% pull(PC4)
      PC5 = t %>% pull(PC5)
      mu <- theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4
      }

      # negative log likelihood
      nLL_mu_and_var <- function(theta) {
      theta_mu = theta[1:5]
      theta_var = theta[6:10]
      x = population_data %>% pull(SCORE1_SUM)
      sum(log(sqrt(f_sigma2(population_data, theta_var))) + (1/2)*(x-f_mu(population_data, theta_mu))^2/f_sigma2(population_data, theta_var))
      }

      # gradient of negative log likelihood function
      grr <- function(theta) {
      theta_mu = theta[1:5]
      theta_var = theta[6:10]
      d_mu_1 <- 1
      d_mu_2 <- population_data %>% pull(PC1)
      d_mu_3 <- population_data %>% pull(PC2)
      d_mu_4 <- population_data %>% pull(PC3)
      d_mu_5 <- population_data %>% pull(PC4)
      d_sig_7 <- 1 * f_sigma2(population_data, theta_var)
      d_sig_8 <- population_data %>% pull(PC1) * f_sigma2(population_data, theta_var)
      d_sig_9 <- population_data %>% pull(PC2) * f_sigma2(population_data, theta_var)
      d_sig_10 <- population_data %>% pull(PC3) * f_sigma2(population_data, theta_var)
      d_sig_11 <- population_data %>% pull(PC4) * f_sigma2(population_data, theta_var)

      x <- population_data %>% pull(SCORE1_SUM)
      mu_coeff <- -(x - f_mu(population_data, theta_mu))/f_sigma2(population_data, theta_var)
      sig_coeff <- 1/(2*f_sigma2(population_data, theta_var)) -(1/2)*(x - f_mu(population_data, theta_mu))^2/(f_sigma2(population_data, theta_var)^2)

      grad <- c(sum(mu_coeff*d_mu_1),
      sum(mu_coeff*d_mu_2),
      sum(mu_coeff*d_mu_3),
      sum(mu_coeff*d_mu_4),
      sum(mu_coeff*d_mu_5),
      sum(sig_coeff*d_sig_7),
      sum(sig_coeff*d_sig_8),
      sum(sig_coeff*d_sig_9),
      sum(sig_coeff*d_sig_10),
      sum(sig_coeff*d_sig_11)
      )
      }

      # use linear model fits as initial parameters for full likelihood fit
      initial_pars <- c(population_model$coefficients, population_var_model$coefficients)
      initial_pars <- setNames(initial_pars, c("Beta0_mu", "Beta1_mu", "Beta2_mu", "Beta3_mu", "Beta4_mu",
                               "Beta0_var", "Beta1_var", "Beta2_var", "Beta3_var", "Beta4_var"))
      fit_mu_and_var <- optim(nLL_mu_and_var, par = initial_pars, gr = grr, method = "BFGS")

      write(ifelse(fit_mu_and_var$convergence == 0, "true", "false"), "fit_converged.txt")

      write_tsv(enframe(fit_mu_and_var$par), "~{output_basename}_fitted_model_params.tsv")

      population_adjusted <- population_data %>% select(-residual_score2) %>% mutate(adjusted_score =
                                                                        (SCORE1_SUM - f_mu(population_data, fit_mu_and_var$par[1:5]))/
                                                                          sqrt(f_sigma2(population_data, fit_mu_and_var$par[6:10])))
      population_adjusted <- population_adjusted %>% mutate(percentile=pnorm(adjusted_score,0))

      write_tsv(population_adjusted, "population_adjusted_scores.tsv")
    EOF
  >>>

  output {
    File    fitted_params              = "~{output_basename}_fitted_model_params.tsv"
    File    adjusted_population_scores = "population_adjusted_scores.tsv"
    Boolean fit_converged              = read_boolean("fit_converged.txt")
  }

  runtime {
    docker: "~{docker_image}"
    disks: "local-disk ~{disk_size} HDD"
    memory: "~{mem_size} GB"
    preemptible: preemptible
  }
}

task AdjustScores {
  input {
    File   fitted_model_params
    File   pcs
    File   scores
    String output_basename
    String docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
    Int    disk_size    = 100
    Int    mem_size     = 2
    Int    preemptible  = 1
  }

  command <<<
    Rscript -<< "EOF"
      library(dplyr)
      library(readr)

      # read in model params
      params_tibble <- read_tsv("~{fitted_model_params}")
      params <- params_tibble %>% pull(value)

      # linear transformation to predict variance
      f_sigma2 <- function(t, theta) {
        PC1 = t %>% pull(PC1)
        PC2 = t %>% pull(PC2)
        PC3 = t %>% pull(PC3)
        PC4 = t %>% pull(PC4)
        PC5 = t %>% pull(PC5)
        sigma2 <- exp(theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4)
      }

      # linear transformation to predict mean
      f_mu <- function(t, theta) {
        PC1 = t %>% pull(PC1)
        PC2 = t %>% pull(PC2)
        PC3 = t %>% pull(PC3)
        PC4 = t %>% pull(PC4)
        PC5 = t %>% pull(PC5)
        mu <- theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4
      }

      scores = inner_join(read_tsv("~{pcs}"),
                                read_tsv("~{scores}"), by=c("IID" = "#IID"))

      adjusted_scores <- scores %>% mutate(adjusted_score =
                                                        (SCORE1_SUM - f_mu(scores, params[1:5]))/
                                                        sqrt(f_sigma2(scores, params[6:10]))
                                                      )
      adjusted_scores <- adjusted_scores %>% mutate(percentile=pnorm(adjusted_score,0))

      # return array scores
      write_tsv(adjusted_scores, "~{output_basename}_adjusted_scores.tsv")
    EOF
  >>>

  output {
    File adjusted_scores = "~{output_basename}_adjusted_scores.tsv"
  }

  runtime {
    docker: "~{docker_image}"
    disks: "local-disk ~{disk_size} HDD"
    memory: "~{mem_size} GB"
    preemptible: preemptible
  }
}

task ExtractIDsPlink {
  input {
    File   vcf
    String docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    Int    disk_size = 2 * ceil(size(vcf, "GB")) + 100
    Int    mem_size  = 8
    Int    preemptible  = 1
  }

  Int base_memory    = if mem_size < 8 then 8 else mem_size
  Int plink_memory   = base_memory * 1000
  Int runtime_memory = base_memory + 2

  command <<<
    /plink2                                     \
        --allow-extra-chr                       \
        --memory                ~{plink_memory} \
        --new-id-max-allele-len 1000 missing    \
        --rm-dup exclude-all                    \
        --set-all-var-ids       '@:#:$r:$a'     \
        --write-snplist         allow-dups      \
        --vcf                   '~{vcf}'
  >>>

  output {
    File ids = "plink2.snplist"
  }

  runtime {
    docker: "~{docker_image}"
    disks: "local-disk ~{disk_size} HDD"
    memory: "~{runtime_memory} GB"
    preemptible: preemptible
  }
}
