version 1.0

import "../steps/PRSStructs.wdl"

task ScoreVcf {
    input {
        File vcf
        String basename
        File weights
        String? extra_args
        File? sites
        File? exclude_sites
        String? chromosome_encoding
        Boolean use_ref_alt_for_ids = false
        Boolean use_dosage_annotation = true
        String docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        Int base_mem = 8
        Int mem_size = base_mem + 2
        Int plink_mem = ceil(base_mem * 0.75 * 1000)
        Int disk_size = 3*ceil(size(vcf, "GB")) + 20
        Int preemptible = 1
    }

    String var_ids_string = "@:#:" + if use_ref_alt_for_ids then "\\$r:\\$a" else "\\$1:\\$2"

    command <<<
        /plink2 --score ~{weights} header ignore-dup-ids list-variants no-mean-imputation \
        cols=maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums --set-all-var-ids ~{var_ids_string} --allow-extra-chr ~{extra_args} -vcf ~{vcf} ~{if use_dosage_annotation then "dosage=DS" else ""} \
        --new-id-max-allele-len 1000 missing ~{"--extract " + sites} ~{"--exclude " + exclude_sites} --out ~{basename} --memory ~{plink_mem} ~{"--output-chr " + chromosome_encoding}
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File score = "~{basename}.sscore"
        File log = "~{basename}.log"
        File sites_scored = "~{basename}.sscore.vars"
    }
}

task AddInteractionTermsToScore {
    input {
        File vcf
        File interaction_weights
        File scores
        File? sites
        String basename
        SelfExclusiveSites? self_exclusive_sites # The interaction term will only be added in no more than selfExclusiveSites.maxAllowed of the
        # effect alleles listed in SelfExclusizeSites.sites is observed
        String docker_image = "us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e"
        Int disk_size = 3*ceil(size(vcf, "GB")) + 20
        Float mem_size = 8
        Int block_buffer = 10000000
        Int preemptible = 1
    }

    command <<<

        tabix ~{vcf}

        python3 << "EOF"
        from cyvcf2 import VCF
        import pandas as pd
        import csv

        vcf = VCF("~{vcf}", lazy=True)
        samples = vcf.samples

        def add_allele_to_count(site, allele, dictionary):
            if site in dictionary:
                dictionary[site][allele]=[0]*len(samples)
            else:
                dictionary[site]={allele:[0]*len(samples)}

        interactions_allele_counts = dict()
        interactions_dict = dict()
        positions = set()
        if ~{if defined(sites) then "True" else "False"}:
            with open("~{sites}") as f_sites:
                sites = {s.strip() for s in f_sites}
        else:
            sites = {}
        with open("~{interaction_weights}") as f:
            for line in csv.DictReader(f, delimiter='\t'):
                site_1 = line['id_1']
                site_2 = line['id_2']
                if len(sites) == 0 or site_1 in sites and site_2 in sites:
                    allele_1 = line['allele_1']
                    allele_2 = line['allele_2']
                    chrom_1 = line['chrom_1']
                    chrom_2 = line['chrom_2']
                    pos_1 = int(line['pos_1'])
                    pos_2 = int(line['pos_2'])
                    weight = float(line['weight'])

                    add_allele_to_count(site_1, allele_1, interactions_allele_counts)
                    add_allele_to_count(site_2, allele_2, interactions_allele_counts)
                    interactions_dict[(site_1, allele_1, site_2, allele_2)] = weight
                    positions.add((chrom_1, pos_1))
                    positions.add((chrom_2, pos_2))

        def add_self_exclusive_site(site, allele, dictionary):
            if site in dictionary:
                dictionary[site].add(allele)
            else:
                dictionary[site]={allele}

        self_exclusive_sites = dict()
        max_self_exclusive_sites = ~{if (defined(self_exclusive_sites)) then select_first([self_exclusive_sites]).maxAllowed else 0}
        self_exclusive_sites_counts = [0]*len(samples)
        if ~{if (defined(self_exclusive_sites)) then "True" else "False"}:
            with open("~{if (defined(self_exclusive_sites)) then select_first([self_exclusive_sites]).sites else ''}") as f_self_exclusive_sites:
                for line in csv.DictReader(f_self_exclusive_sites, delimiter='\t'):
                    id = line['id']
                    if len(sites) == 0 or id in sites:
                        chrom = line['chrom']
                        pos = int(line['pos'])
                        allele = line['allele']
                        add_self_exclusive_site(id, allele, self_exclusive_sites)
                        positions.add((chrom, pos))

        #select blocks to read
        positions = sorted(positions)
        current_chrom=positions[0][0]
        current_start=positions[0][1]
        current_end = current_start+1
        buffer=~{block_buffer}

        blocks_to_read=[]
        for site in positions:
            if site[0] != current_chrom or site[1] - current_end > buffer:
                blocks_to_read.append(current_chrom + ":" + str(current_start) + "-" + str(current_end))
                current_chrom=site[0]
                current_start=site[1]
                current_end = current_start+1
            else:
                current_end = site[1] + 1

        #last block
        blocks_to_read.append(current_chrom + ":" + str(current_start) + "-" + str(current_end))

        #count interaction alleles for each sample
        sites_used_in_score = set()
        for block in blocks_to_read:
            for variant in vcf(block):
                alleles = [a for a_l in [[variant.REF], variant.ALT] for a in a_l]
                vid=":".join(s for s_l in [[variant.CHROM], [str(variant.POS)], sorted(alleles)] for s in s_l)
                if vid in interactions_allele_counts:
                    sites_used_in_score.add(vid)
                    for sample_i,gt in enumerate(variant.genotypes):
                        for gt_allele in gt[:-1]:
                            allele = alleles[gt_allele]
                            if allele in interactions_allele_counts[vid]:
                                interactions_allele_counts[vid][allele][sample_i] += 1
                if vid in self_exclusive_sites:
                    sites_used_in_score.add(vid)
                    for sample_i,gt in enumerate(variant.genotypes):
                        for gt_allele in gt[:-1]:
                            allele = alleles[gt_allele]
                            if allele in self_exclusive_sites[vid]:
                                self_exclusive_sites_counts[sample_i] += 1

        #calculate interaction scores for each sample
        interaction_scores = [0] * len(samples)

        def get_interaction_count(site_and_allele_1, site_and_allele_2, sample_i):
            if site_and_allele_1 == site_and_allele_2:
                return interactions_allele_counts[site_and_allele_1[0]][site_and_allele_1[1]][sample_i]//2
            else:
                return min(interactions_allele_counts[site_and_allele_1[0]][site_and_allele_1[1]][sample_i], interactions_allele_counts[site_and_allele_2[0]][site_and_allele_2[1]][sample_i])

        for interaction in interactions_dict:
            for sample_i in range(len(samples)):
                if self_exclusive_sites_counts[sample_i] <= max_self_exclusive_sites:
                    site_and_allele_1 = (interaction[0], interaction[1])
                    site_and_allele_2 = (interaction[2], interaction[3])
                    interaction_scores[sample_i]+=get_interaction_count(site_and_allele_1, site_and_allele_2, sample_i) * interactions_dict[interaction]

        #add interaction scores to linear scores
        df_interaction_score = pd.DataFrame({"sample_id":samples, "interaction_score":interaction_scores}).set_index("sample_id")
        df_scores=pd.read_csv("~{scores}", sep="\t").astype({'#IID':'string'}).set_index("#IID")
        df_scores = df_scores.join(df_interaction_score)
        df_scores['SCORE1_SUM'] = df_scores['SCORE1_SUM'] + df_scores['interaction_score']
        df_scores.to_csv("~{basename}_scores_with_interactions.tsv", sep="\t")
        with open("~{basename}_sites_used_in_interaction_score.ids", "w") as f_sites_used:
            for site in sites_used_in_score:
                f_sites_used.write("%s\n" % site)
        EOF
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File scores_with_interactions = basename + "_scores_with_interactions.tsv"
        File sites_used_in_interaction_score = basename + "_sites_used_in_interaction_score.ids"
    }
}

task CheckWeightsCoverSitesUsedInTraining {
    input {
        File sites_used_in_training
        WeightSet weight_set
        String docker_image = "python:3.9.10"
        Int disk_size = ceil(size(sites_used_in_training, "GB")) + 10
        Int mem_size = 2
        Int preemptible = 1
    }

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
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }
}

task CompareScoredSitesToSitesUsedInTraining {
    input {
        File sites_used_in_training
        File sites_used_in_scoring
        WeightSet weight_set
        String docker_image = "python:3.9.10"
        Int disk_size = ceil(size(sites_used_in_training, "GB") + size(sites_used_in_scoring, "GB")) + 10
        Int mem_size = 2
        Int preemptible = 1
    }

    command <<<
        python3 << "EOF"
        import csv

        with open("~{sites_used_in_training}") as f_sites_used_in_training:
            sites_used_in_training = {s.strip() for s in f_sites_used_in_training}

        with open("~{sites_used_in_scoring}") as f_sites_used_in_scoring:
            sites_used_in_scoring = {s.strip() for s in f_sites_used_in_scoring}

        missing_sites = sites_used_in_training - sites_used_in_scoring

        with open("missing_sites.txt", "w") as f_missing_sites:
            for site in missing_sites:
                f_missing_sites.write(site)

        with open("n_missing_sites.txt", "w") as f_n_missing_sites:
            f_n_missing_sites.write(f"{len(missing_sites)}")

        max_error_up = 0
        max_error_down = 0

        with open("~{weight_set.linear_weights}") as f_weights:
            weights_reader = csv.reader(f_weights, delimiter = "\t")
            next(weights_reader)
            for line in weights_reader:
                id = line[0]
                if id in missing_sites:
                    missing_site_weight = float(line[2])
                    if missing_site_weight > 0:
                        max_error_up += 2 * missing_site_weight
                    else:
                        max_error_down += 2 * missing_site_weight


        if ~{if defined(weight_set.interaction_weights) then "True" else "False"}:
            with open("~{weight_set.interaction_weights}") as f_interaction_weights:
                for line in csv.DictReader(f_interaction_weights, delimiter='\t'):
                    id_1 = line['id_1']
                    id_2 = line['id_2']
                    if id_1 in missing_sites or id_2 in missing_sites:
                        weight_multiplier = 2 if id_1 != id_2 else 1
                        missing_site_weight = float(line['weight'])
                        if missing_site_weight > 0:
                            max_error_up += weight_multiplier * missing_site_weight
                        else:
                            max_error_down += weight_multiplier * missing_site_weight

        with open("max_error_up.txt", "w") as f_max_error_up:
            f_max_error_up.write(f"{max_error_up}")

        with open("max_error_down.txt", "w") as f_max_error_down:
            f_max_error_down.write(f"{max_error_down}")

        EOF
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File missing_sites = "missing_sites.txt"
        Int n_missing_sites = read_int("n_missing_sites.txt")
        Float max_error_up = read_float("max_error_up.txt")
        Float max_error_down = read_float("max_error_down.txt")
    }
}

task CombineScoringSites {
    input {
        File sites_used_linear_score
        File sites_used_interaction_score
        String basename
        String docker_image = "ubuntu:20.04"
        Int disk_size = ceil(size(sites_used_linear_score, "GB") + 2*size(sites_used_interaction_score, "GB")) + 50
        Int mem_size = 2
        Int preemptible = 1
    }

    command <<<
        cat ~{sites_used_linear_score} ~{sites_used_interaction_score} | sort | uniq > ~{basename}_sites_used_in_score.ids
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File combined_scoring_sites = "~{basename}_sites_used_in_score.ids"
    }
}

task AddShiftToRawScores {
    input {
        File raw_scores
        Float shift
        String basename
        String docker_image = "rocker/tidyverse:4.1.0"
        Int disk_size = 100
        Int mem_size = 2
        Int preemptible = 1
    }

    command <<<
        Rscript -<< "EOF"
        library(dplyr)
        library(readr)

        scores <- read_tsv("~{raw_scores}")
        shifted_scores <- scores %>% mutate(SCORE1_SUM = SCORE1_SUM + ~{shift})

        write_tsv(shifted_scores, "~{basename}.tsv")
        EOF
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File shifted_scores = "~{basename}.tsv"
    }
}

task CombineMissingSitesAdjustedScores {
    input {
        File adjusted_scores_shifted_up
        File adjusted_scores_shifted_down
        File adjusted_scores
        Int n_missing_sites
        String condition_name
        String docker_image = "rocker/tidyverse:4.1.0"
        Int disk_size = 100
        Int mem_size = 2
        Int preemptible = 1
    }

    command <<<
        Rscript -<< "EOF"
        library(dplyr)
        library(readr)

        adjusted_scores <- read_tsv("~{adjusted_scores}") %>% transmute(IID, condition = "~{condition_name}", n_missing_sites = ~{n_missing_sites}, adjusted_score, percentile)
        adjusted_scores_shifted_up <- read_tsv("~{adjusted_scores_shifted_up}") %>% transmute(IID, potential_high_adjusted_score = adjusted_score, potential_high_percentile = percentile)
        adjusted_scores_shifted_down <- read_tsv("~{adjusted_scores_shifted_down}") %>% transmute(IID, potential_low_adjusted_score = adjusted_score, potential_low_percentile = percentile)

        adjusted_scores_shifts <- inner_join(inner_join(adjusted_scores, adjusted_scores_shifted_up), adjusted_scores_shifted_down)
        write_tsv(adjusted_scores_shifts, "missing_sites_shifted_scores.tsv")
        EOF
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File missing_sites_shifted_scores = "missing_sites_shifted_scores.tsv"
    }
}

task TrainAncestryModel {
    input {
        File population_pcs
        File population_scores
        String output_basename
        String docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
        Int disk_size = 100
        Int mem_size = 2
        Int preemptible = 1
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

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File fitted_params = "~{output_basename}_fitted_model_params.tsv"
        File adjusted_population_scores = "population_adjusted_scores.tsv"
        Boolean fit_converged = read_boolean("fit_converged.txt")
    }
}

task AdjustScores {
    input {
        File fitted_model_params
        File pcs
        File scores
        String docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
        Int disk_size = 100
        Int mem_size = 2
        Int preemptible = 1
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
        write_tsv(adjusted_scores, "adjusted_scores.tsv")
        EOF
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File adjusted_scores = "adjusted_scores.tsv"
    }   
}

task MakePCAPlot {
    input {
        File population_pcs
        File target_pcs
        String docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
        Int disk_size = 100
        Int mem_size = 2
        Int preemptible = 1
    }

    command <<<
        Rscript -<< "EOF"
        library(dplyr)
        library(readr)
        library(ggplot2)

        population_pcs <- read_tsv("~{population_pcs}")
        target_pcs <- read_tsv("~{target_pcs}")

        ggplot(population_pcs, aes(x=PC1, y=PC2, color="Population")) +
            geom_point(size=0.1, alpha=0.1) +
            geom_point(data = target_pcs, aes(x=PC1, y=PC2, color="Target")) +
            labs(x="PC1", y="PC2") +
            theme_bw()

        ggsave(filename = "PCA_plot.png", dpi=300, width = 6, height = 6)

        EOF
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File pca_plot = "PCA_plot.png"
    }
}


task ExtractIDsPlink {
    input {
        File vcf
        Boolean use_ref_alt_for_ids = false
        String? chromosome_encoding
        String docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        Int disk_size = 2 * ceil(size(vcf, "GB")) + 100
        Int mem_size = 8
        Int preemptible = 1
        Int plink_mem = ceil(mem_size * 0.75 * 1000)
    }

    String var_ids_string = "@:#:" + if use_ref_alt_for_ids then "\\$r:\\$a" else "\\$1:\\$2"

    command <<<
        /plink2 \
            --vcf ~{vcf} \
            --set-all-var-ids ~{var_ids_string} \
            --new-id-max-allele-len 1000 missing \
            --rm-dup exclude-all \
            --allow-extra-chr \
            --write-snplist allow-dups \
            --memory ~{plink_mem} ~{"--output-chr " + chromosome_encoding}
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File ids = "plink2.snplist"
    }
}

#plink chromosome encoding rules: https://www.cog-genomics.org/plink/2.0/data#irreg_output
task DetermineChromosomeEncoding {
    input {
        File weights
        String docker_image = "python:3.9.10"
        Int disk_size = ceil(size(weights, "GB")) + 10
        Int mem_size = 2
        Int preemptible = 1
    }

    command <<<
        python3 << "EOF"
        code = 'MT'
        with open("~{weights}") as weights_file:
            chroms = {s.split("\t")[0].split(":")[0] for i, s in enumerate(weights_file) if i > 0}
            if any('chr' in c for c in chroms):
                    if 'chrM' in chroms:
                            code = 'chrM'
                    else:
                            code = 'chrMT'
            elif 'M' in chroms:
                    code = 'M'

        with open('chr_encode_out.txt', 'w') as write_code_file:
                write_code_file.write(f'{code}\n')
        EOF
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        String chromosome_encoding = read_string("chr_encode_out.txt")
    }
}

task PerformPCA {
    input {
        File bim
        File bed
        File fam
        String basename
        String docker_image = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
        Int nthreads = 16
        Int disk_size = 400
        Int mem_size = 8
        Int preemptible = 1
    }

    # again, based on Wallace commands
    command <<<
        cp ~{bim} ~{basename}.bim
        cp ~{bed} ~{basename}.bed
        cp ~{fam} ~{basename}.fam

    ~/flashpca/flashpca \
        --bfile ~{basename} \
        -n ~{nthreads} \
        -d 20 \
        --outpc ~{basename}.pc \
        --outpve ~{basename}.pc.variance \
        --outload ~{basename}.pc.loadings \
        --outmeansd ~{basename}.pc.meansd
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File pcs = "~{basename}.pc"
        File pc_variance = "~{basename}.pc.variance"
        File pc_loadings = "~{basename}.pc.loadings"
        File mean_sd = "~{basename}.pc.meansd"
        File eigenvectors = "eigenvectors.txt"
        File eigenvalues = "eigenvalues.txt"
    }
}

# This projects the array dataset using the previously generated PCs, using flashPCA
task ProjectArray {
    input {
        File bim
        File bed
        File fam
        File pc_loadings
        File pc_meansd
        String basename
        String? divisor
        String docker_image = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
        Int nthreads = 16
        Int disk_size = 400
        Int mem_size = 8
        Int preemptible = 1
    }
  
    command <<<
        cp ~{bim} ~{basename}.bim
        cp ~{bed} ~{basename}.bed
        cp ~{fam} ~{basename}.fam
    
        cp ~{pc_loadings} loadings.txt
        cp ~{pc_meansd} meansd.txt
    
        # Check if .bim file, pc loadings, and pc meansd files have the same IDs
        # 1. extract IDs, removing first column of .bim file and first rows of the pc files
        awk '{print $2}' ~{basename}.bim > bim_ids.txt
        awk '{print $1}' loadings.txt | tail -n +2 > pcloadings_ids.txt
        awk '{print $1}' meansd.txt | tail -n +2 > meansd_ids.txt
    
        diff bim_ids.txt pcloadings_ids.txt > diff1.txt
        diff bim_ids.txt meansd_ids.txt > diff2.txt
        diff pcloadings_ids.txt meansd_ids.txt > diff3.txt
    
        if [[ -s diff3.txt ]]; then
            echo "PC loadings file and PC means file do not contain the same IDs; check your input files and run again."
            exit 1
        fi
    
        # check if diff files are not empty
        if [[ -s diff1.txt || -s diff2.txt ]]; then
            echo "IDs in .bim file are not the same as the IDs in the PCA files; check that you have the right files and run again."
            exit 1
        fi
    
        ~/flashpca/flashpca \
            --bfile ~{basename} \
            --numthreads ~{nthreads} \
            --project \
            --inmeansd meansd.txt \
            --outproj projections.txt \
            --inload loadings.txt \
            -v \
            ~{"--div " + divisor}
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File projections = "projections.txt"
    }
}

task ArrayVcfToPlinkDataset {
    input {
        File vcf
        File pruning_sites
        File? subset_to_sites
        String basename
        Boolean use_ref_alt_for_ids = false
        String? chromosome_encoding
        String docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        Int disk_size = 3 * ceil(size(vcf, "GB")) + 20
        Int base_mem = 8
        Int mem_size = base_mem + 2
        Int plink_mem = ceil(base_mem * 0.75 * 1000)
        Int preemptible = 1
    }

    String var_ids_string = "@:#:" + if use_ref_alt_for_ids then "\\$r:\\$a" else "\\$1:\\$2"
  
    command <<<
        /plink2 --vcf ~{vcf} --extract-intersect ~{pruning_sites} ~{subset_to_sites} --allow-extra-chr --set-all-var-ids ~{var_ids_string} \
        --new-id-max-allele-len 1000 missing --out ~{basename} --make-bed --rm-dup force-first ~{"--output-chr " + chromosome_encoding} --memory ~{plink_mem}
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
    }

    output {
        File bed = "~{basename}.bed"
        File bim = "~{basename}.bim"
        File fam = "~{basename}.fam"
    } 
}

task GetBaseMemory {
    # NB: This task computes the memory (in gigabytes) required by vcf,
    # according to the recommendations given in
    # https://www.cog-genomics.org/plink/2.0/other#memory

    input {
        File? vcf
        Int? nvariants
        String docker_image = "python:3.11"
    }

    Int     storage   = 20 + 2 * ceil(size(vcf, "GB"))
    Boolean ERROR     = defined(vcf) == defined(nvariants)
    String  OUTPUTDIR = "OUTPUT"
    String  NVARIANTS = OUTPUTDIR + "/nvariants.txt"
    String  GIGABYTES = OUTPUTDIR + "/gigabytes.txt"

    command <<<
        set -o errexit
        set -o pipefail
        set -o nounset
        # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
        # set -o xtrace

        if ~{if ERROR then "true" else "false"}
        then
            printf -- 'INTERNAL ERROR: too few or too many arguments specified' >&2
            exit 1
        fi

        # --------------------------------------------------------------------------

        mkdir --verbose --parents '~{OUTPUTDIR}'

        NVARIANTS=~{if defined(nvariants)
                    then nvariants
                    else "\"$( zgrep --count --invert-match '^#' '" + vcf + "' | tee '" + NVARIANTS + "' )\""}

        python3 <<EOF > '~{GIGABYTES}'
        import math
        print(8 + max(0, math.ceil((${NVARIANTS} - 50000000)/10000000)))
        EOF
    >>>

    runtime {
        disks : "local-disk ~{storage} HDD"
        docker: "~{docker_image}"
    }

    output {
        Int gigabytes  = read_int(GIGABYTES)
        Int nvariants_ = if defined(nvariants) then nvariants else read_int(NVARIANTS)
    }
}

task RenameChromosomesInTsv {
    input {
        File tsv
        Boolean skipheader
        File lookup = "gs://fc-secure-9ea53c3d-d71a-4f59-92c3-63c75c622a88/reference/etc/rename_chromosomes.tsv"
        String docker_image = "python:3.11"
    }

    Int    storage   = 20 + 2 * ceil(size(tsv, "GB"))
    String OUTPUTDIR = "OUTPUT"
    String RENAMED   = OUTPUTDIR + "/renamed_" + basename(tsv)

    command <<<
        python3 <<EOF
        import sys
        import os
        import re

        def error(message):
            print(f'ERROR: {message}', file=sys.stderr)
            sys.exit(1)


        def read_lookup():
            lookup = dict()
            with open('~{lookup}') as reader:
                for rawline in reader:
                    key, value = rawline.rstrip('\r\n').split('\t')
                    lookup[key] = value
            return lookup


        def main():
            def rename(
                chromosomename,
                _lookup=read_lookup(),
                _parse_re=re.compile(r'^([\da-z]+)(.*)', flags=re.I)):

                match = _parse_re.search(chromosomename)

                if match:
                    core, rest = match.groups()
                    if core in _lookup:
                        return f'{_lookup[core]}{rest}'

                error(f'Unsupported H. sapiens chromosome name: {chromosomename}')

            inputtsv = '~{tsv}'
            outputtsv = '~{RENAMED}'

            os.makedirs(os.path.dirname(outputtsv), exist_ok=True)

            with open(outputtsv, 'w') as writer:
                with open(inputtsv) as reader:
                    skipheader = ~{if skipheader then "True" else "False"}

                    if skipheader:
                        writer.write(next(reader))

                    for rawline in reader:
                        row = rawline.rstrip('\r\n').split('\t')
                        parts = row[0].split(':')
                        newid = ':'.join([rename(parts[0])] + parts[1:])
                        print('\t'.join([newid] + row[1:]), file=writer)
                        continuing = True

        # --------------------------------------------------------------------------

        main()
        EOF
    >>>

    runtime {
        disks : "local-disk ~{storage} HDD"
        docker: "~{docker_image}"
    }

    output {
        File renamed = RENAMED
    }
}

task RenameChromosomesInVcf {
    input {
        File vcf
        File rename = "gs://fc-secure-9ea53c3d-d71a-4f59-92c3-63c75c622a88/reference/etc/rename_chromosomes.tsv"
        String docker_image = "biocontainers/bcftools:v1.9-1-deb_cv1"
    }

    Int    storage   = 20 + 2 * ceil(size(vcf, "GB"))
    String OUTPUTDIR = "OUTPUT"
    String RENAMED   = OUTPUTDIR + "/renamed_" + basename(vcf)

    command <<<
        set -o errexit
        # set -o pipefail
        # set -o nounset
        # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
        # set -o xtrace

        # ---------------------------------------------------------------------------

        mkdir --verbose --parents '~{OUTPUTDIR}'

        WORKDIR="$( mktemp --directory )"
        INPUTVCF="${WORKDIR}/input.vcf.gz"

        ln --symbolic --verbose '~{vcf}' "${INPUTVCF}"

        bcftools index    \
            --force       \
            --tbi         \
            "${INPUTVCF}"

        bcftools annotate            \
            --no-version             \
            --output='~{RENAMED}'    \
            --output-type=z          \
            --rename-chr='~{rename}' \
            "${INPUTVCF}"
    >>>

    runtime {
        disks : "local-disk ~{storage} HDD"
        docker: "~{docker_image}"
    }

    output {
        File renamed = RENAMED
    }
}

task SubsetVcf {
    input {
        File inputvcf
        File regions
        String label = "data"
        Boolean nocleanup = false
        String docker_image = "biocontainers/bcftools:v1.9-1-deb_cv1"
    }

    String OUTPUTDIR = "OUTPUT"
    String OUTPUTVCF = OUTPUTDIR + "/" + label + ".vcf.gz"
    String NREGIONS  = OUTPUTDIR + "/NREGIONS"
    Int    storage   = 20 + 3 * ceil(size(inputvcf, "GB"))

    command <<<
        set -o pipefail
        set -o errexit
        set -o nounset
        export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
        set -o xtrace

        # ---------------------------------------------------------------------------

        printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{storage}'
        printf -- 'INITIAL STORAGE UTILIZATION:\n'
        df --human
        printf -- '\n'

        # ---------------------------------------------------------------------------

        mkdir --verbose --parents '~{OUTPUTDIR}'
        WORKDIR="$( mktemp --directory )"
        INPUTVCF="${WORKDIR}/input.vcf.gz"

        ln --symbolic --verbose '~{inputvcf}' "${INPUTVCF}"

        bcftools index      \
            --force         \
            --tbi           \
            "${INPUTVCF}"

        cleanup() {

            bcftools                             \
                norm                             \
                --multiallelics -any             \
                --no-version                     \
                --output-type v                  \
        | bcftools                             \
                annotate                         \
                --no-version                     \
                --output-type v                  \
                --remove 'INFO,FORMAT'           \
                --set-id '%CHROM:%POS:%REF:%ALT'

        }

        if ~{if nocleanup then "true" else "false"}
        then
            POSTPROCESS=cat
        else
            POSTPROCESS=cleanup
        fi

        bcftools                             \
            view                             \
            --no-version                     \
            --output-type v                  \
            --regions-file '~{regions}'      \
            "${INPUTVCF}"                    \
          | "${POSTPROCESS}"                 \
          | bcftools                         \
                view                         \
                --no-version                 \
                --output-type z              \
                --output-file '~{OUTPUTVCF}'

        wc --lines < '~{regions}' > '~{NREGIONS}'

        # ---------------------------------------------------------------------------

        printf -- 'FINAL STORAGE UTILIZATION:\n'
        df --human

        # ---------------------------------------------------------------------------
    >>>

    runtime {
        docker: "~{docker_image}"
        disks : "local-disk ~{storage} HDD"
    }

    output {
        File result   = OUTPUTVCF
        Int  nregions = read_int(NREGIONS)
    }
}