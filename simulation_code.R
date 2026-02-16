# ==============================================================================
# abm_polarization
# ==============================================================================

library(tidyverse)
library(igraph)
library(jsonlite)

# ============================ GLOBAL PARAMETERS ============================

PARAMS <- list(
  N = 300,
  ISSUES = 7,
  MU = 0.05,
  MAX_CHANGE = 0.15,
  EPSILON = c(0.25, 0.35),
  MAX_TICKS = 6000,          
  SAVE_EVERY = 50,
  ISOLATED_NOISE = 0.01,
  TIECUT_THRESHOLD = 3,
  SEEDS = 1:20
)

BASE_DIR <- "/work/files/simulation/runs"
KANDID_CSV <- "/work/files/data_kandid_with_blocs.csv"  

CONDITIONS <- c("null", "repulsion", "tiecut")

NETWORKS <- list(
  er = function(n, seed) {
    set.seed(seed)
    sample_gnp(n, p = 0.02, directed = FALSE)
  },
  smallworld = function(n, seed) {
    set.seed(seed)
    sample_smallworld(1, n, 6, 0.1)
  },
  ba = function(n, seed) {
    set.seed(seed)
    sample_pa(n, m = 3, directed = FALSE)
  }
)

ISSUE_WEIGHTS <- rep(1 / PARAMS$ISSUES, PARAMS$ISSUES)

# ============================ SKIP LIST ============================

SKIP_CONDITIONS <- c(
  "null_er_eps0.25",
  "null_er_eps0.35",
  "null_smallworld_eps0.25",
  "null_smallworld_eps0.35"
)

# ============================ HELPERS ============================

safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

l1_distance <- function(x, y) mean(abs(x - y))

load_population <- function(seed) {
  set.seed(seed)
  df <- read_csv(KANDID_CSV, show_col_types = FALSE)
  
  opinions <- df %>% select(starts_with("issue_")) %>% as.matrix()
  blocs <- df$bloc
  
  idx <- sample(1:nrow(opinions), PARAMS$N, replace = TRUE)
  pop <- opinions[idx, ]
  agent_blocs <- blocs[idx]
  
  pop <- pop + matrix(rnorm(PARAMS$N * PARAMS$ISSUES, 0, 0.03), PARAMS$N)
  pop <- pmin(pmax(pop, -1), 1)
  
  list(opinions = pop, blocs = agent_blocs)
}

opinion_update <- function(xi, xj, mu, max_change, repulsion = FALSE) {
  delta <- mu * (xj - xi)
  if (repulsion) delta <- -delta
  delta <- pmax(pmin(delta, max_change), -max_change)
  list(
    xi = pmin(pmax(xi + delta, -1), 1),
    xj = pmin(pmax(xj - delta, -1), 1)
  )
}

# ============================ METRICS ============================

compute_all_metrics <- function(pop, G, agent_blocs) {
  
  N <- nrow(pop)
  K <- ncol(pop)
  
  if (is_directed(G)) {
    G <- as.undirected(G)
  }
  
  n_samples <- 1000
  i_sample <- sample(1:N, n_samples, replace = TRUE)
  j_sample <- sample(1:N, n_samples, replace = TRUE)
  valid <- i_sample != j_sample
  
  global_distance <- mean(mapply(function(i, j) {
    l1_distance(pop[i, ], pop[j, ])
  }, i_sample[valid], j_sample[valid]))
  
  issue_distances <- sapply(1:K, function(k) {
    mean(abs(outer(pop[, k], pop[, k], "-")))
  })
  
  issue_sds <- apply(pop, 2, sd)
  
  if (ecount(G) > 0) {
    local_agreements <- sapply(1:N, function(i) {
      neighbors <- neighbors(G, i)
      if (length(neighbors) == 0) return(NA)
      
      similarities <- sapply(neighbors, function(j) {
        1 - l1_distance(pop[i, ], pop[j, ]) / 2
      })
      mean(similarities)
    })
    local_agreement <- mean(local_agreements, na.rm = TRUE)
  } else {
    local_agreement <- NA
  }
  
  comp <- components(G)
  num_components <- comp$no
  largest_component_size <- max(comp$csize)
  largest_component_prop <- largest_component_size / N
  
  density <- edge_density(G)
  
  if (ecount(G) > 0 && num_components == 1) {
    tryCatch({
      communities <- cluster_louvain(G)
      modularity_q <- modularity(communities)
      num_communities <- length(communities)
    }, error = function(e) {
      tryCatch({
        communities <- cluster_fast_greedy(as.undirected(G))
        modularity_q <<- modularity(communities)
        num_communities <<- length(communities)
      }, error = function(e2) {
        modularity_q <<- NA
        num_communities <<- num_components
      })
    })
  } else {
    modularity_q <- NA
    num_communities <- num_components
  }
  
  unique_blocs <- unique(agent_blocs)
  
  if (length(unique_blocs) >= 2) {
    
    within_dists <- sapply(unique_blocs, function(b) {
      members <- which(agent_blocs == b)
      if (length(members) < 2) return(NA)
      
      n_within <- min(200, choose(length(members), 2))
      if (n_within < 1) return(NA)
      
      i_within <- sample(members, n_within, replace = TRUE)
      j_within <- sample(members, n_within, replace = TRUE)
      valid_within <- i_within != j_within
      
      if (sum(valid_within) == 0) return(NA)
      
      mean(mapply(function(i, j) l1_distance(pop[i, ], pop[j, ]),
                  i_within[valid_within], j_within[valid_within]))
    })
    within_bloc_distance <- mean(within_dists, na.rm = TRUE)
    
    bloc_pairs <- combn(unique_blocs, 2)
    between_dists <- apply(bloc_pairs, 2, function(bp) {
      b1_members <- which(agent_blocs == bp[1])
      b2_members <- which(agent_blocs == bp[2])
      
      n_between <- min(200, length(b1_members) * length(b2_members))
      i_between <- sample(b1_members, n_between, replace = TRUE)
      j_between <- sample(b2_members, n_between, replace = TRUE)
      
      mean(mapply(function(i, j) l1_distance(pop[i, ], pop[j, ]),
                  i_between, j_between))
    })
    between_bloc_distance <- mean(between_dists)
    
    polarization_ratio <- between_bloc_distance / within_bloc_distance
    
  } else {
    within_bloc_distance <- NA
    between_bloc_distance <- NA
    polarization_ratio <- NA
  }
  
  tibble(
    global_distance = global_distance,
    issue_dist_1 = issue_distances[1],
    issue_dist_2 = issue_distances[2],
    issue_dist_3 = issue_distances[3],
    issue_dist_4 = issue_distances[4],
    issue_dist_5 = issue_distances[5],
    issue_dist_6 = issue_distances[6],
    issue_dist_7 = issue_distances[7],
    issue_sd_1 = issue_sds[1],
    issue_sd_2 = issue_sds[2],
    issue_sd_3 = issue_sds[3],
    issue_sd_4 = issue_sds[4],
    issue_sd_5 = issue_sds[5],
    issue_sd_6 = issue_sds[6],
    issue_sd_7 = issue_sds[7],
    local_agreement = local_agreement,
    components = num_components,
    largest_component_prop = largest_component_prop,
    density = density,
    modularity = modularity_q,
    num_communities = num_communities,
    edges = ecount(G),
    within_bloc_dist = within_bloc_distance,
    between_bloc_dist = between_bloc_distance,
    polarization_ratio = polarization_ratio
  )
}

# ============================ CORE SIMULATION ============================

run_simulation <- function(condition, network_name, epsilon, seed) {
  
  set.seed(seed)
  
  pop_data <- load_population(seed)
  pop <- pop_data$opinions
  agent_blocs <- pop_data$blocs
  
  G <- NETWORKS[[network_name]](PARAMS$N, seed)
  
  if (is_directed(G)) {
    G <- as.undirected(G, mode = "collapse")
  }
  
  disagreement_counter <- matrix(0, PARAMS$N, PARAMS$N)
  
  history <- tibble()
  
  for (t in 1:PARAMS$MAX_TICKS) {
    
    if (ecount(G) == 0) {
      pop <- pop + matrix(
        rnorm(PARAMS$N * PARAMS$ISSUES, 0, PARAMS$ISOLATED_NOISE),
        PARAMS$N, PARAMS$ISSUES
      )
      pop <- pmin(pmax(pop, -1), 1)
      next
    }
    
    e <- ends(G, sample(ecount(G), 1))
    i <- e[1]; j <- e[2]
    issue <- sample(1:PARAMS$ISSUES, 1, prob = ISSUE_WEIGHTS)
    
    d_issue <- abs(pop[i, issue] - pop[j, issue])
    
    if (d_issue <= epsilon) {
      
      if (condition != "null") {
        upd <- opinion_update(pop[i, ], pop[j, ],
                              PARAMS$MU, PARAMS$MAX_CHANGE,
                              repulsion = FALSE)
        pop[i, ] <- upd$xi
        pop[j, ] <- upd$xj
      }
      
    } else {
      
      disagreement_counter[i, j] <- disagreement_counter[i, j] + 1
      disagreement_counter[j, i] <- disagreement_counter[j, i] + 1
      
      if (condition == "repulsion") {
        upd <- opinion_update(pop[i, ], pop[j, ],
                              PARAMS$MU, PARAMS$MAX_CHANGE,
                              repulsion = TRUE)
        pop[i, ] <- upd$xi
        pop[j, ] <- upd$xj
      }
      
      if (condition == "tiecut" &&
          disagreement_counter[i, j] >= PARAMS$TIECUT_THRESHOLD &&
          are_adjacent(G, i, j)) {
        
        G <- delete_edges(G, get.edge.ids(G, c(i, j)))
      }
    }
    
    if (t %% PARAMS$SAVE_EVERY == 0) {
      metrics <- compute_all_metrics(pop, G, agent_blocs)
      metrics$tick <- t
      history <- bind_rows(history, metrics)
    }
  }
  
  list(history = history, pop = pop, graph = G, blocs = agent_blocs)
}

# ============================ FIND NEXT RUN NUMBER ============================

find_next_run_number <- function(condition_dir) {
  if (!dir.exists(condition_dir)) return(1)
  
  existing_runs <- list.dirs(condition_dir, recursive = FALSE, full.names = FALSE)
  run_numbers <- as.integer(gsub("run_", "", existing_runs))
  run_numbers <- run_numbers[!is.na(run_numbers)]
  
  if (length(run_numbers) == 0) return(1)
  
  return(max(run_numbers) + 1)
}

# ============================ RUN ALL ============================

run_all <- function() {
  
  safe_dir(BASE_DIR)
  
  total_runs <- length(CONDITIONS) * length(NETWORKS) * 
    length(PARAMS$EPSILON) * length(PARAMS$SEEDS)
  
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("STARTING SIMULATION WITH SKIP LIST\n")
  cat(rep("=", 80), "\n")
  cat("Total runs (including skipped):", total_runs, "\n")
  cat("Skipping 4 completed conditions (80 runs)\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  completed <- 0
  skipped <- 0
  run_count <- 0
  start_time <- Sys.time()
  
  for (cond in CONDITIONS) {
    for (net in names(NETWORKS)) {
      for (eps in PARAMS$EPSILON) {
        
        condition_name <- sprintf("%s_%s_eps%.2f", cond, net, eps)
        condition_dir <- file.path(BASE_DIR, condition_name)
        
        # CHECK IF THIS CONDITION SHOULD BE SKIPPED
        if (condition_name %in% SKIP_CONDITIONS) {
          cat(sprintf("\n=== %s | %s | eps=%.2f === [SKIPPING - ALREADY COMPLETE]\n", 
                      cond, net, eps))
          skipped <- skipped + length(PARAMS$SEEDS)
          run_count <- run_count + length(PARAMS$SEEDS)
          next
        }
        
        safe_dir(condition_dir)
        
        cat(sprintf("\n=== %s | %s | eps=%.2f ===\n", cond, net, eps))
        
        for (s in PARAMS$SEEDS) {
          
          run_count <- run_count + 1
          
          # Find next available run number
          next_run_num <- find_next_run_number(condition_dir)
          run_dir <- file.path(condition_dir, sprintf("run_%03d", next_run_num))
          safe_dir(run_dir)
          
          # Check if ALL 4 files exist (in case of partial completion)
          required_files <- c(
            file.path(run_dir, "history.csv"),
            file.path(run_dir, "final_opinions.csv"),
            file.path(run_dir, "final_network.edgelist"),
            file.path(run_dir, "meta.json")
          )
          
          if (all(file.exists(required_files))) {
            cat(sprintf("  [%d/%d]  SKIP run_%03d (seed %d) - already complete\n",
                        run_count, total_runs, next_run_num, s))
            skipped <- skipped + 1
            next
          }
          
          cat(sprintf("  [%d/%d] RUN run_%03d (seed %d)\n",
                      run_count, total_runs, next_run_num, s))
          
          run_start <- Sys.time()
          
          tryCatch({
            res <- run_simulation(cond, net, eps, s)
            
            run_time <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
            
            write_csv(res$history, file.path(run_dir, "history.csv"))
            
            write_csv(
              tibble(
                agent_id = 1:PARAMS$N,
                bloc = res$blocs,
                as_tibble(res$pop, .name_repair = "minimal")
              ),
              file.path(run_dir, "final_opinions.csv")
            )
            
            write_graph(res$graph,
                        file.path(run_dir, "final_network.edgelist"),
                        format = "edgelist")
            
            write_json(
              list(
                condition = cond,
                network = net,
                epsilon = eps,
                seed = s,
                run_number = next_run_num,
                max_ticks = PARAMS$MAX_TICKS,
                runtime_seconds = round(run_time, 2)
              ),
              file.path(run_dir, "meta.json"),
              pretty = TRUE,
              auto_unbox = TRUE
            )
            
            cat(sprintf("    Done in %.1fs\n", run_time))
            completed <- completed + 1
            
          }, error = function(e) {
            cat(sprintf("    ERROR: %s\n", e$message))
            write_lines(
              sprintf("Error at: %s\nMessage: %s\n", Sys.time(), e$message),
              file.path(run_dir, "error.log")
            )
          })
        }
      }
    }
  }
  
  total_time <- difftime(Sys.time(), start_time, units = "hours")
  
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("ALL SIMULATIONS COMPLETE\n")
  cat(rep("=", 80), "\n")
  cat("Total runs:", total_runs, "\n")
  cat("Completed:", completed, "\n")
  cat("Skipped:", skipped, "\n")
  cat("Total time:", round(total_time, 2), "hours\n")
  if (completed > 0) {
    cat("Average per run:", round(as.numeric(total_time) * 60 / completed, 1), "minutes\n")
  }
  cat(rep("=", 80), "\n\n", sep = "")
}

# ============================ MAIN ============================

if (sys.nframe() == 0) {
  run_all()
}
