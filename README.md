# Agent-Based Model of Political Opinion Dynamics and Polarization

This repository contains data collection, candidate clustering, agent-based simulation, and statistical hypothesis testing workflows for studying opinion dynamics and political polarization in Danish local elections using candidate data from Altinget's Kandidattest.

---

## 1. Scripts

R and Python scripts for data collection, dimensional reduction, clustering, simulation modeling, and statistical hypothesis testing.

---

### Data Collection

**`Scraping_script.ipynb`**

- Scrapes Aarhus candidate responses from Altinget Kandidattest (KV25) using Playwright
- Automatically selects Aarhus municipality (ID = 271)
- Scrolls dynamically to load all candidates
- Extracts party–candidate groupings via DOM parsing
- Extracts 26 issue responses from the "Sammenlign svar" section
- Maps visual answer positions to numeric scale `[-2, -1, +1, +2]`
- Implements:
  - Incremental JSON progress saving
  - Three-retry error handling per candidate
  - Timeout protection

**Dependencies**: `pandas`, `playwright` (async)

**Input**  
`https://www.altinget.dk/kandidattest/KV25/valgkort`

**Output**  
`/work/pilot_project/data/raw/altinget_answers_playwright.csv`  
Columns: `candidate_name`, `party`, `Q1`–`Q26`

---

### Candidate Clustering & Bloc Assignment

**`cluster_code.Rmd`**

Identifies ideological blocs via Principal Component Analysis (PCA) and k-means clustering on seven selected policy issues.

#### Data Cleaning

- Loads raw Kandidattest responses (26 questions per candidate)
- Selects seven representative questions:
  - Q1 (Economy)
  - Q2 (Social)
  - Q3 (Transport)
  - Q4 (Environment)
  - Q7 (Education)
  - Q11 (Sports)
  - Q19 (Housing)
- Filters complete cases (no missing responses)
- Rescales from `[-2, 2]` to `[-1, 1]`
- Renames variables to `issue_1`–`issue_7`
- Generates descriptive statistics per issue

Final dataset: approximately 155 candidates.

---

#### Principal Component Analysis (PCA)

- Z-standardizes all seven issues
- Extracts ideological dimensions
- PC1 explains approximately 35% of variance (interpreted as left–right dimension)
- PC2 explains approximately 18% of variance (secondary ideological axis)
- Produces scree plot and variable contribution plots
- Saves candidate scores for PC1 and PC2

---

#### Between-Party Differences

- One-way ANOVA:
  - `PC1 ~ party`
  - `PC2 ~ party`
- Levene’s test for homogeneity of variance
- Tukey HSD pairwise comparisons

---

#### K-Means Clustering

- Clusters candidates using PC1 and PC2
- Silhouette analysis for k = 2–6
- Optimal solution: k = 2 (silhouette ≈ 0.30)
- Uses `nstart = 50` for stability

---

#### Cluster Stability (Bootstrap B = 100)

- Jaccard index:
  - Cluster 1 ≈ 0.75
  - Cluster 2 ≈ 0.77
- Values above 0.75 interpreted as stable clusters

---

#### Bloc Characterization

- Mean PC scores per bloc
- Independent t-tests for bloc separation
- Opinion distance analysis:
  - Within-bloc Manhattan distance (normalized across 7 issues)
  - Between-bloc sampled distance
  - Polarization Ratio (PR) = Between / Within
  - Empirical baseline PR ≈ 1.2–1.4

**Input**  
`/work/files/altinget_answers_playwright_271.csv`

**Output**  
`/work/files/data_kandid_with_blocs.csv`  
Columns: `party`, `issue_1`–`issue_7`, `PC1`, `PC2`, `cluster_cand`, `bloc`

This file initializes agent opinions in the simulation.

---

### Agent-Based Simulation

**`simulation_code.R`**

Implements a bounded confidence opinion dynamics model with disagreement mechanisms.

---

#### Global Parameters

- N = 300 agents
- 7 issues
- μ = 0.05 (convergence rate)
- Maximum opinion shift per interaction = 0.15
- ε ∈ {0.25, 0.35, 0.40, 0.50}
- Simulation length = 9000 ticks
- Metrics saved every 50 ticks
- Seeds = 1–20 (20 replications per condition)

---

#### Network Topologies

- Erdős–Rényi (p = 0.02)
- Small-world (neighbors = 6, rewiring probability = 0.1)
- Barabási–Albert (m = 3)

---

#### Experimental Conditions

1. **Null** – bounded confidence only  
2. **Repulsion** – agents move away when disagreement exceeds ε  
3. **Tie-cutting** – edge removed after three consecutive disagreements  
4. **Repulsion + Tie-cutting** – combined mechanism  

---

#### Opinion Initialization

- Samples 300 candidates with replacement
- Adds Gaussian noise (SD = 0.03)
- Clips values to `[-1, 1]`
- Preserves bloc labels

---

#### Interaction Rule (Per Tick)

1. Select random edge  
2. Select random issue  
3. Compute issue-specific distance  

If distance ≤ ε:
- Assimilation (agents move toward each other by μ)

If distance > ε:
- Repulsion: agents move away by μ
- Tie-cutting: disagreement counter incremented; edge removed after 3 disagreements
- Null: no update

Opinion updates are capped at ±0.15 per interaction.

---

#### Metrics (Recorded Every 50 Ticks)

- Global L1 distance (1000 sampled agent pairs)
- Issue-level dispersion (`issue_dist_1–7`)
- Issue-level standard deviations (`issue_sd_1–7`)
- Network metrics:
  - Edges
  - Components
  - Density
  - Modularity (Louvain)
  - Number of communities
- Bloc polarization:
  - Within-bloc distance
  - Between-bloc distance
  - Polarization Ratio (PR)

---

#### Output Per Run

- `history.csv`
- `final_opinions.csv`
- `final_network.edgelist`
- `meta.json`

Directory structure:

```
/work/files/simulation/runs/
  {condition}_{network}_eps{epsilon}/
    run_{seed}/
      history.csv
      final_opinions.csv
      final_network.edgelist
      meta.json
```

Total runs: approximately 480–960.

---

### Statistical Hypothesis Testing

**`hypothesis_testing.Rmd`**

Conducts mixed-effects modeling and robustness checks.

---

#### Data Preparation

- Loads all `history.csv` files
- Parses condition, network, epsilon, and seed from directory structure
- Averages terminal timesteps (ticks 8000–9000)
- Primary dataset: Small-world networks with ε ∈ {0.40, 0.50}
- Approximately 160 runs in main analysis

---

### RQ1: Global Opinion Dispersion

Model:

```r
lmer(global_distance ~ condition + (1|seed))
```

Tests:
- ANOVA
- Estimated marginal means (EMMs)
- Bonferroni contrasts
- Cohen’s d

Outputs:
- `h1_anova.csv`
- `h1_emmeans.csv`
- `h1_contrasts.csv`
- `h1_cohens_d.csv`

---

### RQ2: Moderation by Tolerance (ε)

Model:

```r
lmer(global_distance ~ condition * epsilon + (1|seed))
```

Tests:
- Interaction F-test
- Simple effects analysis
- Epsilon contrasts within conditions

Outputs:
- `h2_anova.csv`
- `h2_emmeans.csv`
- `h2_epsilon_effects.csv`

---

### RQ3: Bloc Polarization

Model:

```r
lmer(polarization_ratio ~ condition + (1|seed))
```

Tests:
- ANOVA
- EMMs
- Bonferroni contrasts
- One-sample tests against PR = 1.0

Outputs:
- `h3_anova.csv`
- `h3_emmeans.csv`
- `h3_contrasts.csv`
- `h3_test_vs_1.csv`

---

### Supplementary Analyses

Network fragmentation:

```r
lmer(edges ~ condition + (1|seed))
lmer(components ~ condition + (1|seed))
lmer(modularity ~ condition + (1|seed))
```

Robustness across network types:

```r
lmer(global_distance ~ condition * network + (1|seed))
```

Outputs include:
- `edges_emmeans.csv`
- `components_emmeans.csv`
- `modularity_emmeans.csv`
- `robustness_anova.csv`
- `robustness_emmeans.csv`

---

### Model Diagnostics

- Shapiro–Wilk test (residual normality)
- Q-Q plots
- Levene’s test (homogeneity)
- Intraclass Correlation Coefficient (ICC)

Outputs:
- `diagnostics_summary.csv`
- `H1_qqplot.png`
- `H2_qqplot.png`
- `H3_qqplot.png`

---

## 2. Key Output Files

### Candidate Data

- `altinget_answers_playwright_271.csv`
- `data_kandid_with_blocs.csv`

### Simulation Results

Directory:  
`/work/files/simulation/runs/`

Key metrics include:
- `global_distance`
- `polarization_ratio`
- `within_bloc_dist`
- `between_bloc_dist`
- `edges`
- `components`
- `modularity`
- `density`
- `local_agreement`

Time series recorded every 50 ticks.

---

### Statistical Analysis Outputs

Directory:  
`/work/files/simulation/analysis_final/`

Includes:
- ANOVA tables
- Estimated marginal means
- Pairwise contrasts
- Effect sizes
- Diagnostics
- Summary tables (`summary_table.csv`, `hypothesis_summary.csv`)

---

## 3. Experimental Design Summary

| Factor | Levels |
|--------|--------|
| Condition | Null, Repulsion, Tie-cutting, Repulsion+Tiecut |
| Network | Erdős–Rényi, Small-world, Barabási–Albert |
| Epsilon (ε) | 0.40, 0.50 |
| Seed | 1–20 |
 
Total simulation runs: 480.

---

## 4. Kandidattest Questions (Selected Issues)

Scaling: Original responses `[-2, -1, +1, +2]` normalized to `[-1, +1]`.

- Q1 – Economy  
- Q2 – Social  
- Q3 – Transport  
- Q4 – Environment  
- Q7 – Education  
- Q11 – Sports  
- Q19 – Housing  

---

## 5. Dependencies

### Python
- pandas  
- playwright  

### R
- tidyverse  
- lme4  
- lmerTest  
- emmeans  
- effectsize  
- ggplot2  
- patchwork  
- car  
- igraph  
- jsonlite  
- factoextra  
- cluster  
- fpc  



