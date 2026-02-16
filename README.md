# Agent-Based Model of Political Opinion Dynamics and Polarization

This repository contains data collection, candidate clustering, agent-based simulation, and statistical hypothesis testing workflows for studying opinion dynamics and political polarization in Danish local elections using candidate data from Altinget's Kandidattest.

---

## 1. Scripts

R and Python scripts for data collection, dimensional reduction, clustering, simulation modeling, and statistical hypothesis testing.

---

### Data Collection

**`Scraping_script.ipynb`**

- Scrapes Aarhus candidate responses from Altinget Kandidattest (KV25) using Playwright
- Automated browser navigation:
  - Selects Aarhus municipality (ID = 271)
  - Scrolls dynamically to load all candidates
- Extracts party–candidate structure via DOM parsing (JavaScript evaluation)
- Extracts 26 issue responses from "Sammenlign svar" section
- Maps visual answer positions to numeric scale `[-2, -1, +1, +2]`
- Implements:
  - Incremental JSON progress saving
  - 3-retry error handling per candidate
  - Timeout protection

**Dependencies**: pandas, playwright (async)

**Input**:  
`https://www.altinget.dk/kandidattest/KV25/valgkort`

**Output**:  
`/work/pilot_project/data/raw/altinget_answers_playwright.csv`  
Columns: `candidate_name`, `party`, `Q1`–`Q26`

---

### Candidate Clustering & Bloc Assignment

**`cluster_code.Rmd`**

Identifies ideological blocs via PCA and k-means clustering on 7 selected policy issues.

#### Data Cleaning
- Loads raw Kandidattest responses (26 questions)
- Selects 7 representative questions:
  - Q1 (Economy)
  - Q2 (Social)
  - Q3 (Transport)
  - Q4 (Environment)
  - Q7 (Education)
  - Q11 (Sports)
  - Q19 (Housing)
- Filters complete cases
- Rescales `[-2, 2] → [-1, 1]`
- Renames to `issue_1`–`issue_7`
- Generates descriptive statistics per issue

Final N ≈ 155 candidates.

#### Principal Component Analysis (PCA)
- Z-standardizes all 7 issues
- Extracts ideological dimensions
- PC1 ≈ 35% variance (left–right dimension)
- PC2 ≈ 18% variance (secondary axis)
- Outputs scree plot and variable contributions
- Saves candidate scores (PC1, PC2)

#### Between-Party Differences
- One-way ANOVA:
  - `PC1 ~ party`
  - `PC2 ~ party`
- Levene's test for variance homogeneity
- Tukey HSD pairwise comparisons

#### K-Means Clustering
- Clustering on PC1 and PC2
- Silhouette analysis (k = 2–6)
- Optimal k = 2 (silhouette ≈ 0.30)
- `nstart = 50` for stability

#### Jaccard Stability (Bootstrap B = 100)
- Cluster 1 ≈ 0.75
- Cluster 2 ≈ 0.77
- > 0.75 interpreted as stable

#### Bloc Characterization
- Mean PC scores per bloc
- T-tests for bloc separation
- Opinion distance analysis:
  - Within-bloc Manhattan distance
  - Between-bloc sampled distance
  - Polarization Ratio (PR) = Between / Within
  - Baseline PR ≈ 1.2–1.4

**Input**:  
`/work/files/altinget_answers_playwright_271.csv`

**Output**:  
`/work/files/data_kandid_with_blocs.csv`  
Columns: `party`, `issue_1`–`issue_7`, `PC1`, `PC2`, `cluster_cand`, `bloc`

This file initializes agent opinions in the simulation.

---

### Agent-Based Simulation

**`simulation_code.R`**

Implements a bounded confidence opinion dynamics model with disagreement mechanisms.

#### Global Parameters
- N = 300 agents
- 7 issues
- μ = 0.05 (convergence rate)
- Max opinion shift = 0.15
- ε ∈ {0.25, 0.35, 0.40, 0.50}
- Ticks = 6000–10000
- Save interval = 50 ticks
- Seeds = 1–20

#### Network Topologies
- Erdős–Rényi (p = 0.02)
- Small-world (neighbors = 6, rewiring = 0.1)
- Barabási–Albert (m = 3)

#### Experimental Conditions
1. Null (bounded confidence only)
2. Repulsion
3. Tie-cutting (edge removed after 3 disagreements)
4. Repulsion + Tie-cutting

#### Opinion Initialization
- Sample 300 candidates with replacement
- Add Gaussian noise (SD = 0.03)
- Clip to [-1, 1]
- Preserve bloc labels

#### Interaction Rule (Per Tick)
- Select random edge
- Select random issue
- Compute issue-specific distance

If distance ≤ ε:
- Assimilation (shift toward each other by μ)

If distance > ε:
- **Repulsion**: shift away by μ
- **Tie-cutting**: increment disagreement counter, remove edge after 3 disagreements
- **Null**: no update (only interact when distance ≤ ε)

Updates capped at ±0.15 per interaction.

#### Metrics (Every 50 Ticks)
- Global L1 distance (1000 sampled pairs)
- Issue-specific dispersion
- Local neighbor agreement
- Network structure:
  - Edges
  - Components
  - Largest component
  - Density
  - Modularity (Louvain)
- Bloc polarization:
  - Within-bloc distance
  - Between-bloc distance
  - Polarization Ratio (PR)

#### Output Per Run
- `history.csv`
- `final_opinions.csv`
- `final_network.edgelist`
- `meta.json`

Directory structure:
