# bpgmm 1.3.1

* Revised the vignette set to reduce repeated examples across vignettes.
* Made the model-and-sampler, data-preparation, model-selection, variable-prioritization, and diagnostics vignettes more formula-focused and closer to the notation in Lu, Li, and Love (2021).
* Changed the variable-prioritization simulation so it no longer duplicates the larger MFA model-selection example.

# bpgmm 1.3.0

* Standardized the native Rcpp interface to snake_case function names.
* Renamed C++ source files to snake_case and removed stale compiled artifacts from `src/`.
* Cleaned internal native wrapper documentation while preserving sampler behavior.

# bpgmm 1.2.9

* Added data-preparation and posterior-diagnostics vignettes.
* Expanded the vignette collection so each vignette covers a distinct analysis task: input preparation, sampler details, small worked examples, larger model selection, exploratory variable prioritization, and multi-chain diagnostics.

# bpgmm 1.2.8

* Added larger model-selection and exploratory variable-prioritization vignettes with runnable examples and plots.
* Documented how to use posterior allocations and loading matrices for exploratory variable prioritization without presenting it as formal Bayesian variable selection.

# bpgmm 1.2.7

* Moved latent-factor score updates from R loops to native C++.
* Moved prior density evaluation for `lambda` and `psi` to native C++ while preserving the existing internal R API.
* Added strict native tests for latent-score sampling, input validation, and closed-form prior density checks.

# bpgmm 1.2.6

* Added `pgmm_rjmcmc_chains()` for CRAN-safe multi-core execution of independent RJMCMC chains.
* Added deterministic per-chain seeding and unit tests for the independent-chain wrapper.

# bpgmm 1.2.5

* Cleaned native C++ readability without changing sampler behavior.
* Replaced stale include guards and prototypes, removed selected namespace pollution, and clarified covariance-constraint branch checks with model-label comments.
* Regenerated Rcpp exports after tightening native helper signatures.

# bpgmm 1.2.4

* Optimized the native allocation sampler by precomputing component covariance factorizations and sampling from normalized log probabilities directly.
* Removed avoidable dense probability/log-density matrices and repeated per-observation covariance validation in `update_PostZ()`.
* Tightened native helper signatures to avoid unnecessary copies and removed stale debug-only C++ comments from the posterior lambda/psi update.

# bpgmm 1.2.3

* Added CRAN-safe base R plots to the getting-started, worked-examples, and model-and-sampler vignettes.
* Visualized simulated clusters, posterior modal allocations, posterior model counts, and mixture-of-factor-analyzers covariance geometry.

# bpgmm 1.2.2

* Optimized the native `Calculate_Cxy()` helper by accumulating cluster sufficient statistics directly from labels instead of building a dense allocation matrix and multiplying through zero weights.
* Removed repeated temporary vector construction inside the Cxy accumulation loop.
* Preserved the existing C++11 baseline for CRAN portability.

# bpgmm 1.2.1

* Hardened Rcpp entry points with explicit validation for dimensions, finite numeric inputs, covariance constraints, and native hyperparameter vectors.
* Kept native build settings on the existing CRAN-friendly C++11 baseline while modernizing validation and headers.
* Cleaned native headers and removed stale C++ declarations and debug-only includes/comments.
* Expanded native unit tests for invalid Cxy, posterior lambda/psi, and hyperparameter-update inputs.

# bpgmm 1.2.0

* Made the public API consistently snake_case now that the package has no CRAN reverse dependencies.
* Removed legacy exported wrappers `pgmmRJMCMC()`, `summarizePgmmRJMCMC()`, and `summerizePgmmRJMCMC()`.
* Renamed public sampler arguments from camelCase to snake_case: `m_init`, `m_range`, `q_new`, `d_vec`, `s_vec`, `m_step`, `v_step`, and `split_combine`.
* Renamed summary arguments and result fields to snake_case, including `true_cluster`, `allocation`, `n_clusters`, `n_constraints`, `allocation_samples`, and `constraint_samples`.
* Updated README, pkgdown vignettes, examples, and tests for the breaking 1.2.0 API.

# bpgmm 1.1.6

* Added a model-and-sampler vignette that explains the mixture-of-factor-analyzers likelihood, PGMM covariance constraints, conjugate priors, allocation update, RJMCMC moves, and package output fields with formulas and runnable examples.
* Linked the new mathematical vignette from the getting-started guide and pkgdown vignette index.

# bpgmm 1.1.5

* Fixed the allocation prior contribution in RJMCMC acceptance calculations to use the log product of allocated mixture weights, matching the paper's joint posterior.

# bpgmm 1.1.4

* Standardized internal R helper names to snake_case while preserving the exported compatibility wrappers.
* Added snake_case wrappers around generated Rcpp entry points and routed package internals through those wrappers.
* Cleaned the source layout by renaming R files to lowercase hyphenated names and removing the ignored duplicate `R/choosem/` tree.
* Kept result-list names and legacy public arguments stable to avoid breaking existing user code.

# bpgmm 1.1.3

* Fixed `Mstep = 1` by loading the cluster-number proposal helpers as package code rather than leaving them in an ignored nested `R/` directory.
* Fixed split/combine RJMCMC moves by passing `X` explicitly to split-allocation helpers and adding the missing latent-factor update helper.
* Fixed the C++ allocation update to add `log(tao)` to component log densities instead of adding raw mixture weights.
* Added C++ input validation for multivariate normal densities, log-ratio calculations, allocation dimensions, and mixture weights.
* Expanded unit coverage for all eight PGMM covariance constraints, C++ helpers, cluster-number RJMCMC, covariance RJMCMC, and split/combine moves.

# bpgmm 1.1.2

* Modernized README with badges, installation guidance, paper citation, and model-constraint helper examples.
* Added `summarizePgmmRJMCMC()` as the correctly spelled summary function and retained `summerizePgmmRJMCMC()` for backward compatibility.
* Added tidyverse-style `pgmm_rjmcmc()` and `summarize_pgmm_rjmcmc()` as the preferred public API.
* Deprecated `pgmmRJMCMC()`, `summarizePgmmRJMCMC()`, and `summerizePgmmRJMCMC()`; these names remain available as compatibility wrappers.
* Added helpers `model_to_constraint()` and `constraint_to_model()` for translating between paper model labels and legacy constraint vectors.
* Improved package startup citation guidance for users publishing results from `bpgmm`.
* Added unit tests for the public API, covariance-constraint mapping, summary helpers, and native C++ wrappers.
* Fixed zero-iteration handling in `pgmm_rjmcmc()` and added validation for sampler inputs and summary result objects.
* Added a `verbose` argument to suppress per-iteration progress output in examples, tests, and scripted workflows.
* Updated C++ build settings and validation for current R/Rcpp best practices.
