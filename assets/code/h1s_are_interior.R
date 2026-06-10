source("utils.R")

# Ensure the core function is loaded (Included here for completeness)
P_T2_general <- function(h1, eta_star, B, h_eta, h_L, h_H) {
  sigma1 <- sqrt(h1 / (h_eta * (h_eta + h1)))
  tau2_H <- sqrt(h_H / ((h_eta + h1) * (h_eta + h1 + h_H)))
  tau2_L <- sqrt(h_L / ((h_eta + h1) * (h_eta + h1 + h_L)))
  
  z_star <- (eta_star - B) / sigma1
  
  integrand_H <- function(z) {
    hat_eta_1 <- B + z * sigma1
    (1 - pnorm((eta_star - hat_eta_1) / tau2_H)) * dnorm(z)
  }
  
  integrand_L <- function(z) {
    hat_eta_1 <- B + z * sigma1
    (1 - pnorm((eta_star - hat_eta_1) / tau2_L)) * dnorm(z)
  }
  
  val_H <- integrate(integrand_H, lower = -Inf, upper = z_star)$value
  val_L <- integrate(integrand_L, lower = z_star, upper = Inf)$value
  
  return(val_H + val_L)
}

# =========================================================
# Parameters
# =========================================================
B       <- 0.4
h_eta   <- 1.0
h_L     <- 0.5
hH_vals <- c(2.0, 5.0)

# List to store the full dataframe for each h_H run
all_results <- list()

# =========================================================
# SIMULATION LOOP
# =========================================================
for (j in seq_along(hH_vals)) {
  h_H <- hH_vals[j]
  
  cat(sprintf("\nParameters: B = %.1f, h_eta = %.1f, h_L = %.1f, h_H = %.1f\n\n", B, h_eta, h_L, h_H))
  
  # Grid Sweep for eta*
  eta_grid <- seq(B - 0.5, B + 0.5, length.out = 101)
  
  # Note: Changed 'max_P' to 'P_prom' to match the get_vals function later
  results  <- data.frame(eta_star = eta_grid, delta = eta_grid - B, h1_star = NA, P_prom = NA)
  
  cat(sprintf("%10s %10s %10s %10s\n", "eta*", "delta", "h1*", "P_prom"))
  cat(strrep("-", 45), "\n", sep = "")
  
  for (i in seq_along(eta_grid)) {
    e_star <- eta_grid[i]
    
    # Optimize h1 in [h_L, h_H]
    res <- optimize(function(h1) P_T2_general(h1, e_star, B, h_eta, h_L, h_H), 
                    interval = c(h_L, h_H), maximum = TRUE, tol = 1e-6)
    
    results$h1_star[i] <- res$maximum
    results$P_prom[i]  <- res$objective
    
    # Print a subset of the grid to keep output clean
    if (i %% 5 == 1) {
      cat(sprintf("%10.3f %10.3f %10.4f %10.4f\n", 
                  e_star, results$delta[i], res$maximum, res$objective))
    }
  }
  
  # Save the completed dataframe into our list
  all_results[[j]] <- results
}

# =========================================================
# LATEX TABLE GENERATION
# =========================================================
eta_select <- c(0.90, 0.95, 1.00, 1.05, 1.10)

get_vals <- function(df, vals) {
  # Added round() to prevent floating point mismatch errors causing NA
  idx <- match(round(vals, 4), round(df$eta_star, 4))
  df$h1_star[idx]
}

row_h2 <- get_vals(all_results[[1]], eta_select)
row_h5 <- get_vals(all_results[[2]], eta_select)

latex_lines <- c(
  "\\begin{center}",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "$h_H$ & $0.90$ & $0.95$ & $1.00$ & $1.05$ & $1.10$ \\\\",
  "\\midrule",
  paste0(
    "$h_1^*(h_H = 2)$ & ",
    paste(sprintf("%.4f", row_h2), collapse = " & "),
    " \\\\"
  ),
  paste0(
    "$h_1^*(h_H = 5)$ & ",
    paste(sprintf("%.4f", row_h5), collapse = " & "),
    " \\\\"
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{center}"
)

cat("\n\n")
cat("====================================================\n")
cat("LATEX TABLE\n")
cat("====================================================\n\n")

cat(paste(latex_lines, collapse = "\n"))
cat("\n")

# =========================================================
# PERTURBATION ANALYSIS AT eta* = B
#
# Replaces simAprime_grid for the T = 2 case. Uses the closed
# form P_T2 (exact at eta* = B) and a cross-check via
# optimize() on the adaptive-quadrature P_T2_general, so the
# numbers are free of the Gauss-Hermite kink bias that
# contaminates simAprime_grid_results.
# =========================================================

cals <- c(list(baseline = baseline), perturbations)

cat("\n\n====================================================\n")
cat("PERTURBATIONS AT eta* = B  (T = 2)\n")
cat("====================================================\n\n")
cat(sprintf("%-12s %5s %5s %5s %5s   %8s %8s   %10s   %s\n",
            "label", "B", "h_eta", "h_L", "h_H",
            "h1*_cf", "h1*_opt", "P*", "interior"))
cat(strrep("-", 90), "\n", sep = "")

pert_rows <- list()
for (lab in names(cals)) {
  cal <- cals[[lab]]

  # Closed-form maximiser (B-invariant at eta* = B; depends only on h_eta, h_L, h_H).
  cf <- optimize(function(h) P_T2(h, cal$h_L, cal$h_H, cal$h_eta),
                 interval = c(cal$h_L, cal$h_H), maximum = TRUE, tol = 1e-10)

  # Independent check via the integrate-based P_T2_general at eta* = B.
  og <- optimize(function(h) P_T2_general(h, cal$B, cal$B, cal$h_eta, cal$h_L, cal$h_H),
                 interval = c(cal$h_L, cal$h_H), maximum = TRUE, tol = 1e-8)

  # Interior status: flag a corner if the closed-form maximiser is within
  # 1e-4 of either boundary.
  eps <- 1e-4
  interior <- if (cf$maximum < cal$h_L + eps) "corner @ h_L"
              else if (cf$maximum > cal$h_H - eps) "corner @ h_H"
              else "yes"

  cat(sprintf("%-12s %5.2f %5.2f %5.2f %5.2f   %8.5f %8.5f   %10.7f   %s\n",
              lab, cal$B, cal$h_eta, cal$h_L, cal$h_H,
              cf$maximum, og$maximum, cf$objective, interior))

  pert_rows[[lab]] <- list(B = cal$B, h_eta = cal$h_eta,
                           h_L = cal$h_L, h_H = cal$h_H,
                           h1_star = cf$maximum, P_star = cf$objective,
                           interior = interior)
}

# LaTeX table mirroring tab:simA_compstat (built as a single vector,
# printed in one block to match the style of latex_lines above).
fmt_row <- function(lab, r) {
  interior_str <- if (r$interior == "yes") "yes"
                  else if (r$interior == "corner @ h_L") "no (corner at $\\underline h$)"
                  else "no (corner at $\\bar h$)"
  sprintf("%s & %.1f & $(%.2f,\\,%.2f)$ & %.4f & %s \\\\",
          lab, r$h_eta, r$h_L, r$h_H, r$h1_star, interior_str)
}

pert_latex_lines <- c(
  "\\begin{table}",
  "\\centering",
  "\\caption{Optimal first-period precision at $\\eta^*=B$, single-parameter perturbations (continuous solver).}",
  "\\label{tab:simA_compstat}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "Perturbation & $h_\\eta$ & $(\\underline h,\\bar h)$ & $h_1^*$ & Interior \\\\",
  "\\midrule",
  vapply(names(pert_rows), function(lab) fmt_row(lab, pert_rows[[lab]]), character(1)),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

cat("\n\n")
cat("====================================================\n")
cat("LATEX TABLE  (perturbations at eta* = B, T = 2)\n")
cat("====================================================\n\n")
cat(paste(pert_latex_lines, collapse = "\n"))
cat("\n")

# Persist for downstream use.
saveRDS(pert_rows, "h1s_interior_perturbations.rds")
cat("\nWrote h1s_interior_perturbations.rds\n")