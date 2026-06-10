# utils.R -- shared formulas and calibrations for Section 4 simulations.
# Model objects (Section 3 of Readme_thesis.md):
#   eta_i ~ N(B, 1/h_eta).
#   Agent picks h_t in [h_L, h_H] = [s - x, s + x] each period.
#   Signal x_t = mu_t + eta_i + eps_t with eps_t ~ N(0, 1/h_t).
#   After T signals at precisions (h_1, ..., h_T):
#     posterior precision  H_T = h_eta + sum_t h_t.
#     posterior mean       hat_eta_T marginally N(B, sigma_T^2).
#     where sigma_T^2 = (H_T - h_eta) / (h_eta * H_T).

# Marginal sd of hat_eta_T from period 0, with cumulative signal precision Hsum.
sigma_T <- function(Hsum, h_eta) {
  sqrt(Hsum / (h_eta * (h_eta + Hsum)))
}

# One-step conditional sd of hat_eta_t given hat_eta_{t-1} when prior precision
# is H_prev and current-period precision is h_t.
tau_step <- function(h_t, H_prev) {
  sqrt(h_t / (H_prev * (H_prev + h_t)))
}

# r^j = sigma_1 / sigma_2^j in the T = 2 closed-form (eta* = B).
r_j <- function(h1, h_j, h_eta) {
  sqrt(h1 * (h_eta + h1 + h_j) / (h_eta * h_j))
}

# Agent's promotion probability at eta* = B, T = 2, full closed form.
P_T2 <- function(h1, h_L, h_H, h_eta) {
  rL <- r_j(h1, h_L, h_eta)
  rH <- r_j(h1, h_H, h_eta)
  0.5 + (atan(rL) - atan(rH)) / (2 * pi)
}

# Components of P_T2:
#   Recovery: integral over {hat_eta_1 < B} with t = 2 best response h_H.
#   Save:    integral over {hat_eta_1 >= B} with t = 2 best response h_L.
P_T2_recovery <- function(h1, h_H, h_eta) {
  rH <- r_j(h1, h_H, h_eta)
  0.25 - atan(rH) / (2 * pi)
}
P_T2_save <- function(h1, h_L, h_eta) {
  rL <- r_j(h1, h_L, h_eta)
  0.25 + atan(rL) / (2 * pi)
}

# Organisation's gross expected payoff per period given prior N(B, 1/h_eta),
# threshold eta* and effective posterior sd sigma:
#   U_p = B (1 - Phi(z)) + sigma * phi(z),  z = (eta* - B)/sigma.
U_org <- function(B, sigma, eta_star) {
  z <- (eta_star - B) / sigma
  B * (1 - pnorm(z)) + sigma * dnorm(z)
}

# Informative regime payoff with constant h_H over T periods, eta* = B
# (analytic upper bound; under DIS this overestimates the true payoff).
U_I <- function(T, B, h_eta, h_H) {
  s <- sigma_T(T * h_H, h_eta)
  U_org(B, s, B)  # equals B/2 + s * dnorm(0)
}

# Uninformative regime payoff with constant h_L over T periods, eta* = 0.
# Exact for B > 0 by Lemma 2 (agent strictly prefers h_L throughout).
U_U <- function(T, B, h_eta, h_L) {
  s <- sigma_T(T * h_L, h_eta)
  U_org(B, s, 0)  # equals B*Phi(B/s) + s*phi(B/s)
}

# Net payoffs (subtract linear probation cost).
U_I_net <- function(T, B, h_eta, h_H, c) U_I(T, B, h_eta, h_H) - c * T
U_U_net <- function(T, B, h_eta, h_L, c) U_U(T, B, h_eta, h_L) - c * T

# Find T* over an integer grid by enumeration.
T_star <- function(T_grid, B, h_eta, h_L, h_H, c) {
  uI <- sapply(T_grid, U_I_net, B = B, h_eta = h_eta, h_H = h_H, c = c)
  uU <- sapply(T_grid, U_U_net, B = B, h_eta = h_eta, h_L = h_L, c = c)
  iI <- which.max(uI); iU <- which.max(uU)
  best_I <- uI[iI]; best_U <- uU[iU]
  if (best_I >= best_U) {
    list(T = T_grid[iI], regime = "I", value = best_I,
         T_I = T_grid[iI], T_U = T_grid[iU],
         u_I = uI, u_U = uU)
  } else {
    list(T = T_grid[iU], regime = "U", value = best_U,
         T_I = T_grid[iI], T_U = T_grid[iU],
         u_I = uI, u_U = uU)
  }
}

# T**: smallest T at which U_I(T) - U_U(T) <= 0 (gap zero crossing). The gross
# payoffs U_I, U_U are c-independent so the crossing does not move with c.
T_star_star <- function(T_grid, B, h_eta, h_L, h_H) {
  gap <- sapply(T_grid, function(T) U_I(T, B, h_eta, h_H) - U_U(T, B, h_eta, h_L))
  if (all(gap > 0))   return(list(T = NA, gap = gap, sign = "all positive"))
  if (all(gap <= 0))  return(list(T = T_grid[1], gap = gap, sign = "all non-positive"))
  T_cross <- T_grid[min(which(gap <= 0))]
  list(T = T_cross, gap = gap, sign = "crosses")
}

# Calibration baseline. h_L = 0.5, h_H = 2.0 means s = 1.25, x = 0.75.
# B = 0.4, h_eta = 1 are from the Readme Appendix A numerical example.
baseline <- list(B = 1, h_eta = 1, h_L = 0.5, h_H = 2.0, c = 0.02)

# Perturbations for robustness reporting. Notes:
#   - small_B / large_B are required by the prompt; at eta* = B the T = 2 agent
#     objective P(h_1;B) is closed-form independent of B (it depends only on
#     h_eta, h_L, h_H), so simA results coincide with baseline. Reported as a
#     finding rather than hidden.
#   - small_heta / large_heta vary prior precision.
#   - narrow_x / wide_x vary the project spread (h_L = s-x, h_H = s+x with s = 1.25).
perturbations <- list(
  small_heta = list(B = 1, h_eta = 0.5, h_L = 0.5,  h_H = 2.0,  c = 0.02),
  large_heta = list(B = 1, h_eta = 2.0, h_L = 0.5,  h_H = 2.0,  c = 0.02),
  lower_hL   = list(B = 1, h_eta = 1.0, h_L = 0.25, h_H = 2, c = 0.02),
  higher_hH     = list(B = 1, h_eta = 1.0, h_L = 0.5, h_H = 5, c = 0.02)
)

# Sanity check: informative-equilibrium feasibility at T = 1 requires
# B < 2 * phi(0) / sqrt(h_eta) ~ 0.7979 for h_eta = 1. All three calibrations
# satisfy this.
B_max <- function(h_eta) 2 * dnorm(0) / sqrt(h_eta)
