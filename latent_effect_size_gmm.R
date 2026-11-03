## =============================================================================
##  Standardized Latent Mean Differences from Binary Indicators by GMM
##  Complete analysis script: simulation study, empirical study, tables, figures
##
##  Output folders (created in the working directory)
##    output_simulation/   per-cell checkpoints (cell_###.rds) and summary.csv
##    output_empirical/    downloaded data, cleaned data, estimates, item tests
##    output_tables/       CSV and Markdown tables for the Results section
##    output_figures/      PNG figures (300 dpi)
##
##  Sections
##    0  configuration
##    1  population and data generation (simulation)
##    2  model probabilities and analytic gradients
##    3  model specification with equality constraints (partial invariance)
##    4  estimators: iterated GMM, delta-method effect size, robust variance,
##       J test, item-wise C (subset-moment) and D (restriction) tests,
##       sequential exclusion and partial-invariance procedures, comparators
##    5  simulation driver with checkpointing and stage-wise failure recording
##    6  empirical study: Rosenberg Self-Esteem Scale (Open Psychometrics)
##    7  tables
##    8  figures
##    9  main
##
##  HOW TO RUN (one clean run of everything):
##      setwd("<folder containing this file>"); source("latent_effect_size_gmm.R")
##  or from a terminal:  Rscript latent_effect_size_gmm.R > run_log.txt 2>&1
##  With cfg$clean_start = TRUE the four output folders are deleted first, a short
##  consistency check of the three latent estimators is run on one dataset, and then
##  the simulation, the empirical study, the tables and the figures are produced.
##  Expected duration: roughly 10-14 hours on 9 cores. Interrupted runs resume from
##  the last completed cell if clean_start is set to FALSE before re-sourcing.
##
##  Requires: pbivnorm.  Optional: lavaan (WLSMV comparator), mirt (MML comparator),
##  ggplot2 (figures), parallel (multicore on unix).
##  Design: n in {500, 1000, 2000}, J in {5, 10, 20},
##  d in {0, 0.4}, four conditions (invariance, uniform DIF, nonuniform DIF,
##  skewed latent trait), 1000 replications per cell. Replication seeds are
##  unique per (cell, replication) pair.
##  Comparators: observed d (Hedges), WLSMV (lavaan), MML 2PL (mirt).
##  Post-selection inference: sample-splitting versions of both step-2 procedures.
## =============================================================================

need <- function(p) { if (!requireNamespace(p, quietly = TRUE)) install.packages(p); requireNamespace(p, quietly = TRUE) }
stopifnot(need("pbivnorm")); library(pbivnorm)
HAS_LAVAAN <- need("lavaan"); HAS_MIRT <- need("mirt"); HAS_GG <- need("ggplot2")
if (HAS_GG) library(ggplot2)

## ============================ 0. Configuration ===============================
cfg <- list(
  clean_start = TRUE,                  # delete previous outputs and start from scratch
  run_simulation = TRUE, run_empirical = TRUE, make_outputs = TRUE,
  n_rep = 1000, seed = 20261030,
  n_grid = c(500, 1000, 2000), J_grid = c(5, 10, 20), d_grid = c(0, 0.4),
  dif_grid = c("none", "uniform", "nonuniform", "skewlatent"),
  dif_prop = 0.20, dif_shift = 0.5, dif_load = 0.7, latent_skew = 0.5,
  drop_small_n_long_test = TRUE,       # skip n = 500 & J = 20 (documented boundary cells)
  N_pseudo = 200000, N_pseudo_cmp = 50000, use_lavaan = HAS_LAVAAN, use_mirt = HAS_MIRT, alpha = 0.05, min_anchor = 2,
  emp_url = "https://openpsychometrics.org/_rawdata/RSE.zip",
  emp_neg_items = c(3, 5, 8, 9, 10),   # verify against the downloaded codebook
  emp_subsamples = 200, emp_sub_n = 500,
  cores = max(1L, parallel::detectCores() - 1L),
  dir_sim = "output_simulation", dir_emp = "output_empirical",
  dir_tab = "output_tables", dir_fig = "output_figures"
)
if (isTRUE(cfg$clean_start) && !isTRUE(getOption("es.skip_main")))
  unlink(c(cfg$dir_sim, cfg$dir_emp, cfg$dir_tab, cfg$dir_fig), recursive = TRUE)
for (d in c(cfg$dir_sim, cfg$dir_emp, cfg$dir_tab, cfg$dir_fig)) dir.create(d, showWarnings = FALSE)
design <- expand.grid(n = cfg$n_grid, J = cfg$J_grid, d = cfg$d_grid, dif = cfg$dif_grid, stringsAsFactors = FALSE)
if (cfg$drop_small_n_long_test) design <- design[!(design$n == 500 & design$J == 20), ]
rownames(design) <- NULL

## ======================= 1. Population and data ==============================
item_params <- function(J) {
  set.seed(1000 + J)
  list(lam = sample(seq(0.55, 0.85, length.out = J)), tau = seq(-0.8, 0.8, length.out = J),
       dif_items = unique(round(seq(2, J - 1, length.out = max(1, ceiling(cfg$dif_prop * J))))))
}
pair_index <- function(J) { pr <- combn(J, 2); list(pj = pr[1, ], pl = pr[2, ]) }
## Standardized skew-normal latent trait with a given skewness (Azzalini, 1985)
skew_delta <- function(gamma) uniroot(function(dl) { b <- dl * sqrt(2 / pi); (4 - pi) / 2 * b^3 / (1 - b^2)^1.5 - gamma }, c(1e-6, 0.995))$root
SKEW_DELTA <- skew_delta(cfg$latent_skew)
r_latent <- function(m, dif) {
  if (dif != "skewlatent") return(rnorm(m))
  z <- SKEW_DELTA * abs(rnorm(m)) + sqrt(1 - SKEW_DELTA^2) * rnorm(m); mu <- SKEW_DELTA * sqrt(2 / pi)
  (z - mu) / sqrt(1 - mu^2)
}
## Moment structure S from an indicator matrix and group vector (used by both studies)
mom_from_data <- function(X, g) {
  X <- as.matrix(X); J <- ncol(X); P <- pair_index(J)
  I <- cbind(X, X[, P$pj, drop = FALSE] * X[, P$pl, drop = FALSE])
  list(J = J, NM = ncol(I), NPH = 4 * J + 2, pj = P$pj, pl = P$pl, n = length(g),
       w = c(mean(g == 0), mean(g == 1)),
       Ibar = lapply(0:1, function(gg) colMeans(I[g == gg, , drop = FALSE])),
       M2 = lapply(0:1, function(gg) crossprod(I[g == gg, , drop = FALSE]) / sum(g == gg)),
       I = I, g = g, ix = list(tau = list(1:J, 2 * J + 1:J), lam = list(J + 1:J, 3 * J + 1:J), mu = 4 * J + 1, v = 4 * J + 2))
}
## Simulated data, accumulated in chunks (keep_I = FALSE for very large samples)
gen_moments <- function(n_per, J, d, dif, ip, seed, keep_I = TRUE, chunk = 20000) {
  set.seed(seed); P <- pair_index(J); NM <- J + length(P$pj)
  lam_g <- list(ip$lam, ip$lam); shift_g <- list(numeric(J), numeric(J))
  if (dif == "uniform")    shift_g[[1]][ip$dif_items] <- -cfg$dif_shift
  if (dif == "nonuniform") lam_g[[1]][ip$dif_items]   <- cfg$dif_load * ip$lam[ip$dif_items]
  Ibar <- vector("list", 2); M2 <- vector("list", 2); I_all <- NULL; g_all <- NULL
  for (gg in 0:1) {
    s1 <- numeric(NM); s2 <- matrix(0, NM, NM); done <- 0
    while (done < n_per) {
      m <- min(chunk, n_per - done); eta <- (if (gg == 1) d else 0) + r_latent(m, dif)
      lam <- lam_g[[gg + 1]]; sh <- shift_g[[gg + 1]]
      ystar <- outer(eta, lam) + matrix(rnorm(m * J), m, J) %*% diag(sqrt(1 - lam^2), J) + matrix(sh, m, J, byrow = TRUE)
      X <- (ystar > matrix(ip$tau, m, J, byrow = TRUE)) * 1L
      I <- cbind(X, X[, P$pj, drop = FALSE] * X[, P$pl, drop = FALSE])
      s1 <- s1 + colSums(I); s2 <- s2 + crossprod(I); done <- done + m
      if (keep_I) { I_all <- rbind(I_all, I); g_all <- c(g_all, rep(gg, m)) }
    }
    Ibar[[gg + 1]] <- s1 / n_per; M2[[gg + 1]] <- s2 / n_per
  }
  list(J = J, NM = NM, NPH = 4 * J + 2, pj = P$pj, pl = P$pl, n = 2 * n_per, w = c(0.5, 0.5),
       Ibar = Ibar, M2 = M2, I = I_all, g = g_all,
       ix = list(tau = list(1:J, 2 * J + 1:J), lam = list(J + 1:J, 3 * J + 1:J), mu = 4 * J + 1, v = 4 * J + 2))
}

## ================= 2. Model probabilities and analytic gradients ==============
## phi = (tau0[J], lam0[J], tau1[J], lam1[J], mu1, sig1sq); residual var = 1 - lam0^2
grp_par <- function(phi, g, S) list(tau = phi[S$ix$tau[[g + 1]]], lam = phi[S$ix$lam[[g + 1]]], lam0 = phi[S$ix$lam[[1]]],
                                    mu = if (g == 1) phi[S$ix$mu] else 0, v = if (g == 1) phi[S$ix$v] else 1)
probs <- function(phi, g, S) {
  p <- grp_par(phi, g, S); s <- sqrt(p$lam^2 * p$v + 1 - p$lam0^2); a <- (p$lam * p$mu - p$tau) / s
  rho <- p$lam[S$pj] * p$lam[S$pl] * p$v / (s[S$pj] * s[S$pl])
  c(pnorm(a), pbivnorm(a[S$pj], a[S$pl], rho))
}
phi2 <- function(x, y, r) exp(-(x * x - 2 * r * x * y + y * y) / (2 * (1 - r^2))) / (2 * pi * sqrt(1 - r^2))
jac_probs <- function(phi, g, S) {
  J <- S$J; NPH <- S$NPH; p <- grp_par(phi, g, S); lam <- p$lam; lam0 <- p$lam0; tau <- p$tau; mu <- p$mu; v <- p$v
  s2 <- lam^2 * v + 1 - lam0^2; s <- sqrt(s2); a <- (lam * mu - tau) / s
  it <- S$ix$tau[[g + 1]]; il <- S$ix$lam[[g + 1]]; il0 <- S$ix$lam[[1]]
  dA <- matrix(0, J, NPH)
  dA[cbind(1:J, it)] <- -1 / s
  dA[cbind(1:J, il)] <- dA[cbind(1:J, il)] + mu / s - (lam * mu - tau) * lam * v / s^3
  dA[cbind(1:J, il0)] <- dA[cbind(1:J, il0)] + (lam * mu - tau) * lam0 / s^3
  if (g == 1) { dA[, S$ix$mu] <- lam / s; dA[, S$ix$v] <- -(lam * mu - tau) * lam^2 / (2 * s^3) }
  pj <- S$pj; pl <- S$pl; K <- length(pj); rho <- lam[pj] * lam[pl] * v / (s[pj] * s[pl])
  dR <- matrix(0, K, NPH)
  dR[cbind(1:K, il[pj])] <- dR[cbind(1:K, il[pj])] + rho * (1 / lam[pj] - lam[pj] * v / s2[pj])
  dR[cbind(1:K, il[pl])] <- dR[cbind(1:K, il[pl])] + rho * (1 / lam[pl] - lam[pl] * v / s2[pl])
  dR[cbind(1:K, il0[pj])] <- dR[cbind(1:K, il0[pj])] + rho * lam0[pj] / s2[pj]
  dR[cbind(1:K, il0[pl])] <- dR[cbind(1:K, il0[pl])] + rho * lam0[pl] / s2[pl]
  if (g == 1) dR[, S$ix$v] <- rho / v - rho * (lam[pj]^2 / (2 * s2[pj]) + lam[pl]^2 / (2 * s2[pl]))
  G1 <- dnorm(a) * dA; x <- a[pj]; y <- a[pl]; q <- sqrt(1 - rho^2)
  Fx <- dnorm(x) * pnorm((y - rho * x) / q); Fy <- dnorm(y) * pnorm((x - rho * y) / q); Fr <- phi2(x, y, rho)
  rbind(G1, Fx * dA[pj, , drop = FALSE] + Fy * dA[pl, , drop = FALSE] + Fr * dR)
}
mbar <- function(phi, S) c((S$Ibar[[1]] - probs(phi, 0, S)) * S$w[1], (S$Ibar[[2]] - probs(phi, 1, S)) * S$w[2])
Gphi <- function(phi, S) rbind(-S$w[1] * jac_probs(phi, 0, S), -S$w[2] * jac_probs(phi, 1, S))
Wmat <- function(phi, S) {
  NM <- S$NM; W <- matrix(0, 2 * NM, 2 * NM)
  for (gg in 0:1) { p <- probs(phi, gg, S); Ib <- S$Ibar[[gg + 1]]
    W[gg * NM + 1:NM, gg * NM + 1:NM] <- S$w[gg + 1] * (S$M2[[gg + 1]] - outer(p, Ib) - outer(Ib, p) + outer(p, p)) }
  W
}
u_mat <- function(phi, S) {
  U <- matrix(0, S$n, 2 * S$NM)
  for (gg in 0:1) { sel <- S$g == gg; U[sel, gg * S$NM + 1:S$NM] <- sweep(S$I[sel, , drop = FALSE], 2, probs(phi, gg, S)) }
  U
}
numjac <- function(f, x, eps = 1e-5) {
  f0 <- f(x); Jm <- matrix(0, length(f0), length(x))
  for (k in seq_along(x)) { e <- numeric(length(x)); e[k] <- eps; Jm[, k] <- (f(x + e) - f(x - e)) / (2 * eps) }
  Jm
}

## ================== 3. Specification with equality constraints ================
item_moments <- function(j, S) { idx <- c(j, S$J + which(S$pj == j | S$pl == j)); c(idx, S$NM + idx) }
make_spec <- function(S, untied = integer(0), excl = integer(0)) {
  keep <- 1:(2 * S$NM); for (j in excl) keep <- setdiff(keep, item_moments(j, S))
  list(untied = sort(untied), excl = sort(excl), keep = keep)
}
build_map <- function(spec, S, phi_ref) {
  J <- S$J; NPH <- S$NPH; ix <- S$ix; active <- setdiff(1:J, spec$excl)
  free <- c(ix$tau[[1]][active], ix$lam[[1]][active], ix$tau[[2]][intersect(active, spec$untied)],
            ix$lam[[2]][intersect(active, spec$untied)], ix$mu, ix$v)
  M <- matrix(0, NPH, length(free)); M[cbind(free, seq_along(free))] <- 1
  tied <- setdiff(active, spec$untied)
  M[cbind(ix$tau[[2]][tied], match(ix$tau[[1]][tied], free))] <- 1
  M[cbind(ix$lam[[2]][tied], match(ix$lam[[1]][tied], free))] <- 1
  c0 <- phi_ref; c0[rowSums(M) > 0] <- 0
  type <- rep("id", length(free)); type[free %in% c(ix$lam[[1]], ix$lam[[2]])] <- "lam"; type[free == ix$v] <- "var"
  list(M = M, c0 = c0, free = free, type = type)
}
nat_from_z <- function(z, type) ifelse(type == "lam", plogis(z), ifelse(type == "var", exp(z), z))
z_from_nat <- function(pf, type) ifelse(type == "lam", qlogis(pf), ifelse(type == "var", log(pf), pf))
dnat_dz    <- function(pf, type) ifelse(type == "lam", pf * (1 - pf), ifelse(type == "var", pf, 1))
start_phi <- function(S) { p0 <- pmin(pmax(S$Ibar[[1]][1:S$J], .02), .98); tau <- -qnorm(p0); c(tau, rep(.6, S$J), tau, rep(.6, S$J), .2, 1) }

## ============================ 4. Estimators ==================================
gmm_fit <- function(S, spec, phi_init = NULL, Wfix = NULL, maxit = 50) {
  if (is.null(phi_init)) phi_init <- start_phi(S)
  mp <- build_map(spec, S, phi_init); M <- mp$M; c0 <- mp$c0; keep <- spec$keep; nk <- length(keep)
  phi_of <- function(pf) drop(M %*% pf) + c0
  pf <- phi_init[mp$free]; phi <- phi_of(pf)
  for (it in 1:maxit) {
    W <- if (is.null(Wfix)) Wmat(phi, S)[keep, keep] else Wfix; Winv <- solve(W)
    crit <- function(z) { m <- mbar(phi_of(nat_from_z(z, mp$type)), S)[keep]; drop(crossprod(m, Winv %*% m)) }
    grad <- function(z) { pfz <- nat_from_z(z, mp$type); ph <- phi_of(pfz); m <- mbar(ph, S)[keep]
      Gz <- (Gphi(ph, S)[keep, , drop = FALSE] %*% M) * rep(dnat_dz(pfz, mp$type), each = nk)
      drop(2 * crossprod(Gz, Winv %*% m)) }
    res <- optim(z_from_nat(pf, mp$type), crit, grad, method = "BFGS", control = list(maxit = 2000, reltol = 1e-14))
    pf_new <- nat_from_z(res$par, mp$type); conv <- max(abs(pf_new - pf)) < 1e-8; pf <- pf_new; phi <- phi_of(pf)
    if (conv || !is.null(Wfix)) break
  }
  m <- mbar(phi, S)[keep]; Jstat <- S$n * drop(crossprod(m, Winv %*% m)); df <- nk - length(mp$free)
  Gk <- Gphi(phi, S)[keep, , drop = FALSE]; Gf <- Gk %*% M
  Vf <- solve(crossprod(Gf, Winv %*% Gf)) / S$n
  list(phi = phi, J = Jstat, df = df, pJ = 1 - pchisq(Jstat, df), Gk = Gk, M = M, W = W, Winv = Winv,
       V = M %*% Vf %*% t(M), spec = spec, map = mp, iters = it, converged = conv)
}
effect_size <- function(phi, V, S) {
  m1 <- phi[S$ix$mu]; v1 <- phi[S$ix$v]; sd_p <- sqrt((1 + v1) / 2)
  gr <- numeric(S$NPH); gr[S$ix$mu] <- 1 / sd_p; gr[S$ix$v] <- -m1 / (4 * sd_p^3)
  c(d = m1 / sd_p, se = sqrt(drop(crossprod(gr, V %*% gr))), cl = pnorm(m1 / sqrt(1 + v1)))
}
robust_V <- function(fit, S) {                 # Hansen & Lee (2021)
  n <- S$n; keep <- fit$spec$keep; M <- fit$M; c0 <- fit$map$c0; pf <- fit$phi[fit$map$free]
  phi_of <- function(p) drop(M %*% p) + c0
  Ffun <- function(p) { ph <- phi_of(p); m <- mbar(ph, S)[keep]; W <- Wmat(ph, S)[keep, keep]
    G <- Gphi(ph, S)[keep, , drop = FALSE] %*% M; drop(crossprod(G, solve(W, m))) }
  H <- numjac(Ffun, pf, eps = 1e-4)
  U <- u_mat(fit$phi, S)[, keep, drop = FALSE]; W <- Wmat(fit$phi, S)[keep, keep]; Winv <- solve(W); m <- colMeans(U)
  G <- fit$Gk %*% M; GW <- crossprod(G, Winv); Wm <- Winv %*% m
  qg <- sapply(0:1, function(gg) {
    Qi <- rbind(-jac_probs(fit$phi, 0, S) * (gg == 0), -jac_probs(fit$phi, 1, S) * (gg == 1))[keep, , drop = FALSE] %*% M
    drop(crossprod(Qi, Wm)) })
  A <- U %*% t(GW); psi <- A + t(qg[, S$g + 1, drop = FALSE]) - A * drop(U %*% Wm)
  Om <- crossprod(psi) / n; Hinv <- solve(H)
  M %*% (Hinv %*% Om %*% t(Hinv) / n) %*% t(M)
}
c_tests <- function(fit, S, items) {           # subset-moment tests (Eichenbaum, Hansen & Singleton 1988)
  keep <- fit$spec$keep
  out <- sapply(items, function(j) {
    sp2 <- make_spec(S, fit$spec$untied, c(fit$spec$excl, j)); pos <- match(sp2$keep, keep)
    sub <- gmm_fit(S, sp2, phi_init = fit$phi, Wfix = fit$W[pos, pos])
    C <- fit$J - sub$J; dfC <- (length(keep) - length(sp2$keep)) - (length(fit$map$free) - length(sub$map$free))
    c(stat = C, df = dfC, p = 1 - pchisq(C, dfC)) })
  colnames(out) <- items; rbind(out, p_holm = p.adjust(out["p", ], method = "holm"))
}
d_tests <- function(fit, S, items) {           # restriction tests (Newey & West 1987)
  out <- sapply(items, function(j) {
    sp2 <- make_spec(S, c(fit$spec$untied, j), fit$spec$excl)
    un <- gmm_fit(S, sp2, phi_init = fit$phi, Wfix = fit$W)
    D <- fit$J - un$J; c(stat = D, df = 2, p = 1 - pchisq(D, 2)) })
  colnames(out) <- items; rbind(out, p_holm = p.adjust(out["p", ], method = "holm"))
}
sequential_exclusion <- function(fit, S) {
  cur <- fit
  repeat {
    items <- setdiff(1:S$J, cur$spec$excl)
    if (cur$pJ >= cfg$alpha || length(items) < 3) break
    ct <- c_tests(cur, S, items); k <- which.min(ct["p_holm", ])
    if (ct["p_holm", k] >= cfg$alpha) break
    cur <- gmm_fit(S, make_spec(S, cur$spec$untied, c(cur$spec$excl, items[k])), phi_init = cur$phi)
  }
  list(fit = cur, set = cur$spec$excl)
}
sequential_partial <- function(fit, S) {
  cur <- fit
  repeat {
    tied <- setdiff(1:S$J, c(cur$spec$untied, cur$spec$excl))
    if (cur$pJ >= cfg$alpha || length(tied) <= cfg$min_anchor) break
    dt <- d_tests(cur, S, tied); k <- which.min(dt["p_holm", ])
    if (dt["p_holm", k] >= cfg$alpha) break
    cur <- gmm_fit(S, make_spec(S, c(cur$spec$untied, tied[k]), cur$spec$excl), phi_init = cur$phi)
  }
  list(fit = cur, set = cur$spec$untied)
}
observed_d <- function(S) {
  J <- S$J; Ssc <- rowSums(S$I[, 1:J, drop = FALSE]); g <- S$g
  a <- Ssc[g == 1]; b <- Ssc[g == 0]; na <- length(a); nb <- length(b)
  sp <- sqrt(((na - 1) * var(a) + (nb - 1) * var(b)) / (na + nb - 2)); d <- (mean(a) - mean(b)) / sp
  cJ <- 1 - 3 / (4 * (na + nb) - 9); p <- colMeans(S$I[, 1:J, drop = FALSE])
  c(d = cJ * d, se = cJ * sqrt(1 / na + 1 / nb + d^2 / (2 * (na + nb))), kr20 = J / (J - 1) * (1 - sum(p * (1 - p)) / var(Ssc)))
}
wlsmv_d <- function(S) {
  if (!cfg$use_lavaan) return(c(d = NA, se = NA))
  J <- S$J; df <- as.data.frame(S$I[, 1:J, drop = FALSE]); names(df) <- paste0("x", 1:J)
  df$g <- factor(S$g, levels = c(0, 1), labels = c("g0", "g1"))   # explicit order: g0 is the reference group
  items <- names(df)[1:J]
  ## Theta parameterization. All fixed/free settings are given as NUMERIC vectors:
  ## a mixed vector such as c(0, m1) is coerced to character and read by lavaan as labels.
  ## Residual variances and latent-response intercepts are fixed in both groups; with
  ## binary items this is required for the group-2 mean and variance to be identified
  ## (Millsap & Yun-Tein, 2004). Parameters are read back by free-parameter index.
  model <- paste(c(paste0("F =~ NA*x1 + ", paste(items[-1], collapse = " + ")),
                   "F ~ c(0, NA)*1",                        # latent mean: 0 in g0, free in g1
                   "F ~~ c(1, NA)*F",                       # latent variance: 1 in g0, free in g1
                   paste0(items, " ~~ c(1, 1)*", items),    # residual variances = 1 in both groups
                   paste0(items, " ~ c(0, 0)*1")),          # intercepts = 0 in both groups
                 collapse = "\n")
  fit <- lavaan::cfa(model, data = df, group = "g", group.label = c("g0", "g1"),   # fix group order
                     ordered = items, estimator = "WLSMV",
                     parameterization = "theta", group.equal = c("loadings", "thresholds"))
  pt <- lavaan::parTable(fit)
  i_mu <- pt$free[pt$lhs == "F" & pt$op == "~1" & pt$group == 2]
  i_v  <- pt$free[pt$lhs == "F" & pt$op == "~~" & pt$rhs == "F" & pt$group == 2]
  if (length(i_mu) != 1 || length(i_v) != 1 || i_mu == 0 || i_v == 0) stop("wlsmv_d: group-2 mean/variance not free")
  th <- lavaan::coef(fit); V <- lavaan::lavInspect(fit, "vcov")
  mu <- unname(th[i_mu]); v <- unname(th[i_v]); Vs <- V[c(i_mu, i_v), c(i_mu, i_v)]
  sd_p <- sqrt((1 + v) / 2); gr <- c(1 / sd_p, -mu / (4 * sd_p^3))
  c(d = mu / sd_p, se = sqrt(drop(crossprod(gr, Vs %*% gr))))
}
## Full-information MML comparator: multiple-group 2PL (Bock & Aitkin, 1981) with
## invariant slopes and intercepts; group 1 fixed N(0,1), group 2 mean/variance free.
## Note the logistic link; the latent-scale effect size is comparable to the probit model.
mml_d <- function(S) {
  if (!cfg$use_mirt) return(c(d = NA, se = NA))
  J <- S$J; dat <- as.data.frame(S$I[, 1:J, drop = FALSE]); names(dat) <- paste0("x", 1:J)
  grp <- factor(S$g, levels = c(0, 1), labels = c("g0", "g1"))
  fit <- mirt::multipleGroup(dat, 1, group = grp, itemtype = "2PL", SE = TRUE, verbose = FALSE,
                             invariance = c("slopes", "intercepts", "free_means", "free_var"))
  mv <- mirt::mod2values(fit)
  sel <- mv[mv$group == "g1" & mv$item == "GROUP" & mv$name %in% c("MEAN_1", "COV_11"), ]
  sel <- sel[match(c("MEAN_1", "COV_11"), sel$name), ]
  V <- mirt::extract.mirt(fit, "vcov"); nm <- colnames(V)
  i_mu <- grep("^MEAN_1\\.", nm); i_v <- grep("^COV_11\\.", nm)
  if (length(i_mu) != 1 || length(i_v) != 1) stop("mml_d: cannot locate MEAN_1/COV_11 in vcov")
  Vs <- V[c(i_mu, i_v), c(i_mu, i_v)]; mu <- sel$value[1]; v <- sel$value[2]
  sd_p <- sqrt((1 + v) / 2); gr <- c(1 / sd_p, -mu / (4 * sd_p^3))
  c(d = mu / sd_p, se = sqrt(drop(crossprod(gr, Vs %*% gr))))
}
## Sample splitting (Cox, 1975): select items on half A, estimate on half B.
## Inference on half B is conditionally valid given the selection.
split_fit <- function(S, seed) {
  set.seed(seed); i0 <- which(S$g == 0); i1 <- which(S$g == 1)
  A <- c(sample(i0, floor(length(i0) / 2)), sample(i1, floor(length(i1) / 2))); B <- setdiff(seq_len(S$n), A)
  SA <- mom_from_data(S$I[A, 1:S$J, drop = FALSE], S$g[A]); SB <- mom_from_data(S$I[B, 1:S$J, drop = FALSE], S$g[B])
  fA <- gmm_fit(SA, make_spec(SA)); exA <- sequential_exclusion(fA, SA); paA <- sequential_partial(fA, SA)
  fBx <- gmm_fit(SB, make_spec(SB, excl = exA$set)); fBp <- gmm_fit(SB, make_spec(SB, untied = paA$set))
  list(ex = list(es = effect_size(fBx$phi, fBx$V, SB), set = exA$set),
       pa = list(es = effect_size(fBp$phi, fBp$V, SB), set = paA$set))
}
## Complete analysis of one dataset (shared by simulation and empirical study)
analyze_dataset <- function(S, true_set = NULL, split_seed = 1) {
  od <- observed_d(S)
  fit <- gmm_fit(S, make_spec(S)); es <- effect_size(fit$phi, fit$V, S)
  errs <- character(0)
  stg <- function(expr, tag) tryCatch(expr, error = function(e) { errs[tag] <<- conditionMessage(e); NULL })
  Vr <- stg(robust_V(fit, S), "rob"); es_r <- if (is.null(Vr)) c(d = NA, se = NA, cl = NA) else effect_size(fit$phi, Vr, S)
  wl <- stg(wlsmv_d(S), "wlsmv"); if (is.null(wl)) wl <- c(d = NA, se = NA)
  ml <- stg(mml_d(S), "mml"); if (is.null(ml)) ml <- c(d = NA, se = NA)
  ex <- stg(sequential_exclusion(fit, S), "excl"); pa <- stg(sequential_partial(fit, S), "part")
  sp <- stg(split_fit(S, split_seed), "split")
  es_ex <- if (is.null(ex)) c(d = NA, se = NA) else effect_size(ex$fit$phi, ex$fit$V, S)
  es_pa <- if (is.null(pa)) c(d = NA, se = NA) else effect_size(pa$fit$phi, pa$fit$V, S)
  xs <- function(set) if (is.null(sp) || is.null(true_set)) NA else setequal(as.integer(set), true_set)
  list(fit = fit, Vr = Vr, ex = ex, pa = pa, sp = sp, errors = errs,
       row = data.frame(d_obs = od[["d"]], se_obs = od[["se"]], kr20 = od[["kr20"]],
             d_wlsmv = wl[["d"]], se_wlsmv = wl[["se"]], d_mml = ml[["d"]], se_mml = ml[["se"]],
             d_sx = if (is.null(sp)) NA else sp$ex$es[["d"]], se_sx = if (is.null(sp)) NA else sp$ex$es[["se"]], exact_sx = xs(if (is.null(sp)) NA else sp$ex$set),
             d_sp = if (is.null(sp)) NA else sp$pa$es[["d"]], se_sp = if (is.null(sp)) NA else sp$pa$es[["se"]], exact_sp = xs(if (is.null(sp)) NA else sp$pa$set),
             fail_mml = all(is.na(ml)), fail_split = is.null(sp),
             d_gmm = es[["d"]], se_conv = es[["se"]], se_rob = es_r[["se"]], cl = es[["cl"]],
             Jstat = fit$J, df = fit$df, pJ = fit$pJ, converged = fit$converged,
             err_msg = if (length(errs)) paste(names(errs), errs, sep = ": ", collapse = " | ") else "",
             d_excl = es_ex[["d"]], se_excl = es_ex[["se"]], n_excl = if (is.null(ex)) NA else length(ex$set),
             exact_excl = if (is.null(ex) || is.null(true_set)) NA else setequal(as.integer(ex$set), true_set),
             d_part = es_pa[["d"]], se_part = es_pa[["se"]], n_untied = if (is.null(pa)) NA else length(pa$set),
             exact_part = if (is.null(pa) || is.null(true_set)) NA else setequal(as.integer(pa$set), true_set),
             fail_rob = is.null(Vr), fail_wlsmv = all(is.na(wl)), fail_excl = is.null(ex), fail_part = is.null(pa)))
}

## ========================== 5. Simulation driver =============================
## Pseudo-true values: each estimator has its own population limit under misspecification,
## so d* is computed separately for GMM, WLSMV and MML (White, 1982). The GMM limit uses
## N_pseudo per group from moment summaries; the comparators need raw data and use N_pseudo_cmp.
pseudo_true <- function(J, d, dif, ip) {
  if (dif == "none") return(c(d_star = d, mc_se = 0, d_star_wlsmv = d, d_star_mml = d))
  sd0 <- cfg$seed + 99000000L + 1000L * J + 10L * match(dif, cfg$dif_grid) + as.integer(d > 0)
  S <- gen_moments(cfg$N_pseudo, J, d, dif, ip, seed = sd0, keep_I = FALSE)
  fit <- gmm_fit(S, make_spec(S)); es <- effect_size(fit$phi, fit$V, S)
  Sc <- gen_moments(cfg$N_pseudo_cmp, J, d, dif, ip, seed = sd0 + 1L, keep_I = TRUE)
  wl <- tryCatch(wlsmv_d(Sc), error = function(e) c(d = NA, se = NA))
  ml <- tryCatch(mml_d(Sc),   error = function(e) c(d = NA, se = NA))
  c(d_star = es[["d"]], mc_se = es[["se"]], d_star_wlsmv = wl[["d"]], d_star_mml = ml[["d"]])
}
one_rep <- function(r, cell, ip, i) {
  S <- gen_moments(cell$n, cell$J, cell$d, cell$dif, ip, seed = cfg$seed + 100000L * i + r)
  true_set <- if (cell$dif %in% c("none", "skewlatent")) integer(0) else as.integer(ip$dif_items)
  cbind(rep = r, analyze_dataset(S, true_set, split_seed = cfg$seed + 100000L * i + r + 50000L)$row)
}
run_cell <- function(i) {
  cell <- design[i, ]; f_out <- file.path(cfg$dir_sim, sprintf("cell_%03d.rds", i))
  if (file.exists(f_out)) { message(sprintf("cell %d/%d: cached", i, nrow(design))); return(invisible(readRDS(f_out))) }
  ip <- item_params(cell$J); pt <- pseudo_true(cell$J, cell$d, cell$dif, ip)
  message(sprintf("cell %d/%d: n=%d J=%d d=%.1f dif=%s | d*=%.4f (mc se %.4f) | %d reps, %d core(s) [%s]",
                  i, nrow(design), cell$n, cell$J, cell$d, cell$dif, pt[1], pt[2], cfg$n_rep, cfg$cores, format(Sys.time(), "%H:%M:%S")))
  f <- function(r) tryCatch(one_rep(r, cell, ip, i), error = function(e) structure(conditionMessage(e), class = "sim_error"))
  reps <- if (cfg$cores > 1 && .Platform$OS.type == "unix") parallel::mclapply(seq_len(cfg$n_rep), f, mc.cores = cfg$cores, mc.set.seed = FALSE)
          else lapply(seq_len(cfg$n_rep), f)
  ok <- vapply(reps, is.data.frame, logical(1)); errs <- unlist(lapply(reps[!ok], as.character))
  if (!any(ok)) stop("cell ", i, ": all replications failed: ", errs[1])
  if (any(!ok)) message(sprintf("  main-fit failures: %d (%.1f%%)", sum(!ok), 100 * mean(!ok)))
  res <- cbind(cell[rep(1, sum(ok)), ], do.call(rbind, reps[ok]), d_star = pt[["d_star"]], d_star_se = pt[["mc_se"]],
               d_star_wlsmv = pt[["d_star_wlsmv"]], d_star_mml = pt[["d_star_mml"]])
  out <- list(results = res, cell = cell, pseudo_true = pt, item_params = ip, n_failed_main = sum(!ok), n_attempted = cfg$n_rep, errors = errs)
  saveRDS(out, f_out); invisible(out)
}
row_stats <- function(est, se, target) {
  ok <- !is.na(est) & !is.na(se); est <- est[ok]; se <- se[ok]
  if (length(est) < 2) return(c(mean = NA, bias = NA, emp_SD = NA, mean_SE = NA, SE_over_SD = NA, RMSE = NA, cover = NA, cover_mcse = NA, R_eff = length(est)))
  cv <- mean(abs(est - target) <= qnorm(.975) * se)
  c(mean = mean(est), bias = mean(est) - target, emp_SD = sd(est), mean_SE = mean(se), SE_over_SD = mean(se) / sd(est),
    RMSE = sqrt(mean((est - target)^2)), cover = cv, cover_mcse = sqrt(cv * (1 - cv) / length(est)), R_eff = length(est))
}
EST_LABELS <- c(obs = "Observed d (Hedges)", wlsmv = "WLSMV latent d", mml = "MML 2PL latent d",
                conv = "GMM latent d (conventional SE)", rob = "GMM latent d (robust SE)",
                excl = "GMM after exclusion (C test)", part = "GMM partial invariance (D test)",
                sx = "Split-sample exclusion", sp = "Split-sample partial invariance")
summarize_simulation <- function() {
  files <- list.files(cfg$dir_sim, "^cell_\\d+\\.rds$", full.names = TRUE); rows <- list()
  for (f in files) {
    o <- readRDS(f); res <- o$results; cell <- o$cell; d <- cell$d; ds <- o$pseudo_true[["d_star"]]
    dsw <- if ("d_star_wlsmv" %in% names(o$pseudo_true)) o$pseudo_true[["d_star_wlsmv"]] else NA
    dsm <- if ("d_star_mml"   %in% names(o$pseudo_true)) o$pseudo_true[["d_star_mml"]]   else NA
    fail_main <- o$n_failed_main / o$n_attempted
    ests <- list(obs = list(res$d_obs, res$se_obs, d), wlsmv = list(res$d_wlsmv, res$se_wlsmv, dsw),
                 mml = list(res$d_mml, res$se_mml, dsm),
                 conv = list(res$d_gmm, res$se_conv, ds), rob = list(res$d_gmm, res$se_rob, ds),
                 excl = list(res$d_excl, res$se_excl, d), part = list(res$d_part, res$se_part, d),
                 sx = list(res$d_sx, res$se_sx, d), sp = list(res$d_sp, res$se_sp, d))
    for (nm in names(ests)) {
      st <- row_stats(ests[[nm]][[1]], ests[[nm]][[2]], ests[[nm]][[3]])
      rows[[length(rows) + 1]] <- data.frame(cell, est_id = nm, estimator = EST_LABELS[[nm]], t(st), d_star = ds,
        d_star_wlsmv = dsw, d_star_mml = dsm,
        target = ests[[nm]][[3]], J_reject = mean(res$pJ < cfg$alpha), mean_J = mean(res$Jstat), J_df = res$df[1],
        exact_excl = mean(res$exact_excl, na.rm = TRUE), mean_n_excl = mean(res$n_excl, na.rm = TRUE),
        exact_part = mean(res$exact_part, na.rm = TRUE), mean_n_untied = mean(res$n_untied, na.rm = TRUE),
        exact_sx = mean(res$exact_sx, na.rm = TRUE), exact_sp = mean(res$exact_sp, na.rm = TRUE),
        fail_main = fail_main, fail_rob = mean(res$fail_rob), fail_wlsmv = mean(res$fail_wlsmv), fail_mml = mean(res$fail_mml),
        fail_excl = mean(res$fail_excl), fail_part = mean(res$fail_part), fail_split = mean(res$fail_split),
        pred_atten_d = sqrt(mean(res$kr20)) * d, row.names = NULL)
    }
  }
  do.call(rbind, rows)
}

## ======================= 6. Empirical study (RSE) ============================
download_rse <- function() {
  zipf <- file.path(cfg$dir_emp, "RSE.zip")
  if (!file.exists(zipf)) {
    ok <- tryCatch({ download.file(cfg$emp_url, zipf, mode = "wb", quiet = TRUE); TRUE }, error = function(e) FALSE)
    if (!ok) stop("Download failed. Place RSE.zip from ", cfg$emp_url, " into ", cfg$dir_emp)
  }
  unzip(zipf, exdir = file.path(cfg$dir_emp, "raw"), overwrite = TRUE)
  csv <- list.files(file.path(cfg$dir_emp, "raw"), "\\.csv$", recursive = TRUE, full.names = TRUE)[1]
  raw <- read.delim(csv, stringsAsFactors = FALSE)
  if (!all(paste0("Q", 1:10) %in% names(raw))) raw <- read.csv(csv, stringsAsFactors = FALSE)
  raw
}
prepare_rse <- function(raw) {
  items <- paste0("Q", 1:10); stopifnot(all(items %in% names(raw)), "gender" %in% names(raw))
  d <- raw[raw$gender %in% c(1, 2), c(items, "gender", "age")]
  d <- d[apply(d[, items], 1, function(r) all(r %in% 1:4)), ]
  d <- d[!is.na(d$age) & d$age >= 13 & d$age <= 90, ]
  d[, paste0("Q", cfg$emp_neg_items)] <- 5 - d[, paste0("Q", cfg$emp_neg_items)]   # high = high self-esteem
  X <- (d[, items] >= 3) * 1L                      # dichotomize at the scale midpoint (agree vs disagree)
  g <- as.integer(d$gender == 1)                   # 0 = female (reference), 1 = male
  list(data = d, X = X, g = g, ordinal = d[, items])
}
## Moment residuals of the full-invariance GMM fit: observed minus model-implied
## marginal and pairwise endorsement proportions in each group. Standardized by the
## binomial standard error of the observed proportion. Large pairwise residuals that
## cluster among a subset of items indicate residual dependence (multidimensionality)
## rather than group-specific item parameters.
moment_residuals <- function(fit, S) {
  out <- list()
  for (gg in 0:1) {
    obs <- S$Ibar[[gg + 1]]; p <- probs(fit$phi, gg, S); ng <- round(S$n * S$w[gg + 1])
    z <- (obs - p) / sqrt(p * (1 - p) / ng)
    out[[gg + 1]] <- data.frame(group = c("female", "male")[gg + 1],
                                item_i = c(1:S$J, S$pj), item_j = c(rep(NA, S$J), S$pl),
                                type = c(rep("marginal", S$J), rep("pairwise", length(S$pj))),
                                observed = obs, fitted = p, residual = obs - p, z = z)
  }
  do.call(rbind, out)
}
run_empirical <- function() {
  f_out <- file.path(cfg$dir_emp, "empirical_results.rds")
  if (file.exists(f_out)) { message("empirical: cached"); return(invisible(readRDS(f_out))) }
  message("empirical: downloading and preparing RSE data")
  prep <- prepare_rse(download_rse()); X <- prep$X; g <- prep$g
  saveRDS(prep, file.path(cfg$dir_emp, "rse_clean.rds"))
  desc <- data.frame(group = c("female", "male"), n = c(sum(g == 0), sum(g == 1)),
                     rbind(colMeans(X[g == 0, ]), colMeans(X[g == 1, ])))
  names(desc)[-(1:2)] <- paste0("endorse_Q", 1:ncol(X))
  write.csv(desc, file.path(cfg$dir_emp, "descriptives.csv"), row.names = FALSE)
  S <- mom_from_data(X, g)
  message(sprintf("empirical: n = %d (female %d, male %d); fitting", S$n, sum(g == 0), sum(g == 1)))
  an <- analyze_dataset(S, split_seed = cfg$seed + 77)
  ct <- c_tests(an$fit, S, 1:S$J); dt <- d_tests(an$fit, S, 1:S$J)
  resid <- moment_residuals(an$fit, S)
  write.csv(resid, file.path(cfg$dir_emp, "moment_residuals.csv"), row.names = FALSE)
  item_tab <- data.frame(item = 1:S$J, endorse_f = colMeans(X[g == 0, ]), endorse_m = colMeans(X[g == 1, ]),
                         tau = an$fit$phi[S$ix$tau[[1]]], lambda = an$fit$phi[S$ix$lam[[1]]],
                         C = ct["stat", ], C_df = ct["df", ], C_p_holm = ct["p_holm", ],
                         D = dt["stat", ], D_p_holm = dt["p_holm", ])
  ## real-data Monte Carlo: subsamples of emp_sub_n per group, compared with the full-sample value
  set.seed(cfg$seed + 55)
  sub <- do.call(rbind, lapply(seq_len(cfg$emp_subsamples), function(b) {
    idx <- c(sample(which(g == 0), cfg$emp_sub_n), sample(which(g == 1), cfg$emp_sub_n))
    Sb <- mom_from_data(X[idx, ], g[idx])
    tryCatch({ fb <- gmm_fit(Sb, make_spec(Sb)); eb <- effect_size(fb$phi, fb$V, Sb)
               Vb <- robust_V(fb, Sb); er <- effect_size(fb$phi, Vb, Sb); ob <- observed_d(Sb)
               data.frame(b = b, d_obs = ob[["d"]], d_gmm = eb[["d"]], se_conv = eb[["se"]], se_rob = er[["se"]], pJ = fb$pJ) },
             error = function(e) NULL) }))
  if (length(an$errors)) message("empirical: stage errors: ", paste(names(an$errors), an$errors, sep = ": ", collapse = " | "))
  out <- list(S_summary = list(n = S$n, n_f = sum(g == 0), n_m = sum(g == 1), J = S$J), row = an$row,
              phi = an$fit$phi, J = an$fit$J, df = an$fit$df, pJ = an$fit$pJ, item_tab = item_tab, errors = an$errors,
              excl_set = if (is.null(an$ex)) NA else an$ex$set, part_set = if (is.null(an$pa)) NA else an$pa$set,
              split_excl_set = if (is.null(an$sp)) NA else an$sp$ex$set, split_part_set = if (is.null(an$sp)) NA else an$sp$pa$set,
              residuals = resid, n_sub_ok = if (is.null(sub)) 0 else nrow(sub), subsamples = sub)
  saveRDS(out, f_out); write.csv(item_tab, file.path(cfg$dir_emp, "item_tests.csv"), row.names = FALSE)
  write.csv(an$row, file.path(cfg$dir_emp, "effect_sizes.csv"), row.names = FALSE)
  if (!is.null(sub)) write.csv(sub, file.path(cfg$dir_emp, "subsample_estimates.csv"), row.names = FALSE)
  invisible(out)
}

## ================================ 7. Tables ==================================
df_to_md <- function(df, digits = 3) {
  f <- function(x) if (is.numeric(x)) formatC(x, digits = digits, format = "f") else as.character(x)
  m <- as.data.frame(lapply(df, f), stringsAsFactors = FALSE)
  paste(c(paste("|", paste(names(df), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
          apply(m, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))), collapse = "\n")
}
save_table <- function(df, name, digits = 3) {
  write.csv(df, file.path(cfg$dir_tab, paste0(name, ".csv")), row.names = FALSE)
  writeLines(df_to_md(df, digits), file.path(cfg$dir_tab, paste0(name, ".md")))
}
make_tables <- function(summ, emp) {
  ## Table 1: performance under invariance
  t1 <- summ[summ$dif == "none" & summ$est_id %in% c("obs", "wlsmv", "mml", "conv", "rob"), ]
  save_table(t1[order(t1$J, t1$n, t1$d, t1$est_id), c("n", "J", "d", "estimator", "mean", "bias", "emp_SD", "mean_SE", "SE_over_SD", "cover", "cover_mcse", "R_eff")],
             "table1_invariance_performance")
  ## Table 2: full-model estimators under noninvariance, evaluated against d*
  t2 <- summ[summ$dif != "none" & summ$est_id %in% c("obs", "wlsmv", "mml", "conv", "rob"), ]
  save_table(t2[order(t2$dif, t2$J, t2$n, t2$d, t2$est_id), c("dif", "n", "J", "d", "estimator", "target", "mean", "bias", "SE_over_SD", "cover", "cover_mcse", "R_eff")],
             "table2_noninvariance_pseudotrue")
  ## Table 3: detection and localisation
  t3 <- summ[summ$est_id == "excl", c("dif", "n", "J", "d", "J_reject", "mean_J", "J_df", "exact_excl", "mean_n_excl", "exact_part", "mean_n_untied", "exact_sx", "exact_sp")]
  save_table(t3[order(t3$dif, t3$J, t3$n, t3$d), ], "table3_detection_localisation")
  ## Table 4: step-2 estimators against the true d
  t4 <- summ[summ$dif != "none" & summ$est_id %in% c("excl", "part", "sx", "sp"), ]
  save_table(t4[order(t4$dif, t4$J, t4$n, t4$d, t4$est_id), c("dif", "n", "J", "d", "estimator", "mean", "bias", "SE_over_SD", "cover", "cover_mcse", "R_eff")],
             "table4_step2_recovery")
  ## Table 5: failure rates by stage
  t5 <- unique(summ[, c("n", "J", "d", "dif", "fail_main", "fail_rob", "fail_wlsmv", "fail_mml", "fail_excl", "fail_part", "fail_split")])
  save_table(t5[order(t5$J, t5$n, t5$d, t5$dif), ], "table5_failure_rates")
  ## Table 6: predicted versus observed attenuation
  t6 <- summ[summ$est_id == "obs" & summ$dif == "none" & summ$d > 0, c("n", "J", "d", "mean", "pred_atten_d", "bias")]
  names(t6)[4] <- "observed_d_mean"; save_table(t6, "table6_attenuation_check")
  ## Tables 7-8: empirical
  if (!is.null(emp)) {
    r <- emp$row; z <- qnorm(.975)
    t7 <- data.frame(estimator = EST_LABELS[c("obs", "wlsmv", "mml", "conv", "rob", "excl", "part", "sx", "sp")],
                     estimate = c(r$d_obs, r$d_wlsmv, r$d_mml, r$d_gmm, r$d_gmm, r$d_excl, r$d_part, r$d_sx, r$d_sp),
                     se = c(r$se_obs, r$se_wlsmv, r$se_mml, r$se_conv, r$se_rob, r$se_excl, r$se_part, r$se_sx, r$se_sp))
    t7$ci_low <- t7$estimate - z * t7$se; t7$ci_high <- t7$estimate + z * t7$se
    t7$note <- c(sprintf("KR-20 = %.3f", r$kr20), "full invariance", "full invariance, logistic link",
                 sprintf("J = %.1f, df = %d, p = %.4f", r$Jstat, r$df, r$pJ), "Hansen-Lee variance",
                 sprintf("excluded: %s", paste(emp$excl_set, collapse = ",")), sprintf("freed: %s", paste(emp$part_set, collapse = ",")),
                 sprintf("half A excluded: %s", paste(emp$split_excl_set, collapse = ",")), sprintf("half A freed: %s", paste(emp$split_part_set, collapse = ",")))
    save_table(t7, "table7_empirical_effect_sizes")
    save_table(emp$item_tab, "table8_empirical_item_tests")
    if (!is.null(emp$residuals)) {
      rp <- emp$residuals[emp$residuals$type == "pairwise", ]
      rp$neg_pair <- rp$item_i %in% cfg$emp_neg_items & rp$item_j %in% cfg$emp_neg_items
      rp$pos_pair <- !(rp$item_i %in% cfg$emp_neg_items) & !(rp$item_j %in% cfg$emp_neg_items)
      rp$pair_type <- ifelse(rp$neg_pair, "both negatively worded", ifelse(rp$pos_pair, "both positively worded", "mixed"))
      t9 <- aggregate(cbind(mean_residual = residual, mean_abs_z = abs(z)) ~ group + pair_type, data = rp, FUN = mean)
      t9$n_pairs <- aggregate(residual ~ group + pair_type, data = rp, FUN = length)$residual
      save_table(t9[order(t9$group, t9$pair_type), ], "table9_empirical_pairwise_residuals", digits = 4)
    }
  }
  invisible(NULL)
}

## ================================ 8. Figures =================================
save_fig <- function(p, name, w = 9, h = 6) {
  f_png <- file.path(cfg$dir_fig, paste0(name, ".png"))
  tryCatch({ ggsave(f_png, p, width = w, height = h, dpi = 300); message("saved ", f_png); TRUE },
           error = function(e) { message(name, ": PNG write failed: ", conditionMessage(e)); FALSE })
}
make_figures <- function(summ, emp) {
  if (!HAS_GG) { message("ggplot2 not available; figures skipped"); return(invisible(NULL)) }
  summ$Jlab <- factor(paste0("J = ", summ$J), levels = paste0("J = ", sort(unique(summ$J))))
  summ$dlab <- factor(paste0("d = ", summ$d)); summ$nlab <- factor(summ$n)
  summ$estimator <- factor(summ$estimator, levels = EST_LABELS)
  ## Figure 1: coverage of each estimator's own target; rows = condition, columns = J, line type = d
  f1d <- summ[summ$est_id %in% c("obs", "wlsmv", "mml", "conv", "rob") & !is.na(summ$cover), ]
  p1 <- ggplot(f1d, aes(x = nlab, y = cover, colour = estimator, group = interaction(estimator, d), linetype = dlab)) +
    geom_hline(yintercept = .95, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(ymin = cover - 1.96 * cover_mcse, ymax = cover + 1.96 * cover_mcse), width = .1, alpha = .5, linetype = 1) +
    geom_line() + geom_point(size = 1.8) + facet_grid(dif ~ Jlab) +
    labs(x = "n per group", y = "Coverage of the estimator's own target (true d, or its pseudo-true value under misspecification)",
         colour = NULL, linetype = NULL, title = "Coverage of nominal 95% intervals") +
    theme_bw() + theme(legend.position = "bottom", legend.box = "vertical")
  save_fig(p1, "fig1_coverage", 11, 9)
  ## Figure 2: SE/SD ratio, conventional vs robust (slope chart across n)
  f2d <- summ[summ$est_id %in% c("conv", "rob") & !is.na(summ$SE_over_SD), ]
  p2 <- ggplot(f2d, aes(x = nlab, y = SE_over_SD, colour = estimator, group = interaction(estimator, d))) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey40") + geom_line(aes(linetype = dlab)) + geom_point(size = 2) +
    facet_grid(dif ~ Jlab) + labs(x = "n per group", y = "Mean SE / empirical SD", colour = NULL, linetype = NULL,
    title = "Variance estimation: conventional versus misspecification-robust") + theme_bw() + theme(legend.position = "bottom")
  save_fig(p2, "fig2_se_ratio", 11, 7)
  ## Figure 3: bias against target with Monte Carlo error bars (dot plot)
  f3d <- summ[summ$est_id %in% c("obs", "wlsmv", "mml", "conv", "excl", "part", "sx", "sp") & !is.na(summ$bias), ]
  f3d$bias_mcse <- f3d$emp_SD / sqrt(f3d$R_eff)
  p3 <- ggplot(f3d, aes(x = bias, y = estimator, colour = nlab)) + geom_vline(xintercept = 0, colour = "grey40") +
    geom_errorbarh(aes(xmin = bias - 1.96 * bias_mcse, xmax = bias + 1.96 * bias_mcse), height = .25) + geom_point(size = 2) +
    facet_grid(dif + dlab ~ Jlab, scales = "free_x") + labs(x = "Bias relative to target", y = NULL, colour = "n per group",
    title = "Bias of effect-size estimators (full-model GMM/WLSMV measured against d*)") + theme_bw() + theme(legend.position = "bottom")
  save_fig(p3, "fig3_bias", 12, 9)
  ## Figure 4: localisation heat map, exact recovery of the noninvariant item set, four procedures
  f4 <- summ[summ$est_id == "excl" & summ$dif %in% c("uniform", "nonuniform"), c("n", "J", "d", "dif", "exact_excl", "exact_part", "exact_sx", "exact_sp")]
  f4l <- rbind(data.frame(f4[, 1:4], procedure = "Exclusion (C)", exact = f4$exact_excl),
               data.frame(f4[, 1:4], procedure = "Partial inv. (D)", exact = f4$exact_part),
               data.frame(f4[, 1:4], procedure = "Split excl. (C)", exact = f4$exact_sx),
               data.frame(f4[, 1:4], procedure = "Split partial (D)", exact = f4$exact_sp))
  f4l$procedure <- factor(f4l$procedure, levels = c("Exclusion (C)", "Partial inv. (D)", "Split excl. (C)", "Split partial (D)"))
  f4l$dlab <- paste0("d = ", f4l$d)
  p4 <- ggplot(f4l, aes(x = factor(n), y = factor(J), fill = exact)) + geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", exact)), size = 2.8) + facet_grid(dif + dlab ~ procedure) +
    scale_fill_gradient(low = "#fee8c8", high = "#b30000", limits = c(0, 1)) +
    labs(x = "n per group", y = "J", fill = "Exact set\nrecovered", title = "Localisation of noninvariant items") + theme_bw()
  save_fig(p4, "fig4_localisation", 10, 7.5)
  ## Figure 5: J test rejection and failure rates
  f5 <- unique(summ[, c("n", "J", "d", "dif", "J_reject", "fail_main")])
  p5 <- ggplot(f5, aes(x = factor(n), y = J_reject, fill = dif)) + geom_col(position = position_dodge(.8), width = .7) +
    geom_hline(yintercept = .05, linetype = 2) + facet_grid(d ~ J, labeller = label_both) +
    labs(x = "n per group", y = "J rejection rate", fill = "Condition", title = "Overidentification test: size (none) and power (DIF)") +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(p5, "fig5_J_test", 10, 6)
  ## Figure 6: observed vs latent estimates against reliability (attenuation)
  f6 <- summ[summ$est_id %in% c("obs", "conv") & summ$dif == "none" & summ$d > 0, ]
  p6 <- ggplot(f6, aes(x = Jlab, y = mean, colour = estimator, shape = nlab)) + geom_hline(yintercept = unique(f6$d), linetype = 2) +
    geom_point(size = 3, position = position_dodge(.5)) +
    geom_point(aes(y = pred_atten_d), data = f6[f6$est_id == "obs", ], shape = 4, size = 3, colour = "black", position = position_dodge(.5)) +
    labs(x = NULL, y = "Mean estimate", colour = NULL, shape = "n per group",
         title = "Attenuation of the observed effect size (x = sqrt(KR-20) x d predicted)") + theme_bw() + theme(legend.position = "bottom")
  save_fig(p6, "fig6_attenuation", 8, 5)
  ## Empirical figures
  if (!is.null(emp)) {
    r <- emp$row; z <- qnorm(.975)
    lv <- c("Observed d", "WLSMV", "MML 2PL", "GMM conventional", "GMM robust", "GMM exclusion", "GMM partial", "Split exclusion", "Split partial")
    f7 <- data.frame(estimator = factor(lv, levels = rev(lv)),
                     est = c(r$d_obs, r$d_wlsmv, r$d_mml, r$d_gmm, r$d_gmm, r$d_excl, r$d_part, r$d_sx, r$d_sp),
                     se = c(r$se_obs, r$se_wlsmv, r$se_mml, r$se_conv, r$se_rob, r$se_excl, r$se_part, r$se_sx, r$se_sp))
    p7 <- ggplot(f7, aes(x = est, y = estimator)) + geom_vline(xintercept = 0, colour = "grey50") +
      geom_errorbarh(aes(xmin = est - z * se, xmax = est + z * se), height = .2) + geom_point(size = 3) +
      labs(x = "Standardized gender difference in self-esteem (male - female)", y = NULL,
           title = sprintf("Rosenberg Self-Esteem Scale (n = %d): effect sizes with 95%% intervals", emp$S_summary$n)) + theme_bw()
    save_fig(p7, "fig7_empirical_forest", 8, 4.5)
    it <- emp$item_tab
    f8 <- rbind(data.frame(item = it$item, stat = it$C, p = it$C_p_holm, test = "C (subset moment)"),
                data.frame(item = it$item, stat = it$D, p = it$D_p_holm, test = "D (restriction)"))
    f8$flag <- factor(ifelse(f8$p < cfg$alpha, "flagged (Holm p < .05)", "not flagged"), levels = c("not flagged", "flagged (Holm p < .05)"))
    p8 <- ggplot(f8, aes(x = factor(item), y = stat, colour = flag)) + geom_segment(aes(xend = factor(item), y = 0, yend = stat)) +
      geom_point(size = 3) + facet_wrap(~ test, scales = "free_y") +
      scale_colour_manual(values = c("not flagged" = "grey50", "flagged (Holm p < .05)" = "#b30000"), drop = FALSE) +
      labs(x = "Item", y = "Test statistic", colour = NULL, title = "Item-level invariance tests, RSE gender comparison") +
      theme_bw() + theme(legend.position = "bottom")
    save_fig(p8, "fig8_empirical_item_tests", 9, 4.5)
    if (!is.null(emp$residuals)) {
      rp <- emp$residuals[emp$residuals$type == "pairwise", ]
      rp2 <- rbind(rp, transform(rp, item_i = item_j, item_j = rp$item_i))   # symmetric matrix
      p10 <- ggplot(rp2, aes(x = factor(item_i), y = factor(item_j), fill = residual)) + geom_tile(colour = "white") +
        geom_text(aes(label = sprintf("%.3f", residual)), size = 2.4) + facet_wrap(~ group) +
        scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
        labs(x = "Item", y = "Item", fill = "Observed minus\nfitted joint\nendorsement",
             title = "Pairwise moment residuals of the full-invariance model (items 3, 5, 8, 9, 10 are negatively worded)") +
        theme_bw() + theme(legend.position = "right")
      save_fig(p10, "fig10_empirical_residuals", 11, 5.5)
    }
    if (!is.null(emp$subsamples) && nrow(emp$subsamples) > 5) {
      sb <- emp$subsamples
      f9 <- rbind(data.frame(est = sb$d_obs, which = "Observed d"), data.frame(est = sb$d_gmm, which = "GMM latent d"))
      ref <- data.frame(which = c("Observed d", "GMM latent d"), full = c(r$d_obs, r$d_gmm))
      p9 <- ggplot(f9, aes(x = est, fill = which)) + geom_density(alpha = .4) +
        geom_vline(data = ref, aes(xintercept = full, colour = which), linetype = 2) +
        labs(x = sprintf("Estimate from subsamples of %d per group", cfg$emp_sub_n), y = "Density", fill = NULL, colour = NULL,
             title = "Real-data Monte Carlo: sampling distributions versus full-sample values (dashed)") + theme_bw() + theme(legend.position = "bottom")
      save_fig(p9, "fig9_empirical_subsampling", 8, 5)
    }
  }
  invisible(NULL)
}


## ============ Utility: recompute comparators on saved cells (no GMM rerun) ===========
refit_comparators <- function() {
  files <- list.files(cfg$dir_sim, "^cell_\\d+\\.rds$", full.names = TRUE)
  for (f in files) {
    o <- readRDS(f); cell <- o$cell; res <- o$results
    i <- as.integer(sub("cell_(\\d+)\\.rds", "\\1", basename(f)))
    one <- function(r) {
      S <- gen_moments(cell$n, cell$J, cell$d, cell$dif, o$item_params, seed = cfg$seed + 100000L * i + r)
      wl <- tryCatch(wlsmv_d(S), error = function(e) c(d = NA, se = NA))
      ml <- tryCatch(mml_d(S),   error = function(e) c(d = NA, se = NA))
      data.frame(d_wlsmv = wl[["d"]], se_wlsmv = wl[["se"]], d_mml = ml[["d"]], se_mml = ml[["se"]])
    }
    new <- do.call(rbind, if (cfg$cores > 1 && .Platform$OS.type == "unix")
                            parallel::mclapply(res$rep, one, mc.cores = cfg$cores) else lapply(res$rep, one))
    res$d_wlsmv <- new$d_wlsmv; res$se_wlsmv <- new$se_wlsmv; res$fail_wlsmv <- is.na(new$d_wlsmv)
    res$d_mml <- new$d_mml; res$se_mml <- new$se_mml; res$fail_mml <- is.na(new$d_mml)
    o$results <- res; saveRDS(o, f)
    message(sprintf("cell %d: WLSMV ok %.3f | MML ok %.3f  [%s]", i, mean(!res$fail_wlsmv), mean(!res$fail_mml), format(Sys.time(), "%H:%M:%S")))
  }
  invisible(NULL)
}

## ======= Utility: recompute pseudo-true values (GMM, WLSMV, MML) on saved cells =======
refit_pseudo_true <- function() {
  files <- list.files(cfg$dir_sim, "^cell_\\d+\\.rds$", full.names = TRUE)
  for (f in files) {
    o <- readRDS(f); cell <- o$cell
    pt <- pseudo_true(cell$J, cell$d, cell$dif, o$item_params)
    o$pseudo_true <- pt
    o$results$d_star <- pt[["d_star"]]; o$results$d_star_se <- pt[["mc_se"]]
    o$results$d_star_wlsmv <- pt[["d_star_wlsmv"]]; o$results$d_star_mml <- pt[["d_star_mml"]]
    saveRDS(o, f)
    message(sprintf("%s: d*(GMM) = %.4f, d*(WLSMV) = %.4f, d*(MML) = %.4f", basename(f), pt[["d_star"]], pt[["d_star_wlsmv"]], pt[["d_star_mml"]]))
  }
  invisible(NULL)
}

## Consistency check: on one correctly specified dataset the three latent estimators
## estimate the same quantity and should agree to within sampling noise. Aborts the
## long run if a comparator is misconfigured (the symptom seen earlier was a WLSMV
## value of 0.16 against a GMM value of 0.34 on the same data).
consistency_check <- function() {
  ip <- item_params(5); S <- gen_moments(2000, 5, 0.4, "none", ip, seed = 123)
  fit <- gmm_fit(S, make_spec(S)); g <- effect_size(fit$phi, fit$V, S)
  w <- if (cfg$use_lavaan) wlsmv_d(S) else c(d = NA, se = NA)
  m <- if (cfg$use_mirt)   mml_d(S)   else c(d = NA, se = NA)
  message(sprintf("consistency check (n = 2000/group, J = 5, d = 0.4): GMM %.3f (%.3f) | WLSMV %.3f (%.3f) | MML %.3f (%.3f)",
                  g[["d"]], g[["se"]], w[["d"]], w[["se"]], m[["d"]], m[["se"]]))
  tol <- 4 * g[["se"]]
  if (cfg$use_lavaan && (is.na(w[["d"]]) || abs(w[["d"]] - g[["d"]]) > tol)) stop("WLSMV comparator disagrees with GMM; check wlsmv_d()")
  if (cfg$use_mirt   && (is.na(m[["d"]]) || abs(m[["d"]] - g[["d"]]) > tol)) stop("MML comparator disagrees with GMM; check mml_d()")
  invisible(TRUE)
}

## ================================= 9. Main ===================================
if (!isTRUE(getOption("es.skip_main"))) {
  t0 <- Sys.time(); summ <- NULL; emp <- NULL
  consistency_check()
  if (cfg$run_simulation) {
    message("Simulation: ", nrow(design), " cells x ", cfg$n_rep, " reps; lavaan = ", cfg$use_lavaan, "; mirt = ", cfg$use_mirt, " -> ", normalizePath(cfg$dir_sim))
    for (i in seq_len(nrow(design))) run_cell(i)
    summ <- summarize_simulation()
    write.csv(summ, file.path(cfg$dir_sim, "summary.csv"), row.names = FALSE); saveRDS(summ, file.path(cfg$dir_sim, "summary.rds"))
  } else if (file.exists(file.path(cfg$dir_sim, "summary.rds"))) summ <- readRDS(file.path(cfg$dir_sim, "summary.rds"))
  if (cfg$run_empirical) emp <- tryCatch(run_empirical(), error = function(e) { message("empirical study failed: ", conditionMessage(e)); NULL })
  if (cfg$make_outputs) {
    if (!is.null(summ)) { make_tables(summ, emp); make_figures(summ, emp) } else message("no simulation summary; outputs skipped")
  }
  if (!is.null(emp)) {
    r <- emp$row
    cat(sprintf("\nEmpirical (RSE, n = %d): observed d = %.3f (%.3f); GMM d = %.3f (conv %.3f, robust %.3f); J = %.1f, df = %d, p = %.4f\n",
                emp$S_summary$n, r$d_obs, r$se_obs, r$d_gmm, r$se_conv, r$se_rob, r$Jstat, r$df, r$pJ))
    cat("Items flagged: exclusion =", paste(emp$excl_set, collapse = ","), "| partial invariance =", paste(emp$part_set, collapse = ","), "\n")
  }
  message("Done in ", format(round(Sys.time() - t0, 1)), ". Tables -> ", cfg$dir_tab, "; figures -> ", cfg$dir_fig)
}
