
# Funzione Metropolis-Hasting
mh_sample <- function(current, log_target, proposal_sd, max_iter = 10) {
  for (iter in 1:max_iter) {
    candidate <- rnorm(1, current, proposal_sd)
    log_ratio <- log_target(candidate) - log_target(current)
    if (length(log_ratio) > 1) log_ratio <- log_ratio[1]
    if (log(runif(1)) < log_ratio) current <- candidate
  }
  return(current)
}

# Funzione principale MCMC
run_mcmc_base <- function(data, n_iter, burn_in, thin, seed) {
  
  set.seed(seed)
  
  # Preparazione dei dati
  teams <- unique(c(data$Team.1, data$Team.2))
  T <- length(teams)
  
  team_ids <- 1:T
  names(team_ids) <- teams
  
  results <- strsplit(data$FT, "-")
  data$HomeGoals <- as.numeric(sapply(results, function(x) x[1]))
  data$AwayGoals <- as.numeric(sapply(results, function(x) x[2]))
  
  G <- nrow(data)
  
  home_team <- sapply(data$Team.1, function(x) team_ids[x])
  away_team <- sapply(data$Team.2, function(x) team_ids[x])
  
  # Inizializzazione parametri
  att <- rnorm(T, 0, 0.1)
  def <- rnorm(T, 0, 0.1)
  home <- 0.2
  
  # Iperparametri
  mu_att <- 0.1
  sigma2_att <- 0.5
  mu_def <- -0.1
  sigma2_def <- 0.5
  
  # Centratura iniziale
  att <- att - mean(att)
  def <- def - mean(def)
  
  # Strutture per salvare i risultati
  att_all <- matrix(0, nrow = T, ncol = n_iter)
  def_all <- matrix(0, nrow = T, ncol = n_iter)
  home_all <- numeric(n_iter)
  
  mu_att_all <- numeric(n_iter)
  sigma2_att_all <- numeric(n_iter)
  mu_def_all <- numeric(n_iter)
  sigma2_def_all <- numeric(n_iter)
  
  
  for (iter in 1:n_iter) {
    
    # 1. Campiono att_{t}
    for (t in 1:T) {
      home_matches <- which(home_team == t)
      away_matches <- which(away_team == t)
      
      log_target_att <- function(att_val) {
        log_prior <- -0.5 * (att_val - mu_att)^2 / sigma2_att
        
        log_lik <- 0
        # Partite in casa (t = Home)
        if(length(home_matches) > 0) {
          lambda <- exp(home + att_val + def[away_team[home_matches]])
          log_lik <- log_lik + sum(dpois(data$HomeGoals[home_matches], lambda, log = TRUE))
        }
        # Partite in trasferta (t = Away)
        if(length(away_matches) > 0) {
          lambda <- exp(att_val + def[home_team[away_matches]])
          log_lik <- log_lik + sum(dpois(data$AwayGoals[away_matches], lambda, log = TRUE))
        }
        return(log_prior + log_lik)
      }
      
      att[t] <- mh_sample(att[t], log_target_att, proposal_sd = 0.1, max_iter = 1)
    }
    
    # 2. Campiono def_{t} 
    for (t in 1:T) {
      home_matches <- which(home_team == t)
      away_matches <- which(away_team == t)
      
      log_target_def <- function(def_val) {
        log_prior <- -0.5 * (def_val - mu_def)^2 / sigma2_def
        
        log_lik <- 0
        # Partite in casa (t = Home): qui t difende contro away_team (che attacca)
        # Nota: nel modello Poisson standard lambda_away = exp(att_away + def_home)
        if(length(home_matches) > 0) {
          lambda <- exp(att[away_team[home_matches]] + def_val)
          log_lik <- log_lik + sum(dpois(data$AwayGoals[home_matches], lambda, log = TRUE))
        }
        # Partite in trasferta (t = Away): qui t difende contro home_team (che attacca)
        if(length(away_matches) > 0) {
          lambda <- exp(home + att[home_team[away_matches]] + def_val)
          log_lik <- log_lik + sum(dpois(data$HomeGoals[away_matches], lambda, log = TRUE))
        }
        return(log_prior + log_lik)
      }
      
      def[t] <- mh_sample(def[t], log_target_def, proposal_sd = 0.1, max_iter = 1)
    }
    
    # 3. Campiono home
    log_target_home <- function(home_val) {
      log_prior <- -0.5 * home_val^2 / 100
      
      lambda <- exp(home_val + att[home_team] + def[away_team])
      log_lik <- sum(dpois(data$HomeGoals, lambda, log = TRUE))
      
      return(log_prior + log_lik)
    }
    home <- mh_sample(home, log_target_home, proposal_sd = 0.05, max_iter = 1)
    
    
    # 4. Campiono Iperparametri (Gibbs steps)
    
    # mu_att, sigma2_att
    prec_prior <- 1/100
    prec_lik <- T / sigma2_att
    mu_att_post <- (sum(att) / sigma2_att) / (prec_prior + prec_lik)
    sigma2_att_post <- 1 / (prec_prior + prec_lik)
    mu_att <- rnorm(1, mu_att_post, sqrt(sigma2_att_post))
    
    shape_post <- 0.1 + T/2
    rate_post <- 0.1 + 0.5 * sum((att - mu_att)^2)
    sigma2_att <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    
    # mu_def, sigma2_def
    prec_lik <- T / sigma2_def
    mu_def_post <- (sum(def) / sigma2_def) / (prec_prior + prec_lik)
    sigma2_def_post <- 1 / (prec_prior + prec_lik)
    mu_def <- rnorm(1, mu_def_post, sqrt(sigma2_def_post))
    
    rate_post <- 0.1 + 0.5 * sum((def - mu_def)^2)
    sigma2_def <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    
    # Vincoli di identificabilità (Centratura)
    att <- att - mean(att)
    def <- def - mean(def)
    
    # Salvataggio
    att_all[, iter] <- att
    def_all[, iter] <- def
    home_all[iter] <- home
    mu_att_all[iter] <- mu_att
    sigma2_att_all[iter] <- sigma2_att
    mu_def_all[iter] <- mu_def
    sigma2_def_all[iter] <- sigma2_def
    
    if(iter %% 5000 == 0) cat(sprintf("Iterazione %d/%d completata.\n", iter, n_iter))
  }
  
  # Burn-in e thinning
  indices <- seq(burn_in + 1, n_iter, by = thin)
  n_samples <- length(indices)
  
  results <- list(
    att = att_all[, indices],
    def = def_all[, indices],
    home = home_all[indices],
    mu_att = mu_att_all[indices],
    sigma2_att = sigma2_att_all[indices],
    mu_def = mu_def_all[indices],
    sigma2_def = sigma2_def_all[indices],
    teams = teams,
    n_iter = n_iter,
    burn_in = burn_in,
    thin = thin,
    n_samples = n_samples
  )
  
  return(results)
}


###################################################################################
# ESECUZIONE
data <- read.csv("./2015-16/eng.1.csv", stringsAsFactors = FALSE)
# Uso lo stesso seed del mixed model per confrontabilità (dove possibile)
samples <- run_mcmc_base(data, n_iter = 50000, burn_in = 10000, thin = 10, seed = 306256)

teams <- samples$teams
T <- length(teams)
n_samples <- samples$n_samples

# Calcolo medie a posteriori
att_mean <- rowMeans(samples$att)
def_mean <- rowMeans(samples$def)
home_mean <- mean(samples$home)

cat("Effetto casa medio:", home_mean, "\n")
cat('Media Attacco',mean(samples$mu_att))
cat('Media Difesa',mean(samples$mu_def))
cat('Sigma Quadro Difesa',mean(samples$sigma2_def))
cat('Sigma Quadro Attacco',mean(samples$sigma2_att))

# Classifica basata sui parametri
ranking <- data.frame(
  Team = teams,
  Attack = att_mean,
  Defense = def_mean,
  Overall_score = att_mean - def_mean # Attacco alto è bene, Difesa alta è male (gol subiti)
)

ranking <- ranking[order(-ranking$Overall_score), ]

cat("\nCLASSIFICA BASATA SUI PARAMETRI (BASE MODEL)\n")

print(ranking, row.names = FALSE)


# Grafici 

# 1. Barplot Parametri
par(mfrow = c(1, 2))
barplot(att_mean, names.arg = teams, las = 2, cex.names = 0.6,
        main = "Parametri di Attacco (Statico)",
        ylab = "Valore", col = "skyblue")

barplot(def_mean, names.arg = teams, las = 2, cex.names = 0.6,
        main = "Parametri di Difesa (Statico)",
        ylab = "Valore", col = "lightcoral")
par(mfrow = c(1, 1))


# 2. Effetto Casa
par(mfrow = c(1, 2))
hist(samples$home, breaks = 30, main = "Distribuzione Effetto Casa",
     xlab = "Valore", ylab = "Frequenza", col = "lightgreen", border = "white")

plot(samples$home, type = "l", main = "Trace plot effetto casa", 
     xlab = "Campione", ylab = "Valore", col = "blue")
par(mfrow = c(1, 1))


# 3. Trace plots e acf Iperparametri
par(mfrow = c(2, 4))
plot(samples$mu_att, type = "l", main = "Trace mu_att", xlab = "Iter", col="red")
plot(samples$sigma2_att, type = "l", main = "Trace sigma2_att", xlab = "Iter", col="orange")
plot(samples$mu_def, type = "l", main = "Trace mu_def", xlab = "Iter", col="green")
plot(samples$sigma2_def, type = "l", main = "Trace sigma2_def", xlab = "Iter", col="brown")

acf(samples$mu_att, main = "ACF mu_att")
acf(samples$sigma2_att, main = "ACF sigma2_att")
acf(samples$mu_def, main = "ACF mu_def")
acf(samples$sigma2_def, main = "ACF sigma2_def")
par(mfrow = c(1, 1))


# 4. Grafici per una squadra specifica
team_index <- which(teams == "Manchester City")

if (length(team_index) > 0) {
  par(mfrow = c(2, 2))
  
  plot(samples$att[team_index, ], type = "l", main = "Trace att (Man City)",
       xlab = "Campione", ylab = "Valore", col = "blue")
  
  acf(samples$att[team_index, ], main = "ACF att (Man City)")
  
  plot(samples$def[team_index, ], type = "l", main = "Trace def (Man City)",
       xlab = "Campione", ylab = "Valore", col = "red")
  
  acf(samples$def[team_index, ], main = "ACF def (Man City)")
  
  par(mfrow = c(1, 1))
  
  cat("\n=== STATISTICHE MANCHESTER CITY ===\n")
  cat("Attacco medio:", mean(samples$att[team_index, ]), "\n")
  cat("Difesa media:", mean(samples$def[team_index, ]), "\n")
}


# Intervallo di credibilità per l'effetto casa
home_ci <- quantile(samples$home, c(0.025, 0.975))
cat("\nIntervallo di credibilità 95% per l'effetto casa:\n")
print(home_ci)


###############################################################################
# Predizione 

predict_match_base <- function(samples, team_home, team_away, data) {
  
  teams <- samples$teams
  n_samples <- samples$n_samples
  
  # Trovo gli indici delle squadre
  team_home_idx <- which(teams == team_home)
  team_away_idx <- which(teams == team_away)
  
  # Calcolo i parametri per simulare le poisson per ogni campione
  lambda_home <- numeric(n_samples)
  lambda_away <- numeric(n_samples)
  
  for (m in 1:n_samples) {
    # Qui non c'è lo stato z, prendiamo direttamente il parametro  campionato
    att_home <- samples$att[team_home_idx, m]
    def_home <- samples$def[team_home_idx, m]
    
    att_away <- samples$att[team_away_idx, m]
    def_away <- samples$def[team_away_idx, m]
    
    home_effect <- samples$home[m]
    
    lambda_home[m] <- exp(home_effect + att_home + def_away)
    lambda_away[m] <- exp(att_away + def_home)
  }
  
  goals_home_samples <- rpois(n_samples, lambda_home)
  goals_away_samples <- rpois(n_samples, lambda_away)
  
  
  prob_home_win <- mean(goals_home_samples > goals_away_samples)
  prob_draw <- mean(goals_home_samples == goals_away_samples)
  prob_away_win <- mean(goals_home_samples < goals_away_samples)
  
  exp_goals_home <- mean(goals_home_samples)
  exp_goals_away <- mean(goals_away_samples)
  
  median_goals_home <- median(goals_home_samples)
  median_goals_away <- median(goals_away_samples)
  
  ci_home <- quantile(goals_home_samples, c(0.025, 0.975))
  ci_away <- quantile(goals_away_samples, c(0.025, 0.975))
  
  max_goals <- max(max(goals_home_samples), max(goals_away_samples))
  joint_pmf <- matrix(0, nrow = max_goals + 1, ncol = max_goals + 1)
  
  for (i in 0:max_goals) {
    for (j in 0:max_goals) {
      joint_pmf[i + 1, j + 1] <- mean(
        (goals_home_samples == i) & (goals_away_samples == j)
      )
    }
  }
  
  rownames(joint_pmf) <- 0:max_goals
  colnames(joint_pmf) <- 0:max_goals
  
  
  result <- list(
    team_home = team_home,
    team_away = team_away,
    goals_home_samples = goals_home_samples,
    goals_away_samples = goals_away_samples,
    exp_goals_home = exp_goals_home,
    exp_goals_away = exp_goals_away,
    median_goals_home = median_goals_home,
    median_goals_away = median_goals_away,
    ci_home_95 = ci_home,
    ci_away_95 = ci_away,
    prob_home_win = prob_home_win,
    prob_draw = prob_draw,
    prob_away_win = prob_away_win,
    joint_pmf = joint_pmf
  )
  
  return(result)
}


# Funzione di stampa 
print_prediction <- function(prediction) {
  
  cat("\n")
  cat("════════════════════════════════════════════════════════════════\n")
  cat(sprintf("%-60s\n", paste0(prediction$team_home, " vs ", prediction$team_away)))
  cat("════════════════════════════════════════════════════════════════\n")
  
  cat("\n--- VALORE ATTESO GOL ---\n")
  cat(sprintf("%s: %.2f gol\n", prediction$team_home, prediction$exp_goals_home))
  cat(sprintf("%s: %.2f gol\n", prediction$team_away, prediction$exp_goals_away))
  
  cat("\n--- MEDIANA E INTERVALLO CREDIBILE 95% ---\n")
  cat(sprintf("%s: mediana = %d, IC = [%.2f, %.2f]\n", 
              prediction$team_home,
              prediction$median_goals_home,
              prediction$ci_home_95[1],
              prediction$ci_home_95[2]))
  cat(sprintf("%s: mediana = %d, IC = [%.2f, %.2f]\n", 
              prediction$team_away,
              prediction$median_goals_away,
              prediction$ci_away_95[1],
              prediction$ci_away_95[2]))
  
  cat("\n--- PROBABILITÀ RISULTATI ---\n")
  cat(sprintf("Vittoria %s: %.2f%%\n", prediction$team_home, prediction$prob_home_win * 100))
  cat(sprintf("Pareggio:  %.2f%%\n", prediction$prob_draw * 100))
  cat(sprintf("Vittoria %s: %.2f%%\n", prediction$team_away, prediction$prob_away_win * 100))
  
  cat("\n--- TOP 3 RISULTATI PIÙ PROBABILI ---\n")
  
  home_goals <- as.numeric(rownames(prediction$joint_pmf))
  away_goals <- as.numeric(colnames(prediction$joint_pmf))
  
  results_list <- list()
  idx <- 1
  
  for (i in 1:length(home_goals)) {
    for (j in 1:length(away_goals)) {
      results_list[[idx]] <- data.frame(
        Home_Goals = home_goals[i],
        Away_Goals = away_goals[j],
        Probability = prediction$joint_pmf[i, j]
      )
      idx <- idx + 1
    }
  }
  
  results_df <- do.call(rbind, results_list)
  results_df <- results_df[order(results_df$Probability, decreasing = TRUE), ]
  top_results <- head(results_df, 3)
  
  print(top_results, row.names = FALSE)
  
  cat("\n")
}

# Esempio di predizione
pred <- predict_match_base(samples, team_home = "Leicester City", team_away = "Aston Villa", data = data)

print_prediction(pred)
