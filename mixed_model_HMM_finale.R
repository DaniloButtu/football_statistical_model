rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  x <- matrix(rgamma(n * k, shape = alpha, rate = 1), ncol = k, byrow = TRUE)
  row_sums <- rowSums(x)
  x / row_sums
}

log_lik_match <- function(home_goals, away_goals, home_att, away_att, home_def, away_def, home_eff) {
  lambda1 <- exp(home_eff + home_att + away_def)
  lambda2 <- exp(away_att + home_def)
  dpois(home_goals, lambda1, log = TRUE) + dpois(away_goals, lambda2, log = TRUE)
}

mh_sample <- function(current, log_target, proposal_sd, max_iter = 10) {
  for (iter in 1:max_iter) {
    candidate <- rnorm(1, current, proposal_sd)
    log_ratio <- log_target(candidate) - log_target(current)
    if (length(log_ratio) > 1) log_ratio <- log_ratio[1]
    if (log(runif(1)) < log_ratio) current <- candidate
  }
  return(current)
}

#Funzione principale MCMC
run_mcmc <- function(data, n_iter, burn_in, thin, K, seed) {
  
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
  
  
  #Indici temporali
  match_counter <- rep(1, T)
  home_match_idx <- rep(0, G)
  away_match_idx <- rep(0, G)
  
  for (g in 1:G) {
    ht <- home_team[g]
    at <- away_team[g]
    
    home_match_idx[g] <- match_counter[ht]
    away_match_idx[g] <- match_counter[at]
    
    match_counter[ht] <- match_counter[ht] + 1
    match_counter[at] <- match_counter[at] + 1
  }
  
  S_t <- match_counter - 1
  max_matches <- max(S_t)
  
 
  #Inizializzazione dei parametri 
  
  # Calcolo prestazioni iniziali per ogni squadra (prime 5 partite)
  team_performance <- numeric(T)
  for (t in 1:T) {
    team_name <- teams[t]
    
    # Trovo le prime 5 partite della squadra per iniziallizzare in modo sensato i parametri
    home_games <- which(data$Team.1 == team_name)[1:min(5, sum(data$Team.1 == team_name))]
    away_games <- which(data$Team.2 == team_name)[1:min(5, sum(data$Team.2 == team_name))]
    
    #gol fatti e subiti
    home_goals_for <- sum(data$HomeGoals[home_games], na.rm = TRUE)
    away_goals_for <- sum(data$AwayGoals[away_games], na.rm = TRUE)
    home_goals_against <- sum(data$AwayGoals[home_games], na.rm = TRUE)
    away_goals_against <- sum(data$HomeGoals[away_games], na.rm = TRUE)
    
    team_performance[t] <- (home_goals_for + away_goals_for) - (home_goals_against + away_goals_against)
  }
  
  # Inizializzazione degli z
  z <- matrix(NA, nrow = T, ncol = max_matches)
  
  for (t in 1:T) {
    if (team_performance[t] > median(team_performance)) {
      initial_state <- 1 #se la performance del team nelle prime partite è maggiore della mediana 
                         #della performance di tutto il campionato assegno stato 1=in forma
    } else {
      initial_state <- 2
    }
    
    for (s in 1:S_t[t]) {
      if (s == 1) {
        z[t, s] <- initial_state
      } else {
        #80% probabilità di rimanere nello stesso stato in perchè voglio modellare il fatto che le squadre rimangono per un tot di partite in uno stato
        if (runif(1) < 0.8) {
          z[t, s] <- z[t, s-1]
        } else {
          z[t, s] <- sample(setdiff(1:K, z[t, s-1]), 1) #con setdiff(1:K, z[t, s-1]) tolgo dall'insieme degli stati lo stato assunto alla partita precedente
        }
      }
    }
  }
  
  # Inizializzazione fattori attacco/difesa
  att <- matrix(0, nrow = T, ncol = K)
  def <- matrix(0, nrow = T, ncol = K)
  
  good_teams <- which(team_performance > median(team_performance))
  poor_teams <- which(team_performance <= median(team_performance))
  
  #Stato 1 (in forma), allora do dei valori un pochino piu' alti
  att[good_teams, 1] <- rnorm(length(good_teams), 0.2, 0.1) 
  def[good_teams, 1] <- rnorm(length(good_teams), -0.1, 0.1)
  
  #Stato 2 (fuori forma)
  att[poor_teams, 2] <- rnorm(length(poor_teams), -0.2, 0.1)
  def[poor_teams, 2] <- rnorm(length(poor_teams), 0.1, 0.1)
  
  #Vincoli di identificabilità
  for (k in 1:K) {
    att[, k] <- att[, k] - mean(att[, k])
    def[, k] <- def[, k] - mean(def[, k])
  }
  
  home <- rnorm(1, 0.2, 0.1)
  
  mu_att <- c(0.1, -0.1)
  sigma2_att <- c(0.5, 0.5)
  mu_def <- c(-0.1, 0.1)
  sigma2_def <- c(0.5, 0.5)
  
  pi <- rep(0.5, K)
  P <- matrix(0.5, nrow = K, ncol = K, byrow = TRUE)
  
  #Uso array per salvare tutte le iterazioni 
  att_all <- array(0, dim = c(T, K, n_iter))
  def_all <- array(0, dim = c(T, K, n_iter))
  home_all <- numeric(n_iter)
  z_all <- array(0, dim = c(T, max_matches, n_iter))
  mu_att_all <- matrix(0, nrow = n_iter, ncol = K)
  sigma2_att_all <- matrix(0, nrow = n_iter, ncol = K)
  mu_def_all <- matrix(0, nrow = n_iter, ncol = K)
  sigma2_def_all <- matrix(0, nrow = n_iter, ncol = K)
  pi_all <- matrix(0, nrow = n_iter, ncol = K)
  P_all <- array(0, dim = c(K, K, n_iter))
  
  for (iter in 1:n_iter) {
    
    # Campiono z_{t,s} | ...
    for (t in 1:T) {
      for (s in 1:S_t[t]) {
        match_g <- which((home_team == t & home_match_idx == s) | 
                           (away_team == t & away_match_idx == s))
        g <- match_g[1]
        
        log_probs <- numeric(K)
        
        for (k in 1:K) {
          if (s == 1) {
            if (S_t[t] > 1) {
              log_prior <- log(pi[k]) + log(P[k, z[t, 2]])
            } else {
              log_prior <- log(pi[k])
            }
          } else if (s < S_t[t]) {
            log_prior <- log(P[z[t, s-1], k]) + log(P[k, z[t, s+1]])
          } else {
            log_prior <- log(P[z[t, s-1], k])
          }
          
          if (home_team[g] == t) {
            home_att_val <- att[t, k]
            home_def_val <- def[t, k]
            away_att_val <- att[away_team[g], z[away_team[g], away_match_idx[g]]]
            away_def_val <- def[away_team[g], z[away_team[g], away_match_idx[g]]]
          } else {
            home_att_val <- att[home_team[g], z[home_team[g], home_match_idx[g]]]
            home_def_val <- def[home_team[g], z[home_team[g], home_match_idx[g]]]
            away_att_val <- att[t, k]
            away_def_val <- def[t, k]
          }
          
          log_lik <- log_lik_match(data$HomeGoals[g], data$AwayGoals[g],
                                   home_att_val, away_att_val,
                                   home_def_val, away_def_val, home)
          
          log_probs[k] <- log_prior + log_lik
        }
        
        max_log_prob <- max(log_probs)
        probs <- exp(log_probs - max_log_prob)
        probs <- probs / sum(probs)
        
        z[t, s] <- sample(1:K, 1, prob = probs)
      }
    }
    
    # Campiono att_{t,k} | ...
    for (t in 1:T) {
      for (k in 1:K) {
        home_matches <- which(home_team == t & z[t, home_match_idx] == k)
        away_matches <- which(away_team == t & z[t, away_match_idx] == k)
        
        log_target_att <- function(att_val) {
          log_prior <- -0.5 * (att_val - mu_att[k])^2 / sigma2_att[k]
          
          log_lik <- 0
          for (g in home_matches) {
            lambda <- exp(home + att_val + def[away_team[g], z[away_team[g], away_match_idx[g]]])
            log_lik <- log_lik + dpois(data$HomeGoals[g], lambda, log = TRUE)
          }
          for (g in away_matches) {
            lambda <- exp(att_val + def[home_team[g], z[home_team[g], home_match_idx[g]]])
            log_lik <- log_lik + dpois(data$AwayGoals[g], lambda, log = TRUE)
          }
          
          return(log_prior + log_lik)
        }
        
        att[t, k] <- mh_sample(att[t, k], log_target_att, proposal_sd = 0.25, max_iter = 5)
      }
    }
    
    # Campiono def_{t,k} | ...
    for (t in 1:T) {
      for (k in 1:K) {
        home_matches <- which(home_team == t & z[t, home_match_idx] == k)
        away_matches <- which(away_team == t & z[t, away_match_idx] == k)
        
        log_target_def <- function(def_val) {
          log_prior <- -0.5 * (def_val - mu_def[k])^2 / sigma2_def[k]
          
          log_lik <- 0
          for (g in home_matches) {
            lambda <- exp(att[away_team[g], z[away_team[g], away_match_idx[g]]] + def_val)
            log_lik <- log_lik + dpois(data$AwayGoals[g], lambda, log = TRUE)
          }
          for (g in away_matches) {
            lambda <- exp(home + att[home_team[g], z[home_team[g], home_match_idx[g]]] + def_val)
            log_lik <- log_lik + dpois(data$HomeGoals[g], lambda, log = TRUE)
          }
          
          return(log_prior + log_lik)
        }
        
        def[t, k] <- mh_sample(def[t, k], log_target_def, proposal_sd = 0.25, max_iter = 5)
      }
    }
    
    # Campiono home | ...
    log_target_home <- function(home_val) {
      log_prior <- -0.5 * home_val^2 / 10000
      
      log_lik <- 0
      for (g in 1:G) {
        lambda <- exp(home_val + 
                        att[home_team[g], z[home_team[g], home_match_idx[g]]] + 
                        def[away_team[g], z[away_team[g], away_match_idx[g]]])
        log_lik <- log_lik + dpois(data$HomeGoals[g], lambda, log = TRUE)
      }
      
      return(log_prior + log_lik)
    }
    
    home <- mh_sample(home, log_target_home, proposal_sd = 0.2, max_iter = 5)
    
    # Campiono mu_att, sigma2_att | ...
    for (k in 1:K) {
      prec_prior <- 1/10000
      prec_lik <- T / sigma2_att[k]
      mu_att_post <- (sum(att[, k]) / sigma2_att[k]) / (prec_prior + prec_lik)
      sigma2_att_post <- 1 / (prec_prior + prec_lik)
      mu_att[k] <- rnorm(1, mu_att_post, sqrt(sigma2_att_post))
      
      shape_post <- 0.01 + T/2
      rate_post <- 0.01 + 0.5 * sum((att[, k] - mu_att[k])^2)
      sigma2_att[k] <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    }
    
    # Campiono mu_def, sigma2_def | ...
    for (k in 1:K) {
      prec_lik <- T / sigma2_def[k]
      mu_def_post <- (sum(def[, k]) / sigma2_def[k]) / (prec_prior + prec_lik)
      sigma2_def_post <- 1 / (prec_prior + prec_lik)
      mu_def[k] <- rnorm(1, mu_def_post, sqrt(sigma2_def_post))
      
      rate_post <- 0.01 + 0.5 * sum((def[, k] - mu_def[k])^2)
      sigma2_def[k] <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    }
    
    # Campiono P| ...
    n_transitions <- matrix(0, nrow = K, ncol = K)
    
    for (t in 1:T) {
      for (s in 2:S_t[t]) {
        from_state <- z[t, s-1]
        to_state <- z[t, s]
        n_transitions[from_state, to_state] <- n_transitions[from_state, to_state] + 1
      }
    }
    
    for (i in 1:K) {
      alpha_post <- c(1, 1) + n_transitions[i, ]
      P[i, ] <- rdirichlet(1, alpha_post)
    }
    
    # Campiono pi | ...
    n_initial <- numeric(K) + 1
    for (t in 1:T) {
      initial_state <- z[t, 1]
      n_initial[initial_state] <- n_initial[initial_state] + 1
    }
    
    pi <- rdirichlet(1, n_initial)
    
    # Vincoli di identificabilità
    for (k in 1:K) {
      att[, k] <- att[, k] - mean(att[, k])
      def[, k] <- def[, k] - mean(def[, k])
      mu_att[k] <- mu_att[k] - mean(att[, k])
      mu_def[k] <- mu_def[k] - mean(def[, k])
    }
    
    att_all[, , iter] <- att
    def_all[, , iter] <- def
    home_all[iter] <- home
    z_all[, , iter] <- z
    mu_att_all[iter, ] <- mu_att
    sigma2_att_all[iter, ] <- sigma2_att
    mu_def_all[iter, ] <- mu_def
    sigma2_def_all[iter, ] <- sigma2_def
    pi_all[iter, ] <- pi
    P_all[, , iter] <- P
  }
  
 
  #burn-in e thinning
  indices <- seq(burn_in + 1, n_iter, by = thin)
  n_samples <- length(indices)
  
  results <- list(
    att = att_all[, , indices],
    def = def_all[, , indices],
    home = home_all[indices],
    z = z_all[, , indices],
    mu_att = mu_att_all[indices, ],
    sigma2_att = sigma2_att_all[indices, ],
    mu_def = mu_def_all[indices, ],
    sigma2_def = sigma2_def_all[indices, ],
    pi = pi_all[indices, ],
    P = P_all[, , indices],
    teams = teams,
    n_iter = n_iter,
    burn_in = burn_in,
    thin = thin,
    K = K,
    n_samples = n_samples
  )
  
  return(results)
}



###################################################################################
data <- read.csv("./2015-16/eng.1.csv", stringsAsFactors = FALSE)
samples <- run_mcmc(data, n_iter = 50000, burn_in = 10000, thin = 10, K = 2, seed = min(353244,306256))


teams <- samples$teams
T <- length(teams)
K <- samples$K
n_samples <- samples$n_samples
S_t <- apply(samples$z[, , 1], 1, function(x) sum(!is.na(x)))  # Conta i valori non-NA per ogni squadra

att_mean <- apply(samples$att, c(1, 2), mean)
def_mean <- apply(samples$def, c(1, 2), mean)
home_mean <- mean(samples$home)

mu_att_mean <- colMeans(samples$mu_att)
sigma2_att_mean <- colMeans(samples$sigma2_att)
mu_def_mean <- colMeans(samples$mu_def)
sigma2_def_mean <- colMeans(samples$sigma2_def)

pi_mean <- colMeans(samples$pi)
P_mean <- apply(samples$P, c(1, 2), mean)

cat("Effetto casa medio:", home_mean, "\n")
cat("\nMatrice di transizione media:\n")
print(P_mean)


#Proporzione di tempo che ogni squadra passa in ciascuno stato
state_proportions <- matrix(0, nrow = T, ncol = K)
for (t in 1:T) {
  # Usiamo tutti i campioni per una stima più robusta
  states <- matrix(samples$z[t, 1:S_t[t], ], nrow = S_t[t], ncol = n_samples)
  for (k in 1:K) {
    state_proportions[t, k] <- mean(states == k, na.rm = TRUE)
  }
}

rownames(state_proportions) <- teams
colnames(state_proportions) <- c("In forma", "Fuori forma")

cat("\nPROPORZIONE DI TEMPO IN OGNI STATO\n")
print(round(state_proportions, 3))


#Classifica basata sui parametri di attacco e difesa 
ranking <- data.frame(
  Team = teams,
  Attack_in_forma = att_mean[, 1],
  Defense_in_forma = def_mean[, 1],
  Attack_fuori_forma = att_mean[, 2],
  Defense_fuori_forma = def_mean[, 2],
  Overall_score = (att_mean[, 1]-def_mean[, 1])*state_proportions[, 1] + (att_mean[, 2]-def_mean[, 2])*state_proportions[, 2],  
  Prop_in_forma = state_proportions[, 1]
)

ranking <- ranking[order(-ranking$Overall_score), ]

cat("\nCLASSIFICA BASATA SUI PARAMETRI\n")
print(ranking, row.names = FALSE)


#Grafici
par(mfrow = c(2, 2))
barplot(att_mean[, 1], names.arg = teams, las = 2, cex.names = 0.6,
        main = "Parametri di Attacco - Stato 1 (in forma)",
        ylab = "Valore", col = "skyblue")

barplot(def_mean[, 1], names.arg = teams, las = 2, cex.names = 0.6,
        main = "Parametri di Difesa - Stato 1 (in forma)",
        ylab = "Valore", col = "lightcoral")

barplot(att_mean[, 2], names.arg = teams, las = 2, cex.names = 0.6,
        main = "Parametri di Attacco - Stato 2",
        ylab = "Valore", col = "purple")

barplot(def_mean[, 2], names.arg = teams, las = 2, cex.names = 0.6,
        main = "Parametri di Difesa - Stato 2",
        ylab = "Valore", col = "orange")
par(mfrow = c(1, 1))


par(mfrow = c(1, 2))
hist(samples$home, breaks = 30, main = "Distribuzione Effetto Casa",
     xlab = "Valore", ylab = "Frequenza", col = "lightgreen", border = "white")

plot(samples$home, type = "l", main = "Trace plot effetto casa", 
     xlab = "Campione", ylab = "Valore", col = "blue")
par(mfrow = c(1, 1))


par(mfrow = c(1, 1))
image(1:K, 1:K, P_mean, col = heat.colors(12), 
      main = "Matrice di Transizione Media",
      xlab = "Stato precedente", ylab = "Stato successivo", axes = FALSE)
axis(1, at = 1:K, labels = c("In forma", "Fuori forma"))
axis(2, at = 1:K, labels = c("In forma", "Fuori forma"))
text(expand.grid(1:K, 1:K), sprintf("%.3f", P_mean), cex = 1.2)
par(mfrow = c(1, 1))


#Trace plots e acf 
par(mfrow = c(2, 4))
plot(samples$mu_att[, 1], type = "l", main = "Trace plot mu_att (stato 1)", 
     xlab = "Campione", ylab = "Valore", col = "red")

plot(samples$mu_att[, 2], type = "l", main = "Trace plot mu_att (stato 2)", 
     xlab = "Campione", ylab = "Valore", col = "blue")

plot(samples$mu_def[, 1], type = "l", main = "Trace plot mu_def (stato 1)", 
     xlab = "Campione", ylab = "Valore", col = "green")

plot(samples$mu_def[, 2], type = "l", main = "Trace plot mu_def (stato 2)", 
     xlab = "Campione", ylab = "Valore", col = "grey")

acf(samples$mu_att[, 1], main = "acf mu_att (stato 1)")
acf(samples$mu_att[, 2], main = "acf mu_att (stato 2)")
acf(samples$mu_def[, 1], main = "acf mu_def (stato 1)")
acf(samples$mu_def[, 2], main = "acf mu_def (stato 2)")
par(mfrow = c(1, 1))


par(mfrow = c(2, 4))
plot(samples$sigma2_att[,1], type = "l", main = "Trace plot sigma2_att (stato 1)", 
     xlab = "Campione", ylab = "Valore", col = "orange")

plot(samples$sigma2_att[, 2], type = "l", main = "Trace plot sigma2_att (stato 2)", 
     xlab = "Campione", ylab = "Valore", col = "pink")

plot(samples$sigma2_def[, 1], type = "l", main = "Trace plot sigma2_def (stato 1)", 
     xlab = "Campione", ylab = "Valore", col = "brown")

plot(samples$sigma2_def[, 2], type = "l", main = "Trace plot sigma2_def (stato 2)", 
     xlab = "Campione", ylab = "Valore", col = "black")

acf(samples$sigma2_att[, 1], main = "acf sigma2_att (stato 1)")
acf(samples$sigma2_att[, 2], main = "acf sigma2_att (stato 2)")
acf(samples$sigma2_def[, 1], main = "acf sigma2_def (stato 1)")
acf(samples$sigma2_def[, 2], main = "acf sigma2_def (stato 2)")
par(mfrow = c(1, 1))



#Grafici per una squadra specifica
team_index <- which(teams == "Manchester City")

if (length(team_index) > 0) {
  par(mfrow = c(2, 4))
  
  plot(samples$att[team_index, 1, ], type = "l", main = "Trace plot att (Manchester City, stato 1)",
       xlab = "Campione", ylab = "Valore", col = "blue")
  
  plot(samples$att[team_index, 2, ], type = "l", main = "Trace plot att (Manchester City, stato 2)",
       xlab = "Campione", ylab = "Valore", col = "red")
  
  acf(samples$att[team_index, 1, ], main = "ACF att (Manchester City, stato 1)")
  acf(samples$att[team_index, 2, ], main = "ACF att (Manchester City, stato 2)")
  
  plot(samples$def[team_index, 1, ], type = "l", main = "Trace plot def (Manchester City, stato 1)",
       xlab = "Campione", ylab = "Valore", col = "blue")
  
  plot(samples$def[team_index, 2, ], type = "l", main = "Trace plot def (Manchester City, stato 2)",
       xlab = "Campione", ylab = "Valore", col = "red")
  
  acf(samples$def[team_index, 1, ], main = "ACF def (Manchester City, stato 1)")
  acf(samples$def[team_index, 2, ], main = "ACF def (Manchester City, stato 2)")
  
  par(mfrow = c(1, 1))
  
  cat("\n=== STATISTICHE MANCHESTER CITY ===\n")
  cat("Attacco stato 1:", mean(samples$att[team_index, 1, ]), "\n")
  cat("Difesa stato 1:", mean(samples$def[team_index, 1, ]), "\n")
  cat("Attacco stato 2:", mean(samples$att[team_index, 2, ]), "\n")
  cat("Difesa stato 2:", mean(samples$def[team_index, 2, ]), "\n")
  cat("Proporzione in forma:", state_proportions[team_index, 1], "\n")
}


# Intervallo di credibilità per l'effetto casa
home_ci <- quantile(samples$home, c(0.025, 0.975))
cat("\nIntervallo di credibilità 95% per l'effetto casa:\n")
print(home_ci)









###############################################################################
# Predizione

predict_match <- function(samples, team_home, team_away, data) {
  
  teams <- samples$teams
  n_samples <- samples$n_samples
  K <- samples$K
  
  #Trovo gli indici delle squadre
  team_home_idx <- which(teams == team_home)
  team_away_idx <- which(teams == team_away)
  
  # Preparazione dei dati
  results <- strsplit(data$FT, "-")
  data$HomeGoals <- as.numeric(sapply(results, function(x) x[1]))
  data$AwayGoals <- as.numeric(sapply(results, function(x) x[2]))
  
  G <- nrow(data)
  T_tot <- length(teams)
  
  home_team <- sapply(data$Team.1, function(x) which(teams == x))
  away_team <- sapply(data$Team.2, function(x) which(teams == x))
  
  match_counter <- rep(1, T_tot)
  home_match_idx <- rep(0, G)
  away_match_idx <- rep(0, G)
  
  for (g in 1:G) {
    ht <- home_team[g]
    at <- away_team[g]
    
    home_match_idx[g] <- match_counter[ht]
    away_match_idx[g] <- match_counter[at]
    
    match_counter[ht] <- match_counter[ht] + 1
    match_counter[at] <- match_counter[at] + 1
  }
  
  S_t <- match_counter - 1
  
  
  # Estraggo l'ulitmo stato per ogni squadra, da ogni campione
  last_state_home <- numeric(n_samples)
  last_state_away <- numeric(n_samples)
  
  for (m in 1:n_samples) {
    last_s_home <- S_t[team_home_idx]
    last_s_away <- S_t[team_away_idx]
    
    last_state_home[m] <- samples$z[team_home_idx, last_s_home, m]
    last_state_away[m] <- samples$z[team_away_idx, last_s_away, m]
  }
  
  
  #Calcolo lo stato in cui dovrà essere la squadra per la partita che devo predirre usando 
  #la matrice di transizone associata al campione
  new_state_home <- numeric(n_samples)
  new_state_away <- numeric(n_samples)
  
  for (m in 1:n_samples) {
    prob_home <- samples$P[last_state_home[m], , m]
    new_state_home[m] <- sample(1:K, size = 1, prob = prob_home)
    
    prob_away <- samples$P[last_state_away[m], , m]
    new_state_away[m] <- sample(1:K, size = 1, prob = prob_away)
  }
  
  
  #Calcolo i parametri per simulare le poisson per ogni campione
  lambda_home <- numeric(n_samples)
  lambda_away <- numeric(n_samples)
  
  for (m in 1:n_samples) {
    att_home <- samples$att[team_home_idx, new_state_home[m], m]
    def_home <- samples$def[team_home_idx, new_state_home[m], m]
    
    att_away <- samples$att[team_away_idx, new_state_away[m], m]
    def_away <- samples$def[team_away_idx, new_state_away[m], m]
    
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
  
  dist_goals_home <- table(goals_home_samples) / n_samples
  dist_goals_away <- table(goals_away_samples) / n_samples
  
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
    dist_goals_home = dist_goals_home,
    dist_goals_away = dist_goals_away,
    prob_home_win = prob_home_win,
    prob_draw = prob_draw,
    prob_away_win = prob_away_win,
    joint_pmf = joint_pmf,
    n_samples = n_samples,
    lambda_home = lambda_home,
    lambda_away = lambda_away,
    new_state_home = new_state_home,
    new_state_away = new_state_away,
    last_state_home = last_state_home,
    last_state_away = last_state_away
  )
  
  return(result)
}


#Nota: la funzione per avere un output raccolto e pulito in formato ASCII della funzione è stata creata con AI
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
  
  # Converti la PMF in data frame correttamente
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



pred <- predict_match(samples, team_home = "Arsenal", team_away = "Chelsea", data = data)
print_prediction(pred)

