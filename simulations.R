#!/usr/bin/env Rscript

library(parallel) 
library(purrr) # install.packages("purrr")
library(tibble) # install.packages("tibble")
library(glue) # install.packages("glue")
library(dplyr) # install.packages("dplyr")
library(vroom) # install.packages("vroom")
library(pbmcapply) # install.packages("pbmcapply")
library(Rcpp) # install.packages("Rcpp")
library(openxlsx) # install.packages("openxlsx")
library(tictoc) # install.packages("tictoc")

rm(list = ls(all = TRUE))

Rcpp::sourceCpp("src/log_likelihood_rcpp.cpp")

s_zero_adjusted <- function(S){
  function(t, p0, ...){
    if(t == 0) return(1 - p0)
    (1 - p0) * S(t, ...)
  }
}

pdf_zero_adjusted <- function(f){
  function(t, p0, ...){
    if(t == 0) return(p0)
    (1 - p0) * f(t, ...)
  }
}

f_gompertz <- function(t, a, b, theta) {
  b * exp(a * t) * (1 + (theta * b / a) * (exp(a * t) - 1))^(-1 - 1 / theta)
}

S_gompertz <- function(t, a, b, theta) {
  (1 + theta * b / a * (exp(a * t) - 1))^(-1/theta)
}

pdf <- pdf_zero_adjusted(f = f_gompertz) # Densidade da Zero adjusted Gompertz
s <- s_zero_adjusted(S = S_gompertz) # Sobrevivencia da Zero adjusted Gompertz

geracao <- function(n, a, beta00, beta01, beta10, beta11, theta, tau) {
  # n: número de amostras
  # a: parâmetro a da Gompertz
  # theta: parâmetro da gama
  # beta00, beta01: parâmetros de regressão para p0x
  # beta10, beta11: parâmetros de regressão para b(x)
  # tau: parâmetro da censura (taxa exponencial)
  
  # 1. Gerar covariáveis X e ui ~ Uniform(0,1)
  x <- rbinom(n, 1, 0.5)
  u <- runif(n)
  
  # 2. Calcular p0x e b(x)
  p0x <- exp(beta00 + beta01 * x) / (1 + exp(beta00 + beta01 * x))  # p0(x)
  b_x <- exp(beta10 + beta11 * x)                                   # b(x) com base na covariável
  
  # 3. Calcular p1x para a fração de cura
  p1x <- (1 - p0x) * (1 - theta * b_x / a)^(-1 / theta)
  
  # 4. Calcular tempos de sobrevivência s
  quantile_gompertz <- function(u, a, b, p0, theta) {
    return((1 / a) * log(1 + a / (b * theta) * (((1 - u) / (1 - p0))^(-theta) - 1)))
  }
  
  # Vetorização
  s <- ifelse(
    u <= p0x, 
    0,  # Evento imediato (t = 0)
    ifelse(
      u > (1 - p1x),
      Inf,  # Curado (t = Inf)
      quantile_gompertz(runif(n, min = p0x, max = 1 - p1x), a, b_x, p0x, theta)  # Tempo gerado
    )
  )
  
  # 5. Gerar censura wi ~ Exp(tau) e tempos observados t_i = min(s_i, wi)
  w <- rexp(n, rate = tau)
  t <- pmin(s, w)  # Observamos o menor entre o tempo gerado e o tempo censurado
  
  # 6. Indicador de censura: delta = 1 se o evento foi observado, 0 se censurado
  delta <- as.integer(t < w)
  
  # 7. Indicador de se o tempo é maior que zero: delta_star = 1 se t > 0, 0 se t = 0
  delta_star <- as.integer(t > 0)
  
  # Retornar um data frame com os dados de sobrevivência
  return(data.frame(t = t, x = x, delta = delta, delta_star = delta_star))
}


log_likelihood <- function(par, dados) {
  
  beta00 <- par[1L]
  beta01 <- par[2L]
  beta10 <- par[3L]
  beta11 <- par[4L]
  a <- par[5L]
  theta<- par[6L]
  
  t <- dados$t
  delta <- dados$delta
  x <- dados$x
  
  p0x <- exp(beta00 + beta01 * x) / (1 + exp(beta00 + beta01 * x))
  b_x <- exp(beta10 + beta11 * x)
  
  # Log-verossimilhança
  log_like <- numeric(length(t))  # Vetor para armazenar a log-verossimilhança
  
  one_step <- function(i){
    if (delta[i] == 1) {  # Evento observado
      if (t[i] == 0) {
        log_like[i] <- log(p0x[i])  # Caso t = 0 (evento imediato)
      } else {
        log_like[i] <- log(pdf(t = t[i], a = a, b = b_x[i], p0 = p0x[i], theta = theta))  # Para t > 0
      }
    } else {  # Censura
      log_like[i] <- log(s(t = t[i], a = a, b = b_x[i], p0 = p0x[i], theta = theta))  # Censura
    }
  }
  sum(sapply(X = 1L:length(t), FUN = one_step))
}

# Parametros reais --------------------------------------------------------
n <- 500          # Tamanho da amostra
a <- -1          # Parâmetro a da Gompertz
beta00 <- -1     # Intercepto para p0x
beta01 <- 1      # Efeito da covariável em p0x
beta10 <- -1      # Intercepto para b(x)
beta11 <- 1      # Efeito da covariável em b(x)
tau <- 0.2         # Taxa de censura (para a distribuição exponencial)
theta<- 0.75 
p00 <- exp(beta00 + beta01 * 0) / (1 + exp(beta00 + beta01 * 0))  # p0(x)
p01 <- exp(beta00 + beta01 * 1) / (1 + exp(beta00 + beta01 * 1))  # p0(x)
b_0 <- exp(beta10 + beta11 * 0)
b_1 <- exp(beta10 + beta11 * 1)
cura_0 <- (1-p00)*(1-theta*b_0/a)^(-1/theta)
cura_1 <- (1-p01)*(1-theta*b_1/a)^(-1/theta)

mc <- 
  function(
    n_mc = 10L,
    n_boot = 5L,
    sig = 0.05,
    n = 500,
    tau = tau,
    beta00 = beta00, 
    beta01 = beta01, 
    beta10 = beta10, 
    beta11 = beta11, 
    a = a,
    theta=theta,
    paralelo = TRUE
  ) {
    
    try_optim <- function(...){
      tryCatch(
        expr = optim(...),
        error = function(e) NULL
      )
    }
    
    one_step <- function(m){
      # Amostra original da replica de Monte-Carlo - MC
      repeat {
        # Gerar dados
        dados <- geracao(
          n = n,
          a = a,
          beta00 = beta00,
          beta01 = beta01,
          beta10 = beta10,
          beta11 = beta11,
          theta = theta,
          tau = tau
        )
        
        # Executar otimização
        result_mc <- suppressWarnings(try_optim(
          par = c(0.5, 0.5, 0.5, 0.5, -0.5, 0.5),
          fn = \(par) log_likelihood_rcpp(par = par, dados = dados),
          control = list(fnscale = -1),
          method = "BFGS"
        ))
        
        # Verificar se a otimização convergiu
        if (!is.null(result_mc) && result_mc$convergence == 0L) {
          break
        }
      }
      
      est_mc <- result_mc$par
      
      hat_beta00 <- est_mc[1L]
      hat_beta01 <- est_mc[2L]
      hat_beta10 <- est_mc[3L]
      hat_beta11 <- est_mc[4L]
      hat_a <- est_mc[5L]
      hat_theta <- est_mc[6L]
      
      hat_p00 <- exp(hat_beta00 + hat_beta01 * 0)/(1 + exp(hat_beta00 + hat_beta01 * 0))  # p0(x)
      hat_p01 <- exp(hat_beta00 + hat_beta01 * 1)/(1 + exp(hat_beta00 + hat_beta01 * 1))  # p0(x)
      hat_b_0 <- exp(hat_beta10 + hat_beta11 * 0)
      hat_b_1 <- exp(hat_beta10 + hat_beta11 * 1)
      hat_cura_0 <- (1 - hat_p00) * (1 - hat_theta*hat_b_0/hat_a)^(-1/hat_theta)
      hat_cura_1 <- (1 - hat_p01) * (1 - hat_theta*hat_b_1/hat_a)^(-1/hat_theta)
      
      est_mc <- c(est_mc, hat_p00, hat_p01, hat_b_0, hat_b_1, hat_cura_0, hat_cura_1)
      
      names(est_mc) <- c(
        "beta00",
        "beta01",
        "beta10",
        "beta11",
        "a",
        "theta",
        "p00",
        "p01",
        "b_0",
        "b_1",
        "cura_0",
        "cura_1"
      )
      
      # Realizando bootstrap
      one_boot <- function(b){
        repeat{
          id <- sample(1L:n, size = n, replace = TRUE)
          dados_boot <- dados[id, ]
          
          result_boot <-
            try_optim(
              par = c(0.5, 0.5, 0.5, 0.5, -0.5, 0.5),
              fn = \(par) log_likelihood_rcpp(par = par, dados = dados_boot),
              control = list(fnscale = -1),
              method = "BFGS"
            )
          
          if (!is.null(result_boot) && result_boot$convergence == 0L) {
            break
          }
        }
        
        est_boot <- result_boot$par
        
        hat_beta00 <- est_boot[1L]
        hat_beta01 <- est_boot[2L]
        hat_beta10 <- est_boot[3L]
        hat_beta11 <- est_boot[4L]
        hat_a <- est_boot[5L]
        hat_theta <- est_boot[6L]
        
        hat_p00 <- exp(hat_beta00 + hat_beta01 * 0)/(1 + exp(hat_beta00 + hat_beta01 * 0))  # p0(x)
        hat_p01 <- exp(hat_beta00 + hat_beta01 * 1)/(1 + exp(hat_beta00 + hat_beta01 * 1))  # p0(x)
        hat_b_0 <- exp(hat_beta10 + hat_beta11 * 0)
        hat_b_1 <- exp(hat_beta10 + hat_beta11 * 1)
        hat_cura_0 <- (1 - hat_p00) * (1 - hat_theta * hat_b_0/hat_a) ^ (-1/hat_theta)
        hat_cura_1 <- (1 - hat_p01) * (1 - hat_theta * hat_b_1/hat_a) ^ (-1/hat_theta)
        
        est_boot <- c(est_boot, hat_p00, hat_p01, hat_b_0, hat_b_1, hat_cura_0, hat_cura_1)
        
        names(est_boot) <- c(
          "beta00",
          "beta01",
          "beta10",
          "beta11",
          "a",
          "theta",
          "p00",
          "p01",
          "b_0",
          "b_1",
          "cura_0",
          "cura_1"
        )
        est_boot
      }
      
      # Tibble com as estimativas bootstrap (cada columa um parametro)
      est_boot <- purrr::map_dfr(1L:n_boot, one_boot)
      
      # Estimador coorigido por vies via bootstrap
      est_coorigida_boot <- 2 * est_mc - apply(est_boot, MARGIN = 2L, FUN = mean)
      
      ic_boot <- 
        apply(est_boot, MARGIN = 2L, FUN = quantile, probs = c(sig/2, 1 - sig/2)) |>
        as_tibble() # Primeira linha é o limite inferior e a segunda é o limite superior
      
      # Checando cobertura dos intervalos
      cobertura_beta00 <- if(beta00 > ic_boot$beta00[1L] & beta00 < ic_boot$beta00[2L]) TRUE else FALSE
      cobertura_beta01 <- if(beta01 > ic_boot$beta01[1L] & beta01 < ic_boot$beta01[2L]) TRUE else FALSE
      cobertura_beta10 <- if(beta10 > ic_boot$beta10[1L] & beta10 < ic_boot$beta10[2L]) TRUE else FALSE
      cobertura_beta11 <- if(beta11 > ic_boot$beta11[1L] & beta11 < ic_boot$beta11[2L]) TRUE else FALSE
      cobertura_a <- if(a > ic_boot$a[1L] & a < ic_boot$a[2L]) TRUE else FALSE
      cobertura_theta <- if(theta > ic_boot$theta[1L] & theta < ic_boot$theta[2L]) TRUE else FALSE
      cobertura_p00 <- if(p00 > ic_boot$p00[1L] & p00 < ic_boot$p00[2L]) TRUE else FALSE
      cobertura_p01 <- if(p01 > ic_boot$p01[1L] & p01 < ic_boot$p01[2L]) TRUE else FALSE
      cobertura_b_0 <- if(b_0 > ic_boot$b_0[1L] & b_0 < ic_boot$b_0[2L]) TRUE else FALSE
      cobertura_b_1 <- if(b_1 > ic_boot$b_1[1L] & b_1 < ic_boot$b_1[2L]) TRUE else FALSE
      cobertura_cura_0 <- if(cura_0 > ic_boot$cura_0[1L] & cura_0 < ic_boot$cura_0[2L]) TRUE else FALSE
      cobertura_cura_1 <- if(cura_1 > ic_boot$cura_1[1L] & cura_1 < ic_boot$cura_1[2L]) TRUE else FALSE
      
      ic_boot <- ic_boot |> t() |> as.data.frame()
      
      colnames(ic_boot) <- c("li", "ls")
      ic_boot <- 
        ic_boot |>
        dplyr::mutate(amplitudade = ls - li) 
      
      tibble(
        id_mc = m,
        n = n,
        nomes_parametros = c(
          "beta00",
          "beta01",
          "beta10",
          "beta11",
          "a",
          "theta",
          "p00",
          "p01",
          "b_0",
          "b_1",
          "cura_0",
          "cura_1"
        ),
        parametros_verdadeiros = c(
          beta00,
          beta01,
          beta10,
          beta11,
          a,
          theta,
          p00,
          p01,
          b_0,
          b_1,
          cura_0,
          cura_1
        ),
        estimativas_mc = est_mc,
        estimativas_corrigida_boot = est_coorigida_boot,
        cobertura_ic_boot = c(
          cobertura_beta00,
          cobertura_beta01,
          cobertura_beta10,
          cobertura_beta11,
          cobertura_a,
          cobertura_theta,
          cobertura_p00,
          cobertura_p01,
          cobertura_b_0,
          cobertura_b_1,
          cobertura_cura_0,
          cobertura_cura_1
        ),
        amplitude = ic_boot$amplitudade,
        li_boot = ic_boot$li,
        ls_boot = ic_boot$ls
      )
      
    } # Fim da funcao one_step
    
    try_one_step <- function(...){
      tryCatch(one_step(...), error = function(e) NULL)
    }
    
    replicas_mc <- function(m){
      repeat{
        r <- try_one_step(m)
        if(!is.null(r)){
          break
        }
      }
      return(r)
    }
    
    if(paralelo){
      cores <- parallel::detectCores()
    } else { 
      cores <- 1L
    }
    
    t0 <- Sys.time()
    r <- pbmclapply(
      X = 1L:n_mc,         
      FUN = replicas_mc,    # Função a ser executada
      mc.cores = cores      # Número de núcleos
    )
    t1 <- Sys.time()
    
    time_secs <- difftime(t1, t0, units = "secs")[[1L]]
    time_mins <- difftime(t1, t0, units = "mins")[[1L]]
    time_hours <- difftime(t1, t0, units = "hours")[[1L]]
    
    r <- dplyr::bind_rows(r)
    r$tau <- tau
    r$tempo_segundos <- time_secs
    r$tempo_minutos <- time_mins
    r$tempo_horas <- time_hours
    
    
    r <-
      r |>
      dplyr::group_by(nomes_parametros) |>
        dplyr::mutate(
          vies_quadrado_mc = (parametros_verdadeiros - estimativas_mc) ^ 2,
          vies_quadrado_boot = (parametros_verdadeiros - estimativas_corrigida_boot) ^ 2,
          var_mc = var(estimativas_mc),
          var_boot = var(estimativas_corrigida_boot),
          eqm_quadrado_mc = vies_quadrado_mc + var_mc,
          eqm_estimativa_corrigida_boot = vies_quadrado_boot + var_boot,
          probabilidade_cobertura_boot = sum(cobertura_ic_boot) / n_mc
        ) |>
      summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
    
    openxlsx::write.xlsx(
      r,
      glue::glue("data/simulacoes_n_{n}_tau_{tau}.xlsx")
    )
    
    return(r)
    
  } # Fim da funcao MC

simulacao <- 
  function(
    n = c(50, 100, 250, 500, 1000, 2500, 5000),
    n_mc = 2L,
    n_boot = 10L,
    sig = 0.05,
    tau = tau,
    beta00 = beta00, 
    beta01 = beta01, 
    beta10 = beta10, 
    beta11 = beta11, 
    a = a,
    theta=theta,
    paralelo = TRUE
  ){
    
    r <- purrr::map_dfr(tau, function(current_tau) {
      purrr::map_dfr(n, function(current_n) {
        mc(
          n_mc = n_mc,
          n_boot = n_boot,
          sig = sig,
          n = current_n,
          tau = current_tau,
          beta00 = beta00,
          beta01 = beta01,
          beta10 = beta10,
          beta11 = beta11,
          a = a,
          theta = theta,
          paralelo = paralelo
        )
      })
    })
    r$tau <- rep(x = tau, each = nrow(r)/length(tau))
    r
  }

RNGkind("L'Ecuyer-CMRG")
set.seed(123)
mc.reset.stream()

tic()
resultados <- simulacao(
  n = c(50, 100, 250, 500, 1000, 2500, 5000),
  n_mc = 5000,
  n_boot = 500,
  sig = 0.05,
  tau = c(0.2, 0.4, 0.6),
  beta00 = beta00,
  beta01 = beta01,
  beta10 = beta10,
  beta11 = beta11,
  a = a,
  theta = theta,
  paralelo = TRUE
) 
toc()

# Salvando os resultados da simulacao -------------------------------------
openxlsx::write.xlsx(resultados, "data/resultados_simulacao.xlsx")