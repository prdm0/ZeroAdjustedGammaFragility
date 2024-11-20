#include <Rcpp.h>
using namespace Rcpp;

// Função Gompertz para a densidade de probabilidade (pdf)
double f_gompertz(double t, double a, double b, double theta) {
  return b * exp(a * t) * pow(1 + (theta * b / a) * (exp(a * t) - 1), -1 - 1 / theta);
}

// Função de sobrevivência Gompertz
double S_gompertz(double t, double a, double b, double theta) {
  return pow(1 + theta * b / a * (exp(a * t) - 1), -1 / theta);
}

// Função pdf ajustada para t = 0
double pdf_zero_adjusted(double t, double p0, double a, double b, double theta) {
  if (t == 0) {
    return p0;  // Retorna p0 se t = 0
  } else {
    return (1 - p0) * f_gompertz(t, a, b, theta);  // Para t > 0
  }
}

// Função s ajustada para t = 0
double s_zero_adjusted(double t, double p0, double a, double b, double theta) {
  if (t == 0) {
    return 1 - p0;  // Retorna (1 - p0) se t = 0
  } else {
    return (1 - p0) * S_gompertz(t, a, b, theta);  // Para t > 0
  }
}

// [[Rcpp::export]]
double log_likelihood_rcpp(NumericVector par, DataFrame dados) {
  // Extração de parâmetros
  double beta00 = par[0];
  double beta01 = par[1];
  double beta10 = par[2];
  double beta11 = par[3];
  double a = par[4];
  double theta = par[5];
  
  // Extração dos dados
  NumericVector t = dados["t"];
  NumericVector delta = dados["delta"];
  NumericVector x = dados["x"];
  
  int n = t.size();
  
  // Cálculo de p0x e b(x)
  NumericVector p0x(n);
  NumericVector b_x(n);
  
  for (int i = 0; i < n; i++) {
    double exp_p0x = exp(beta00 + beta01 * x[i]);
    p0x[i] = exp_p0x / (1 + exp_p0x);
    b_x[i] = exp(beta10 + beta11 * x[i]);
  }
  
  // Cálculo da log-verossimilhança
  double log_like = 0.0;
  
  for (int i = 0; i < n; i++) {
    if (delta[i] == 1) {  // Evento observado
      if (t[i] == 0) {
        log_like += log(p0x[i]);  // Caso t = 0
      } else {
        log_like += log(pdf_zero_adjusted(t[i], p0x[i], a, b_x[i], theta));  // Para t > 0
      }
    } else {  // Censura
      log_like += log(s_zero_adjusted(t[i], p0x[i], a, b_x[i], theta));
    }
  }
  
  return log_like;
}
