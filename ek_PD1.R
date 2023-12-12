library(mvtnorm)
library(manipulate)

rawData <- read.csv("heartRisk.csv", header = TRUE, sep = ",", dec = ".")
data2 <- rawData[, c("Risk","Age","Systolic","HDL")]
data2

#sampling
set.seed(123)
randomIndex <- sample(nrow(data2), 30)
data <- data2[randomIndex, ]
row.names(data) <- NULL
summary(data)

#hist(data$Age, main="Age")
#hist(data$Systolic, main="Systolic")
#hist(data$HDL, main="HDL")
#hist(data$Risk, main="Risk")

data$Age = log(data$Age)
data$Systolic = log(data$Systolic)
data$HDL = log(data$HDL)
data$Risk = log(data$Risk)
summary(data)

#model i statystyki dostateczne
OLS_results <- lm(Risk ~ Age + Systolic + HDL, data = data)
summary(OLS_results)

y <- as.matrix(data$Risk)
y
N.data <- length(y)
N.data
X <- cbind(as.matrix(rep(1, N.data)), 
           as.matrix(data[, c("Age","Systolic","HDL")]))
X
Beta.ols.data <- OLS_results$coefficients
Beta.ols.data
v.data <- OLS_results$df.residual
v.data
XTX.data <- t(X) %*% X
XTX.data
s2.data <- sum((OLS_results$residuals) ^ 2) / v.data
s2.data


#parametry a priori
#wspolczynniki regresji
Beta.prior <- c(-50, 0.8, 0.5, -0.7) #wartosci oczekiwane
#wariancja bledu
sm2.prior <- 4 #Wartość oczekiwana a priori odchylenia standardowego reszt (s): 0.5, zatem jego wariancji (s2) 0.25, a precyzji (s−2) 4
#macierz kowariancji
U.prior <- 0.4 * diag(4) # zerowe kowariancje i wariancja (co do warto?ci oczekiwanej) 0.1 * 4 * (1/4) = 0.1 (tzn. odchylenie standardowe (0.1)^0.5 = 0.31)
U.prior[1, 1] <- 50 #dla stalej wysoko bo nic o niej nie wiemy
v.prior <- 100 #Traktujemy informację a priori jako ekwiwalent modelu oszacowanego na 100 obserwacjach – podczas gdy nasza próba zawiera 30.
vs2.prior <- v.prior / sm2.prior
k <- ncol(X)
N <- length(y)


#parametry a posteriori
Beta.posterior <- solve(solve(U.prior) + XTX.data) %*% (solve(U.prior) %*% Beta.prior + XTX.data %*% Beta.ols.data)
U.posterior <- solve(solve(U.prior) + XTX.data)
v.posterior <- v.prior + N.data
vs2.posterior <- v.prior / sm2.prior + v.data * s2.data + t(Beta.ols.data-Beta.prior) %*% solve(U.prior + solve(XTX.data)) %*% (Beta.ols.data - Beta.prior)
sm2.posterior <- 1 / (vs2.posterior / v.posterior)

Beta.ols.data
Beta.prior
Beta.posterior

# gestosci
beta.space <- seq(from = -2, to = 4, by = 0.01)
n_eval_points <- length(beta.space)
n_parameters <- length(Beta.posterior)
prior.marg.dens.beta <- matrix(NA, nrow = n_parameters, ncol = n_eval_points)
posterior.marg.dens.beta <- matrix(NA, nrow = n_parameters, ncol = n_eval_points)
for(ii in 1:length(Beta.posterior)) {
  prior.marg.dens.beta[ii, ] <- apply(as.matrix(beta.space), 1, dmvt,
                                      delta = Beta.prior[ii], sigma = as.matrix(U.prior[ii, ii] / sm2.prior), df = v.prior, log = FALSE)
  
  posterior.marg.dens.beta[ii, ] <- apply(as.matrix(beta.space), 1, dmvt,
                                          delta = Beta.posterior[ii], sigma = as.matrix(U.posterior[ii, ii] / sm2.posterior), df = v.posterior, log = FALSE)
}

grey_area <- rgb(160, 160, 160, 80, names = NULL, maxColorValue = 255)
grey_line <- rgb(80, 80, 80, 160, names = NULL, maxColorValue = 255)
green_area <- rgb(24, 121, 104, 80, names = NULL, maxColorValue = 255)
green_line <- rgb(13, 85, 72, 160, names = NULL, maxColorValue = 255)

#wykresy
par(mfrow = c(2, 2))
for(ii in 1:length(Beta.posterior)) {
  title <- ifelse(ii==1, "intercept", colnames(data)[ii])
  plot(beta.space, prior.marg.dens.beta[ii, ], las = 1, lwd = 2, bty = "n", col = grey_area,
       ylim = c(0, max(c(max(prior.marg.dens.beta[ii, ]),max(posterior.marg.dens.beta[ii, ]))) + 1), type = "l", ylab = "gestosc", main = title)
  
  polygon(c(beta.space, rev(beta.space)), c(prior.marg.dens.beta[ii, ], 
                                            rep(0, length(beta.space))), col = grey_area, border = NA)
  abline(v = Beta.prior[ii], col = grey_line, lwd = 3)
  text(Beta.prior[ii], max(prior.marg.dens.beta[ii, ]) + 0.4, paste("E(beta) a priori = ", Beta.prior[ii]), col = grey_line)
  abline(v = Beta.ols.data[ii], col = rgb(0, 0, 0, 1), lwd = 3)
  text(Beta.prior[ii], max(posterior.marg.dens.beta[ii, ]) -0.5, paste("oszacowanie OLS = ", round(Beta.ols.data[ii], 4)), col = rgb(0, 0, 0, 1))
  lines(beta.space, posterior.marg.dens.beta[ii, ], lwd = 2, col = green_line)
  
  polygon(c(beta.space, rev(beta.space)), c(posterior.marg.dens.beta[ii, ], 
                                            rep(0, length(beta.space))), col = green_area, border = NA)
  abline(v = Beta.posterior[ii], col = green_line, lwd = 3)
  text(Beta.posterior[ii], max(posterior.marg.dens.beta[ii, ]) + 0.6, paste("E(beta) a posteriori = ", round(Beta.posterior[ii], digits = 4)), col = green_line)
}

#HPDI i czynniki Bayesa
ii <- 4

grey_area <- rgb(160, 160, 160, 80, names = NULL, maxColorValue = 255)
grey_line <- rgb(80, 80, 80, 160, names = NULL, maxColorValue = 255)
green_area <- rgb(24, 121, 104, 80, names = NULL, maxColorValue = 255)
green_line <- rgb(13, 85, 72, 160, names = NULL, maxColorValue = 255)
red_area <- rgb(255, 100, 123, 80, names = NULL, maxColorValue = 255)
red_line <- rgb(200, 0, 30, 160, names = NULL, maxColorValue = 255)

par(mfrow = c(1, 1))
manipulate( 
  {#zmienna binarna wskazuje, gdzie bedzie HPDI - tzn. o najwyzszej gestosci a posteriori ponad zadany poziom
    credible_set_indicator <- as.vector(as.integer(posterior.marg.dens.beta[ii, ] > line_level))
    credible_set_begin <- match(1, credible_set_indicator)
    credible_set_end <- length(credible_set_indicator) - match(1, rev(credible_set_indicator))
    #Lewy i prawy brzeg HPDI
    x1 <- beta.space[credible_set_begin]
    x2 <- beta.space[credible_set_end]
    #Na potrzeby wykresu tworzymy wektor, ktory przyjmuje wartosc gestosci a posteriori w HPDI i zero poza nim
    posterior.cs <- posterior.marg.dens.beta[ii, ] * credible_set_indicator
    #Poziom ufnosci
    HPDI_probab <- sum(posterior.cs) * 0.01
    #Wykres gestosci a posteriori
    plot(beta.space, posterior.marg.dens.beta[ii, ], las = 1, lwd = 2, bty = "n", col = green_line,
         ylim = c(0, max(posterior.marg.dens.beta[ii, ] + 1)), type = "l", ylab = "gestosc", main = colnames(data)[ii])
    polygon(c(beta.space, rev(beta.space)), 
            c(posterior.marg.dens.beta[ii, ], rep(0, length(beta.space))), 
            col = green_area, border = NA)
    text(Beta.posterior[ii], max(posterior.marg.dens.beta[ii, ]) + 0.6, paste("E(beta) a posteriori = ", round(Beta.posterior[ii], digits = 4)), col = green_line)
    abline(v = Beta.posterior[ii], col = green_line, lwd = 3)
    #Linia pozioma odcinajaca "najwyzsze" gestosci a posteriori (HPD)
    abline(h = line_level, col = red_line, lwd = 3)
    #Pole oznaczajace gestosc a posteriori w przedziale ufnosci HPDI
    polygon(c(beta.space, rev(beta.space)), 
            c(posterior.cs, rep(0, length(beta.space))), 
            col = red_area, border = NA)
    
    #Wyswietl poziom ufnosci i granice przedzialu
    text(0, max(posterior.marg.dens.beta[ii, ]) + 0.2, paste(round(HPDI_probab * 100, digits = 1), "% przedział HPDI: (", round(x1, digits = 2), " , ", round(x2, digits = 2), ")"), col = red_line)
  },
  line_level = slider(0, max(posterior.marg.dens.beta[ii, ]) + 0.002, step = 0.001, initial = max(posterior.marg.dens.beta[ii, ]) + 0.001)
)

#Wyznaczanie wiarygodnosci brzegowe 2 modeli: ze wszystkimi zmiennymi i bez zmiennej Age

## Model 1: ze wszystkimi 3 zmiennymi
P_y_M1 <- ((det(U.posterior) ^ 0.5) * gamma(v.posterior / 2) * ((vs2.posterior) ^ (- v.posterior / 2))) / ((pi ^ (N.data / 2)) * (det(U.prior) ^ 0.5) * gamma(v.prior / 2) * ((vs2.prior) ^ (- v.prior / 2)))

## Model 2: bez jednej zmiennej (3 modele)
P_y_M2 <- rep(NA, 3)

for (ii in 2:4) {
  X_2 <- X[, -c(ii)]
  eval(parse(text = paste("OLS_results_2 <- lm(Risk ~ ", paste(colnames(X_2), collapse = "+"), ", data = data)", sep = "")))
  Beta.ols.data_2 <- OLS_results_2$coefficients
  v.data_2 <- OLS_results_2$df.residual
  XTX.data_2 <- t(X_2) %*% X_2
  s2.data_2 <- sum((OLS_results_2$residuals) ^ 2) / v.data_2
  
  Beta.prior_2 <- Beta.prior[-c(ii)]
  U.prior_2 <- U.prior[- c(ii), -c(ii)]

  sm2.prior_2 <- sm2.prior
  v.prior_2 <- v.prior
  vs2.prior_2 <- vs2.prior
  
  Beta.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2) %*% (solve(U.prior_2) %*% Beta.prior_2 + XTX.data_2 %*% Beta.ols.data_2)
  U.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2)
  v.posterior_2 <- v.prior_2 + N.data
  vs2.posterior_2 <- v.prior_2 / sm2.prior_2 + v.data_2 * s2.data_2 + t(Beta.ols.data_2 - Beta.prior_2) %*% solve(U.prior_2 + solve(XTX.data_2)) %*% (Beta.ols.data_2 - Beta.prior_2)
  sm2.posterior_2 <- 1 / (vs2.posterior_2 / v.posterior_2)
  
  #gestosci brzegowe
  P_y_M2[ii - 1] <- ((det(U.posterior_2) ^ 0.5) * gamma(v.posterior_2 / 2) * ((vs2.posterior_2) ^ (- v.posterior_2 / 2))) / ((pi ^ (N.data / 2)) * (det(U.prior_2) ^ 0.5) * gamma(v.prior_2 / 2) * ((vs2.prior_2) ^ (- v.prior_2 / 2)))
}

#Czynniki Bayesa dla poszczegolnych zmiennych
BF_1_2 <- rep(P_y_M1, 3) / P_y_M2
BF_1_2_table <- data.frame(names(Beta.ols.data[2:4]), BF_1_2)
colnames(BF_1_2_table) <- c("zmienna", "czynnik Bayesa (analitycznie)")
(BF_1_2_table)