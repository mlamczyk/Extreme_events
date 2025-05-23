library(copula)
library(VineCopula)
library(kdecopula)
library(ggplot2)
library(ggExtra)
library(cowplot)
library(dplyr)
library(gamlss)
library(fitdistrplus)
library(MASS)


# ANALIZA II


# Zadanie 1

## Preprocessing

dane_meteo <- read.csv("C:/Users/magda/OneDrive/Pulpit/Zdarzenia_ekstremalne/Data_BB/dane_meteo.csv", encoding = "UTF-8", sep = "")

dane_chojnice <- dane_meteo %>%
  filter(Nazwa == "CHOJNICE")

dane_lato <- dane_chojnice %>%
  filter(Miesiac %in% c(6, 7, 8))

X1 <- dane_lato$temperatura
X2 <- dane_lato$cisnienie
X3 <- dane_lato$predkosc
length(X1)
length(X2)
length(X3)

# Wybrane zmienne: (X2 - ciśnienie, X3 - wiatry)

# Zmienna X3 jest dyskretna, dlatego
# generujemy próbę z rozkładu jednostajnego (szum)
x_unif <- runif(length(X3), -1 / 2, 1 / 2)
# Dodajemy do danych wiatrów
X3 <- X3 + x_unif

# Dopasowanie rozkładów brzegowych (z projektu Analiza I)
X123 <- data.frame(X1 = X1, X2 = X2, X3 = X3) # temperatura, ciśnienie, wiatry
X123 <- X123[complete.cases(X123), ]

Fits123 <- lapply(1:3, function(i) fitDist(X123[, i], type="realline")) # modele zdefiniowane na linii rzeczywistej
sapply(1:3, function(i) Fits123[[i]]$family)

# temperatura "SHASHo2" "Sinh-Arcsinh"
# ciśnienie "ST1" "Skew t (Azzalini type 1)"
# wiatry "JSU" "Johnson SU"

# Histogramy rozkładów brzegowych z gęstościami
png("histogramy-gestosc.png", width=600, height=500)
par(mfrow = c(3, 1))
hist(X123$X1, prob = TRUE, main = paste("Rozkład brzegowy temperatury"), xlab = "Temperatura [°C]")
curve(dSHASHo2(x,
               mu = Fits123[[1]]$mu,
               sigma = Fits123[[1]]$sigma,
               nu = Fits123[[1]]$nu,
               tau = Fits123[[1]]$tau
), col = 2, add = TRUE)
hist(X123$X2, prob = TRUE, main = paste("Rozkład brzegowy ciśnienia"), xlab = "Ciśnienie [hPa]")
curve(dST1(x,
            mu = Fits123[[2]]$mu,
            sigma = Fits123[[2]]$sigma,
            nu = Fits123[[2]]$nu,
            tau = Fits123[[2]]$tau
), col = 2, add = TRUE)
hist(X123$X3, prob = TRUE, main = paste("Rozkład brzegowy prędkości wiatru"), xlab = "Wiatr [m/s]")
curve(dJSU(x,
             mu = Fits123[[3]]$mu,
             sigma = Fits123[[3]]$sigma,
             nu = Fits123[[3]]$nu,
             tau = Fits123[[3]]$tau
), col = 2, add = TRUE)
dev.off()


## Analiza kopułowa zmiennych X1 i X2 (temperatury, ciśnienie)


## Analiza kopułowa zmeinnych X3 i X1 (wiatry, temperatura)


## Analiza kopułowa zmiennych X2 i X3 (ciśnienie, wiatry)

# Tworzymy dataframe dla zmiennych X2 i X3
X23 <- data.frame(X2 = X2, X3 = X3) # ciśnienie, wiatry
dim(X23)
X23 <- X23[complete.cases(X23), ]
dim(X23)

# Wykres rozrzutu z histogramami rozkładów brzegowych
png("rozrzut-histogramy.png", width=400, height=400)
p <- ggplot(X23, aes(X2, X3)) + # ciśnienie, wiatry
  geom_point()
ggMarginal(p, type = "histogram")
dev.off()


### Dobór kopuły - metoda  nieparametryczna

# Tworzymy pseudo-obserwacje nieparametryczne
U <- pobs(X23)
colnames(U) <- c("u2", "u3")
Y <- qnorm(U) # 'normalna' normalizacja
colnames(Y) <- c("y2", "y3")

# Wykresy rozrzutu
df <- data.frame(X23, U, Y)
head(df)

p1 <- ggplot(df, aes(X2, X3)) +
  geom_point()
p2 <- ggplot(df, aes(u2, u3)) +
  geom_point()
p3 <- ggplot(df, aes(y2, y3)) +
  geom_point()
p1.hist <- ggMarginal(p1, type = "histogram")
p2.hist <- ggMarginal(p2, type = "histogram")
p3.hist <- ggMarginal(p3, type = "histogram")

png("rozrzut1.png", width=600, height=500)
cowplot::plot_grid(p1.hist, p2.hist, p3.hist, ncol = 1, nrow = 3)
dev.off()


# Dobór kopuły
t1 <- Sys.time()
cop.npar <- BiCopSelect(U[, 1], U[, 2])
t2 <- Sys.time()
t2 - t1 # ok. 48 s

cop.npar # Bivariate copula: Rotated BB1 270 degrees (par = -0.28, par2 = -1.02, tau = -0.14)
cop.npar$family # 37

# Tau teoretyczny: -0.13
# Tau empiryczny:
cor(U[, 1], U[, 2], method = "kendall") # -0.145692

# Inne kopuły - porównanie, sortujemy wzgledem AIC, BIC
t1 <- Sys.time()
comp.npar <- BiCopEstList(U[, 1], U[, 2])
t2 <- Sys.time()
t2 - t1 # 1.18 min

# Pierwsze trzy ,,najlepsze kopuly''
AIC3.npar <- head(comp.npar$summary[order(comp.npar$summary$AIC), ], 3)
BIC3.npar <- head(comp.npar$summary[order(comp.npar$summary$BIC), ], 3)
logLik3.npar <- head(comp.npar$summary[order(comp.npar$summary$logLik, decreasing = TRUE), ], 3)
AIC3.npar # 37, 39, 33
BIC3.npar # 37, 39, 33
logLik3.npar # 37, 39, 33
# Dla każdego kryterium trzy najlepsze kopuły są takie same, w kryterium BIC druga i trzecia kopuła
# czasem różnią się kolejnością (inne ziarno).


### Dobór kopuły - metoda parametryczna

# Korzystamy z dopasowanych wcześniej rozkładów brzegowych
# ciśnienie "ST1" "Skew t (Azzalini type 1)"
# wiatry "JSU" "Johnson SU"

# Histogramy i QQ-ploty
png("histogramy-qq.png", width=400, height=400)
par(mfrow = c(2, 2))
hist(X23$X2, prob = TRUE, xlab = "Ciśnienie [hPa]", main = "Rozkład ciśnienia")
curve(dST1(x,
            mu = Fits123[[2]]$mu,
            sigma = Fits123[[2]]$sigma,
            nu = Fits123[[2]]$nu,
            tau = Fits123[[2]]$tau
), col = 2, add = TRUE)
hist(X23$X3, prob = TRUE, xlab = "Wiatr [m/s]", main = "Rozkład prędkości wiatru")
curve(dJSU(x,
             mu = Fits123[[3]]$mu,
             sigma = Fits123[[3]]$sigma,
             nu = Fits123[[3]]$nu,
             tau = Fits123[[3]]$tau
), col = 2, add = TRUE)

alpha <- ppoints(100)
X2emp <- quantile(X23$X2, alpha)
X2teo <- qST1(alpha,
               mu = Fits123[[2]]$mu,
               sigma = Fits123[[2]]$sigma,
               nu = Fits123[[2]]$nu,
               tau = Fits123[[2]]$tau
)
X3emp <- quantile(X23$X3, alpha)
X3teo <- qJSU(alpha,
                mu = Fits123[[3]]$mu,
                sigma = Fits123[[3]]$sigma,
                nu = Fits123[[3]]$nu,
                tau = Fits123[[3]]$tau
)

plot(X2emp, X2teo, xlab = "X2 empiryczne", ylab = "X2 teorotyczne")
abline(a = 0, b = 1, col = 2)
plot(X3emp, X3teo, xlab = "X3 empiryczne", ylab = "X3 teorotyczne")
abline(a = 0, b = 1, col = 2)
dev.off()

# Tworzymy pseudo-obserwacje parametryczne
V <- cbind(
  pST1(X23[, 1],
        mu = Fits123[[2]]$mu,
        sigma = Fits123[[2]]$sigma,
        nu = Fits123[[2]]$nu,
        tau = Fits123[[2]]$tau
  ),
  pJSU(X23[, 2],
         mu = Fits123[[3]]$mu,
         sigma = Fits123[[3]]$sigma,
         nu = Fits123[[3]]$nu,
         tau = Fits123[[3]]$tau
  )
)

colnames(V) <- c("v2", "v3")
# V zawiera NA, ale jeśli usuniemy je teraz to nie pójdzie dataframe do wykresów rozrzutu, bo mamy różną liczbę wierszy
# V <- V[complete.cases(V), ]
Y <- qnorm(V)
colnames(Y) <- c("y2", "y3")

# Wykresy rozrzutu
df <- data.frame(X23, V, Y)
head(df)

# p1 <- ggplot(df, aes(X1,X2)) + geom_point()
p4 <- ggplot(df, aes(v2, v3)) +
  geom_point()
p5 <- ggplot(df, aes(y2, y3)) +
  geom_point()
p4.hist <- ggMarginal(p4, type = "histogram")
p5.hist <- ggMarginal(p5, type = "histogram")

png("rozrzut2.png", width=600, height=500)
cowplot::plot_grid(p1.hist, p4.hist, p5.hist, ncol = 1, nrow = 3)
dev.off()

# V zawiera NA, usuwamy je
V <- V[complete.cases(V), ]
Y <- qnorm(V)
colnames(Y) <- c("y2", "y3")

# Dobór kopuły
t1 <- Sys.time()
cop.par <- BiCopSelect(V[, 1], V[, 2])
t2 <- Sys.time()
t2 - t1 # ok. 50 s

cop.par # Bivariate copula: Rotated BB1 270 degrees (par = -0.28, par2 = -1.02, tau = -0.14)
cop.par$family # 37

# Tau teoretyczny: -0.13
# Tau empiryczny:
cor(V[, 1], V[, 2], method = "kendall") # -0.145692

# Inne kopuły - porównanie, sortujemy względem AIC, BIC
t1 <- Sys.time()
comp.par <- BiCopEstList(V[, 1], V[, 2])
t2 <- Sys.time()
t2 - t1 # 1.17 min

# Pierwsze trzy ,,najlepsze kopuły''
AIC3.par <- head(comp.par$summary[order(comp.par$summary$AIC), ], 3)
BIC3.par <- head(comp.par$summary[order(comp.par$summary$BIC), ], 3)
logLik3.par <- head(comp.par$summary[order(comp.par$summary$logLik, decreasing = TRUE), ], 3)
AIC3.par # 37, 39, 33
BIC3.par # 37, 39, 33
logLik3.par # 37, 39, 33
# Dla każdego kryterium trzy najlepsze kopuły są takie same, jednak dla BIC kolejność kopuł czasem znacznie się różni
# (np. wyszło 33, 37, 39 dla jakiegoś innego ziarna).
# Wybrane kopuły są takie same jak w metodzie nieparametrycznej.

# Porównanie wykresów rozrzutu X, pseudo-obserwacji nieparametrycznych U i pseudo-obserwacji
# parametrycznych V
png("porownanie.png", width=800, height=400)
cowplot::plot_grid(p1.hist, p2.hist, p4.hist, ncol = 3, nrow = 1)
dev.off()


### Diagnostyka ###
# Porównanie gęstosci empirycznej i teoretycznej kopuły
p6 <- BiCopKDE(U[, 1], U[, 2], type = "surface")
p7 <- plot(cop.npar)
p8 <- plot(cop.par)

png("gestosci-kopula.png", width = 400, height = 800)
cowplot::plot_grid(p6, p7, p8, ncol = 1, nrow = 3)
dev.off()

# Porównanie konturu 'empirycznego' i teoretycznego (kopuły)
png("kontury-porownanie.png", width = 600, height = 400)
par(mfrow = c(2, 2))
# Kontur empiryczny
BiCopKDE(U[, 1], U[, 2], type = "contour", main = "Empiryczny kontur (KDE U)")
BiCopKDE(V[, 1], V[, 2], type = "contour", main = "Empiryczny kontur (KDE V)")
# Kontur teoretyczny
contour(cop.npar, main = "Teoretyczny kontur (kopuła npar)")
contour(cop.par, main = "Teoretyczny kontur (kopuła par)")
dev.off()

# Kontur 'empiryczny' inaczej
UU <- as.copuladata(U)
VV <- as.copuladata(V)
pairs(UU)
pairs(VV)


# Próby wygenerowane z kopuł
N <- dim(X23)[1]
N

sem.npar <- BiCopSim(N, cop.npar)
sem.par <- BiCopSim(N, cop.par)

sem <- data.frame(sem.npar = sem.npar, sem.par = sem.par)

p9 <- ggplot(sem, aes(sem.npar.1, sem.npar.2)) +
  geom_point()
p10 <- ggplot(sem, aes(sem.par.1, sem.par.2)) +
  geom_point()
p9.hist <- ggMarginal(p9, type = "histogram")
p10.hist <- ggMarginal(p10, type = "histogram")

cowplot::plot_grid(p9.hist, p10.hist, ncol = 1, nrow = 2)


# Próby z rozkładów F=C(F1,F2)
u2 <- sem.npar[, 1]
u3 <- sem.npar[, 2]
v2 <- sem.par[, 1]
v3 <- sem.par[, 2]

x2.npar <- quantile(X23$X2, u2)
x3.npar <- quantile(X23$X3, u3)
x2.par <- qST1(v2,
                mu = Fits123[[2]]$mu,
                sigma = Fits123[[2]]$sigma,
                nu = Fits123[[2]]$nu,
                tau = Fits123[[2]]$tau
)
x3.par <- qJSU(v3,
                 mu = Fits123[[3]]$mu,
                 sigma = Fits123[[3]]$sigma,
                 nu = Fits123[[3]]$nu,
                 tau = Fits123[[3]]$tau
)

sem.F <- data.frame(
  x2.npar = x2.npar, x3.npar = x3.npar,
  x2.par = x2.par, x3.par = x3.par
)
head(sem.F)

p11 <- ggplot(sem.F, aes(x2.npar, x3.npar)) +
  geom_point()
p12 <- ggplot(sem.F, aes(x2.par, x3.par)) +
  geom_point()
p11.hist <- ggMarginal(p11, type = "histogram")
p12.hist <- ggMarginal(p12, type = "histogram")

cowplot::plot_grid(p11.hist, p12.hist, ncol = 1, nrow = 2)


# Porównanie współczynnikow ekstremalnych empirycznych  i teoretycznych
# Estymacja dolnego i górnego współczynnika

# Empiryczne
p <- 0.01 # cut-off
(lam.C <- c(
  lower = fitLambda(U, p = p)[2, 1],
  upper = fitLambda(U, p = p, lower.tail = FALSE)[2, 1]
)) # lower 0, upper 0
(lam.C <- c(
  lower = fitLambda(V, p = p)[2, 1],
  upper = fitLambda(V, p = p, lower.tail = FALSE)[2, 1]
)) # lower 0, upper 0

# Współczynniki dla wyestymowanej kopuły (teoretyczne)
BiCopPar2TailDep(cop.npar) # lower 0, upper 0
BiCopPar2TailDep(cop.par) # lower 0, upper 0


# Współczynniki zależnosci ekstremalnych dla pierwszych
# trzech kopuł wybranych przez AIC
nAIC3.npar <- as.numeric(rownames(AIC3.npar))
lu.coeff.AIC3 <- t(sapply(nAIC3.npar, function(i) BiCopPar2TailDep(comp.npar$models[[i]])))
rownames(lu.coeff.AIC3) <- BiCopName(AIC3.npar[, 1])
lu.coeff.AIC3
#         lower upper
# BB1_270 0     0
# BB7_270 0     0
# C270    0     0

# 37: rotated BB1 270 degrees
# 39: rotated BB7 270 degrees
# 33: rotated Clayton 270 degrees


# Zadanie 2

# Metoda parametryczna
# Pseudo-obserwacje parametryczne
V123 <- cbind(
  pSHASHo2(X123[, 1],
           mu = Fits123[[1]]$mu,
           sigma = Fits123[[1]]$sigma,
           nu = Fits123[[1]]$nu,
           tau = Fits123[[1]]$tau
  ),
  pST1(X123[, 2],
        mu = Fits123[[2]]$mu,
        sigma = Fits123[[2]]$sigma,
        nu = Fits123[[2]]$nu,
        tau = Fits123[[2]]$tau
  ),
  pJSU(X123[, 3],
         mu = Fits123[[3]]$mu,
         sigma = Fits123[[3]]$sigma,
         nu = Fits123[[3]]$nu,
         tau = Fits123[[3]]$tau
  )
)

colnames(V123) <- c("v1", "v2", "v3")
# V123 zawiera NA, usuwamy je
V123 <- V123[complete.cases(V123), ]
Y123 <- qnorm(V123)
colnames(Y123) <- c("y1", "y2", "y3")

# Dobór kopuł
cop.par12 <- BiCopSelect(V123[, 1], V123[, 2]) # temperatura, ciśnienie
cop.par13 <- BiCopSelect(V123[, 1], V123[, 3]) # temperatura, wiatry
cop.par12 # Rotated Tawn type 1 180 degrees (par = 1.29, par2 = 0.06, tau = 0.03)
cop.par13 # Survival BB8 (par = 1.75, par2 = 0.75, tau = 0.13)

# Metoda nieparametryczna (pobs)
# Tworzymy pseudo-obserwacje nieparametryczne
U123 <- pobs(X123)
colnames(U123) <- c("u1", "u2", "u3")
Y123 <- qnorm(U123) # 'normalna' normalizacja
colnames(Y123) <- c("y1", "y2", "y3")

# Dobór kopuł
cop.npar12 <- BiCopSelect(U123[, 1], U123[, 2]) # temperatura, ciśnienie
cop.npar13 <- BiCopSelect(U123[, 1], U123[, 3]) # temperatura, wiatry
cop.npar12 # Rotated Tawn type 1 180 degrees (par = 1.31, par2 = 0.06, tau = 0.04)
cop.npar13 # Survival BB8 (par = 1.76, par2 = 0.74, tau = 0.13)

### Rozkłady warunkowe ###

# Rozkład warunkowy temperatur (X1), przy warunku, że zmienna X2 (odpowiednio X3)
# ma wartość na poziomie 20-letniego i odpowiednio 50-letniego poziomu zwotu

N <- length(X1) # liczność danych

# Poziomy zwrotu (metoda pierwsza z Analiza I) dla X2 (ciśnienie) i X3 (wiatry)
ny <- 30 + 31 + 31 # liczba dni w sezonie letnim
k20 <- 20 * ny * 24
k50 <- 50 * ny * 24

# Obliczenie poziomów zwrotu r20 i r50 (parametryczne?)
r20.X2 <- qST1(1 - 1 / k20,
                mu = Fits123[[2]]$mu,
                sigma = Fits123[[2]]$sigma,
                nu = Fits123[[2]]$nu,
                tau = Fits123[[2]]$tau
)
r50.X2 <- qST1(1 - 1 / k50,
                mu = Fits123[[2]]$mu,
                sigma = Fits123[[2]]$sigma,
                nu = Fits123[[2]]$nu,
                tau = Fits123[[2]]$tau
)

r20.X3 <- qJSU(1 - 1 / k20,
                 mu = Fits123[[3]]$mu,
                 sigma = Fits123[[3]]$sigma,
                 nu = Fits123[[3]]$nu,
                 tau = Fits123[[3]]$tau
)
r50.X3 <- qJSU(1 - 1 / k50,
                 mu = Fits123[[3]]$mu,
                 sigma = Fits123[[3]]$sigma,
                 nu = Fits123[[3]]$nu,
                 tau = Fits123[[3]]$tau
)

# ciśnienie
r20.X2 # 1017.222
r50.X2 # 1018.563
# wiatry
r20.X3 # 13.61127
r50.X3 # 14.37286


### Temperatura przy założeniu X2 = r20.X2 (ciśnienie) ###

# Kroki:
# 1. Wyznaczamy warunek $\{U_2 = u_2\}$ z $F_2(x_2)=u_2$.
# 2. Generujemy $N$-elementową próbę $(u_i)_{i=1}^N$ z rozkładu jednostajnego $U(0,1)$.
# 3. Wyznaczamy wartości funkcji odwrotnej do $h_{1|2}(\cdot \mid \nu)$ dla wygenerowanej próby:
# $$w_i = h_{1|2}^{-1}(u_i \mid \nu), \quad i = 1, \dots, N.$$
# 4. Składamy z funkcją kwantylową rozkładu $F_1$:
# $$x_i = F_1^{-1}(w_i), \quad i = 1, \dots, N.$$



# Model nieparametryczny
# --- 1.
# Wykorzystujemy dystrybuantę empiryczną ciśnienia
F2.npar <- ecdf(X123$X2) # empirical cumulative distribution function
u2_r20_npar <- F2.npar(r20.X2) # jaki to kwantyl?
u2_r20_npar # 1

# --- 2.
U2_r20_cond_npar <- rep(u2_r20_npar, N) # warunek powtórzony N razy
# Losujemy U1 ~ U(0,1) do wykorzystania w kopule odwrotnej
set.seed(123)
U1 <- runif(N)

# --- 3.
# H-inv: z kopuły nieparametrycznej
# odwrotna do C(u1|u2), czyli do h1|2
w1_cond_r20_X2_npar <- BiCopHinv2(U1, U2_r20_cond_npar, cop.npar12)

# --- 4.
# Transformacja odwrotna z dystrybuanty empirycznej temperatury
# (odwrotna do F1)
F1.npar <- ecdf(X123$X1)
x1_cond_r20_X2_npar <- quantile(X123$X1, probs = w1_cond_r20_X2_npar)



# Model parametryczny
# --- 1.
# Dystrybuanta parametryczna ciśnienia
u2_r20_par <- pST1(r20.X2,
                    mu = Fits123[[2]]$mu,
                    sigma = Fits123[[2]]$sigma,
                    nu = Fits123[[2]]$nu,
                    tau = Fits123[[2]]$tau
)
u2_r20_par # 0.9999774

# --- 2.
U2_r20_cond_par <- rep(u2_r20_par, N)
# U1 tak samo jak w modelu npar
set.seed(123)
U1 <- runif(N)

# --- 3.
# H-inv: z kopuły parametrycznej
# odwrotna do C(u1|u2), czyli do h1|2
w1_cond_r20_X2_par <- BiCopHinv2(U1, U2_r20_cond_par, cop.par12)

# --- 4.
# Odwrotność dystrybuanty SHASHo2 temperatury
x1_cond_r20_X2_par <- qSHASHo2(w1_cond_r20_X2_par,
                               mu = Fits123[[1]]$mu,
                               sigma = Fits123[[1]]$sigma,
                               nu = Fits123[[1]]$nu,
                               tau = Fits123[[1]]$tau
)


# Dane i próby warunkowe na histogramach
png("warunkowe-x2-20.png", width=600, height=500)
par(mfrow = c(3, 1))
hist(X123$X1,
     prob = TRUE, main = "Temperatura latem (dane empiryczne)",
     xlab = "Temperatura [°C]"
)
curve(dSHASHo2(x,
               mu = Fits123[[1]]$mu,
               sigma = Fits123[[1]]$sigma,
               nu = Fits123[[1]]$nu,
               tau = Fits123[[1]]$tau
), col = 2, add = TRUE)
hist(x1_cond_r20_X2_npar,
     prob = TRUE, main = "Temperatura | X2=r20 (nieparam.)",
     xlab = "Temperatura [°C]"
)
hist(x1_cond_r20_X2_par,
     prob = TRUE, main = "Temperatura | X2=r20 (param.)",
     xlab = "Temperatura [°C]"
)
dev.off()

# Prognoza jako średnia z rozkładu warunkowego
mean(x1_cond_r20_X2_npar)
mean(x1_cond_r20_X2_par) 
# średnie: npar 17.76016; par 17.75569

# Przedziały 5% 95%
quantile(x1_cond_r20_X2_npar, c(0.05, 0.95)) # 10.5 26.6
quantile(x1_cond_r20_X2_par, c(0.05, 0.95)) # 10.47292 26.43860



### Temperatura przy założeniu X2 = r50.X2 (ciśnienie) ###

# --- 1.
u2_r50_npar <- F2.npar(r50.X2)
u2_r50_npar # 1
u2_r50_npar.par <- pST1(r50.X2,
                         mu = Fits123[[2]]$mu,
                         sigma = Fits123[[2]]$sigma,
                         nu = Fits123[[2]]$nu,
                         tau = Fits123[[2]]$tau
)
u2_r50_npar.par # 0.9999909

# --- 2.
U2_r50_cond_npar <- rep(u2_r50_npar, N)
U2_r50_cond_par <- rep(u2_r50_npar.par, N)

# --- 3.
w1_cond_r50_X2_npar <- BiCopHinv2(U1, U2_r50_cond_npar, cop.npar12)
w1_cond_r50_X2_par <- BiCopHinv2(U1, U2_r50_cond_par, cop.par12)

# --- 4.
x1_cond_r50_X2_npar <- quantile(X123$X1, probs = w1_cond_r50_X2_npar)
x1_cond_r50_X2_par <- qSHASHo2(w1_cond_r50_X2_par,
                               mu = Fits123[[1]]$mu,
                               sigma = Fits123[[1]]$sigma,
                               nu = Fits123[[1]]$nu,
                               tau = Fits123[[1]]$tau
)


# Dane i próby warunkowe na histogramach
png("warunkowe-x2-50.png", width=600, height=500)
par(mfrow = c(3, 1))
hist(X123$X1,
     prob = TRUE, main = "Temperatura latem (dane empiryczne)",
     xlab = "Temperatura [°C]"
)
curve(dSHASHo2(x,
               mu = Fits123[[1]]$mu,
               sigma = Fits123[[1]]$sigma,
               nu = Fits123[[1]]$nu,
               tau = Fits123[[1]]$tau
), col = 2, add = TRUE)
hist(x1_cond_r50_X2_npar,
     prob = TRUE, main = "Temperatura | X2=r50 (nieparam.)",
     xlab = "Temperatura [°C]"
)
hist(x1_cond_r50_X2_par,
     prob = TRUE, main = "Temperatura | X2=r50 (param.)",
     xlab = "Temperatura [°C]"
)
dev.off()

# Prognoza jako średnia z rozkładu warunkowego
mean(x1_cond_r50_X2_npar)
mean(x1_cond_r50_X2_par)
# średnie: npar 17.76016; par 17.75724

# Przedziały 5% 95%
quantile(x1_cond_r50_X2_npar, c(0.05, 0.95)) # 10.5 26.6
quantile(x1_cond_r50_X2_par, c(0.05, 0.95)) # 10.47317 26.44238



### Temperatura przy założeniu X3 = r20.X3 (wiatry) ###

# --- 1.
F3.npar <- ecdf(X123$X3)
u3_r20_npar <- F3.npar(r20.X3)
u3_r20_npar # 1
u3_r20_npar.par <- pJSU(r20.X3,
                          mu = Fits123[[3]]$mu,
                          sigma = Fits123[[3]]$sigma,
                          nu = Fits123[[3]]$nu,
                          tau = Fits123[[3]]$tau
)
u3_r20_npar.par # 0.9999774

# --- 2.
U3_r20_cond_npar <- rep(u3_r20_npar, N)
U3_r20_cond_par <- rep(u3_r20_npar.par, N)

# --- 3.
w1_cond_r20_X3_npar <- BiCopHinv2(U1, U3_r20_cond_npar, cop.npar13)
w1_cond_r20_X3_par <- BiCopHinv2(U1, U3_r20_cond_par, cop.par13)

# --- 4.
x1_cond_r20_X3_npar <- quantile(X123$X1, probs = w1_cond_r20_X3_npar)
x1_cond_r20_X3_par <- qSHASHo2(w1_cond_r20_X3_par,
                               mu = Fits123[[1]]$mu,
                               sigma = Fits123[[1]]$sigma,
                               nu = Fits123[[1]]$nu,
                               tau = Fits123[[1]]$tau
)

# Dane i próby warunkowe na histogramach
png("warunkowe-x3-20.png", width=600, height=500)
par(mfrow = c(3, 1))
hist(X123$X1,
     prob = TRUE, main = "Temperatura latem (dane empiryczne)",
     xlab = "Temperatura [°C]"
)
curve(dSHASHo2(x,
               mu = Fits123[[1]]$mu,
               sigma = Fits123[[1]]$sigma,
               nu = Fits123[[1]]$nu,
               tau = Fits123[[1]]$tau
), col = 2, add = TRUE)
hist(x1_cond_r20_X3_npar,
     prob = TRUE, main = "Temperatura | X3=r20 (nieparam.)",
     xlab = "Temperatura [°C]"
)
hist(x1_cond_r20_X3_par,
     prob = TRUE, main = "Temperatura | X3=r20 (param.)",
     xlab = "Temperatura [°C]"
)
dev.off()

# Prognoza jako średnia z rozkładu warunkowego
mean(x1_cond_r20_X3_npar)
mean(x1_cond_r20_X3_par)
# średnie: npar 18.7197; par 18.71697

# Przedziały 5% 95%
quantile(x1_cond_r20_X3_npar, c(0.05, 0.95)) # 11.5 27.0
quantile(x1_cond_r20_X3_par, c(0.05, 0.95)) # 11.43369 26.96084



### Temperatura przy założeniu X3 = r50.X3 (wiatry) ###

# --- 1.
u3_r50_npar <- F3.npar(r50.X3)
u3_r50_npar # 1
u3_r50_npar.par <- pJSU(r50.X3,
                          mu = Fits123[[3]]$mu,
                          sigma = Fits123[[3]]$sigma,
                          nu = Fits123[[3]]$nu,
                          tau = Fits123[[3]]$tau
)

# --- 2.
U3_r50_cond_npar <- rep(u3_r50_npar, N)
U3_r50_cond_par <- rep(u3_r50_npar.par, N)

# --- 3.
w1_cond_r50_X3_npar <- BiCopHinv2(U1, U3_r50_cond_npar, cop.npar13)
w1_cond_r50_X3_par <- BiCopHinv2(U1, U3_r50_cond_par, cop.par13)

# --- 4.
x1_cond_r50_X3_npar <- quantile(X123$X1, probs = w1_cond_r50_X3_npar)
x1_cond_r50_X3_par <- qSHASHo2(w1_cond_r50_X3_par,
                               mu = Fits123[[1]]$mu,
                               sigma = Fits123[[1]]$sigma,
                               nu = Fits123[[1]]$nu,
                               tau = Fits123[[1]]$tau
)

# Dane i próby warunkowe na histogramach
png("warunkowe-x3-50.png", width=600, height=500)
par(mfrow = c(3, 1))
hist(X123$X1,
     prob = TRUE, main = "Temperatura latem (dane empiryczne)",
     xlab = "Temperatura [°C]"
)
curve(dSHASHo2(x,
               mu = Fits123[[1]]$mu,
               sigma = Fits123[[1]]$sigma,
               nu = Fits123[[1]]$nu,
               tau = Fits123[[1]]$tau
), col = 2, add = TRUE)
hist(x1_cond_r50_X3_npar,
     prob = TRUE, main = "Temperatura | X3=r50 (nieparam.)",
     xlab = "Temperatura [°C]"
)
hist(x1_cond_r50_X3_par,
     prob = TRUE, main = "Temperatura | X3=r50 (param.)",
     xlab = "Temperatura [°C]"
)
dev.off()

# Prognoza jako średnia z rozkładu warunkowego
mean(x1_cond_r50_X3_npar)
mean(x1_cond_r50_X3_par)
# średnie: npar 18.7197; par 18.717

# Przedziały 5% 95%
quantile(x1_cond_r50_X3_npar, c(0.05, 0.95)) # 11.5 27.0
quantile(x1_cond_r50_X3_par, c(0.05, 0.95)) # 11.43370 26.96086



# Wyniki zapisujemy do tabel
tabela_poziomy_zwrotu <- data.frame(
  r20 = c(r20.X2, r20.X3),
  r50 = c(r50.X2, r50.X3),
  row.names = c("Ciśnienie (X2)", "Wiatr (X3)")
)
tabela_poziomy_zwrotu

tabela_srednie_temp <- data.frame(
  npar = c(
    mean(x1_cond_r20_X2_npar),
    mean(x1_cond_r50_X2_npar),
    mean(x1_cond_r20_X3_npar),
    mean(x1_cond_r50_X3_npar)
  ),
  par = c(
    mean(x1_cond_r20_X2_par),
    mean(x1_cond_r50_X2_par),
    mean(x1_cond_r20_X3_par),
    mean(x1_cond_r50_X3_par)
  ),
  row.names = c(
    "X2 = r20", "X2 = r50",
    "X3 = r20", "X3 = r50"
  )
)
tabela_srednie_temp

# Średnia temperatura
mean(X1) # 17.54933

tabela_przedzial_npar <- data.frame(
  dolny = c(
    quantile(x1_cond_r20_X2_npar, 0.05),
    quantile(x1_cond_r50_X2_npar, 0.05),
    quantile(x1_cond_r20_X3_npar, 0.05),
    quantile(x1_cond_r50_X3_npar, 0.05)
  ),
  gorny = c(
    quantile(x1_cond_r20_X2_npar, 0.95),
    quantile(x1_cond_r50_X2_npar, 0.95),
    quantile(x1_cond_r20_X3_npar, 0.95),
    quantile(x1_cond_r50_X3_npar, 0.95)
  ),
  row.names = c(
    "X2 = r20", "X2 = r50",
    "X3 = r20", "X3 = r50"
  )
)
tabela_przedzial_npar

tabela_przedzial_par <- data.frame(
  dolny = c(
    quantile(x1_cond_r20_X2_par, 0.05),
    quantile(x1_cond_r50_X2_par, 0.05),
    quantile(x1_cond_r20_X3_par, 0.05),
    quantile(x1_cond_r50_X3_par, 0.05)
  ),
  gorny = c(
    quantile(x1_cond_r20_X2_par, 0.95),
    quantile(x1_cond_r50_X2_par, 0.95),
    quantile(x1_cond_r20_X3_par, 0.95),
    quantile(x1_cond_r50_X3_par, 0.95)
  ),
  row.names = c(
    "X2 = r20", "X2 = r50",
    "X3 = r20", "X3 = r50"
  )
)
tabela_przedzial_par

# Zapisywanie wyników w pliku Wyniki_2.RData
save(
  tabela_poziomy_zwrotu,
  tabela_srednie_temp,
  tabela_przedzial_npar,
  tabela_przedzial_par,
  file = "Wyniki_2.RData"
)

# load("Wyniki_2.RData")