library(dplyr)
library(gamlss)
library(ggplot2)
library(fitdistrplus)
library(ismev)
library(evir)

### WCZYTANIE DANYCH ###

dane_meteo <- read.csv("C:/Users/magda/OneDrive/Pulpit/Zdarzenia_ekstremalne/Data_BB/dane_meteo.csv", encoding="UTF-8", sep="")
unique(dane_meteo$Nazwa)


### INFORMACJE O WYBRANEJ STACJI METEOROLOGICZNEJ ###

# Miejscowość: Chojnice
# Województwo: Pomorskie, Polska
# Szerokość geograficzna: 53.6978°N
# Długość geograficzna: 17.5584°E

# Chojnice znajdują się w północnej Polsce, niedaleko Borów Tucholskich.
# Jest to obszar o klimacie umiarkowanym z wyraźnymi sezonowymi różnicami temperatur i ciśnienia.

library(leaflet)
library(htmlwidgets)

# Tworzymy mapę i dodajemy marker w Chojnicach
leaflet() %>%
  addTiles() %>%  # Domyślna mapa OpenStreetMap
  setView(lng = 17.5584, lat = 53.6978, zoom = 10) %>%
  addMarkers(lng = 17.5584, lat = 53.6978, popup = "Stacja meteorologiczna Chojnice")

# Zapisujemy mapę do pliku HTML
# saveWidget(mapa, file = "mapa_chojnice.html", selfcontained = TRUE)


### PRZYGOTOWANIE DANYCH ###

# Wybieramy stację CHOJNICE
dane_chojnice <- dane_meteo %>% 
  filter(Nazwa == "CHOJNICE")

# Wybieramy sezon: lato (czerwiec, lipiec, sierpień)
dane_lato <- dane_chojnice %>%
  filter(Miesiac %in% c(6, 7, 8))

# Sprawdzamy, czy pomiary ciśnienia są istotne statystycznie, tzn. cisnienie_status = 8
length(dane_chojnice$cisnienie_status==8)
length(dane_lato$cisnienie_status==8)
length(dane_lato$cisnienie)
# Wszystkie pomiary ciśnienia latem mają status 8

# x <- dane_meteo[dane_meteo$Nazwa=="CHOJNICE", 'cisnienie_status']
# length(x==8)

# Obliczamy maksima roczne cisnienia (z okresu letniego)
dane_max <- aggregate(cisnienie ~ Rok, data = dane_lato, FUN = max) # 15 obserwacji
m_dane <- as.numeric(dane_max$cisnienie)

# Wybieramy kolumnę ciśnienia
dane <- dane_lato$cisnienie
length(dane)


### DOPASOWANIE NAJLEPSZEGO ROKŁADU GAMLSS ###

# Nasze dane to maksima godzinowe z każdego dnia lata (czerwiec, lipiec, sierpień) przez 15 lat (2006-2020)
hist(dane, prob=T, main="Rozkład ciśnienia", xlab=NA)

# Dopasowujemy rozkłady z `gamlss`
#?fitdist

t1 <- Sys.time()
fit <- fitDist(dane, type="realline")
t2 <- Sys.time()
t2-t1

# Wyświetlamy rozkłady posortowane według malejącej wartości AIC
fit$fits

# Rozkład o najmniejszej wartości AIC, z pełną nazwą rozkładu i skrótem tej nazwy
fit$family # "ST1" "Skew t (Azzalini type 1)"

# Rozkład Skew t (Azzalini Type 1) (oznaczany jako ST1) to rozszerzona wersja rozkładu t-Studenta,
# która dodatkowo uwzględnia skośność. ST1 może być asymetryczny.

# Nazwy  parametrów najlepszego rozkładu i ich wartości
fit$parameters # "mu"    "sigma" "nu"    "tau"

mu <- fit$mu # parametr lokalizacji
sigma <- fit$sigma # parametr skali (określa rozrzut danych)
nu <- fit$nu # parametr stopni swobody (kontroluje grubość ogonów)
tau <- fit$tau # parametr skośności (> 0 prawoskośny, = 0 symetryczny, < 0 leweoskośny)

mu; sigma; nu; tau

# Porównanie ST1 i t-Studenta
par(mfrow=c(1,1))
x <- seq(-5, 5, length.out=100)
y1 <- dST1(x, mu=0, sigma=1, nu=5, tau=3) # ST1
y2 <- dt(x, df=5) # t-Student

plot(x, y2, type="l", col="black", lwd=2, ylim=c(0, 0.7), ylab="gęstość", main="Porównanie rozkładu ST1 i t-Studenta")
lines(x, y1, col="red", lwd=2)
legend("topright", legend=c("t-Student (symetryczny)", "ST1 (prawoskośny)"),
       col=c("black", "red"), lty=1, lwd=2)


### WYKRESY DIAGNOSTYCZNE ###

par(mfrow=c(2,2))

# Histogram z gęstością teoretyczną
hist(dane, prob=TRUE, main="Rozkład ciśnienia", xlab="Ciśnienie")
curve(dST1(x, mu, sigma, nu, tau), add=T, col="red")

# QQ-plot (kawantyl-kwantyl)
alpha <- ppoints(100)
kwantyle_teo <- qST1(alpha, mu, sigma, nu, tau)
kwantyle_emp <- quantile(dane, alpha, na.rm=TRUE)

plot(kwantyle_emp, kwantyle_teo, main="QQ-plot")
abline(a=0, b=1, col="red")

# Dystrybuanta empiryczna vs teoretyczna
plot(ecdf(dane), main="Dystrybuanta empiryczna vs teoretyczna")
curve(pST1(x, mu, sigma, nu, tau), col="red", add=TRUE)


# -----
# Wykresy diagnostyczne z biblioteki `fitdistrplus` - nie działają;
# trzeba ponownie wyestymowac parametry rozkladu ST1

#X <- as.numeric(na.omit(dane))
#fST1 <- fitdist(X, "ST1", start=list(mu=mu, sigma=sigma, nu=nu, tau=tau))
#plot(fST1)

# -----


### OBLICZENIE POZIOMÓW ZWROTU x20 i x50 ###

# Dopasowujemy rozkład, obliczamy kwantyl x_20 przeliczamy z godzin na lata ->
# 1 na 20 lat = 20 * (30+31+31) * 24 = k

# Obliczenie k dla 20 i 50 lat
ny <- 30 + 31 + 31  # liczba dni w sezonie letnim
k20 <- 20 * ny * 24
k50 <- 50 * ny * 24

# Obliczenie poziomów zwrotu x20 i x50
x20.1 <- qST1(1 - 1/k20, mu, sigma, nu, tau)
x50.1 <- qST1(1 - 1/k50, mu, sigma, nu, tau)

x20.1; x50.1 # 1017.384; 1018.751

# Histogram z oznaczonymi poziomami zwrotu
par(mfrow=c(1,1))
hist(dane, prob=TRUE, main="Rozkład ciśnienia z poziomami zwrotu", xlab="Ciśnienie",
     col="lightgray", border="black", xlim=c(965,1020))
curve(dST1(x, mu, sigma, nu, tau), col="red", lwd=2, add=TRUE)
abline(v=x20.1, col="blue", lwd=2, lty=2)
abline(v=x50.1, col="purple", lwd=2, lty=2)
legend("topleft", legend=c("ST1 fit", "x20", "x50"),
       col=c("red", "blue", "purple"), lty=c(1,2,2), lwd=2)

#============  nowe dane, wyniki estymacji warto zapisac w pliku TwojaNazwa.Rdata
#save(maxday,fGT,file="C:/Users/.../Dane/TwojaNazwa.Rdata")

#ladowanie danych zapisanych w formacie RData
#dane <- lode(file="C:/Users/.../Dane/TwojaNazwa.Rdata")


### METODA MAKSIMÓW BLOKOWYCH (BMM) ###

# Histogram maksimów rocznych
hist(m_dane, prob=TRUE, 
     main="Rozkład maksimów rocznych (letnich) ciśnienia",
     xlab="Ciśnienie [hPa]", 
     border="black")

# Estymacja parametrów GEV w ismev
fit.m <- ismev::gev.fit(m_dane)
fit.m$mle
#         mu       sigma        ksi 
# 1007.425704    4.454460   -1.448937

# Wykresy diagnostyczne ismev
ismev::gev.diag(fit.m)

# Estymacja parametrów GEV w fExtremes
fit.f <- fExtremes::gevFit(m_dane)

# Wykresy i podsumowanie fExtremes
par(mfrow=c(2,2))
summary(fit.f)
#         xi      mu        beta 
# -1.110833 1007.347541    3.501854


### OBLICZENIE POZIOMÓW ZWROTU x20 i x50 ###

parametryGEV <- fit.m$mle

# qgev(p, xi = 1, mu = 0, sigma = 1) z biblioteki evir
x20.2 <- qgev(0.95, parametryGEV[3], parametryGEV[1], parametryGEV[2])[1]
x50.2 <- qgev(0.98, parametryGEV[3], parametryGEV[1], parametryGEV[2])[1]

# fExtremes, daje podobne wyniki do evir
#x20 <- qgev(0.95, -1.110833, 1007.347541, 3.501854)[1]
#x50 <- qgev(0.98, -1.110833, 1007.347541, 3.501854)[1]

x20.2; x50.2 # 1010.458, 1010.489

# Sprawdzamy przekroczenia poziomów zwrotu
przekroczenia_x20 <- dane_lato %>% filter(cisnienie > x20.2)
przekroczenia_x50 <- dane_lato %>% filter(cisnienie > x50.2)

# Liczba przekroczeń
liczba_x20 <- nrow(przekroczenia_x20)
liczba_x50 <- nrow(przekroczenia_x50)

cat("Liczba przekroczeń x20:", liczba_x20, "\n") # 4
cat("Liczba przekroczeń x50:", liczba_x50, "\n") # 4

# Lata i wielkości przekroczeń
przekroczenia_x20 <- przekroczenia_x20[, c("Rok", "Miesiac", "Dzien", "Godzina", "cisnienie")]
przekroczenia_x50 <- przekroczenia_x50[, c("Rok", "Miesiac", "Dzien", "Godzina", "cisnienie")]

print(przekroczenia_x20)
print(przekroczenia_x50)
# Mamy 4 przekroczenia o ciśnieniu 1010.5: 2.06.2011, 3.06.2011, 2 pomiary z 7.07.2013

# Histogram z oznaczonymi poziomami zwrotu
par(mfrow=c(1,1))
hist(dane, prob=TRUE, main="Rozkład ciśnienia z poziomami zwrotu", xlab="Ciśnienie",
     col="lightgray", border="black", xlim=c(965,1015))
abline(v=x20.2, col="red", lwd=2, lty=2)
abline(v=x50.2, col="blue", lwd=2, lty=2)
legend("topleft", legend=c("x20", "x50"),
       col=c("red", "blue"), lty=c(2,2), lwd=2)

# Okienko na ogon :D
par(mfrow=c(1,1))
hist(dane, prob=TRUE, main="Rozkład ciśnienia z poziomami zwrotu", xlab="Ciśnienie",
     col="lightgray", border="black", xlim=c(1005,1015))
abline(v=x20.2, col="red", lwd=2, lty=2)
abline(v=x50.2, col="blue", lwd=2, lty=2)
legend("topleft", legend=c("x20", "x50"),
       col=c("red", "blue"), lty=c(2,2), lwd=2)


### METODA PRZEKROCZEŃ PROGU (POT)###

# Dobór odpowiednio wysokiego progu
# PPatrzymy na wykres kwantylowy i oceniamy, czy dane dopasowują się do funkcji liniowej.
# Jeśli tak, to możemy spróbować zmniejszyć próg. (80% za mały próg, 95% wykres jest okej)

# Wybieramy próg na poziomie kwantyla 95%
u <- quantile(dane, 0.95)
u

# Wykresy rozrzutu z zaznaczonym progiem u 
par(mfrow=c(2,1))
plot(dane, type="h")
abline(h=u, lwd=2, col='red')

# Nadwyżki nad progiem u
Y <- dane[dane > u] - u
plot(Y, type='h')

# Liczba obserwacji przekraczających próg u = kwantyl 95%
przekroczenia_u <- sum(dane > u)
przekroczenia_u # 1624 obserwacji

### DOPASOWANIE ROZKŁADU GPD ###

# Estymujemy parametry rozkladu GPD
fitGPD <- ismev::gpd.fit(dane, u)

# Wyestymowane parametry rozkladu GPD
xi <- fitGPD$mle[[2]]
beta<- fitGPD$mle[[1]]
xi; beta # -0.535749, 4.010741

# Ocena dobroci dopasowania
par(mfrow=c(2,1))
hist(Y, prob=TRUE) # histogram nadwyżek
curve(evir::dgpd(x,xi,0,beta), col='red', lwd=2, add=T) # gęstość rozkładu GPD

qqplot(Y, evir::qgpd(ppoints(1000),xi,0,beta), main="QQ-plot")
abline(a=0,b=1,col=2)

# Wykresy dopasowania
ismev::gpd.diag(fitGPD)

### OBLICZENIE POZIOMÓW ZWROTU x20 i x50 ###

n <- length(dane)

# Ręcznie
x20.3 <- as.numeric(u + (beta/xi) * ((20 * (przekroczenia_u/n))^xi - 1))
x50.3 <- as.numeric(u + (beta/xi) * ((50 * (przekroczenia_u/n))^xi - 1))
x20.3; x50.3

# Z wbudowanej funkcji
fitGPD <- gpd(dane, u)
xi_est <- fitGPD$par.est[[1]]
beta_est <- fitGPD$par.est[[2]]
xi_est; beta_est
x_20 <- evir::riskmeasures(fitGPD,0.95)[2]
x_50 <- evir::riskmeasures(fitGPD,0.98)[2]
x_20; x_50


# Wyniki poziomów zwrotu z trzech modeli
cat("Model 1, x20:", x20.1, ", x50:", x50.1)
cat("Model 2, x20:", x20.2, ", x50:", x50.2)
cat("Model 3, x20:", x20.3, ", x50:", x50.3)
cat("Model 3 (evir), x20:", x_20, ", x50:", x_50)

