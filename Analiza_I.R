library(dplyr)
library(gamlss)
library(ggplot2)
library(fitdistrplus)
library(ismev)
library(evir)
library(leaflet)
library(htmlwidgets)
library(webshot2)

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


# Tworzymy mapę i dodajemy marker w Chojnicach
mymap <- leaflet() %>%
  addTiles() %>%  # Domyślna mapa OpenStreetMap
  setView(lng = 19.5, lat = 52.3, zoom = 6) %>%
  addMarkers(lng = 17.5584, lat = 53.6978, popup = "Stacja meteorologiczna Chojnice")

mymap

# Zapisujemy mapę jako plik HTML
#saveWidget(mymap, "mapa.html", selfcontained = TRUE)

# Zapisujemy zrzut ekranu do pliku PNG
#webshot("mapa.html", file = "mapa.png", vwidth = 800, vheight = 600)


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
sum(is.na(dane))
# nie ma braków w pomiarach

# Podstawowe statystyki
statystyki <- summary(dane)
odchylenie_standardowe <- sd(dane)
tabela_statystyk <- data.frame(
  Statystyka = c("Wartość minimalna", "Mediana", "Średnia", "Odchylenie standardowe", "Wartość maksymalna"),
  Wartość = c(statystyki["Min."],statystyki["Median"],statystyki["Mean"],odchylenie_standardowe,statystyki["Max."])
)
tabela_statystyk


### DOPASOWANIE NAJLEPSZEGO ROKŁADU GAMLSS ###

# Nasze dane to maksima godzinowe z każdego dnia lata (czerwiec, lipiec, sierpień) przez 15 lat (2006-2020)
# Histogram rozkładu pomiarów ciśnienia
#png("histogram-rozrzut.png", width=800, height=300)
par(mfrow=c(1,2))
hist(dane, prob=T, main="(a) Rozkład ciśnienia", xlab="Ciśnienie [hPa]")

# Wykres rozrzutu maksimów roczynch ciśnienia
plot(dane_max$Rok, dane_max$cisnienie, main="(b) Maksima roczne ciśnienia latem (2006-2020)",
     xlab="Rok", ylab="Maksymalne ciśnienie [hPa]", pch=19)
#dev.off()

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
tau <- fit$tau # parametr skośności (> 0 prawoskośny, = 0 symetryczny, < 0 lewoskośny)

mu; sigma; nu; tau

tabela_st1 <- data.frame(
  mu = mu,
  sigma = sigma,
  nu = nu,
  tau = tau
)
tabela_st1

### WYKRESY DIAGNOSTYCZNE ###

# plotdist(dane, histo = TRUE, demp = TRUE)

#png("diagnostyczne-st1.png", width=800, height=500)
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
#dev.off()


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

# Sprawdzamy przekroczenia poziomów zwrotu
przekroczenia_x20.1 <- dane_lato %>% filter(cisnienie > x20.1)
przekroczenia_x50.1 <- dane_lato %>% filter(cisnienie > x50.1)

# Liczba przekroczeń
liczba_x20.1 <- nrow(przekroczenia_x20.1)
liczba_x50.1 <- nrow(przekroczenia_x50.1)

cat("Liczba przekroczeń x20:", liczba_x20.1, "\n") # 0
cat("Liczba przekroczeń x50:", liczba_x50.1, "\n") # 0

# Histogram z oznaczonymi poziomami zwrotu
#png("histogram-x.png", width=800, height=400)
par(mfrow=c(1,1))
hist(dane, prob=TRUE, main="Rozkład ciśnienia z poziomami zwrotu", xlab="Ciśnienie",
     col="lightgray", border="black", xlim=c(965,1020))
curve(dST1(x, mu, sigma, nu, tau), col="red", lwd=2, add=TRUE)
abline(v=x20.1, col="blue", lwd=2, lty=2)
abline(v=x50.1, col="purple", lwd=2, lty=2)
legend("topleft", legend=c("ST1 fit", "x20", "x50"),
       col=c("red", "blue", "purple"), lty=c(1,2,2), lwd=2)
#dev.off()


### METODA MAKSIMÓW BLOKOWYCH (BMM) ###

# Estymacja parametrów GEV w ismev
fit.m <- ismev::gev.fit(m_dane)
fit.m$mle
#         mu       sigma        ksi 
# 1007.425704    4.454460   -1.448937

parametryGEV <- fit.m$mle

tabela_gev <- data.frame(
  mu = parametryGEV[1],
  sigma = parametryGEV[2],
  xi = parametryGEV[3]
)
tabela_gev

# Wykresy diagnostyczne ismev
#png("diagnostyczne-gev.png", width=800, height=500)
ismev::gev.diag(fit.m)
#dev.off()

# Estymacja parametrów GEV w fExtremes
#fit.f <- fExtremes::gevFit(m_dane)

# Wykresy i podsumowanie fExtremes
#par(mfrow=c(2,2))
#summary(fit.f)
#         xi      mu        beta 
# -1.110833 1007.347541    3.501854


### OBLICZENIE POZIOMÓW ZWROTU x20 i x50 ###

# qgev(p, xi = 1, mu = 0, sigma = 1) z biblioteki evir
x20.2 <- qgev(0.95, parametryGEV[3], parametryGEV[1], parametryGEV[2])[1]
x50.2 <- qgev(0.98, parametryGEV[3], parametryGEV[1], parametryGEV[2])[1]

# fExtremes, daje podobne wyniki do evir
#x20 <- qgev(0.95, -1.110833, 1007.347541, 3.501854)[1]
#x50 <- qgev(0.98, -1.110833, 1007.347541, 3.501854)[1]

x20.2; x50.2 # 1010.458, 1010.489

# Sprawdzamy przekroczenia poziomów zwrotu
przekroczenia_x20.2 <- dane_lato %>% filter(cisnienie > x20.2)
przekroczenia_x50.2 <- dane_lato %>% filter(cisnienie > x50.2)

# Liczba przekroczeń
liczba_x20.2 <- nrow(przekroczenia_x20.2)
liczba_x50.2 <- nrow(przekroczenia_x50.2)

cat("Liczba przekroczeń x20:", liczba_x20.2, "\n") # 4
cat("Liczba przekroczeń x50:", liczba_x50.2, "\n") # 4

# Lata i wielkości przekroczeń
przekroczenia_x20.2 <- przekroczenia_x20.2[, c("Rok", "Miesiac", "Dzien", "Godzina", "cisnienie")]
przekroczenia_x50.2 <- przekroczenia_x50.2[, c("Rok", "Miesiac", "Dzien", "Godzina", "cisnienie")]

print(przekroczenia_x20.2)
print(przekroczenia_x50.2)
# Mamy 4 przekroczenia o ciśnieniu 1010.5: 2.06.2011, 3.06.2011, 2 pomiary z 7.07.2013
print(identical(przekroczenia_x20.2,przekroczenia_x50.2))
print(all.equal(przekroczenia_x20.2,przekroczenia_x50.2))
# Przekroczenia dla obu poziomów zwrotu są takie same

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


### METODA PRZEKROCZEŃ PROGU (POT) ###

# Dobór odpowiednio wysokiego progu
# Patrzymy na wykres kwantylowy i oceniamy, czy dane dopasowują się do funkcji liniowej.
# Jeśli tak, to możemy spróbować zmniejszyć próg. (80% za mały próg, 95% wykres jest okej)

# Wybieramy próg na poziomie kwantyla 95%
u <- quantile(dane, 0.95)
u

# Wykresy rozrzutu z zaznaczonym progiem u
#png("rozrzut-nadwyzki.png", width=800, height=500)
par(mfrow=c(2,1))
plot(dane, type="h", main="(a) Wykres rozrzutu z zaznaczonym progiem u")
abline(h=u, lwd=2, col='red')

# Nadwyżki nad progiem u
Y <- dane[dane > u] - u
plot(Y, type='h', main="(b) Nadwyżki nad progiem u")
#dev.off()

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

tabela_gpd <- data.frame(
  xi = xi,
  beta = beta
)
tabela_gpd

# Ocena dobroci dopasowania
par(mfrow=c(2,1))
hist(Y, prob=TRUE, main="Histogram nadwyżek")
curve(evir::dgpd(x,xi,0,beta), col='red', lwd=2, add=T) # gęstość rozkładu GPD

qqplot(Y, evir::qgpd(ppoints(1000),xi,0,beta), main="QQ-plot")
abline(a=0,b=1,col=2)

# Wykresy dopasowania
#png("diagnostyczne-gpd.png", width=800, height=500)
ismev::gpd.diag(fitGPD)
#dev.off()

### OBLICZENIE POZIOMÓW ZWROTU x20 i x50 ###

n <- length(dane)

# Ręcznie
x20.3 <- as.numeric(u + (beta/xi) * ((20 * (przekroczenia_u/n))^xi - 1))
x50.3 <- as.numeric(u + (beta/xi) * ((50 * (przekroczenia_u/n))^xi - 1))
x20.3; x50.3 # 1003.121, 1006.056

# Sprawdzamy przekroczenia poziomów zwrotu
przekroczenia_x20.3 <- dane_lato %>% filter(cisnienie > x20.3)
przekroczenia_x50.3 <- dane_lato %>% filter(cisnienie > x50.3)

# Liczba przekroczeń
liczba_x20.3 <- nrow(przekroczenia_x20.3)
liczba_x50.3 <- nrow(przekroczenia_x50.3)

cat("Liczba przekroczeń x20:", liczba_x20.3, "\n") # 1671
cat("Liczba przekroczeń x50:", liczba_x50.3, "\n") # 658


# Z wbudowanej funkcji
#fitGPD <- gpd(dane, u)
#xi_est <- fitGPD$par.est[[1]]
#beta_est <- fitGPD$par.est[[2]]
#xi_est; beta_est
#x_20 <- evir::riskmeasures(fitGPD,0.95)[2]
#x_50 <- evir::riskmeasures(fitGPD,0.98)[2]
#x_20; x_50


# Wyniki poziomów zwrotu z trzech modeli
cat("Model 1, x20:", x20.1, ", x50:", x50.1)
cat("Model 2, x20:", x20.2, ", x50:", x50.2)
cat("Model 3, x20:", x20.3, ", x50:", x50.3)
#cat("Model 3 (evir), x20:", x_20, ", x50:", x_50) # wychodzi tak samo jak ręcznie

tabela_wynikow <- data.frame(
  x20 = c(x20.1, x20.2, x20.3),
  x50 = c(x50.1, x50.2, x50.3),
  row.names = c("Model 1", "Model 2", "Model 3")
)
tabela_wynikow



# Nowe dane, wyniki estymacji zapisujemy w pliku Wyniki.Rdata
#save(tabela_statystyk, tabela_wynikow, tabela_st1, tabela_gev,
#     przekroczenia_x20.2, tabela_gpd,
#     file="Wyniki.Rdata")

# Ładowanie danych zapisanych w formacie RData
# dane <- load(file="Wyniki.Rdata")