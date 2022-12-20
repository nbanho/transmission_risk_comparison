library(tidyverse)
library(dplyr)
library(ggpubr)
theme_set(theme_pubr())
options(scipen = 999) #damit Zahlen in Dezimalschreibweise angezeigt werden

#Assumptions
##f_bar
C_a <- 0.038 # Annahme aus Rudnick (Verteilung brauchen?)
C_out <- 325 #p.p.m (Verteilung brauchen)? MW https://www.kane.co.uk/knowledge-centre/what-are-safe-levels-of-co-and-co2-in-rooms

##f Wert hinzufügen, durch 1mio teilen, um p.p.m wegzurechnen (vergleiche mit Rechnung in Rudnick Ende 238)
ch1 <- ch %>% 
  mutate(
    f = ((co2-C_out)/C_a)/1000000) # f ist nun in Dezimalschreibweise

satz1 <- satz %>% 
  mutate(
    f = ((co2-C_out)/C_a)/1000000)

##f_bar pro Schule, da wir keine stetigen Daten haben, können wir einfach Durchschnitt rechnen
f_bar_ch <- ch1 %>% 
  group_by(school) %>% 
  summarize(f_bar = mean(f, na.rm = TRUE))

f_bar_af <- satz1 %>% 
  group_by(country) %>% 
  summarize(f_bar = mean(f, na.rm = TRUE)) ## nur Südafrika im Datensatz?

#Vergleich f-Werte zwischen CH und SA (entweder erst hier Verteilung benötigen, anstatt für alle vorherigen, die darin verrechnet sind)
distribution_f_compare <- ggplot(ch1, aes(x=f)) +
  geom_histogram(binwidth=.001, alpha=.9, position="identity") +
  geom_histogram(data =satz1, binwidth =.001, alpha = .03, colour = "red") +
  ggtitle("Comparison rebreathed fraction CH [black] vs. SA [red]")

distribution_f_compare

#quantum pro Stunde [gemittelt über versch. Studien]
#Ich verwende folgende Studien für den Vergleich:
  #Riley (1962): 130 Patienten, q-Wert: 1.25 
  #Escombe (2008): 117 Patienten, q-Wert: 8.2 
  #Nardell (1991) : 1 Patient, q-Wert: 12.5 
  #Andrews (2014) : 571 "Patienten", q-Wert: 0.89
  #Dhamadhakari (2012) : 17 Patienten, q-Wert: 138/34 (ohne Maske/Maske)

q <- (1.25*130+8.2*117+12.5+0.89*571+138*17)/(130+117+1+571+138)

ch1 %>% count(school, n) #check welche klassen wieviele Schüler

class_school1 <- (14*732 + 24*1954)/(732+1954) #gewichtetes Mittel der Klassengrösse (simpler approach)
class_school2 <- 20
class_school_sa <- 30 #Powerpoint
class_school_tz <- 50 #Powerpoint

prev_y_ch <- 0.006 # Prävalenzannahme Powerpoint in Prozent
prev_y_sa <- 0.7 # Prävalenzannahme Powerpoint oder 0.00432 wenn wir nur die relevante Altersgruppe beachten https://www.thelancet.com/action/showPdf?pii=S1473-3099%2822%2900149-9
prev_y_tan <- 0.4 # Prävalenzannahme Powerpoint

I_ch_1_y <- prev_y_ch/100*class_school1 #Prävalenz pro Klasse (jährlich)
I_ch_2_y <- prev_y_ch/100*class_school2
I_sa_y <- prev_y_sa/100*class_school_sa
I_tz  #Daten fehlen

I_ch_1_s <- I_ch_1_y/360/24 #Prävalenz pro Klasse (stündlich)
I_ch_2_s <- I_ch_2_y/360/24
I_sa_s <- I_sa_y/360/24
I_tz  #Daten fehlen

t = 8 #ein Tag (8 Stunden)
t_y = 1600 #ein Jahr in Stunden 8*5*40

## P_Wert für fixen q Wert {jährliches Ansteckungsrisiko}
P_1 = 1 - exp(-(f_bar_ch[1,2, drop=TRUE]*I_ch_1_s*q*t_y)/class_school1)
P_2 = 1 - exp(-(f_bar_ch[2,2, drop=TRUE]*I_ch_2_s*q*t_y)/class_school2)
P_sa = 1 - exp(-(f_bar_af[1,2, drop=TRUE]*I_sa_s*q*t_y)/class_school_sa)

#Simulation
## Nun nehme ich Verteilung an für Parameter q und simuliere die verschiedenen P-Werte [nehme noch alte Annahme LogNormal als Platzhalter]
##tägliches Infektionsrisiko 
mu2 <- log(q) - 1/2*log((21/q)^2 + 1) #21 aus anderer Ausrechnung (nicht relevant, da Verteilungsannahme noch überschrieben wird)
sigma2 <- 1.7
runs <- 100 #simuliere P-Werte für verschiedene q-Werte, welche ich aus angenommener Verteilung annehme, nehme für alle Parameter stündliche Angabe
sim_P_1 <- 1 - exp(-(f_bar_ch[1,2, drop=TRUE]*I_ch_1_s*rlnorm(runs, meanlog=mu2, sdlog=sigma2)*t)/class_school1)
sim_P_2 <- 1 - exp(-(f_bar_ch[2,2, drop=TRUE]*I_ch_2_s*rlnorm(runs, meanlog=mu2, sdlog=sigma2)*t)/class_school2)
sim_P_sa <- 1 - exp(-(f_bar_af[1,2, drop=TRUE]*I_sa_s*rlnorm(runs, meanlog=mu2, sdlog=sigma2)*t)/class_school_sa)

df_sim <- data.frame(sim_P_1, sim_P_2, sim_P_sa) #in Dataframe damit ich mit ggplot abbilden kann, zudem in tidy format
df_sim_tidy <- tibble(school = c(rep("ch_1", length(sim_P_1)), rep("ch_2", length(sim_P_2)), rep("sa", length(sim_P_sa))), P_sim = c(sim_P_1,sim_P_2,sim_P_sa)) %>% 
  mutate(school = as.factor(school))

## Plot der Simulation der P-Werte (pro Tag)

plot_sim_boxplot <- ggplot(df_sim_tidy, aes(x = school, y=P_sim)) + geom_boxplot()
plot_sim_boxplot

## Auf ein Jahr hochgerechnet {jährliches Ansteckungsrisiko}
sim_P_1_year <- 1 - exp(-(f_bar_ch[1,2, drop=TRUE]*I_ch_1_y*rlnorm(runs, meanlog=mu2, sdlog=sigma2)*t*300)/class_school1)
sim_P_2_year <- 1 - exp(-(f_bar_ch[2,2, drop=TRUE]*I_ch_2_y*rlnorm(runs, meanlog=mu2, sdlog=sigma2)*t*300)/class_school2)
sim_P_sa_year <- 1 - exp(-(f_bar_af[1,2, drop=TRUE]*I_sa_y*rlnorm(runs, meanlog=mu2, sdlog=sigma2)*t*300)/class_school_sa)

df_sim_year <- data.frame(sim_P_1_year, sim_P_2_year, sim_P_sa_year) #in Dataframe damit ich mit ggplot abbilden kann, zudem in tidy format
df_sim_tidy_year <- tibble(school = c(rep("ch_1", length(sim_P_1_year)), rep("ch_2", length(sim_P_2_year)), rep("sa", length(sim_P_sa_year))), P_sim = c(sim_P_1_year,sim_P_2_year,sim_P_sa_year)) %>% 
  mutate(school = as.factor(school))

## Plot der Simulation der P-Werte (jähriches Ansteckungsrisiko)

plot_sim_boxplot_year <- ggplot(df_sim_tidy_year, aes(x = school, y=P_sim, colour = school)) + 
  geom_boxplot() +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  xlab("School") + 
  ylab("Annual risk of infection")
plot_sim_boxplot_year

###### folgende Dinge sind offen / unbeachtet ######
# neue Annahme von Verteilung (Student) für Simulation von Nico übernehmen [Code anfragen]
# evtl. noch genauere Aufteilung der Klassengrössen? Mache weighted mean für Schule_1
# Tanzania Daten?
# über welchen Zeitraum? Wie definiere ich t?
# Wollen wir Ansteckungsrisiko vergleichen für gleichgrosse Klassen oder Unterschiede in den Klassengrössen beibehalten?
# Für Co2 Verteilung annehmen? Rechne momentan für jeden Wert den f-Wert aus und damit dann f_bar
# q Wert von Andrews in PP [3.3] (finde dieses Paper nicht)
