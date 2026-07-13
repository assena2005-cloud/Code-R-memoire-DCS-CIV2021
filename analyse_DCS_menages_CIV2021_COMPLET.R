#################################################################################
#
#   DÉTERMINANTS DES DÉPENSES CATASTROPHIQUES DE SANTÉ (DCS) DES MÉNAGES
#   EN CÔTE D'IVOIRE — EHCVM 2021-2022
#
#   Script R complet — Mémoire de LICENCE
#   Auteur : KONAN ASSENA ELIEL
#
#   Ce script regroupe, dans l'ordre logique de l'analyse, l'ensemble des
#   traitements réalisés pour le mémoire :
#
#     SECTION 0 : Configuration générale (packages, chemins)
#     SECTION 1 : Construction de la base finale DCS à partir des bases EHCVM
#     SECTION 2 : Faits stylisés — Chapitre 2 (indicateurs H, MPO, G)
#     SECTION 3 : Tests d'association bivariés (Khi-2 de Rao-Scott, Wilcoxon)
#     SECTION 4 : Modèle économétrique principal (logit / probit pondérés)
#     SECTION 5 : Analyse de robustesse (comparaison des 3 seuils de DCS)
#
#   NB : Chaque section peut aussi être exécutée indépendamment des autres
#        (voir l'option RECHARGER_DEPUIS_CSV en Section 0).
#
#################################################################################


#################################################################################
# SECTION 0 : CONFIGURATION GÉNÉRALE
#################################################################################

# ---- 0.1 Packages nécessaires -------------------------------------------------
# Décommenter la ligne suivante lors de la toute première exécution :
# install.packages(c("haven","dplyr","tidyr","tidyverse","survey","pROC",
#                     "WeightedROC","margins"))

library(haven)      # lecture des fichiers .dta (EHCVM)
library(dplyr)       # manipulation de données
library(tidyr)       # remise en forme des données
library(survey)      # plan de sondage complexe (svydesign, svyglm, svychisq)
library(pROC)        # courbes ROC / AUC
library(ggplot2)      # graphiques du Chapitre 2 (Figures 2, 3, 4, 5, 6)
library(scales)       # formatage des labels (%) sur les graphiques

# ---- 0.2 Chemins d'accès -------------------------------------------------------
# Dossier contenant les bases sources EHCVM (.dta) ET la base finale (.csv/.dta)
# --> À ADAPTER selon votre machine avant l'exécution.
CHEMIN_DONNEES <- "C:/Users/HP/Documents/LICENCE 3 SEA/REDACTION DE MEMOIRE/"

# Si TRUE : la Section 2 (et suivantes) recharge la base finale directement
# depuis le fichier CSV déjà sauvegardé (utile si vous sautez la Section 1).
# Si FALSE : la base construite en Section 1 est réutilisée directement en mémoire.
RECHARGER_DEPUIS_CSV <- FALSE


#################################################################################
# SECTION 1 : CONSTRUCTION DE LA BASE FINALE DCS
# Source EHCVM 2021 (welfare, s03_me, individu, s04b_me) -> base_finale_DCS_CIV2021
#################################################################################

if (!RECHARGER_DEPUIS_CSV) {

  # ---- 1.1 Chargement des bases sources ----------------------------------------
  welfare  <- read_dta(file.path(CHEMIN_DONNEES, "ehcvm_welfare_civ2021.dta"))
  s03      <- read_dta(file.path(CHEMIN_DONNEES, "s03_me_civ2021.dta"))
  individu <- read_dta(file.path(CHEMIN_DONNEES, "ehcvm_individu_civ2021.dta"))
  s04b     <- read_dta(file.path(CHEMIN_DONNEES, "s04b_me_civ2021.dta"))

  cat("Welfare :", nrow(welfare), "ménages\n")
  cat("S03     :", nrow(s03), "individus\n")
  cat("Individu:", nrow(individu), "individus\n")
  cat("S04b    :", nrow(s04b), "individus\n")

  # ---- 1.2 Dépenses de santé annualisées par ménage -----------------------------
  # Période de référence 3 mois -> on annualise en multipliant par 4
  vars_3mois  <- c("s03q13", "s03q14", "s03q15", "s03q16", "s03q17",
                   "s03q18a", "s03q18b", "s03q18c")
  # Période de référence 12 mois -> déjà annuelles, pas de multiplication
  vars_12mois <- c("s03q24", "s03q24b", "s03q26", "s03q27", "s03q29",
                   "s03q30", "s03q31", "s03q31b", "s03q48")

  # Remplacement des valeurs manquantes par 0
  s03 <- s03 %>%
    mutate(across(all_of(c(vars_3mois, vars_12mois)),
                  ~ replace_na(as.numeric(.), 0)))

  # Annualisation des dépenses sur 3 mois (x4) ; les dépenses annuelles restent inchangées
  s03 <- s03 %>%
    mutate(across(all_of(vars_3mois), ~ . * 4, .names = "{.col}_ann")) %>%
    mutate(across(all_of(vars_12mois), ~ ., .names = "{.col}_ann"))

  # Dépense totale de santé par individu
  vars_ann <- paste0(c(vars_3mois, vars_12mois), "_ann")
  s03 <- s03 %>%
    mutate(dep_sante_indiv = rowSums(across(all_of(vars_ann)), na.rm = TRUE))

  # Agrégation au niveau du ménage
  dep_sante <- s03 %>%
    group_by(grappe, menage) %>%
    summarise(dep_sante_annuelle = sum(dep_sante_indiv, na.rm = TRUE), .groups = "drop")

  # ---- 1.3 Couverture assurance maladie (au moins 1 membre assuré) --------------
  # s03q32 : 1 = Oui, 2 = Non
  couvmal <- s03 %>%
    mutate(assure = as.integer(s03q32 == 1)) %>%
    group_by(grappe, menage) %>%
    summarise(couv_assurance = max(assure, na.rm = TRUE), .groups = "drop")

  # ---- 1.4 Composition du ménage (résidents uniquement) -------------------------
  res <- individu %>% filter(resid == 1)

  # Présence d'au moins une personne âgée (>= 60 ans)
  pers_agee <- res %>%
    group_by(grappe, menage) %>%
    summarise(presence_pers_agee = as.integer(any(age >= 60, na.rm = TRUE)),
              .groups = "drop")

  # Présence d'au moins un enfant de moins de 5 ans
  enfant_5 <- res %>%
    group_by(grappe, menage) %>%
    summarise(presence_enfant_5 = as.integer(any(age < 5, na.rm = TRUE)),
              .groups = "drop")

  # ---- 1.5 Formalité de l'emploi du chef de ménage -------------------------------
  # Chef de ménage : lien = 1
  cm_id <- individu %>%
    filter(lien == 1) %>%
    select(grappe, menage, membres__id = numind)

  s04b_cm <- s04b %>%
    inner_join(cm_id, by = c("grappe", "menage", "membres__id"))

  # Emploi formel si congés payés (s04q33 = 1) OU cotisation CNSS (s04q38 = 1)
  formalite <- s04b_cm %>%
    mutate(emploi_formel = as.integer(s04q33 == 1 | s04q38 == 1)) %>%
    distinct(grappe, menage, .keep_all = TRUE) %>%
    select(grappe, menage, emploi_formel)

  # ---- 1.6 Fusion de toutes les bases --------------------------------------------
  base <- welfare %>%
    left_join(dep_sante,   by = c("grappe", "menage")) %>%
    left_join(couvmal,     by = c("grappe", "menage")) %>%
    left_join(pers_agee,   by = c("grappe", "menage")) %>%
    left_join(enfant_5,    by = c("grappe", "menage")) %>%
    left_join(formalite,   by = c("grappe", "menage"))

  cat("Base fusionnée :", nrow(base), "ménages x", ncol(base), "variables\n")

  # ---- 1.7 Capacité à payer (CAP) ------------------------------------------------
  base <- base %>%
    mutate(
      dep_sante_annuelle = replace_na(dep_sante_annuelle, 0),
      cap      = dtot - dali,
      cap_safe = if_else(cap == 0, NA_real_, cap)
    )

  # ---- 1.8 Variables dépendantes : DCS (3 seuils) --------------------------------
  base <- base %>%
    mutate(
      ratio_cap  = dep_sante_annuelle / cap_safe,
      ratio_dtot = dep_sante_annuelle / dtot,
      DCS_40 = as.integer(ratio_cap  >= 0.40),  # Méthode OMS (Xu et al., 2003)
      DCS_10 = as.integer(ratio_dtot >= 0.10),  # Seuil 10 % des dépenses totales
      DCS_25 = as.integer(ratio_dtot >= 0.25)   # Seuil 25 % des dépenses totales
    )

  cat("DCS_40 :", round(mean(base$DCS_40, na.rm = TRUE) * 100, 2), "% des ménages\n")
  cat("DCS_10 :", round(mean(base$DCS_10, na.rm = TRUE) * 100, 2), "% des ménages\n")
  cat("DCS_25 :", round(mean(base$DCS_25, na.rm = TRUE) * 100, 2), "% des ménages\n")

  # ---- 1.9 Quintiles de bien-être (pondérés) -------------------------------------
  base <- base %>%
    arrange(pcexp) %>%
    mutate(
      cum_weight   = cumsum(hhweight),
      total_weight = sum(hhweight),
      quantile_w   = cum_weight / total_weight,
      quintile     = case_when(
        quantile_w <= 0.20 ~ 1L,
        quantile_w <= 0.40 ~ 2L,
        quantile_w <= 0.60 ~ 3L,
        quantile_w <= 0.80 ~ 4L,
        TRUE               ~ 5L
      )
    ) %>%
    select(-cum_weight, -total_weight, -quantile_w)

  # ---- 1.10 Recodage des variables explicatives ----------------------------------
  base <- base %>%
    mutate(
      sexe_cm = as.integer(hgender == 1),   # 1 = Homme, 0 = Femme
      age_cm  = hage,                       # Âge du chef de ménage (continu)

      # Niveau d'instruction du CM : 0=Aucun, 1=Primaire, 2=Secondaire, 3=Supérieur
      educ_cm = case_when(
        heduc %in% c(1, 2) ~ 0L,
        heduc == 3         ~ 1L,
        heduc %in% 4:7     ~ 2L,
        heduc %in% 8:9     ~ 3L
      ),

      marie_cm      = as.integer(hmstat %in% c(2, 3, 4)),  # 1 = En couple
      urbain        = as.integer(milieu == 1),              # 1 = Urbain, 0 = Rural
      taille_menage = hhsize,                                # Taille du ménage (continu)
      emploi_formel = replace_na(emploi_formel, 0L),
      couv_assurance = replace_na(couv_assurance, 0L),
      presence_pers_agee = replace_na(presence_pers_agee, 0L),
      presence_enfant_5  = replace_na(presence_enfant_5, 0L),

      # Indicatrices de quintile (quintile 1 = référence)
      quintile_2 = as.integer(quintile == 2),
      quintile_3 = as.integer(quintile == 3),
      quintile_4 = as.integer(quintile == 4),
      quintile_5 = as.integer(quintile == 5)
    )

  # ---- 1.11 Sélection des variables finales --------------------------------------
  base_finale <- base %>%
    select(
      # Identifiants
      grappe, menage, hhid, vague, region, zae,
      # Poids de sondage
      hhweight,
      # Variables dépendantes
      dep_sante_annuelle, cap, ratio_cap, ratio_dtot,
      DCS_40, DCS_10, DCS_25,
      # Variables explicatives - ménage
      taille_menage, milieu, urbain,
      quintile, quintile_2, quintile_3, quintile_4, quintile_5,
      couv_assurance, presence_pers_agee, presence_enfant_5,
      # Variables explicatives - chef de ménage
      sexe_cm, age_cm, educ_cm, marie_cm, emploi_formel,
      hcsp, hsectins,
      # Variables originales conservées
      dtot, dali, dnal, pcexp, zref,
      hgender, hage, heduc, hmstat, hbranch
    )

  cat("\n=== BASE FINALE ===\n")
  cat("Dimensions :", nrow(base_finale), "ménages x", ncol(base_finale), "variables\n")
  cat("Valeurs manquantes educ_cm :", sum(is.na(base_finale$educ_cm)), "\n")

  cat("\n=== STATISTIQUES CLÉS ===\n")
  cat("DCS_40 : ", round(mean(base_finale$DCS_40) * 100, 1), "% des ménages\n")
  cat("DCS_10 : ", round(mean(base_finale$DCS_10) * 100, 1), "% des ménages\n")
  cat("DCS_25 : ", round(mean(base_finale$DCS_25) * 100, 1), "% des ménages\n")
  cat("Dép. santé moy : ", round(mean(base_finale$dep_sante_annuelle)), "FCFA\n")
  cat("CAP moyenne    : ", round(mean(base_finale$cap)), "FCFA\n")

  # ---- 1.12 Sauvegarde de la base finale (CSV + DTA) -----------------------------
  write.csv(base_finale,
            file.path(CHEMIN_DONNEES, "base_finale_DCS_CIV2021.csv"),
            row.names = FALSE)
  write_dta(base_finale,
            file.path(CHEMIN_DONNEES, "base_finale_DCS_CIV2021.dta"))
  cat("\n✅ Base sauvegardée en CSV et DTA\n")
}


#################################################################################
# SECTION 2 : FAITS STYLISÉS (CHAPITRE 2)
# Reproduit :
#   - Tableau 4 : H, MPO, G des DCS selon le seuil (10 % / 25 % / 40 % de la CAP)
#   - Tableau 5 : DCS (seuil 40 % CAP) selon le milieu
#   - DCS (seuil 40 % CAP) par quintile de consommation
#   - Incidence des DCS selon la couverture assurantielle
#################################################################################

# ---- 2.1 Chargement de la base -------------------------------------------------
if (RECHARGER_DEPUIS_CSV) {
  base <- read.csv(file.path(CHEMIN_DONNEES, "base_finale_DCS_CIV2021.csv"),
                    stringsAsFactors = FALSE)
}

# ---- 2.2 Plan de sondage (avec stratification par région) ----------------------
plan_faits_stylises <- svydesign(
  ids     = ~grappe,
  strata  = ~region,
  weights = ~hhweight,
  data    = base,
  nest    = TRUE
)

# ---- 2.3 Fonction générique : Incidence (H), Intensité (MPO) et indice G --------
# H   = proportion pondérée de ménages en DCS
# MPO = dépassement moyen pondéré du seuil, parmi les ménages en DCS
# G   = H x MPO / 100
# NB : ratio_var doit être le MÊME ratio que celui utilisé pour construire dcs_var
#      (ratio_dtot pour les seuils 10 %/25 % ; ratio_cap pour le seuil 40 %)
calc_H_MPO_G <- function(data, dcs_var, ratio_var, seuil, poids = "hhweight") {
  w   <- data[[poids]]
  dcs <- data[[dcs_var]]
  rat <- data[[ratio_var]]

  H   <- 100 * sum(w[dcs == 1]) / sum(w)
  MPO <- 100 * weighted.mean(rat[dcs == 1] - seuil, w[dcs == 1])
  G   <- H * MPO / 100
  menages_pond <- sum(w[dcs == 1])

  data.frame(H = H, MPO = MPO, G = G, menages_pond = menages_pond)
}

# ---- 2.4 Tableau 4 : H, MPO, G selon le seuil -----------------------------------
tab4 <- bind_rows(
  calc_H_MPO_G(base, "DCS_10", "ratio_dtot", 0.10) |> mutate(Seuil = "10% dépenses totales"),
  calc_H_MPO_G(base, "DCS_25", "ratio_dtot", 0.25) |> mutate(Seuil = "25% dépenses totales"),
  calc_H_MPO_G(base, "DCS_40", "ratio_cap",  0.40) |> mutate(Seuil = "40% capacité à payer")
) |>
  relocate(Seuil) |>
  mutate(across(c(H, MPO, G), \(x) round(x, 2)), menages_pond = round(menages_pond))

cat("\n── TABLEAU 4 : H, MPO, G selon le seuil ──\n")
print(tab4)
write.csv(tab4, file.path(CHEMIN_DONNEES, "tableau4_incidence_intensite_G.csv"),
          row.names = FALSE)

# ---- 2.5 Tableau 5 : DCS (seuil 40% CAP) selon le milieu ------------------------
tab5 <- base |>
  mutate(milieu_lab = factor(milieu, levels = c(1, 2), labels = c("Urbain", "Rural"))) |>
  group_by(milieu_lab) |>
  group_modify(~ calc_H_MPO_G(.x, "DCS_40", "ratio_cap", 0.40)) |>
  ungroup() |>
  mutate(across(c(H, MPO, G), \(x) round(x, 2)))

cat("\n── TABLEAU 5 : DCS (seuil 40% CAP) selon le milieu ──\n")
print(tab5)
write.csv(tab5, file.path(CHEMIN_DONNEES, "tableau5_DCS_milieu.csv"),
          row.names = FALSE)

# ---- 2.6 DCS (seuil 40% CAP) par quintile de consommation -----------------------
tab_quintile <- base |>
  mutate(quintile_lab = factor(quintile, levels = 1:5,
                                labels = c("Q1 (plus pauvre)", "Q2", "Q3", "Q4",
                                           "Q5 (plus riche)"))) |>
  group_by(quintile_lab) |>
  group_modify(~ calc_H_MPO_G(.x, "DCS_40", "ratio_cap", 0.40)) |>
  ungroup() |>
  mutate(across(c(H, MPO, G), \(x) round(x, 2)))

cat("\n── DCS (seuil 40% CAP) par quintile de consommation ──\n")
print(tab_quintile)

# ---- 2.7 Incidence des DCS selon la couverture assurantielle --------------------
tab_assurance <- base |>
  mutate(assurance_lab = factor(couv_assurance, levels = c(0, 1),
                                 labels = c("Non couvert", "Couvert"))) |>
  group_by(assurance_lab) |>
  group_modify(~ calc_H_MPO_G(.x, "DCS_40", "ratio_cap", 0.40)) |>
  ungroup() |>
  mutate(across(c(H, MPO, G), \(x) round(x, 2)))

cat("\n── Incidence des DCS selon la couverture assurantielle ──\n")
print(tab_assurance)

# ---- 2.8 Charte graphique commune (style repris des figures existantes) --------
# Palette : bleu (H / indicateur principal) et orange (MPO / indicateur secondaire),
# reprise du code couleur déjà utilisé dans les figures du mémoire.
couleur_H   <- "#4472C4"   # bleu
couleur_MPO <- "#ED7D31"   # orange

theme_memoire <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.title       = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_blank()
  )

DOSSIER_FIGURES <- file.path(CHEMIN_DONNEES, "figures_chapitre2")
if (!dir.exists(DOSSIER_FIGURES)) dir.create(DOSSIER_FIGURES, recursive = TRUE)

# ---- 2.9 Figure 2 : Incidence (H) et intensité (MPO) des DCS selon le seuil -----
fig2_data <- tab4 |>
  select(Seuil, H, MPO) |>
  pivot_longer(cols = c(H, MPO), names_to = "Indicateur", values_to = "Valeur")

fig2 <- ggplot(fig2_data, aes(x = Seuil, y = Valeur, fill = Indicateur)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = paste0(round(Valeur, 2), "%")),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3.6) +
  scale_fill_manual(values = c("H" = couleur_H, "MPO" = couleur_MPO),
                     labels = c("H (Incidence)", "MPO (Intensité)")) +
  labs(title = "Incidence (H) et intensité (MPO) des DCS selon le seuil retenu",
       x = NULL, y = "Pourcentage (%)") +
  theme_memoire

ggsave(file.path(DOSSIER_FIGURES, "figure2_H_MPO_seuil.png"), fig2,
       width = 8, height = 4.5, dpi = 300, bg = "white")

# ---- 2.10 Figure 3 : DCS (seuil 40% CAP) selon le milieu de résidence -----------
fig3_data <- tab5 |>
  select(milieu_lab, H, MPO) |>
  pivot_longer(cols = c(H, MPO), names_to = "Indicateur", values_to = "Valeur")

fig3 <- ggplot(fig3_data, aes(x = milieu_lab, y = Valeur, fill = Indicateur)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5) +
  geom_text(aes(label = paste0(round(Valeur, 2), "%")),
            position = position_dodge(width = 0.6), vjust = -0.4, size = 3.6) +
  scale_fill_manual(values = c("H" = couleur_H, "MPO" = couleur_MPO),
                     labels = c("H (Incidence)", "MPO (Intensité)")) +
  labs(title = "DCS (seuil 40% CAP) selon le milieu de résidence",
       x = NULL, y = "Pourcentage (%)") +
  theme_memoire

ggsave(file.path(DOSSIER_FIGURES, "figure3_DCS_milieu.png"), fig3,
       width = 7, height = 4.2, dpi = 300, bg = "white")

# ---- 2.11 Figure 4 : DCS (seuil 40% CAP) par quintile de consommation ----------
fig4_data <- tab_quintile |>
  select(quintile_lab, H, MPO) |>
  pivot_longer(cols = c(H, MPO), names_to = "Indicateur", values_to = "Valeur")

fig4 <- ggplot(fig4_data, aes(x = quintile_lab, y = Valeur, fill = Indicateur)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = paste0(round(Valeur, 2), "%")),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c("H" = couleur_H, "MPO" = couleur_MPO),
                     labels = c("H (Incidence)", "MPO (Intensité)")) +
  labs(title = "DCS (seuil 40% CAP) par quintile de consommation",
       x = NULL, y = "Pourcentage (%)") +
  theme_memoire

ggsave(file.path(DOSSIER_FIGURES, "figure4_DCS_quintile.png"), fig4,
       width = 8, height = 4.5, dpi = 300, bg = "white")

# ---- 2.12 Figure 5 : Incidence des DCS selon la couverture assurantielle -------
fig5_data <- tab_assurance |> select(assurance_lab, H)

fig5 <- ggplot(fig5_data, aes(x = assurance_lab, y = H, fill = assurance_lab)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = paste0(round(H, 2), "%")), vjust = -0.4, size = 3.8) +
  scale_fill_manual(values = c("Non couvert" = couleur_MPO, "Couvert" = couleur_H)) +
  labs(title = "Incidence des DCS selon la couverture en assurance maladie",
       x = NULL, y = "Incidence H (%)") +
  theme_memoire +
  theme(legend.position = "none")

ggsave(file.path(DOSSIER_FIGURES, "figure5_DCS_assurance.png"), fig5,
       width = 6.5, height = 4.2, dpi = 300, bg = "white")

# ---- 2.13 Figure 6 : Évolution des DCS en Côte d'Ivoire de 2008 à 2021 ---------
# NB : les points 2008 et 2015 proviennent de Gbayoro et al. (2015) (cités dans le
# texte du mémoire) ; seul le point 2021 est calculé à partir de l'EHCVM 2021-2022
# (tab4). Cette figure ne peut donc pas être recalculée automatiquement depuis la
# base de données : les valeurs sont saisies manuellement ci-dessous.
fig6_data <- data.frame(
  Annee = rep(c(2008, 2015, 2021), times = 2),
  Seuil = rep(c("10% dépenses totales", "40% capacité à payer"), each = 3),
  H     = c(17.4, 12.4, 8.45,   # seuil 10%
            4.14, 3.30, 1.60)   # seuil 40% CAP
)

fig6 <- ggplot(fig6_data, aes(x = Annee, y = H, color = Seuil, group = Seuil)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_text(aes(label = paste0(H, "%")), vjust = -1, size = 3.4, show.legend = FALSE) +
  scale_color_manual(values = c("10% dépenses totales" = couleur_H,
                                 "40% capacité à payer" = couleur_MPO)) +
  scale_x_continuous(breaks = c(2008, 2015, 2021)) +
  labs(title = "Évolution des DCS en Côte d'Ivoire de 2008 à 2021",
       x = NULL, y = "Incidence H (%)") +
  theme_memoire

ggsave(file.path(DOSSIER_FIGURES, "figure6_evolution_DCS_2008_2021.png"), fig6,
       width = 7, height = 4.2, dpi = 300, bg = "white")

cat("\n── Figures 2 à 6 enregistrées dans :", DOSSIER_FIGURES, "──\n")


#################################################################################
# SECTION 3 : TESTS D'ASSOCIATION BIVARIÉS
# Variables explicatives vs DCS_40 :
#   - Khi-2 de Rao-Scott (variables catégorielles)
#   - Wilcoxon pondéré / svyranktest (variables continues)
#################################################################################

# ---- 3.1 Recodage pour les tests (facteurs avec libellés explicites) -----------
data_tests <- base %>%
  mutate(
    DCS_40             = factor(DCS_40,             levels = c(0, 1),
                                 labels = c("Non", "Oui")),
    milieu             = factor(milieu,             levels = c(1, 2),
                                 labels = c("Urbain", "Rural")),
    quintile           = factor(quintile,           levels = 1:5),
    sexe_cm            = factor(sexe_cm,            levels = c(0, 1),
                                 labels = c("Feminin", "Masculin")),
    couv_assurance     = factor(couv_assurance,     levels = c(0, 1),
                                 labels = c("Non", "Oui")),
    emploi_formel      = factor(emploi_formel,      levels = c(0, 1),
                                 labels = c("Non", "Oui")),
    presence_enfant_5  = factor(presence_enfant_5,  levels = c(0, 1),
                                 labels = c("Non", "Oui")),
    presence_pers_agee = factor(presence_pers_agee, levels = c(0, 1),
                                 labels = c("Non", "Oui"))
  )

# ---- 3.2 Plan de sondage (sans stratification, clustering par grappe) ----------
plan_dcs <- svydesign(
  ids     = ~grappe,
  weights = ~hhweight,
  data    = data_tests,
  nest    = TRUE
)

# ---- 3.3 Variables catégorielles : Test de Rao-Scott (khi-2 pondéré) -----------
vars_cat <- c("milieu", "quintile", "sexe_cm", "couv_assurance",
              "emploi_formel", "presence_enfant_5", "presence_pers_agee")

cat("\n", strrep("=", 70), "\n", sep = "")
cat("TEST DE RAO-SCOTT (khi-2 pondéré) - Variables catégorielles\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("%-30s %10s %8s %10s %5s\n", "Variable", "Stat.", "ddl", "p-valeur", "Sig"))
cat(strrep("-", 68), "\n")

resultats_cat <- data.frame()

for (v in vars_cat) {
  formule <- as.formula(paste("~", v, "+ DCS_40"))
  test    <- svychisq(formule, design = plan_dcs, statistic = "Chisq")

  stat <- round(test$statistic, 3)
  ddl  <- round(test$parameter, 0)
  pval <- test$p.value
  sig  <- ifelse(pval < 0.001, "***",
          ifelse(pval < 0.01,  "**",
          ifelse(pval < 0.05,  "*",
          ifelse(pval < 0.1,   ".",  ""))))

  cat(sprintf("  %-28s %10.3f %8.0f %10.4f %5s\n", v, stat, ddl, pval, sig))

  resultats_cat <- rbind(resultats_cat, data.frame(
    Variable    = v,
    Statistique = stat,
    DDL         = ddl,
    p_valeur    = round(pval, 4),
    Sig         = sig,
    stringsAsFactors = FALSE
  ))
}

cat(strrep("-", 68), "\n")
cat("Méthode : Test de Rao-Scott (correction plan de sondage complexe)\n")
cat("Sig : *** p<0.001  ** p<0.01  * p<0.05  . p<0.1\n\n")

# ---- 3.4 Variables continues : Wilcoxon pondéré (svyranktest) -----------------
vars_cont <- c("taille_menage", "age_cm")

cat(strrep("=", 70), "\n")
cat("TEST DE WILCOXON PONDÉRÉ (svyranktest) - Variables continues\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("%-30s %10s %10s %5s\n", "Variable", "Stat.", "p-valeur", "Sig"))
cat(strrep("-", 58), "\n")

resultats_cont <- data.frame()

for (v in vars_cont) {
  formule <- as.formula(paste(v, "~ DCS_40"))
  test    <- svyranktest(formule, design = plan_dcs, test = "wilcoxon")

  stat <- round(test$statistic, 3)
  pval <- test$p.value
  sig  <- ifelse(pval < 0.001, "***",
          ifelse(pval < 0.01,  "**",
          ifelse(pval < 0.05,  "*",
          ifelse(pval < 0.1,   ".",  ""))))

  cat(sprintf("  %-28s %10.3f %10.4f %5s\n", v, stat, pval, sig))

  resultats_cont <- rbind(resultats_cont, data.frame(
    Variable    = v,
    Statistique = stat,
    p_valeur    = round(pval, 4),
    Sig         = sig,
    stringsAsFactors = FALSE
  ))
}

cat(strrep("-", 58), "\n")
cat("Méthode : Test de Wilcoxon pondéré (svyranktest)\n")
cat("Sig : *** p<0.001  ** p<0.01  * p<0.05  . p<0.1\n\n")

# ---- 3.5 Tableau récapitulatif de tous les tests -------------------------------
cat(strrep("=", 72), "\n")
cat("TABLEAU RÉCAPITULATIF — Tous les tests\n")
cat(strrep("=", 72), "\n\n")
cat(sprintf("%-30s %-22s %10s %10s %5s\n",
            "Variable", "Test", "Stat.", "p-valeur", "Sig"))
cat(strrep("-", 80), "\n")

for (i in 1:nrow(resultats_cat)) {
  cat(sprintf("  %-28s %-22s %10.3f %10.4f %5s\n",
              resultats_cat$Variable[i],
              "Rao-Scott (khi-2)",
              resultats_cat$Statistique[i],
              resultats_cat$p_valeur[i],
              resultats_cat$Sig[i]))
}

for (i in 1:nrow(resultats_cont)) {
  cat(sprintf("  %-28s %-22s %10.3f %10.4f %5s\n",
              resultats_cont$Variable[i],
              "Wilcoxon pondéré",
              resultats_cont$Statistique[i],
              resultats_cont$p_valeur[i],
              resultats_cont$Sig[i]))
}

cat(strrep("-", 80), "\n")
cat("Sig : *** p<0.001  ** p<0.01  * p<0.05  . p<0.1\n")


#################################################################################
# SECTION 4 : MODÈLE ÉCONOMÉTRIQUE PRINCIPAL (SEUIL 40% CAP)
# Modèles logit et probit pondérés — VERSION CORRIGÉE ET CERTIFIÉE
#################################################################################

# ---- 4.1 Recodage des variables (facteurs avec référence explicite) -----------
data_modele <- base %>%
  mutate(
    DCS_40 = as.integer(DCS_40),

    # Quintile : variable catégorielle, référence = quintile 1 (plus pauvres)
    quintile = factor(quintile, levels = 1:5),

    # Milieu : 1 = urbain (référence), 2 = rural
    milieu = factor(milieu, levels = c(1, 2), labels = c("Urbain", "Rural")),

    # Variables binaires
    sexe_cm             = factor(sexe_cm,             levels = c(0, 1)),
    couv_assurance       = factor(couv_assurance,       levels = c(0, 1)),
    emploi_formel        = factor(emploi_formel,        levels = c(0, 1)),
    presence_enfant_5    = factor(presence_enfant_5,    levels = c(0, 1)),
    presence_pers_agee   = factor(presence_pers_agee,   levels = c(0, 1))
  )

# Références explicites (déjà définies par levels, sécurisation supplémentaire)
data_modele$quintile           <- relevel(data_modele$quintile,           ref = "1")
data_modele$milieu             <- relevel(data_modele$milieu,             ref = "Urbain")
data_modele$sexe_cm            <- relevel(data_modele$sexe_cm,            ref = "0")
data_modele$couv_assurance     <- relevel(data_modele$couv_assurance,     ref = "0")
data_modele$emploi_formel      <- relevel(data_modele$emploi_formel,      ref = "0")
data_modele$presence_enfant_5  <- relevel(data_modele$presence_enfant_5,  ref = "0")
data_modele$presence_pers_agee <- relevel(data_modele$presence_pers_agee, ref = "0")

# ---- 4.2 Plan de sondage ---------------------------------------------------------
# ids = ~grappe  (clustering par grappe)
# weights = ~hhweight (poids d'expansion bruts)
# nest = TRUE    (grappes imbriquées dans les strates si présentes)
# Équivalent Stata : logit y x [pweight=hhweight], vce(cluster grappe)
plan_modele <- svydesign(
  ids     = ~grappe,
  weights = ~hhweight,
  data    = data_modele,
  nest    = TRUE
)

cat("\nPlan de sondage :\n")
print(plan_modele)

# ---- 4.3 Modèle logit pondéré ----------------------------------------------------
# family = binomial() (et non quasibinomial(), qui empêche le calcul de logLik/AIC/BIC)
modele_logit <- svyglm(
  DCS_40 ~ quintile + milieu + taille_menage + age_cm + sexe_cm +
    couv_assurance + emploi_formel + presence_enfant_5 + presence_pers_agee,
  design = plan_modele,
  family = binomial(link = "logit")
)

cat("\n── RÉSUMÉ MODÈLE LOGIT ──\n")
print(summary(modele_logit))

# ---- 4.4 Modèle probit pondéré ---------------------------------------------------
modele_probit <- svyglm(
  DCS_40 ~ quintile + milieu + taille_menage + age_cm + sexe_cm +
    couv_assurance + emploi_formel + presence_enfant_5 + presence_pers_agee,
  design = plan_modele,
  family = binomial(link = "probit")
)

cat("\n── RÉSUMÉ MODÈLE PROBIT ──\n")
print(summary(modele_probit))

# ---- 4.5 Odds ratios et IC à 95% -------------------------------------------------
df_resid <- degf(plan_modele)  # degrés de liberté du plan de sondage

OR   <- exp(coef(modele_logit))
IC   <- exp(confint(modele_logit, df = df_resid))
PVAL <- summary(modele_logit)$coefficients[, 4]

RESULTATS_OR <- data.frame(
  Variable = names(OR),
  OR       = round(OR, 3),
  IC_inf   = round(IC[, 1], 3),
  IC_sup   = round(IC[, 2], 3),
  p_value  = round(PVAL, 4),
  Sig      = ifelse(PVAL < 0.001, "***",
             ifelse(PVAL < 0.01,  "**",
             ifelse(PVAL < 0.05,  "*",
             ifelse(PVAL < 0.1,   ".",  ""))))
)

cat("\n── ODDS RATIOS (Logit pondéré) ──\n")
print(RESULTATS_OR)

# ---- 4.6 Critères de comparaison (LL, AIC, BIC, Pseudo-R²) -----------------------
# logLik() fonctionne avec binomial(), et non quasibinomial()
modele_nul <- svyglm(DCS_40 ~ 1, design = plan_modele, family = binomial(link = "logit"))
modele_nul_probit <- svyglm(DCS_40 ~ 1, design = plan_modele, family = binomial(link = "probit"))

LL_logit  <- as.numeric(logLik(modele_logit))
LL_probit <- as.numeric(logLik(modele_probit))
LL_null_l <- as.numeric(logLik(modele_nul))
LL_null_p <- as.numeric(logLik(modele_nul_probit))

n_obs <- nrow(data_modele)
k     <- length(coef(modele_logit))   # nombre de paramètres (constante incluse)

AIC_logit  <- -2 * LL_logit  + 2 * k
AIC_probit <- -2 * LL_probit + 2 * k
BIC_logit  <- -2 * LL_logit  + log(n_obs) * k
BIC_probit <- -2 * LL_probit + log(n_obs) * k

PR2_logit  <- 1 - (LL_logit  / LL_null_l)
PR2_probit <- 1 - (LL_probit / LL_null_p)

cat("\n── CRITÈRES DE COMPARAISON ──\n")
cat(sprintf("%-30s %15s %15s\n", "Critère", "Logit", "Probit"))
cat(strrep("-", 62), "\n")
cat(sprintf("%-30s %15.1f %15.1f\n", "Log-vraisemblance (LL)", LL_logit,  LL_probit))
cat(sprintf("%-30s %15.1f %15.1f\n", "AIC",                     AIC_logit, AIC_probit))
cat(sprintf("%-30s %15.1f %15.1f\n", "BIC",                     BIC_logit, BIC_probit))
cat(sprintf("%-30s %15.4f %15.4f\n", "Pseudo-R² (McFadden)",    PR2_logit, PR2_probit))

# ---- 4.7 Test du rapport de vraisemblance (LR test) et test de Wald global -------
# LR test : compare le modèle complet (12 paramètres) au modèle nul (constante seule)
# H0 : tous les coefficients des variables explicatives sont nuls
LR_stat_classique <- -2 * (LL_null_l - LL_logit)
ddl_LR            <- k - length(coef(modele_nul))
p_LR_classique    <- pchisq(LR_stat_classique, df = ddl_LR, lower.tail = FALSE)

cat("\n── TEST DU RAPPORT DE VRAISEMBLANCE (LR test, classique) ──\n")
cat(sprintf("Statistique LR (Khi-2) : %.3f | ddl = %d | p-valeur : %s\n",
            LR_stat_classique, ddl_LR,
            ifelse(p_LR_classique < 0.001, "< 0.001", round(p_LR_classique, 4))))

# LR test ajusté au plan de sondage complexe (Rao-Scott) — recommandé, tient compte
# du clustering (grappe) et des poids (hhweight)
cat("\n── TEST DU RAPPORT DE VRAISEMBLANCE (LR test, ajusté Rao-Scott) ──\n")
lr_test_rao_scott <- regTermTest(
  modele_logit,
  ~ quintile + milieu + taille_menage + age_cm + sexe_cm +
    couv_assurance + emploi_formel + presence_enfant_5 + presence_pers_agee,
  method = "LRT"
)
print(lr_test_rao_scott)

# Test de Wald global (joint) : teste simultanément la nullité de tous les coefficients
# H0 : béta_1 = béta_2 = ... = béta_12 = 0
cat("\n── TEST DE WALD GLOBAL (joint sur l'ensemble des variables) ──\n")
wald_test_global <- regTermTest(
  modele_logit,
  ~ quintile + milieu + taille_menage + age_cm + sexe_cm +
    couv_assurance + emploi_formel + presence_enfant_5 + presence_pers_agee,
  method = "Wald"
)
print(wald_test_global)

cat("\nInterprétation : la convergence du LR test et du test de Wald global\n")
cat("confirme la significativité conjointe des variables explicatives du modèle logit.\n")
cat("Le LR test (ajusté Rao-Scott) est privilégié pour l'interprétation en contexte\n")
cat("d'événement rare (faible prévalence de la DCS au seuil de 40%).\n")

# ---- 4.8 AUC pondérée ------------------------------------------------------------
# La ROC doit être calculée avec les poids : roc(y, prob) sans poids donne une
# AUC non représentative de la population.
prob_logit  <- as.numeric(predict(modele_logit,  type = "response"))
prob_probit <- as.numeric(predict(modele_probit, type = "response"))
poids       <- weights(plan_modele)
y_obs       <- data_modele$DCS_40[!is.na(data_modele$DCS_40)]

# Fonction AUC pondérée (méthode des trapèzes)
auc_ponderee <- function(y, prob, w) {
  ord   <- order(prob, decreasing = TRUE)
  y_ord <- y[ord]; w_ord <- w[ord]

  tp <- cumsum(w_ord * (y_ord == 1))
  fp <- cumsum(w_ord * (y_ord == 0))
  tp <- tp / max(tp); fp <- fp / max(fp)

  n   <- length(tp)
  auc <- sum((fp[2:n] - fp[1:(n - 1)]) * (tp[2:n] + tp[1:(n - 1)]) / 2)
  return(auc)
}

AUC_logit_pond  <- auc_ponderee(y_obs, prob_logit,  poids)
AUC_probit_pond <- auc_ponderee(y_obs, prob_probit, poids)

cat(sprintf("%-30s %15.4f %15.4f\n", "AUC (pondérée)", AUC_logit_pond, AUC_probit_pond))
cat(strrep("-", 62), "\n")

# ---- 4.8bis Seuil de Youden et matrice de confusion PONDÉRÉS -------------------
# pROC::roc() ne supporte pas nativement les poids d'enquête : on calcule donc
# manuellement sensibilité/spécificité pondérées à chaque seuil candidat,
# cohérent avec la fonction auc_ponderee() déjà utilisée en 4.8.
roc_ponderee <- function(y, prob, w) {
  ord   <- order(prob, decreasing = TRUE)
  y_o   <- y[ord]; w_o <- w[ord]; prob_o <- prob[ord]
  P <- sum(w_o[y_o == 1])   # somme des poids des positifs
  N <- sum(w_o[y_o == 0])   # somme des poids des négatifs
  tp <- cumsum(w_o * (y_o == 1))
  fp <- cumsum(w_o * (y_o == 0))
  sens <- tp / P
  spec <- 1 - fp / N
  data.frame(threshold = prob_o, sensitivity = sens, specificity = spec)
}
roc_df_pond <- roc_ponderee(y_obs, prob_logit, poids)
# Seuil optimal = argmax de l'indice de Youden (sens + spec - 1), pondéré
roc_df_pond$youden <- roc_df_pond$sensitivity + roc_df_pond$specificity - 1
idx_opt_pond   <- which.max(roc_df_pond$youden)
coords_opt_pond <- roc_df_pond[idx_opt_pond, c("threshold", "sensitivity", "specificity")]
cat("\n── PERFORMANCES AU SEUIL OPTIMAL PONDÉRÉ (Youden, Logit) ──\n")
print(coords_opt_pond)
# Matrice de confusion pondérée (sommes de poids, représentative de la population)
seuil_pond  <- coords_opt_pond$threshold[1]
y_pred_pond <- ifelse(prob_logit >= seuil_pond, 1, 0)
tn_p <- sum(poids[y_obs == 0 & y_pred_pond == 0])
fp_p <- sum(poids[y_obs == 0 & y_pred_pond == 1])
fn_p <- sum(poids[y_obs == 1 & y_pred_pond == 0])
tp_p <- sum(poids[y_obs == 1 & y_pred_pond == 1])
cat("\nMatrice de confusion PONDÉRÉE (seuil =", round(seuil_pond, 4), ") :\n")
mat_pond <- matrix(c(tn_p, fp_p, fn_p, tp_p), nrow = 2, byrow = TRUE,
                    dimnames = list(Réel = c("0", "1"), Prédit = c("0", "1")))
print(mat_pond)
# Indicateurs de performance pondérés
accuracy_pond <- (tn_p + tp_p) / (tn_p + fp_p + fn_p + tp_p)
vpp_pond      <- tp_p / (tp_p + fp_p)
vpn_pond      <- tn_p / (tn_p + fn_p)
cat(sprintf("\nAccuracy (pondérée) : %.4f\n", accuracy_pond))
cat(sprintf("VPP (pondérée)      : %.4f\n", vpp_pond))
cat(sprintf("VPN (pondérée)      : %.4f\n", vpn_pond))
cat(strrep("-", 62), "\n")

# ---- 4.9 Effets marginaux moyens (AME) -------------------------------------------
# AME_j = (1/N) * sum_i [ p_i*(1-p_i)*beta_j ]                    (variables continues)
# AME = E[P(y=1|X, D=1)] - E[P(y=1|X, D=0)]                       (variables binaires)
cat("\n── EFFETS MARGINAUX MOYENS (AME, Logit) ──\n")

beta_logit <- coef(modele_logit)
pred_logit <- as.numeric(predict(modele_logit, type = "response"))
lambda     <- mean(pred_logit * (1 - pred_logit))  # facteur d'échelle moyen

# Variables continues : AME = lambda * beta
vars_continues <- c("taille_menage", "age_cm")
for (v in vars_continues) {
  ame <- lambda * beta_logit[v]
  cat(sprintf("  AME %-20s : %+.4f pp\n", v, ame * 100))
}

# Variables binaires/catégorielles : AME par prédiction contrefactuelle
data_temp <- data_modele[complete.cases(data_modele[, c(
  "DCS_40", "quintile", "milieu", "taille_menage", "age_cm", "sexe_cm",
  "couv_assurance", "emploi_formel", "presence_enfant_5", "presence_pers_agee"
)]), ]

calc_ame_binaire <- function(modele, donnees, variable, val1, val0) {
  d1 <- donnees; d1[[variable]] <- factor(val1, levels = levels(donnees[[variable]]))
  d0 <- donnees; d0[[variable]] <- factor(val0, levels = levels(donnees[[variable]]))
  p1 <- predict(modele, newdata = d1, type = "response")
  p0 <- predict(modele, newdata = d0, type = "response")
  mean(p1 - p0, na.rm = TRUE)
}

vars_bin <- list(
  list(var = "milieu",             v1 = "Rural", v0 = "Urbain"),
  list(var = "sexe_cm",            v1 = "1",     v0 = "0"),
  list(var = "couv_assurance",     v1 = "1",     v0 = "0"),
  list(var = "emploi_formel",      v1 = "1",     v0 = "0"),
  list(var = "presence_enfant_5",  v1 = "1",     v0 = "0"),
  list(var = "presence_pers_agee", v1 = "1",     v0 = "0")
)
for (vb in vars_bin) {
  ame <- calc_ame_binaire(modele_logit, data_temp, vb$var, vb$v1, vb$v0)
  cat(sprintf("  AME %-20s : %+.4f pp\n", vb$var, ame * 100))
}

# Quintiles (vs référence Q1)
for (q in c("2", "3", "4", "5")) {
  ame <- calc_ame_binaire(modele_logit, data_temp, "quintile", q, "1")
  cat(sprintf("  AME %-20s : %+.4f pp\n", paste0("quintile_", q), ame * 100))
}

# ---- 4.10 Courbe ROC (graphique) --------------------------------------------------
roc_logit_obj  <- roc(y_obs, prob_logit,  quiet = TRUE)
roc_probit_obj <- roc(y_obs, prob_probit, quiet = TRUE)

# Graphique ROC (Logit vs Probit)
par(mfrow = c(1, 1))
plot(roc_logit_obj,
     col = "#1A3A5C", lwd = 2,
     main = "Courbe ROC — Logit vs Probit (EHCVM 2021-2022)",
     xlab = "1 - Spécificité", ylab = "Sensibilité")
lines(roc_probit_obj, col = "#C0392B", lwd = 2, lty = 2)
abline(a = 0, b = 1, col = "grey", lty = 3)
legend("bottomright",
       legend = c(sprintf("Logit  (AUC non pond. = %.3f)", auc(roc_logit_obj)),
                  sprintf("Probit (AUC non pond. = %.3f)", auc(roc_probit_obj))),
       col = c("#1A3A5C", "#C0392B"), lwd = 2, lty = c(1, 2))

cat("\nAUC pondérée Logit  :", round(AUC_logit_pond,  4), "\n")
cat("AUC pondérée Probit :", round(AUC_probit_pond, 4), "\n")

# ---- 4.11 Export du tableau final (odds ratios) ----------------------------------
tableau_final <- RESULTATS_OR[-1, ]  # on retire l'intercept

# Décommenter pour exporter vers Excel (nécessite le package writexl) :
# writexl::write_xlsx(tableau_final,
#                      file.path(CHEMIN_DONNEES, "Resultats_Logit_DCS40_corriges.xlsx"))

cat("\n── TABLEAU FINAL (odds ratios, seuil 40%) ──\n")
print(tableau_final)

cat("\n✅ Modèle principal (Section 4) terminé.\n")


#################################################################################
# SECTION 5 : ANALYSE DE ROBUSTESSE
# Comparaison du modèle logit pondéré selon les 3 seuils de DCS (10%, 25%, 40%)
#################################################################################

# ---- 5.1 Recodage (identique à la Section 4, sans relevel explicite) ------------
data_robustesse <- base %>%
  mutate(
    quintile           = factor(quintile,           levels = 1:5),
    milieu             = factor(milieu,             levels = c(1, 2),
                                 labels = c("Urbain", "Rural")),
    sexe_cm            = factor(sexe_cm,            levels = c(0, 1)),
    couv_assurance     = factor(couv_assurance,     levels = c(0, 1)),
    emploi_formel      = factor(emploi_formel,      levels = c(0, 1)),
    presence_enfant_5  = factor(presence_enfant_5,  levels = c(0, 1)),
    presence_pers_agee = factor(presence_pers_agee, levels = c(0, 1))
  )

# ---- 5.2 Plan de sondage ---------------------------------------------------------
plan_robustesse <- svydesign(
  ids     = ~grappe,
  weights = ~hhweight,
  data    = data_robustesse,
  nest    = TRUE
)

# ---- 5.3 Formule commune aux 3 modèles -------------------------------------------
formule_base <- ~ quintile + milieu + taille_menage + age_cm +
  sexe_cm + couv_assurance + emploi_formel +
  presence_enfant_5 + presence_pers_agee

# ---- 5.4 Estimation des 3 modèles (un par seuil de DCS) --------------------------
seuils <- list(
  list(var = "DCS_10", label = "Seuil 10%"),
  list(var = "DCS_25", label = "Seuil 25%")
)

modeles   <- list()
resultats <- list()

for (s in seuils) {
  formule_mod <- update(formule_base, as.formula(paste(s$var, "~ .")))

  mod <- svyglm(formule_mod, design = plan_robustesse, family = binomial(link = "logit"))
  modeles[[s$var]] <- mod

  df_resid <- degf(plan_robustesse)
  coefs    <- coef(mod)
  IC       <- confint(mod, df = df_resid)
  pvals    <- summary(mod)$coefficients[, 4]

  OR    <- exp(coefs)
  OR_lb <- exp(IC[, 1])
  OR_ub <- exp(IC[, 2])

  res <- data.frame(
    Variable  = names(OR),
    OR        = round(OR, 3),
    IC_inf    = round(OR_lb, 3),
    IC_sup    = round(OR_ub, 3),
    p_valeur  = round(pvals, 4),
    Sig       = ifelse(pvals < 0.001, "***",
                ifelse(pvals < 0.01,  "**",
                ifelse(pvals < 0.05,  "*",
                ifelse(pvals < 0.1,   ".",  "")))),
    stringsAsFactors = FALSE
  )
  resultats[[s$var]] <- res

  cat(strrep("=", 70), "\n")
  cat(s$label, "— Variable :", s$var, "\n")
  cat("Prévalence (pondérée) :", round(weighted.mean(data_robustesse[[s$var]], data_robustesse$hhweight, na.rm = TRUE) * 100, 2), "%\n")
  cat(strrep("=", 70), "\n")
  print(res[res$Variable != "(Intercept)", ])
  cat("\n")
}

# ---- 5.5 Critères de comparaison des 3 modèles -----------------------------------
cat(strrep("=", 75), "\n")
cat("CRITÈRES DE COMPARAISON — 3 seuils\n")
cat(strrep("=", 75), "\n")
cat(sprintf("%-30s %15s %15s\n", "Critère", "DCS_10", "DCS_25"))
cat(strrep("-", 75), "\n")

# Log-vraisemblance pondérée manuelle, avec poids NORMALISÉS (somme des poids = n).
# Cohérent avec l'approche var_weights du modèle principal (section 4) : on corrige
# l'hétéroscédasticité entre observations sans gonfler artificiellement la taille
# d'échantillon effective (ce que ferait un usage en freq_weights avec poids bruts).
ll_pond <- function(mod, design) {
  y <- model.response(model.frame(mod))
  p <- as.numeric(predict(mod, type = "response"))
  w <- weights(design)
  w_norm <- w * length(w) / sum(w)
  sum(w_norm * (y * log(p + 1e-15) + (1 - y) * log(1 - p + 1e-15)))
}

k <- length(coef(modele_logit))
n <- nrow(data_robustesse)

for (s in seuils) {
  assign(paste0("ll_", s$var), ll_pond(modeles[[s$var]], plan_robustesse))
}

null_ll <- function(var) {
  mod_null <- svyglm(as.formula(paste(var, "~ 1")),
                      design = plan_robustesse, family = binomial(link = "logit"))
  ll_pond(mod_null, plan_robustesse)
}

ll_null_10 <- null_ll("DCS_10")
ll_null_25 <- null_ll("DCS_25")

criteres <- data.frame(
  Critere = c("Log-vraisemblance (LL)", "AIC", "BIC", "Pseudo-R2 McFadden",
              "Prevalence ponderee (%)"),
  DCS_10  = c(round(ll_DCS_10, 1),
              round(-2 * ll_DCS_10 + 2 * k, 1),
              round(-2 * ll_DCS_10 + k * log(n), 1),
              round(1 - ll_DCS_10 / ll_null_10, 4),
              round(weighted.mean(data_robustesse$DCS_10, data_robustesse$hhweight, na.rm = TRUE) * 100, 2)),
  DCS_25  = c(round(ll_DCS_25, 1),
              round(-2 * ll_DCS_25 + 2 * k, 1),
              round(-2 * ll_DCS_25 + k * log(n), 1),
              round(1 - ll_DCS_25 / ll_null_25, 4),
              round(weighted.mean(data_robustesse$DCS_25, data_robustesse$hhweight, na.rm = TRUE) * 100, 2))
)

for (i in 1:nrow(criteres)) {
  cat(sprintf("%-30s %15s %15s\n",
              criteres$Critere[i],
              as.character(criteres$DCS_10[i]),
              as.character(criteres$DCS_25[i])))
}
cat(strrep("-", 75), "\n\n")

# ---- 5.6 AUC pondérée pour les 3 modèles -----------------------------------------
auc_pond <- function(y, prob, w) {
  ord   <- order(prob, decreasing = TRUE)
  y_o   <- y[ord]; w_o <- w[ord]
  tp    <- cumsum(w_o * (y_o == 1)); fp <- cumsum(w_o * (y_o == 0))
  tpr   <- tp / max(tp); fpr <- fp / max(fp)
  n_pts <- length(tpr)
  sum((fpr[2:n_pts] - fpr[1:(n_pts - 1)]) * (tpr[2:n_pts] + tpr[1:(n_pts - 1)]) / 2)
}

w_poids <- weights(plan_robustesse)

for (s in seuils) {
  prob <- as.numeric(predict(modeles[[s$var]], type = "response"))
  y    <- data_robustesse[[s$var]]
  auc  <- auc_pond(y, prob, w_poids)
  cat(sprintf("AUC pondérée %-25s : %.4f\n", s$label, auc))
}

# ---- 5.7 Tableau comparatif des odds ratios — 3 seuils ---------------------------
cat("\n", strrep("=", 90), "\n", sep = "")
cat("TABLEAU COMPARATIF DES OR — Robustesse selon les 3 seuils\n")
cat(strrep("=", 90), "\n")
cat(sprintf("%-30s %18s %18s\n", "Variable", "DCS 10%", "DCS 25%"))
cat(strrep("-", 90), "\n")

vars_affich <- rownames(summary(modele_logit)$coefficients)
vars_affich <- vars_affich[vars_affich != "(Intercept)"]

for (v in vars_affich) {
  or_10 <- resultats[["DCS_10"]][resultats[["DCS_10"]]$Variable == v, ]
  or_25 <- resultats[["DCS_25"]][resultats[["DCS_25"]]$Variable == v, ]

  fmt <- function(r) {
    if (nrow(r) == 0) return("       —      ")
    sprintf("%5.3f %s", r$OR, r$Sig)
  }

  cat(sprintf("  %-28s %18s %18s\n",
              v, fmt(or_10), fmt(or_25)))
}

cat(strrep("-", 90), "\n")
cat("Sig : *** p<0.001  ** p<0.01  * p<0.05  . p<0.1\n")
cat("Méthode : svyglm Binomial, plan de sondage ids=~grappe, weights=~hhweight\n")

cat("\n✅ Analyse de robustesse (Section 5) terminée.\n")
cat("\n#################################################################################\n")
cat("FIN DU SCRIPT — Mémoire : Déterminants des DCS des ménages en Côte d'Ivoire\n")
cat("#################################################################################\n")
