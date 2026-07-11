# Déterminants des Dépenses Catastrophiques de Santé (DCS) des ménages en Côte d'Ivoire

Code source R utilisé dans le cadre d'un mémoire de Master portant sur les déterminants des dépenses catastrophiques de santé (DCS) des ménages ivoiriens, à partir des données de l'**Enquête Harmonisée sur les Conditions de Vie des Ménages (EHCVM) 2021-2022**.

**Auteur :** Konan Assena Eliel
**Cadre :** Mémoire de Licence

---

## 📋 Contenu du dépôt

| Fichier | Description |
|---|---|
| `analyse_DCS_menages_CIV2021_COMPLET.R` | Script R complet regroupant l'ensemble des traitements et analyses |

## 🔍 Structure du script

Le script est organisé en 5 sections exécutées séquentiellement :

1. **Construction de la base finale** — Fusion des bases sources EHCVM 2021 (`welfare`, `s03_me`, `individu`, `s04b_me`), calcul des dépenses de santé annualisées, de la capacité à payer (CAP) et des indicatrices de DCS selon 3 seuils.
2. **Faits stylisés** — Calcul des indicateurs d'incidence (H), d'intensité (MPO) et de l'indice G des DCS, selon le seuil, le milieu de résidence, le quintile de bien-être et la couverture assurantielle.
3. **Tests d'association bivariés** — Tests du Khi-2 de Rao-Scott (variables catégorielles) et de Wilcoxon pondéré (variables continues) entre les variables explicatives et la DCS.
4. **Modèle économétrique principal** — Modèles logit et probit pondérés (seuil de 40 % de la capacité à payer), avec odds ratios, critères de comparaison (AIC, BIC, pseudo-R² de McFadden), AUC pondérée et effets marginaux moyens (AME).
5. **Analyse de robustesse** — Comparaison des résultats du modèle logit selon les 3 seuils de DCS retenus (10 %, 25 %, 40 %).

## 📊 Source des données

Enquête Harmonisée sur les Conditions de Vie des Ménages (EHCVM), Côte d'Ivoire, 2021-2022.
Les fichiers de données sources (`.dta`) ne sont **pas inclus** dans ce dépôt pour des raisons de confidentialité et de droits d'accès aux microdonnées ; ils peuvent être obtenus auprès de l'institut statistique national compétent.

## ⚙️ Prérequis

```r
install.packages(c("haven", "dplyr", "tidyr", "survey", "pROC"))
```

## ▶️ Utilisation

1. Ouvrir le script `analyse_DCS_menages_CIV2021_COMPLET.R`
2. Modifier la variable `CHEMIN_DONNEES` (Section 0) pour pointer vers le dossier contenant vos fichiers de données EHCVM
3. Exécuter le script dans l'ordre des sections

Si vous disposez déjà de la base finale construite (`base_finale_DCS_CIV2021.csv`), vous pouvez passer directement aux sections 2 à 5 en réglant `RECHARGER_DEPUIS_CSV <- TRUE`.

## 📖 Méthodologie et références

- Grossman, M. (1972). *On the Concept of Health Capital and the Demand for Health.*
- Arrow, K. J. (1963). *Uncertainty and the Welfare Economics of Medical Care.*
- Xu, K., Evans, D. B., Kawabata, K., Zeramdini, R., Klavus, J., & Murray, C. J. L. (2003). *Household catastrophic health expenditure: a multicountry analysis.* The Lancet.

## Citation

Si vous souhaitez faire référence à ce code :

> Konan Assena Eliel (2026). *Code source — Déterminants des dépenses catastrophiques de santé des ménages en Côte d'Ivoire (EHCVM 2021-2022).* Mémoire de Licence. Disponible sur GitHub : [lien du dépôt]

---

*Ce dépôt accompagne un mémoire de Master et est mis à disposition à des fins de transparence méthodologique et de reproductibilité.*
