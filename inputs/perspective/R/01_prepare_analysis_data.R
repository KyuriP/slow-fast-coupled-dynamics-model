# ───────────────────────────────────────────────────────────────────────────────
# 01_prepare_analysis_data.R
# Load HELIUS data and construct the main analytic dataframe
# ───────────────────────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────────────────────
# 1) LOAD DATA
# ───────────────────────────────────────────────────────────────────────────────
# HELIUS data are in SPSS format; read_sav preserves labelled vectors.
helius <- haven::read_sav("data/helius.sav")


# ───────────────────────────────────────────────────────────────────────────────
# 2) BUILD SLOW-LAYER INPUTS
# ───────────────────────────────────────────────────────────────────────────────
# Candidate slow-layer risk indicators.
# Construction is intentionally transparent: z-score inputs, domain means,
# and then a composite mean.
slow_vars <- helius %>%
  select(
    gender      = H1_geslacht,
    age         = H1_lft,
    inc_dif     = H1_InkHhMoeite,
    work_sit    = H1_WerkSit,
    edu_cat     = H1_Opleid,
    arbeid      = H1_Arbeidsparticipatie,
    occ_level   = H1_BeroepsNiveau,
    pcs         = H1_PCS12,
    activity    = H1_Squash_totmwk,
    smoking     = H1_PackYears,
    alcohol     = H1_AlcoholConsumption,
    health_lit  = H1_SBSQ_meanscore,   # SBSQ = health literacy
    soc_support = H1_SSQT,             # SSQT = social support
    discrim     = H1_Discr_sumscore,
    stress_work = H1_StressWerk,
    stress_home = H1_StressHuish,
    lang_diff   = H1_MoeiteNLtaal,
    cult_dist   = H1_CultDistMeanScore6_corrected,
    ethn        = H1_etniciteit
  ) %>%
  # HELIUS missing codes → NA
  mutate(
    across(everything(), ~ na_if(.x, -9)),
    across(everything(), ~ na_if(.x, -1))
  ) %>%
  # Harmonise fields + cohort-specific recodes
  mutate(
    # Culture distance: ensure numeric, discard invalid values
    cult_dist = as.numeric(cult_dist),
    cult_dist = ifelse(cult_dist < 0, NA, cult_dist),

    # Stress items are intended to be 1–4; everything else becomes NA
    stress_home = ifelse(stress_home %in% 1:4, stress_home, NA),
    stress_work = ifelse(stress_work %in% 1:4, stress_work, NA),

    # Ethnicity code: 6 treated as missing here
    ethn = as.numeric(ethn),
    ethn = ifelse(ethn == 6, NA, ethn),

    # Dutch participants: language difficulty / cultural distance are
    # conceptually not applicable. Here missing is coded as 0 = no difficulty.
    cult_dist = ifelse(ethn == 1 & is.na(cult_dist), 0, cult_dist),
    lang_diff = ifelse(ethn == 1 & is.na(lang_diff), 0, lang_diff)
  ) %>%
  # Z-score all slow-layer indicators that enter the composite.
  mutate(
    across(
      c(
        inc_dif, work_sit, edu_cat, arbeid, occ_level,
        pcs, activity, smoking, alcohol,
        health_lit, soc_support, discrim,
        stress_work, stress_home, lang_diff, cult_dist
      ),
      \(x) as.numeric(scale(x))
    )
  ) %>%
  # Reverse protective indicators so higher = higher risk.
  mutate(
    edu_rev      = -edu_cat,
    occ_rev      = -occ_level,
    pcs_rev      = -pcs,
    activity_rev = -activity,
    lit_rev      = -health_lit,
    support_rev  = -soc_support
  )


# ───────────────────────────────────────────────────────────────────────────────
# 3) SYMPTOMS
# ───────────────────────────────────────────────────────────────────────────────
# PHQ-9 items are ordinal 0–3 in HELIUS.
# They are dichotomized later for the Ising analyses.
symptoms_df <- helius %>%
  select(
    anh    = H1_WlbvRecent1,
    dep    = H1_WlbvRecent2,
    slp    = H1_WlbvRecent3,
    ene    = H1_WlbvRecent4,
    app    = H1_WlbvRecent5,
    glt    = H1_WlbvRecent6,
    con    = H1_WlbvRecent7,
    mot    = H1_WlbvRecent_89,
    sui    = H1_WlbvRecent10,
    PHQsum = H1_PHQ9_sumscore
  ) %>%
  mutate(
    across(everything(), ~ na_if(.x, -9)),
    across(everything(), ~ na_if(.x, -1))
  )


# ───────────────────────────────────────────────────────────────────────────────
# 4) CONSTRUCT SLOW RISK LOAD
# ───────────────────────────────────────────────────────────────────────────────
# Domain scores are row means of z-scored indicators.
# The composite is the mean of the three domain scores.
slow_vars <- slow_vars %>%
  mutate(
    # Socioeconomic/resources domain
    SES_risk = rowMeans(
      cbind(inc_dif, work_sit, arbeid, edu_rev, occ_rev, lit_rev, support_rev),
      na.rm = TRUE
    ),

    # Structural / psychosocial stress domain
    structural_stress = rowMeans(
      cbind(discrim, lang_diff, cult_dist, stress_work, stress_home),
      na.rm = TRUE
    ),

    # Health / lifestyle domain
    health_risk = rowMeans(
      cbind(pcs_rev, activity_rev, smoking, alcohol),
      na.rm = TRUE
    ),

    # Composite slow-layer risk level
    slow_risk_load = rowMeans(
      cbind(SES_risk, structural_stress, health_risk),
      na.rm = TRUE
    )
  )


# ───────────────────────────────────────────────────────────────────────────────
# 5) BUILD ANALYSIS DATAFRAME
# ───────────────────────────────────────────────────────────────────────────────
analysis_df <- bind_cols(slow_vars, symptoms_df) %>%
  filter(!is.na(slow_risk_load), !is.na(PHQsum)) %>%
  mutate(
    # Ethnicity grouping: Dutch vs Non-Dutch
    dutch_grp = ifelse(ethn == 1, "Dutch", "Non-Dutch"),

    # Gender coding for covariate adjustment / summaries
    gender_grp = dplyr::case_when(
      gender == 1 ~ "Male",
      gender == 2 ~ "Female",
      TRUE        ~ NA_character_
    ),

    # Age split: median split
    age_grp = ifelse(
      age >= median(age, na.rm = TRUE),
      "Older",
      "Younger"
    ),

    # Slow-risk split: median split
    slow_group = ifelse(
      slow_risk_load > median(slow_risk_load, na.rm = TRUE),
      "High Slow Risk",
      "Low Slow Risk"
    ),

    # Z-scored Slow Risk Load for figures and per-SD interpretations
    slow_risk_z = as.numeric(scale(slow_risk_load))
  ) %>%
  mutate(
    slow_group = factor(
      slow_group,
      levels = c("Low Slow Risk", "High Slow Risk")
    ),
    dutch_grp = factor(
      as.character(dutch_grp),
      levels = c("Dutch", "Non-Dutch")
    ),
    age_grp = factor(
      as.character(age_grp),
      levels = c("Younger", "Older")
    )
  )
