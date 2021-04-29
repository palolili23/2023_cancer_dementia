w.out1 <-
  weightit(
    cancer_v ~ bs(age_0, 3) + sex + education + apoe4 + as.factor(smoke1) +
      bs(sbp1, 3) + bs(bmi1, 3) + ht1 + as.factor(diabetes_prev) + cohort,
    data = data_wide,
    stabilize = TRUE,
    estimand = "ATE",
    method = "ps"
  )
w.out1

# 1.2. IPCW for death, based on baseline covariates ----------------------------

w.out_death <-
  weightit(
    competing_plr ~ bs(age_0, 3) + sex + education + apoe4 + 
      as.factor(smoke1) + ht1 + bs(sbp1, 3) + bs(bmi1,3) +
      as.factor(diabetes_prev) +  cancer_v + cohort,
    data = data_wide,
    stabilize = TRUE,
    estimand = "ATE",
    method = "ps"
  )


# 2. Time-varying cancer -----------------------------------------------------

cancer_den <-
  glm(
    cancer_v ~ bs(age_0, 3) + sex + education + apoe4 + as.factor(smoke1) + ht1 + 
      bs(sbp1, 3) + bs(bmi1, 3) + as.factor(diabetes_prev) + cohort,
    data = data_long,
    family = binomial
  )


# 2.2. IPCW. Weights on death - dementia ---------------------------------------------

death_den <- glm(
  competing_plr ~ cancer_v + bs(time, 3) + bs(age_0, 3) + sex + education + apoe4 +
    as.factor(smoke1) + ht1 + bs(sbp1, 3) + bs(bmi1, 3) + as.factor(diabetes_prev) +
    as.factor(smoke) + bs(sbp, 3) + bs(bmi, 3) + ht + ht_drug + hd_v + stroke_v + diab_v + cohort,
  data = data_long,
  family = quasibinomial
)

death_num <-
  glm(
    competing_plr ~ cancer_v + bs(time, 3) + bs(age_0, 3) + sex + education + apoe4 +
      as.factor(smoke1) + ht1 + bs(sbp1, 3) + bs(bmi1, 3) + as.factor(diabetes_prev) + cohort,
    data = data_long,
    family = binomial
  )

