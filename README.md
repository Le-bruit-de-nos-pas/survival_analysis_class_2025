# Survival Analysis: A Mathematical Perspective

Survival analysis is a statistical methodology for analyzing and interpreting time-to-event data, where the event of interest could be death, failure, or another specified outcome. This report provides a mathematical foundation for key concepts in survival analysis, including Kaplan-Meier estimation, Cox proportional hazards regression, parametric survival models, and competing risks.

---

## 1. Survival Analysis Basics

### 1.1 Survival Function

The survival function $S(t)$ represents the probability that a subject survives beyond time $t$:

$$
S(t) = P(T > t),
$$

where $T$ is the random variable denoting the time-to-event.

Key properties of $S(t)$ include:

- $S(0) = 1$: Survival probability at time zero is always 1.
- $\lim_{t \to \infty} S(t) = 0$: Survival probability approaches zero as time increases indefinitely.
- $S(t)$ is a non-increasing function.

### 1.2 Probability Density and Cumulative Distribution Functions

The relationship between survival, cumulative distribution, and probability density functions is:

- Cumulative distribution function (CDF): $F(t) = 1 - S(t)$.
- Probability density function (PDF): $f(t) = \frac{dF(t)}{dt} = -\frac{dS(t)}{dt}$.

### 1.3 Hazard Function

The hazard function $h(t)$ describes the instantaneous rate of occurrence of the event at time $t$, given survival up to that time:

$$
h(t) = \lim_{\Delta t \to 0} \frac{P(t \leq T < t + \Delta t \mid T \geq t)}{\Delta t}.
$$

It can also be expressed in terms of the survival function:

$$
h(t) = -\frac{d}{dt} \log S(t).
$$

### 1.4 Cumulative Hazard Function

The cumulative hazard function $H(t)$ is defined as:

$$
H(t) = \int_0^t h(u) \, du.
$$

It relates to the survival function through:

$$
S(t) = e^{-H(t)}.
$$

---

## 2. Kaplan-Meier Estimation

The Kaplan-Meier estimator provides a non-parametric estimate of the survival function based on observed survival times and censoring. Let $t_1, t_2, \ldots, t_k$ denote the observed event times, and $d_i$ and $n_i$ the number of events and individuals at risk at time $t_i$, respectively. The Kaplan-Meier estimator is given by:

$$
\hat{S}(t) = \prod_{t_i \leq t} \left( 1 - \frac{d_i}{n_i} \right).
$$

### Variance of the Kaplan-Meier Estimator

Using Greenwood's formula, the variance of $\hat{S}(t)$ is:

$$
\text{Var}(\hat{S}(t)) = \hat{S}(t)^2 \sum_{t_i \leq t} \frac{d_i}{n_i (n_i - d_i)}.
$$

### Confidence Intervals for Survival Estimates

Confidence intervals for $\hat{S}(t)$ can be derived using transformations, such as the log-minus-log transformation, to stabilize variances:

$$
CI: \left[ \hat{S}(t)^u, \hat{S}(t)^l \right],
$$

where $l$ and $u$ are derived from the standard error and critical values of the normal distribution.

---

## 3. Cox Proportional Hazards Regression

The Cox proportional hazards model is a semi-parametric approach to relate covariates $\mathbf{x}$ to the hazard function. The model assumes:

$$
h(t \mid \mathbf{x}) = h_0(t) \exp(\mathbf{x}^T \boldsymbol{\beta}),
$$

where:

- $h_0(t)$: Baseline hazard function (unspecified).
- $\mathbf{x}$: Vector of covariates.
- $\boldsymbol{\beta}$: Vector of regression coefficients.

### Partial Likelihood

To estimate $\boldsymbol{\beta}$, Cox proposed the partial likelihood:

$$
L(\boldsymbol{\beta}) = \prod_{i=1}^k \frac{\exp(\mathbf{x}_i^T \boldsymbol{\beta})}{\sum_{j \in R(t_i)} \exp(\mathbf{x}_j^T \boldsymbol{\beta})},
$$

$$
L(\boldsymbol{\beta}) = \prod_{i=1}^k \frac{\exp(\mathbf{x}_i^T \boldsymbol{\beta})}{\sum_{j \in R(t_i)} \exp(\mathbf{x}_j^T \boldsymbol{\beta})}
$$

where $R(t_i)$ is the risk set at time $t_i$.

The log-partial likelihood is maximized to obtain estimates of $\boldsymbol{\beta}$.

### Schoenfeld Residuals

Schoenfeld residuals assess the proportional hazards assumption by checking if covariate effects are constant over time. Deviations indicate a potential violation.

---

## 4. Parametric Survival Models

Parametric survival models specify a distribution for the survival times, such as exponential, Weibull, or log-normal.

### 4.1 Exponential Model

The exponential survival model assumes a constant hazard rate $\lambda$:

$$
h(t) = \lambda, \quad S(t) = e^{-\lambda t}.
$$

The likelihood function is:

$$
L(\lambda) = \prod_{i=1}^n \lambda^{d_i} e^{-\lambda t_i},
$$

where $d_i$ indicates whether the event was observed (1) or censored (0).

### 4.2 Weibull Model

The Weibull model generalizes the exponential by allowing a time-varying hazard:

$$
h(t) = \lambda p t^{p-1}, \quad S(t) = e^{-\lambda t^p}.
$$

- $\lambda$: Scale parameter.
- $p$: Shape parameter.

The likelihood is:

$$
L(\lambda, p) = \prod_{i=1}^n \lambda p t_i^{p-1} e^{-\lambda t_i^p}.
$$

### 4.3 Log-Normal Model

The log-normal model assumes the logarithm of survival times follows a normal distribution. The survival function is:

$$
S(t) = 1 - \Phi\left( \frac{\log t - \mu}{\sigma} \right),
$$

where $\Phi(\cdot)$ is the standard normal CDF, and $\mu$, $\sigma$ are the mean and standard deviation of $\log T$.

---

## 5. Competing Risks

Competing risks occur when multiple types of events can prevent the occurrence of the primary event of interest. Let $k$ denote the event type.

### 5.1 Cause-Specific Hazard

The cause-specific hazard for event type $k$ is:

$$
h_k(t) = \lim_{\Delta t \to 0} \frac{P(t \leq T < t + \Delta t, \text{event } k \mid T \geq t)}{\Delta t}.
$$

### 5.2 Cumulative Incidence Function

The cumulative incidence function (CIF) for event type $k$ is:

$$
F_k(t) = P(T \leq t, \text{event } k).
$$

It relates to the cause-specific hazard:

$$
F_k(t) = \int_0^t S(u) h_k(u) \, du,
$$

where $S(u)$ is the overall survival function.

### Subdistribution Hazard

Fine and Gray proposed the subdistribution hazard to directly model the CIF:

$$
h_k^{\text{sub}}(t) = -\frac{d}{dt} \log\left( 1 - F_k(t) \right).
$$

This approach accounts for competing risks directly in regression modeling.

---

## 6. Advanced Topics in Survival Analysis

### 6.1 Frailty Models

Frailty models incorporate random effects to account for unobserved heterogeneity. The hazard function becomes:

$$
h(t \mid \mathbf{x}, z) = z h_0(t) \exp(\mathbf{x}^T \boldsymbol{\beta}),
$$

where $z$ is a random frailty term, often assumed to follow a gamma distribution.

### 6.2 Time-Varying Covariates

In some studies, covariates change over time. The hazard function is modified as:

$$
h(t \mid \mathbf{x}(t)) = h_0(t) \exp(\mathbf{x}(t)^T \boldsymbol{\beta}).
$$

### 6.3 Joint Models

Joint models simultaneously analyze survival data and longitudinal data (e.g., repeated biomarker measurements) to assess their interplay.

### 6.4 Multiple Imputation for Censored Data

Survival datasets often include missing or censored data. Multiple imputation can address these issues:

- Methods include predictive mean matching or Bayesian approaches.
- Adjust survival estimates using imputed datasets and Rubin’s rules for combining results.

### 6.5 Cure Models

Cure models account for scenarios where a fraction of subjects may never experience the event:

- Mixture cure models: Separate populations into “cured” and “susceptible.”
- Non-mixture models: Directly model the survival curve with asymptotic behavior reflecting cure.

### 6.6 Interval Censoring

Interval censoring occurs when exact event times are unknown but lie within an interval:

- Interval censoring likelihoods and extensions of Kaplan-Meier estimation.
- Methods like the Turnbull estimator for nonparametric inference.

### 6.7 Model Diagnostics and Goodness-of-Fit

Diagnostic techniques for survival models:

- Residual analyses (e.g., Martingale, Deviance, and Schoenfeld residuals).
- Visual checks using Cox-Snell residual plots.
- Likelihood ratio tests, AIC/BIC for model comparison.

### 6.8 Accelerated Failure Time Models

The Accelerated Failure Time (AFT) model is a parametric alternative to Cox regression:

$$
T = \exp(\mathbf{x}^T \boldsymbol{\beta} + \epsilon),
$$

where $\epsilon$ is a random error term. The interpretation of $\boldsymbol{\beta}$ relates to time acceleration or deceleration factors.

### 6.9 Simulation Studies in Survival Analysis

- Generating synthetic survival data under various distributions (e.g., Weibull, Log-normal).
- Techniques for assessing power, bias, and variability in survival studies.

### 6.10 Bayesian Survival Analysis

- Use of Bayesian priors for survival modeling, particularly for small or sparse datasets.
- Markov Chain Monte Carlo (MCMC) techniques for posterior estimation.
- Application to parametric, semi-parametric, and competing risks models.

---

## References

1. Kalbfleisch, J. D., & Prentice, R. L. (2002). _The Statistical Analysis of Failure Time Data_.
2. Cox, D. R. (1972). Regression Models and Life Tables. _Journal of the Royal Statistical Society. Series B_, 34(2), 187–220.
3. Klein, J. P., & Moeschberger, M. L. (2003). _Survival Analysis: Techniques for Censored and Truncated Data_.
4. Fine, J. P., & Gray, R. J. (1999). A Proportional Hazards Model for the Subdistribution of a Competing Risk. _Journal of the American Statistical Association_.
