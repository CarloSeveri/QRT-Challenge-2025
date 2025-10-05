### Methodology and Findings

Our study aimed to predict patient survival outcomes by integrating clinical and molecular data, leveraging advanced feature engineering, rigorous preprocessing, and robust survival modeling techniques. The pipeline combined domain knowledge with machine learning to capture both linear and non-linear patterns in the data. 

#### Feature Engineering and Preprocessing

We began with a comprehensive feature engineering phase, focusing on both clinical and molecular predictors. For clinical biochemical features, we considered key blood markers including $$BM \textunderscore BLAST$$ (bone marrow blast percentage), $$HB$$ (hemoglobin), $$PLT$$ (platelet count), $$WBC$$ (white blood cell count), $$ANC$$ (absolute neutrophil count), and $$MONOCYTES$$. These markers are directly related to hematopoietic function and disease burden in patients. To capture non-linear effects and interactions between these features, we created polynomial transformations of degree 3, considering only interaction terms. This allows the model to account for synergistic relationships between multiple biomarkers without introducing unnecessary complexity from higher-degree monomials.  

In addition to polynomial interactions, we applied logarithmic, square root, and squared transformations to the raw values of biochemical features:

$$
X \textunderscore log = \log(1 + X), \quad X \textunderscore sqrt = \sqrt{X}, \quad X \textunderscore square = X^2
$$

The logarithmic transformation mitigates the influence of extreme outliers and skewed distributions, while square and square-root transformations help capture potential curvature in the relationships between biomarkers and survival outcomes.

We further derived biologically meaningful ratios between cell types to capture systemic relationships and relative proportions in hematopoiesis:

$$
WBC \textunderscore to \textunderscore HB = \frac{WBC}{HB + \epsilon}, \quad PLT \textunderscore to \textunderscore ANC = \frac{PLT}{ANC + \epsilon}, \quad BLAST \textunderscore to \textunderscore WBC = \frac{BM \textunderscore BLAST}{WBC + \epsilon}
$$

$$
ANC \textunderscore to \textunderscore MONOCYTES = \frac{ANC}{MONOCYTES + \epsilon}, \quad PLT \textunderscore to \textunderscore HB = \frac{PLT}{HB + \epsilon}, \quad BLAST \textunderscore to \textunderscore ANC = \frac{BM \textunderscore BLAST}{ANC + \epsilon}, \quad BLAST \textunderscore to \textunderscore PLT = \frac{BM \textunderscore BLAST}{PLT + \epsilon}, \quad MONO \textunderscore to \textunderscore WBC = \frac{MONOCYTES}{WBC + \epsilon}
$$

Here, $$\epsilon = 1e-5$$ prevents division by zero. These ratios provide insights into the relative burden of different blood cell populations, reflecting both disease progression and patient physiology. We then applied log-transformations to all ratio features to normalize their distributions and reduce heteroscedasticity.

On the molecular side, we focused on mutations in high-risk genes such as $$TP53$$, $$RUNX1$$, and $$ASXL1$$, which are known to strongly influence prognosis. We encoded the presence of mutations as binary variables, and also calculated aggregate features, such as:

1. **High-risk mutation indicator:**

$$
HIGH \textunderscore RISK \textunderscore MUT = \max(TP53 \textunderscore MUT, RUNX1 \textunderscore MUT, ASXL1 \textunderscore MUT)
$$

This captures whether a patient carries any of the mutations known to significantly impact survival.

2. **Epigenetic mutation load:** 

$$
EPIGENETIC \textunderscore MUT = ASXL1 \textunderscore MUT + TET2 \textunderscore MUT + DNMT3A \textunderscore MUT
$$

Mutations in these genes reflect clonal hematopoiesis and epigenetic dysregulation, which are associated with disease progression.  

3. **Interaction terms among high-risk genes:**  

$$
TP53 \textunderscore RUNX1 = TP53 \textunderscore MUT \times RUNX1 \textunderscore MUT, \quad TP53 \textunderscore ASXL1 = TP53 \textunderscore MUT \times ASXL1 \textunderscore MUT, \quad RUNX1 \textunderscore ASXL1 = RUNX1 \textunderscore MUT \times ASXL1 \textunderscore MUT
$$

These interactions capture potential synergistic effects between mutations that may exacerbate disease severity.

Cytogenetic abnormalities were also included, as they provide crucial prognostic information. We encoded presence of specific chromosomal aberrations and created a composite cytogenetic score:

$$
CYTO \textunderscore SCORE = MONOSOMY \textunderscore 7 + TRISOMY \textunderscore 8 + COMPLEX \textunderscore KARYO
$$

where $$COMPLEX \textunderscore KARYO$$ indicates three or more abnormalities. This score captures the overall chromosomal instability of each patient, which is strongly predictive of poor outcomes.

Finally, we incorporated variant allele frequency (VAF) data to quantify the clonal prevalence of mutations:

$$
VAF \textunderscore MAX \textunderscore log = \log(1 + VAF \textunderscore MAX), \quad VAF \textunderscore MEAN \textunderscore log = \log(1 + VAF \textunderscore MEAN), \quad VAF \textunderscore SUM \textunderscore log = \log(1 + VAF \textunderscore SUM)
$$

VAF statistics provide information about the size of the clone carrying the mutation, offering additional prognostic power beyond binary mutation presence.

#### Correlation Filtering

To prevent multicollinearity and enhance model stability, we computed the absolute correlation matrix among all features and removed one feature from each highly correlated pair (|corr| > 0.9). This step ensures that redundant information does not inflate the variance of model coefficients or bias feature importance estimates.

#### Modeling Rationale

We applied two complementary survival modeling approaches:

1. **Coxnet (penalized Cox regression):**  

The Coxnet model maximizes the penalized partial likelihood:

$$
\mathcal{L}(\beta) = \sum_{i=1}^n \delta_i \left( x_i^\top \beta - \log \sum_{j \in R_i} e^{x_j^\top \beta} \right) - \lambda \left( \alpha \|\beta\|_1 + \frac{1-\alpha}{2} \|\beta\|_2^2 \right)
$$

where $$\delta_i$$ is the event indicator, $$R_i$$ the risk set, $$\lambda$$ the regularization parameter, and $$\alpha$$ the L1/L2 mixing ratio. The L1 penalty encourages sparsity, facilitating feature selection, while the L2 penalty stabilizes coefficient estimates. Hyperparameters were optimized via 3-fold cross-validation using the concordance index as the scoring metric:

$$
C\text{-index} = \frac{\text{# concordant pairs}}{\text{# comparable pairs}}
$$

2. **Random Survival Forests (RSF):**  

RSF is a non-parametric ensemble method that constructs multiple survival trees on bootstrapped samples and aggregates predictions to estimate cumulative hazard functions. This method captures non-linear relationships and complex interactions that may not be well-modeled by linear regression. Each tree splits the data to maximize differences in survival times, making RSF highly flexible for heterogeneous patient populations.

#### Findings

The Coxnet model achieved a concordance index of approximately 0.72 on the training set and 0.70 on the test set, indicating strong predictive performance. RSF produced comparable results with slightly higher concordance, confirming the robustness of the selected features. Feature importance analysis consistently highlighted high-risk mutations (particularly $$TP53$$$ and $$RUNX1$$$), cytogenetic abnormalities, and derived ratio features (e.g., $$WBC \textunderscore to \textunderscore HB$$, $$BLAST \textunderscore to \textunderscore WBC as the strongest predictors. Polynomial and interaction terms improved model flexibility, enabling the capture of non-linear effects between biochemical markers and survival. Overall, the study demonstrates the value of integrating clinical, molecular, and derived features in survival prediction using both parametric and non-parametric machine learning models.
