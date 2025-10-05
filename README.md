### Methodology and Findings

In this study, we aimed to predict patient survival outcomes by integrating clinical and molecular data through a combination of feature engineering, preprocessing, and advanced survival modeling techniques. Our methodology consisted of several stages: data preprocessing, feature transformation, correlation filtering, and model training using both penalized Cox regression (Coxnet) and Random Survival Forests (RSF).

#### Feature Engineering and Preprocessing

We began by preprocessing the clinical dataset to create meaningful derived variables. For key biochemical features such as $$BM\_BLAST$$ (bone marrow blasts percentage), $$HB$$ (hemoglobin), $$PLT$$ (platelets), $$WBC$$ (white blood cell count), $$ANC$$ (absolute neutrophil count), and $$MONOCYTES$$, we generated polynomial features of degree 3 (interaction-only, without bias) to capture non-linear relationships between these biomarkers and survival. Additionally, we created logarithmic, square root, and squared transformations to normalize skewed distributions and highlight potential non-linear effects:

$$
x_{log} = \log(1 + x), \quad x_{sqrt} = \sqrt{x}, \quad x_{square} = x^2
$$

We also derived biologically relevant ratios that capture interactions between cell populations, such as:

$$
WBC\_to\_HB = \frac{WBC}{HB + \epsilon}, \quad PLT\_to\_ANC = \frac{PLT}{ANC + \epsilon}, \quad BLAST\_to\_WBC = \frac{BM\_BLAST}{WBC + \epsilon}
$$

$$
ANC\_to\_MONOCYTES = \frac{ANC}{MONOCYTES + \epsilon}, \quad PLT\_to\_HB = \frac{PLT}{HB + \epsilon}, \quad BLAST\_to\_ANC = \frac{BM\_BLAST}{ANC + \epsilon}, \quad BLAST\_to\_PLT = \frac{BM\_BLAST}{PLT + \epsilon}, \quad MONO\_to\_WBC = \frac{MONOCYTES}{WBC + \epsilon}
$$

where $$\epsilon = 1e-5$$ prevents division by zero. These ratios were then log-transformed to reduce skewness and improve model stability.

From the molecular perspective, we incorporated mutation data from high-risk genes (e.g., $$TP53, RUNX1, ASXL1$$) and epigenetic regulators ($$ASXL1, TET2, DNMT3A$$). For each patient, we calculated:

1. **Mutation counts:** presence of mutations per gene.
2. **High-risk mutations:** a binary feature indicating whether any high-risk mutation was present:

$$
HIGH\_RISK\_MUT = \max(TP53\_MUT, RUNX1\_MUT, ASXL1\_MUT)
$$

3. **Epigenetic mutation load:** sum of mutations in key epigenetic regulators:

$$
EPIGENETIC\_MUT = ASXL1\_MUT + TET2\_MUT + DNMT3A\_MUT
$$

We also created interaction features among high-risk mutations to capture potential synergistic effects, such as:


Cytogenetic abnormalities were incorporated via indicator variables for key lesions, including $$-7$$, $$+8$$, $$del5q$$, $$t(8;21)$$, and $$inv(16)$$. We also computed a composite cytogenetic score:

$$
CYTO\_SCORE = MONOSOMY\_7 + TRISOMY\_8 + COMPLEX\_KARYO
$$

to summarize overall chromosomal risk.

Finally, we included variant allele frequency (VAF) statistics (max, mean, sum) per patient to quantify the burden of mutations and their clonal prevalence. Log transformations were applied to reduce skewness:

$$
VAF\_MAX\_log = \log(1 + VAF\_MAX), \quad VAF\_MEAN\_log = \log(1 + VAF\_MEAN), \quad VAF\_SUM\_log = \log(1 + VAF\_SUM)
$$

#### Correlation Filtering

To reduce multicollinearity and improve model stability, we computed the absolute correlation matrix and removed one feature from each highly correlated pair (|corr| > 0.9). This step preserved interpretability and ensured that features contributing redundant information were eliminated, reducing noise in the penalized regression model.

#### Modeling

For survival prediction, we employed two complementary methods:

1. **Coxnet (penalized Cox proportional hazards model):**  

The Coxnet model maximizes the penalized partial likelihood:

$$
\mathcal{L}(\beta) = \sum_{i=1}^n \delta_i \left( x_i^\top \beta - \log \sum_{j \in R_i} e^{x_j^\top \beta} \right) - \lambda \left( \alpha \|\beta\|_1 + \frac{1-\alpha}{2} \|\beta\|_2^2 \right)
$$

where $$\delta_i$$ is the event indicator, $$R_i$$ the risk set, $$\lambda$$ the penalty strength, and $$\alpha$$ the L1/L2 mixing ratio. Hyperparameters were optimized via 3-fold cross-validation, using the **concordance index** as the scoring metric:

$$
C\text{-index} = \frac{\text{# concordant pairs}}{\text{# comparable pairs}}
$$

2. **Random Survival Forests (RSF):**  

RSF, a non-parametric ensemble method, constructs multiple survival trees using bootstrapped samples and aggregates predictions to estimate cumulative hazard functions. This approach captures complex, non-linear interactions and higher-order feature effects that Coxnet may miss.

#### Findings

The Coxnet model achieved a concordance index of approximately 0.72 on the training set and 0.70 on the test set, while RSF showed slightly higher concordance, demonstrating robustness of our features across parametric and non-parametric methods. Feature importance analysis consistently highlighted high-risk mutations (especially $$TP53$$ and $$RUNX1$$), cytogenetic abnormalities, and derived ratios like $$WBC\_to\_HB$$ and $$BLAST\_to\_WBC$$ as key predictors. Polynomial and interaction terms improved the capture of non-linear relationships between biochemical markers and survival, illustrating the value of integrating clinical, molecular, and derived features in survival prediction models.
