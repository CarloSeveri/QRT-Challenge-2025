### Methodology and Findings

In this study, we aimed to predict survival outcomes using a combination of classical preprocessing, feature engineering, and machine learning models tailored for survival analysis. The workflow started with extensive preprocessing of clinical and molecular data, including the creation of polynomial features for key biochemical variables ($$BM\_BLAST, HB, PLT, WBC, ANC, MONOCYTES$$), as well as derived ratios such as  

$$
WBC\_to\_HB = \frac{WBC}{HB + \epsilon}, \quad PLT\_to\_ANC = \frac{PLT}{ANC + \epsilon}, \quad BLAST\_to\_WBC = \frac{BM\_BLAST}{WBC + \epsilon}
$$  

where $$\epsilon = 1e-5$$ to prevent division by zero. These transformations aimed to capture non-linear relationships and interactions that might be predictive of survival outcomes.

Next, we removed highly correlated features (|corr| > 0.9) to reduce redundancy and multicollinearity, ensuring model stability. For survival modeling, we employed **Coxnet** (a penalized Cox proportional hazards model) and **Random Survival Forests (RSF)**. The Coxnet model solves the following penalized partial likelihood:  

$$
\mathcal{L}(\beta) = \sum_{i=1}^n \delta_i \left( x_i^\top \beta - \log \sum_{j \in R_i} e^{x_j^\top \beta} \right) - \lambda \left( \alpha \|\beta\|_1 + \frac{1-\alpha}{2} \|\beta\|_2^2 \right)
$$  

where $$\delta_i$$ is the event indicator, $$R_i$$ is the risk set at time $$t_i$$, $$\lambda$$ is the regularization strength, and $$\alpha$$ controls the balance between L1 (lasso) and L2 (ridge) penalties. Grid search with cross-validation was used to select optimal hyperparameters ($$\alpha, \lambda_{min}$$) by maximizing the concordance index:  

$$
C\text{-index} = \frac{\text{# concordant pairs}}{\text{# comparable pairs}}
$$  

RSF, a non-parametric ensemble tree-based method, was used in parallel to capture non-linear interactions and complex survival patterns. Trees were built using bootstrapped samples, and predictions were aggregated to estimate the cumulative hazard function.  

The rationale behind these choices was twofold: Coxnet provides interpretability and sparsity, allowing us to identify key molecular and clinical drivers of survival, whereas RSF captures complex non-linear interactions that the Cox model may miss. Our findings showed that penalized Cox models achieved a concordance index of approximately 0.72 on the training set and 0.70 on the test set, while RSF achieved similar or slightly higher concordance. Feature importance analysis revealed that high-risk mutations (e.g., $$TP53, RUNX1$$) and derived ratios (e.g., $$WBC\_to\_HB$$) were consistently among the most predictive features, highlighting the interplay between clinical and molecular factors in determining patient survival.
