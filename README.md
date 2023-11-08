# Multinomial Logistic Growth for SARS-CoV-2 lineages

To estimate growth rates, we used a multinomial logistic regression model of global lineage frequency data. GISAID sequences were obtained from the bulk .fasta download (dated 2023-10-02) and processed with Nextclade (v2.14) using the BA.2 reference (`sars-cov-2-21L`). Using Nextclade quality metrics, only sequences with >90% coverage and an overall QC status of "good" were retained. Since outlier dates could distort model estimates, we required sequences to have a fully specified deposition and collection date (ie. year, month, and day), and to have a collection date within 150 days of deposition (primarily to avoid collection dates where the year was incorrectly annotated). Finally, for each lineage, we excluded sequences with dates that were extreme outliers for that specific lineage, falling outside of 3.5 times the interquartile range of the median. For the analysis, we retained countries with > 500 sequences, and lineages with >200 sequences. Counts for each lineage were aggregated per country, per day.

Lineage counts were modeled using multinomial logistic regression, with global per-lineage growth rates (ie. shared between all countries), and per-lineage per-country intercepts. The model was implemented in the Julia language and the likelihood maximized using gradient descent, using Flux.jl and CUDA.jl to allow for GPU computation. Visualizations of the model fit are available [here](https://github.com/MurrellGroup/MultinomialLogisticGrowth/tree/main/plots)

Growth rates are interpreted such that the ratio of the frequency of two lineages $l_i$ and $l_j$ changes with time $t$:
$$\frac{l_i(t)}{l_j(t)}=K \cdot e^{t(g_i-g_j)}$$
where $g_i$ and $g_j$ are the growth rates for each lineage, $t$ is time in years, and $K$ is a constant determined by their intercepts.

