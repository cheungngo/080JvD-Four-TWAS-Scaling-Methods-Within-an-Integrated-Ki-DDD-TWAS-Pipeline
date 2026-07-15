# **Systematic Comparison of Four TWAS Scaling Methods Within an Integrated Ki–DDD–TWAS Pipeline for Antipsychotic Metabolic Risk Assessment**

## **Abstract**

Antipsychotic-associated metabolic abnormalities, including dyslipidemia, weight gain, hyperglycemia, and metabolic syndrome, vary substantially between drugs. We developed a reproducible computational pipeline that integrates receptor-binding affinity, defined daily dose, and precomputed transcriptome-wide association study data to produce drug-level metabolic-risk scores for antipsychotics. The pipeline retrieves drugs classified under ATC group N05A, normalizes defined daily doses to milligrams, fuzzy-matches drug names to a Ki database, aggregates receptor-binding measurements, maps receptors to genes, and integrates five lipid-related traits: low-density lipoprotein cholesterol, high-density lipoprotein cholesterol, log-transformed triglycerides, non-high-density lipoprotein cholesterol, and total cholesterol. Four transcriptome-wide association study scaling strategies were compared: a baseline without per-receptor scaling, linear scaling by absolute z score, exponential scaling by absolute z score, and exponential scaling by negative log10 p value. The final comparison included 52 drugs, five traits, and 1,040 drug–trait–method risk-score observations. Drug rankings were highly concordant across methods, with a mean Spearman correlation of 0.997. Only 38 of 260 drug–trait combinations showed a rank shift of at least three positions. Clozapine ranked first for every trait and method, while HRH1 and HTR2A were the two most consistent receptor-level contributors. The results indicate that the principal effect of changing the transcriptome-wide association study transformation was to alter score magnitude rather than drug prioritization. The negative-log10-p method is suitable as a primary analysis when evidence-weighted scaling is desired, while the baseline method is useful as a sensitivity analysis.

**Keywords:** Antipsychotics; metabolic risk; transcriptome-wide association study; receptor pharmacology; Ki affinity; dose normalization; reproducible computing

## **1\. Introduction**

### **1.1 Clinical and scientific need**

Antipsychotic drugs remain essential for the treatment of schizophrenia and other severe mental disorders, but their use is associated with clinically important metabolic complications. These complications include weight gain, glucose dysregulation, dyslipidemia, obesity, metabolic syndrome, and increased cardiovascular risk \[1-4\]. The magnitude and pattern of these adverse effects differ across individual drugs, and the risk may be influenced by treatment duration, baseline characteristics, dose, pharmacological profile, and patient susceptibility \[1-3,5\].

Comparative evidence also indicates that antipsychotic metabolic effects are not uniform. In a short-term network meta-analysis, clozapine and olanzapine were among the drugs with the least favorable metabolic profiles, whereas aripiprazole, brexpiprazole, cariprazine, lurasidone, and ziprasidone were comparatively more benign \[6\]. In a separate mid- to long-term network meta-analysis, chlorpromazine, clozapine, olanzapine, and zotepine were among the drugs associated with substantial weight gain, although the confidence in some comparisons was limited and estimates varied across outcomes \[7\]. These differences illustrate why drug-level prioritization requires a method that can integrate several types of evidence without treating a single pharmacological or clinical measure as a complete representation of metabolic risk.

Most existing approaches to antipsychotic metabolic liability are based either on observed clinical outcomes or on receptor-binding profiles. Receptor-based approaches are biologically informative because several neurotransmitter systems have been implicated in appetite regulation, glucose metabolism, lipid metabolism, and weight gain. However, receptor affinity alone does not capture genetically regulated variation in the expression of receptor genes or related biological pathways. Conversely, genome-wide association studies and transcriptome-wide association studies provide information about genetic associations but do not directly describe the binding profile, dose, or pharmacological exposure of a drug.

The present work combines these complementary data types into a single prioritization framework. The resulting score should be interpreted as a computational metabolic-risk index rather than a clinical risk estimate, incidence rate, odds ratio, or probability of an adverse event. The analysis is designed to compare drugs within a common framework and to test whether the principal conclusions remain stable when the way in which transcriptome-wide association evidence is incorporated is changed.

### **1.2 Limitations of current methods**

Pure receptor-affinity scores generally treat receptor binding as the primary determinant of metabolic liability. Such scores can be useful for mechanistic comparison, but they may omit genetically regulated expression evidence and trait-specific molecular associations. They may also be sensitive to the choice of Ki aggregation method, the representation of dose, the set of included receptors, and the weights assigned to different receptor systems.

Transcriptome-wide association studies address a different question. They use genetic information and reference transcriptome data to estimate the association between genetically regulated expression and complex traits. PrediXcan, TWAS, and summary-statistic extensions have made it possible to prioritize genes associated with disease or quantitative traits using individual-level or GWAS summary data \[8-10\]. These methods have increased the biological interpretability of GWAS results, but they do not, by themselves, provide a drug-level pharmacological score. They also do not establish causality automatically, because linkage disequilibrium, co-regulation, pleiotropy, tissue mismatch, and prediction uncertainty can generate associations that require further investigation \[11-13\].

A further limitation is that transcriptome-wide association evidence is often integrated using a single predefined transformation. This makes it difficult to determine whether the final drug ranking is robust to alternative ways of weighting evidence. Absolute z-score weighting emphasizes standardized association magnitude, whereas p-value weighting emphasizes statistical evidence. These measures are related but are not identical. A p value is influenced by both effect size and study precision, while a z score retains direction but does not directly encode the same evidence scale. The relative importance of these properties depends on the scientific question.

There is also a practical reproducibility gap. An end-to-end drug-level analysis requires several data transformations: drug-name matching, receptor mapping, dose conversion, Ki parsing, missing-value handling, transcriptome-wide association lookup, score normalization, and ranking. If these steps are not documented explicitly, it can be difficult to determine whether observed differences reflect biological evidence, a preprocessing decision, or an implementation artifact.

### **1.3 Rationale for the current pipeline**

The pipeline was designed as a modular integration framework with four main components. First, it extracts antipsychotic drugs from ATC group N05A and standardizes available defined daily dose information. Second, it integrates receptor-binding measurements from a Ki database and converts them into affinity–dose contributions. Third, it maps receptor labels to official gene symbols and retrieves precomputed transcriptome-wide association results for five lipid-related traits. Fourth, it compares four transcriptome-wide association scaling strategies using the same drug, receptor, and trait data.

The lipid traits were selected because lipid abnormalities are clinically relevant components of antipsychotic-associated metabolic dysfunction and because LDL cholesterol, HDL cholesterol, triglycerides, non-HDL cholesterol, and total cholesterol represent related but nonidentical aspects of lipid biology. Major lipid traits are also associated with vascular disease and cardiovascular risk \[14,15\].

The downstream comparison framework was designed to distinguish changes in score magnitude from changes in drug ordering \[19\]. It therefore includes rank correlations, Pearson correlations, top-five and top-ten Jaccard overlap, rank-shift analysis, per-drug amplification ratios, receptor-contribution recomputation, consensus ranking, leader-separation metrics, and Friedman tests. This combination is important because a method can produce significantly different numerical scores while preserving essentially the same ranking.

### **1.4 Objectives**

The first objective was to describe an end-to-end Ki–defined daily dose–transcriptome-wide association study pipeline in sufficient detail for independent implementation. The second objective was to compare four transcriptome-wide association scaling strategies using identical input data and downstream processing. The third objective was to determine whether the principal drug rankings and receptor-level contributions were stable across methods. The final objective was to identify a practical primary method and a complementary sensitivity analysis while clearly separating computational prioritization from clinical causality.

## **2\. Materials and Methods**

### **2.1 Data sources**

The antipsychotic drug list was generated from the World Health Organization ATC/DDD Index using the N05A parent code and its child subgroups \[16\]. The scraper identified the N05AA, N05AB, N05AC, N05AD, N05AE, N05AF, N05AG, N05AH, N05AL, N05AN, and N05AX subgroups, then parsed drug-level rows containing ATC code, name, defined daily dose, unit, route, and note fields. The code preserves N05AN, which includes lithium in the ATC hierarchy, and does not expect reserpine because it is classified outside N05A in the source taxonomy. Defined daily dose data were used as a standardized dose proxy rather than as a measurement of an individual patient’s prescribed dose or plasma exposure.

The Ki data were read from an internal CSV file identified in the supplied code as KiDatabase\_2026-06-22.csv. The Ki database was not accompanied by a public DOI or a complete data citation in the supplied materials. Accordingly, the database version should be cited as an internal or deposited dataset in the final submission, subject to licensing and redistribution permissions. Ki values were interpreted as nanomolar affinity measurements. The parser accepted numeric values and strings containing qualifiers such as greater-than, less-than, approximate, or equal-to signs. Commas used as thousands separators were removed before numeric conversion.

Five trait-specific transcriptome-wide association directories were used: LDL, HDL, logTG, non-HDL, and TC. The pipeline consumed precomputed transcriptome-wide association files rather than training expression prediction models or recalculating gene–trait associations from genotype data. This distinction is important: the present method is a downstream integration and scoring framework built on transcriptome-wide association outputs, not a new transcriptome-wide association model.

The transcriptome-wide association reader was designed to tolerate common file-format differences. It attempted comma-separated, tab-separated, and automatically detected delimiters and tried UTF-8, UTF-8 with byte-order mark, Latin-1, and CP1252 encodings. Candidate columns were defined for gene symbol, Ensembl identifier, z score, and p value. When an Ensembl identifier was used, version suffixes were removed. When multiple records for the same gene were found, the record with the smallest available p value was retained. Files with missing gene or z-score columns were skipped and reported.

Receptor labels were mapped to official gene symbols using a manually specified dictionary. The dictionary included serotonin receptors such as 5-HT1A, 5-HT2A, 5-HT2C, and 5-HT7; dopamine receptors D1–D5; adrenergic receptors; histamine receptors H1–H4; muscarinic receptors M1–M5; sigma receptor 1; and the serotonin, dopamine, and noradrenaline transporters. Looser matching removed spaces and hyphens to accommodate minor label differences.

The predefined metabolic-weight dictionary assigned the largest weights to HRH1, HTR2C, and CHRM3. The complete weight set was HRH1, 1.00; HTR2C, 0.92; CHRM3, 0.85; HTR2A, 0.55; ADRA1A and ADRA1B, 0.48 each; HTR6, 0.42; ADRA2A, 0.32; ADRA2B and ADRA2C, 0.30 each; CHRM1, 0.28; HTR7, 0.25; CHRM4, 0.22; CHRM5, 0.20; DRD3, 0.18; DRD2, 0.12; DRD4, 0.10; HTR1A and ADRB1, 0.08 each; ADRB2, 0.07; SLC6A4, 0.06; SIGMAR1 and SLC6A2, 0.05 each; and SLC6A3, 0.04. These are model weights and should not be interpreted as experimentally estimated causal coefficients.

The transcriptome reference context was based on the Genotype-Tissue Expression Project \[17\]. Population reference panels such as 1000 Genomes were considered for linkage-disequilibrium calculations \[18\].

### **2.2 Pipeline architecture**

The pipeline first fetched the N05A index page and discovered immediate child subgroup codes. Each subgroup page was then retrieved with a short delay between requests. The resulting rows were written to a CSV file containing ATC code, name, defined daily dose, unit, route, note, and subgroup. A conventional Python implementation would need to replace the notebook shell commands and parameterize the data paths before redistribution.

Defined daily doses were converted to milligrams using a simple unit-conversion function. Values in grams were multiplied by 1,000, values in milligrams were left unchanged, and microgram values were multiplied by 0.001. If the unit was missing or unrecognized, the code assumed milligrams. This assumption increases robustness to incomplete source files but may introduce error if an unrecognized unit is not actually milligrams.

Drug names were normalized by lowercasing, trimming leading and trailing whitespace, and collapsing repeated spaces. Each normalized name was matched to a normalized Ki-database ligand vocabulary using the RapidFuzz WRatio scorer. Matches with a score below 80 were rejected. For rejected drugs, the pipeline retained a row containing the input name and matching score but no receptor-level affinity data. The use of a fuzzy threshold prevents silent loss of unmatched drugs but does not replace manual inspection of borderline matches.

For accepted matches, all Ki-database records associated with the selected ligand were grouped by receptor. The default aggregation rule was the minimum observed Ki value, representing the strongest measured binding for a drug–receptor pair. Alternative aggregation rules, including the median and arithmetic mean, were supported by the code but were not used in the reported comparison. The inverse affinity was calculated as inv\_Ki \= 1 divided by Ki in nanomolar units.

The affinity–dose contribution was calculated as affinity\_dose\_score \= inv\_Ki multiplied by log of DDD\_mg plus 1\. If defined daily dose was missing, the code used inverse Ki alone as a fallback. Before receptor weighting, the pipeline applied a log1p transformation to the nonnegative affinity–dose score. For recognized metabolic genes, the transformed affinity contribution was multiplied by the corresponding receptor weight. Contributions were then summed within each drug.

Only drugs with at least one valid defined daily dose were retained for downstream scoring. Rows with an affinity–dose score equal to zero were removed. Missing Ki values could be imputed using a receptor-specific mean or median, with a global mean or median as fallback when a receptor-specific estimate was unavailable. The supplied runs used mean imputation. After Ki imputation, inverse Ki, affinity–dose score, and strong-binder status were recalculated.

Transcriptome-wide association values were optionally imputed. The preprocessing configuration used a global mean for missing z scores and p values when the receptor mapped to a recognized gene. Imputed rows were retained but were deliberately not classified as transcriptome-wide association significant. During risk-score calculation, a second fairness procedure was applied to missing metabolic-gene z scores or p values. Missing metabolic-gene z scores were assigned the mean absolute z score observed among available metabolic-receptor rows. For p-value-based scaling, missing p values were represented using the mean negative log10 p value among available metabolic-receptor rows. These two imputation stages should be distinguished because they affect both the stored data and the method-specific scaling calculation.

The strong-binder indicator was defined as Ki below 10 nM. Transcriptome-wide association significance was defined in the code as a nominal p value below 0.05. This threshold was used for the significance bonus and was not a genome-wide Bonferroni- or false-discovery-rate-adjusted threshold. Consequently, the twas\_significant variable should not be interpreted as a definitive genome-wide association claim.

### **2.3 Four TWAS scaling methods**

The baseline method, labelled Original, applied no per-receptor transcriptome-wide association scaling. However, the actual version 0 code retained a drug-level transcriptome-wide association multiplier based on the maximum absolute z score and a nominal significance bonus. Therefore, Original should be described precisely as “no per-receptor TWAS scaling” rather than as complete absence of TWAS influence.

The Linear method multiplied each receptor contribution by 1 plus 0.24 multiplied by the absolute transcriptome-wide association z score. This transformation is monotonic in absolute z score and gives increasing influence to stronger standardized associations. The method also retained a secondary drug-level multiplier of 1 plus 0.04 multiplied by the maximum absolute z score.

The Expo\_z method multiplied each receptor contribution by the exponential of 0.08 multiplied by the absolute z score. This gives stronger amplification to high absolute z scores than the Linear method while remaining controlled for the z-score ranges present in the data. A secondary drug-level multiplier of 1 plus 0.02 multiplied by the maximum absolute z score was retained.

The Logp method multiplied each receptor contribution by the exponential of 0.045 multiplied by negative log10 of the transcriptome-wide association p value. P values were lower-clipped at 1 × 10 to the minus 300 power before transformation, and negative log10 p values were capped at 25 to prevent extreme numerical amplification. The secondary drug-level multiplier was 1 plus 0.02 multiplied by the maximum absolute z score.

The p-value-based transformation is an evidence-weighting choice rather than a canonical transcriptome-wide association procedure. It emphasizes statistical strength but discards direction and is influenced by sample size and precision. In contrast, absolute z-score scaling retains the magnitude of the standardized association but also discards direction because the absolute value is used. Both choices therefore encode an assumption that positive and negative transcriptome-wide association signals should contribute similarly to the risk score. The underlying TWAS literature provides the gene-level association framework but does not establish that one of these downstream transformations is universally optimal \[8-11\].

All four scoring implementations used the log1p-transformed affinity–dose contribution, the predefined receptor weights, normalization to chlorpromazine, and a 0.10 multiplier per nominally significant metabolic receptor. After calculation, the chlorpromazine score was set to 100\. If chlorpromazine was absent or had a zero score, the code normalized to the maximum available drug score instead. The resulting values are relative model indices, not clinical relative risks.

### **2.4 Downstream comparison framework**

The downstream comparison used only the generated output files. For each method and trait, the drug-level risk-score file was read and converted into long and wide formats containing drug, trait, score, rank, and method. The comparison included the four methods and five common traits.

Between-method agreement was assessed using Pearson correlation and Spearman rank correlation \[19,20\]. Spearman correlation was emphasized because the primary scientific question concerns whether drugs retain their ordering when the scaling transformation changes \[20,21\]. Top-five and top-ten sets were compared with the Jaccard index, calculated as the size of the intersection divided by the size of the union of two method-specific drug sets \[22\].

Rank sensitivity was quantified for each drug–trait combination as the maximum method-specific rank minus the minimum method-specific rank. A combination was flagged when this shift was at least three positions. This threshold was an operational sensitivity threshold rather than a clinical definition of instability.

Per-drug amplification was calculated as the score under a transformed method divided by the corresponding Original score. Ratios were summarized using the mean, median, and maximum. Because ratios can become unstable when the Original score is close to zero, amplification factors were interpreted together with absolute score differences.

Receptor contributions were recomputed from the long-format receptor-level output. The calculation applied the log1p affinity transformation, receptor weights, and the method-specific per-receptor scaling rule. Receptor totals were summed within gene and ranked within each trait. Mean receptor rank across traits was then used to identify stable receptor-level contributors. This analysis describes the receptor contribution component and does not reproduce every drug-level multiplier used in the final risk score.

Consensus ranking was generated by averaging each drug’s risk score and rank across all available method–trait combinations. Each of the 52 final drugs had 20 observations, corresponding to four methods and five traits. The consensus score was therefore a descriptive average across method and trait conditions.

The Friedman test was performed separately for each trait using drug as the repeated-measures block and method as the within-drug condition. This nonparametric test evaluates whether at least one method differs systematically in score distribution. It does not establish that the drug ranking has changed. For comparisons across multiple methods, rank-based tests and post hoc procedures are commonly recommended as robust alternatives to parametric analysis when distributional assumptions are uncertain \[23,24\].

Leader separation was defined as the highest drug score divided by the median score for a method and trait. The mean of this ratio across traits was used as a summary of how strongly a method separated the highest-ranked drug from the middle of the distribution. This measure was interpreted as a distributional property and not as evidence that a method was more clinically valid.

### **2.5 Implementation and reproducibility**

The analysis was implemented in a Google Colab-oriented Python workflow using pandas, NumPy, RapidFuzz, tqdm, SciPy, Matplotlib, Seaborn, and openpyxl. The supplied artifact contains several sequential versions of the scoring function: the baseline version, Linear scaling, exponential z-score scaling, and negative-log10-p scaling. It also contains notebook shell commands, repeated function definitions, hard-coded /content paths, and Google Drive copy commands. Consequently, the supplied artifact should be regarded as an analysis notebook export rather than as a directly distributable Python package.

For publication, the workflow should be separated into independent modules for data acquisition, preprocessing, scoring, method comparison, and reporting. Each method should be run from a clean environment with an explicit configuration file. The final repository should include a locked environment specification, input-file checksums, versioned parameter files, unit tests for score calculations, and a small example dataset. These changes would reduce the possibility that function redefinition, execution order, notebook state, or path configuration alters the result.

The original pipeline created one output directory for each trait and method. The downstream comparison created a method\_comparison\_v4 directory containing wide score and rank files, between-method correlations, top-N overlap, score-magnitude summaries, rank-shift files, leader-gap files, receptor-contribution summaries, amplification factors, consensus rankings, Friedman-test results, and figures. The supplied figures included method-correlation heatmaps, a leader-separation plot, and a consensus top-15 heatmap.

## **3\. Results**

### **3.1 Data coverage and processing summary**

The supplied downstream outputs contained five traits for each of the four methods: HDL, LDL, logTG, non-HDL, and TC. Each trait contained 52 drugs, producing 1,040 drug–trait–method risk-score observations and 260 unique drug–trait combinations. Each drug in the consensus file had 20 observations, corresponding to four methods across five traits.

The final risk-score matrix therefore had complete downstream coverage across the planned method and trait combinations. The materials supplied for this manuscript do not include the detailed per-run counts of initial N05A entries, accepted fuzzy matches, rejected fuzzy matches, Ki imputations, or TWAS imputations. Those values should be extracted from the individual pipeline\_summary.txt files and reported in the final version of the article rather than inferred from the final ranking files. The principal analysis configuration and summary results are shown in Table 1\.

**Table 1\. Overall analysis configuration and summary results**

| Item | Result |
| :---- | :---- |
| TWAS-scaling methods | 4: Original, Linear, Expo\_z, and Logp |
| Lipid traits | 5: HDL, LDL, logTG, non-HDL, and TC |
| Drugs in final comparison | 52 |
| Unique drug–trait combinations | 260 |
| Drug–trait–method observations | 1,040 |
| Reference drug | Chlorpromazine, normalized to 100 |
| Leader drug | Clozapine |
| Mean Spearman correlation across method pairs and traits | 0.997 |
| Drug–trait combinations with rank shift of at least 3 positions | 38 of 260 |
| Percentage with rank shift of at least 3 positions | 14.6% |
| Maximum reported rank shift | 6 positions |
| Maximum-shift example | Cariprazine for TC |
| Most consistent receptor-level contributors | HRH1 and HTR2A |

*Note:* The chlorpromazine value of 100 is a model normalization and is not a clinical relative risk, incidence rate, or probability.

### **3.2 Between-method agreement**

The four methods produced highly concordant drug rankings. The mean Spearman correlations across traits were 1.000 for Expo\_z versus Logp, 0.999 for Linear versus Expo\_z, 0.998 for Linear versus Logp, 0.995 for Original versus Linear, 0.995 for Original versus Logp, and 0.994 for Original versus Expo\_z. Averaged across all method pairs and traits, the mean Spearman correlation was 0.997.

The lowest individual Spearman correlation was 0.990 for Original versus Expo\_z in HDL. The remaining correlations were also very high, and Pearson correlations ranged from approximately 0.997 to 1.000. The results indicate that the methods differed slightly more when the baseline method was compared with a TWAS-scaled method than when the three TWAS-scaled methods were compared with one another.

Expo\_z and Logp were the most similar pair. Their mean Spearman correlation was 1.000 when rounded to three decimal places, and their trait-specific correlations were between 0.999 and 1.000. These findings support the conclusion that the two evidence-weighting transformations generated nearly interchangeable drug rankings in this dataset. The pairwise agreement results are summarized in Table 2\.

**Table 2\. Between-method rank agreement**

| Method pair | Mean Spearman correlation across traits | Lowest or reported trait-specific value |
| :---- | :---- | :---- |
| Original vs Linear | 0.995 | Not separately reported |
| Original vs Expo\_z | 0.994 | 0.990 for HDL |
| Original vs Logp | 0.995 | Not separately reported |
| Linear vs Expo\_z | 0.999 | Not separately reported |
| Linear vs Logp | 0.998 | Not separately reported |
| Expo\_z vs Logp | 1.000 | 0.999–1.000 across traits |

*Note:* Pearson correlations across comparisons ranged from approximately 0.997 to 1.000.

### **3.3 Ranking stability**

The top-five and top-ten analyses showed substantial overlap. For the top five drugs, the mean Jaccard index was 1.000 for Expo\_z versus Logp, 0.933 for Linear versus Expo\_z, 0.933 for Linear versus Logp, 0.933 for Original versus Expo\_z, 0.933 for Original versus Logp, and 0.867 for Original versus Linear. For the top ten drugs, the mean Jaccard index was 1.000 for Expo\_z versus Logp, Linear versus Expo\_z, and Linear versus Logp. Original shared approximately 0.927 of the top-ten set with each of the other methods.

Only 38 of 260 drug–trait combinations showed a rank shift of at least three positions. The largest shift was six positions for cariprazine in total cholesterol, whose rank changed from 21 under Original to 27 under Expo\_z. Other larger shifts occurred for acepromazine in HDL, bromperidol in HDL and total cholesterol, fluspirilene in LDL and non-HDL, and levomepromazine in HDL.

The most stable drugs had a maximum rank shift of zero across the five traits. These were brexpiprazole, chlorpromazine, clozapine, levosulpiride, and triflupromazine. Asenapine, iloperidone, fluphenazine, molindone, and olanzapine had maximum shifts of one. The drugs with larger shifts were generally located in the middle or lower part of the ranking, where relatively small score differences can change rank positions without changing the broad risk category. The top-N overlap results are summarized in Table 3\.

**Table 3\. Mean Jaccard overlap of top-ranked drug sets**

| Method pair | Mean top-five Jaccard index | Mean top-ten Jaccard index |
| :---- | :---- | :---- |
| Original vs Linear | 0.867 | 0.927 |
| Original vs Expo\_z | 0.933 | 0.927 |
| Original vs Logp | 0.933 | 0.927 |
| Linear vs Expo\_z | 0.933 | 1.000 |
| Linear vs Logp | 0.933 | 1.000 |
| Expo\_z vs Logp | 1.000 | 1.000 |

*Note:* Values are means across the five lipid traits.

### **3.4 Amplification effects**

Relative to Original, the Linear method had a mean amplification factor of 1.178, a median of 1.149, and a maximum of 1.677. Expo\_z had a mean of 1.209, a median of 1.189, and a maximum of 1.802. Logp had a mean of 1.203, a median of 1.191, and a maximum of 1.812.

These results indicate that the transformed methods generally increased score magnitudes by approximately 18% to 21% relative to the baseline. The largest amplification was observed under Logp, closely followed by Expo\_z. The most amplified rows included pipotiazine, sultopride, flupentixol, triflupromazine, levosulpiride, periciazine, levomepromazine, acepromazine, promazine, and acetophenazine, particularly for HDL.

Amplification ratios were highest for drugs with low Original scores. For example, a small absolute increase from a very low score can produce a large fold-change. This was evident for drugs such as levosulpiride and flupentixol. Therefore, amplification factors should not be interpreted in isolation, especially for scores near zero. The rank and absolute-score analyses provide a more stable representation of the impact of method choice. Amplification and leader-separation results are summarized in Table 4\.

**Table 4\. Score amplification relative to Original and leader separation**

| Method | Mean amplification factor | Median amplification factor | Maximum amplification factor | Mean top-to-median leader ratio |
| :---- | :---- | :---- | :---- | :---- |
| Original | 1.000 | 1.000 | 1.000 | 15.992 |
| Linear | 1.178 | 1.149 | 1.677 | 13.995 |
| Expo\_z | 1.209 | 1.189 | 1.802 | 13.337 |
| Logp | 1.203 | 1.191 | 1.812 | 13.272 |

*Note:* Amplification factors are calculated relative to the Original score for the same drug and trait. The Original values equal 1.000 by definition. Ratios can be large when the Original score is close to zero and should be interpreted together with absolute scores.

### **3.5 Receptor-level contributions**

The receptor-contribution analysis identified HRH1 and HTR2A as the first- and second-ranked contributors across all four methods. In the Original method, the mean receptor ranks were HRH1, 1.0; HTR2A, 2.0; HTR2C, 3.0; DRD2, 4.0; DRD3, 5.0; ADRA1A, 6.0; HTR7, 7.0; and ADRA1B, 8.0.

The Linear method produced a small reordering, with DRD2 at mean rank 3.4 and HTR2C at 3.8. Expo\_z placed HTR2C at 3.2 and DRD2 at 3.8. Logp returned a pattern close to Original, with HTR2C at 3.0 and DRD2 at 4.0. HRH1 and HTR2A were present in the top three under every method.

The stability of HRH1 is biologically plausible. Histamine H1 receptor affinity has been associated with short-term antipsychotic-related weight gain, and receptor-occupancy analyses have identified H1 blockade as an important contributor to weight gain and diabetes-related outcomes \[25,26\]. Broader pharmacological reviews also implicate H1, 5-HT2C, muscarinic, adrenergic, and dopaminergic mechanisms in antipsychotic-associated metabolic effects \[27-29\].

The model-derived position of HTR2A should be interpreted carefully. It reflects the combination of the predefined weight, drug-specific binding affinity, dose proxy, transcriptome-wide association values, and aggregation procedure. It does not establish that HTR2A is the second most important causal mechanism for metabolic disease across all antipsychotics.

Selected receptor-contribution ranks and the reported rank-shift summary are presented in Table 5\.

**Table 5\. Selected receptor-contribution ranks and rank-shift sensitivity**

**Panel A. Selected mean receptor ranks**

| Method | HRH1 | HTR2A | HTR2C | DRD2 |
| :---- | :---- | :---- | :---- | :---- |
| Original | 1.0 | 2.0 | 3.0 | 4.0 |
| Linear | 1.0 | 2.0 | 3.8 | 3.4 |
| Expo\_z | 1.0 | 2.0 | 3.2 | 3.8 |
| Logp | 1.0 | 2.0 | 3.0 | 4.0 |

**Panel B. Rank-shift summary**

| Metric | Result |
| :---- | :---- |
| Total drug–trait combinations | 260 |
| Combinations with rank shift of at least 3 positions | 38 |
| Largest reported shift | 6 positions |
| Largest-shift example | Cariprazine for TC; rank 21 under Original and 27 under Expo\_z |
| Drugs with maximum shift of 0 | Brexpiprazole, chlorpromazine, clozapine, levosulpiride, and triflupromazine |
| Drugs with maximum shift of 1 | Asenapine, iloperidone, fluphenazine, molindone, and olanzapine |
| Other reported larger-shift examples | Acepromazine for HDL; bromperidol for HDL and TC; fluspirilene for LDL and non-HDL; levomepromazine for HDL |

*Note:* Panel A reports selected receptor ranks explicitly available in the supplied summary. Panel B reports the complete count and the larger-shift examples identified in the supplied results narrative; the full underlying per-drug rank-shift file was not included in the supplied appendix text.

### **3.6 Clozapine and consensus high-risk drugs**

Clozapine ranked first for every trait and every method. Relative to chlorpromazine, which was fixed at 100, clozapine scores under Original were 128.8 for HDL, 115.9 for LDL, 121.7 for logTG, 115.9 for non-HDL, and 123.2 for total cholesterol. The corresponding Linear scores were 131.5, 114.2, 120.0, 114.2, and 123.6. Expo\_z scores were 130.2, 115.3, 120.9, 115.3, and 123.5. Logp scores were 130.1, 115.9, 121.2, 115.9, and 123.9.

The method-specific differences for clozapine were modest and were not uniformly positive. Relative to Original, the mean score difference was approximately −0.4 for Linear, −0.1 for Expo\_z, and \+0.3 for Logp. The largest positive difference was 2.6 points for HDL under Linear. Thus, the high clozapine score was not created primarily by the choice of TWAS transformation.

The consensus ranking averaged scores and ranks across all four methods and five traits. Clozapine had a mean risk score of 121.1 and a mean rank of 1.00. Chlorpromazine had a mean risk score of 100.0 and a mean rank of 2.00. The next six drugs were asenapine, mean risk 83.9 and mean rank 3.05; zotepine, 75.8 and 4.45; risperidone, 74.6 and 4.70; ziprasidone, 71.8 and 5.80; olanzapine, 66.8 and 7.10; and thioridazine, 62.9 and 7.90.

The remainder of the first 20 consensus positions were chlorprothixene, sertindole, iloperidone, brexpiprazole, fluphenazine, quetiapine, mesoridazine, perphenazine, haloperidol, loxapine, aripiprazole, and levomepromazine. Their mean ranks ranged from 9.25 for chlorprothixene to 20.25 for levomepromazine.

The consensus placement of clozapine is consistent with extensive clinical evidence that clozapine carries substantial metabolic liability, although the model should not be interpreted as a clinical outcome validation \[2,1,30\]. The reference value of 100 for chlorpromazine is a normalization choice and does not represent a clinical incidence rate or relative risk estimate. The strongest leader separation was observed under Original, with an average top-to-median ratio of 15.992, compared with 13.995 for Linear, 13.337 for Expo\_z, and 13.272 for Logp.

Clozapine scores, consensus leaders, and the corresponding statistical comparison results are shown in Table 6\.

**Table 6\. Clozapine scores, consensus ranking, and Friedman statistics**

**Panel A. Clozapine score by method and trait**

| Method | HDL | LDL | logTG | non-HDL | TC |
| :---- | :---- | :---- | :---- | :---- | :---- |
| Original | 128.8 | 115.9 | 121.7 | 115.9 | 123.2 |
| Linear | 131.5 | 114.2 | 120.0 | 114.2 | 123.6 |
| Expo\_z | 130.2 | 115.3 | 120.9 | 115.3 | 123.5 |
| Logp | 130.1 | 115.9 | 121.2 | 115.9 | 123.9 |

**Panel B. Leading consensus ranking**

| Consensus rank | Drug | Mean risk score | Mean rank |
| :---- | :---- | :---- | :---- |
| 1 | Clozapine | 121.1 | 1.00 |
| 2 | Chlorpromazine | 100.0 | 2.00 |
| 3 | Asenapine | 83.9 | 3.05 |
| 4 | Zotepine | 75.8 | 4.45 |
| 5 | Risperidone | 74.6 | 4.70 |
| 6 | Ziprasidone | 71.8 | 5.80 |
| 7 | Olanzapine | 66.8 | 7.10 |
| 8 | Thioridazine | 62.9 | 7.90 |

**Panel C. Friedman tests across methods**

| Trait | Friedman statistic | Reported p value |
| :---- | :---- | :---- |
| HDL | 51.424 | 0.0000 |
| LDL | 48.953 | 0.0000 |
| logTG | 54.365 | 0.0000 |
| non-HDL | 50.882 | 0.0000 |
| TC | 44.576 | 0.0000 |

*Note:* The displayed p values were rounded to four decimal places in the supplied output; exact p values cannot be recovered from the supplied summary.

### **3.7 Statistical comparison**

The Friedman test was significant for every trait. The reported statistics were 51.424 for HDL, 48.953 for LDL, 54.365 for logTG, 50.882 for non-HDL, and 44.576 for total cholesterol. The output displayed p values of 0.0000 for all five tests; because the summary reports four decimal places, the exact p values cannot be recovered from the supplied file. The statistics are included in Table 6\.

These results show that the methods produced systematically different paired score magnitudes. They do not, by themselves, demonstrate substantial ranking changes. The rank correlations, top-N overlap, and rank-shift analysis provide the more direct evidence for the stability of drug prioritization. In this dataset, the combination of significant Friedman tests and very high rank concordance indicates that the transformations changed the scale and distribution of scores more than they changed the ordering of drugs.

## **4\. Discussion**

### **4.1 Methodological strengths**

The main strength of this work is that it makes the full sequence of transformations explicit. Drug names are obtained from an established classification system, dose units are normalized, Ki measurements are parsed and aggregated, receptor labels are mapped to genes, receptor contributions are weighted, transcriptome-wide association data are imported, and method-specific risk scores are normalized to a common reference. Each stage can therefore be inspected independently.

A second strength is the explicit comparison of multiple transcriptome-wide association transformations. Rather than selecting one arbitrary scaling rule and presenting a single ranking, the analysis evaluates a baseline, a linear absolute-z transformation, an exponential absolute-z transformation, and a negative-log10-p transformation. This allows the user to identify conclusions that persist across methods and to distinguish stable findings from method-sensitive findings.

A third strength is the use of several complementary robustness measures. A high Spearman correlation indicates stable ordering, but it does not show whether the highest-risk drugs are retained. Jaccard overlap addresses top-ranked set stability, rank shifts identify individual drug–trait combinations that are sensitive to method choice, amplification factors describe score rescaling, and consensus rankings summarize the average result. The receptor-contribution recomputation also tests whether the dominant mechanistic components persist when the scaling rule changes.

The results provide evidence of strong computational robustness. The mean Spearman correlation was 0.997, Expo\_z and Logp were nearly identical, and only 38 of 260 drug–trait combinations shifted by at least three positions. Clozapine and chlorpromazine were completely stable, and the top receptor contributors remained largely unchanged. These findings support the use of the pipeline for prioritization and comparative analysis.

### **4.2 Choice of scaling method**

The Logp transformation is a reasonable primary option when the objective is to emphasize transcriptome-wide association evidence according to statistical strength. The transformation is transparent, monotonic, and bounded by an explicit cap. It also avoids directly using 1 divided by p, which would be numerically unstable and would produce extreme amplification for very small p values.

However, Logp should not be presented as an established or universally superior TWAS weighting method. Negative-log10 p values do not preserve effect direction and are influenced by sample size, prediction accuracy, and study precision. A highly significant association may reflect a small effect estimated precisely, while a less significant association may reflect a larger effect estimated with greater uncertainty. The Logp transformation therefore represents a deliberate evidence-weighting choice rather than a direct estimate of biological effect.

Expo\_z and Logp produced nearly identical rankings in the present data. This indicates that, over the range of z scores and p values represented in the files, both transformations emphasized similar drug–receptor patterns. Either can therefore be used when a TWAS-enhanced primary score is desired. The Original method remains useful as a sensitivity analysis because it produced the largest top-to-median separation. Nevertheless, the Original method should be described accurately as having no per-receptor scaling, not as having no TWAS influence, because the version 0 implementation retained a drug-level z-score multiplier and a nominal significance bonus.

### **4.3 Biological interpretation**

The persistence of HRH1 among the leading receptor contributors is consistent with pharmacological evidence linking H1 receptor affinity to antipsychotic-associated weight gain. The role of HTR2C is also supported by experimental and pharmacological literature, while muscarinic M3, adrenergic, dopaminergic, and other serotonin receptors may contribute through overlapping mechanisms affecting appetite, insulin secretion, glucose metabolism, lipid metabolism, and energy balance \[25-29\].

The model placed HTR2A second across all methods, but this should not be interpreted as a direct causal ranking of receptor mechanisms. The result emerges from the interaction of the manually specified receptor weights, empirical Ki values, dose normalization, transcriptome-wide association evidence, and aggregation rules. In particular, HTR2C received a larger predefined weight than HTR2A, yet HTR2A ranked higher in the final contribution analysis. This demonstrates that the contribution ranking is a composite output rather than a direct readout of the prior weights.

Clozapine’s stable position is broadly compatible with its well-documented clinical metabolic liability \[2,1,30\]. However, the analysis does not demonstrate that clozapine’s risk is caused primarily by receptor binding. Its high score was preserved across scaling methods, which indicates that the ranking was not strongly dependent on the additional TWAS weighting. It does not establish that receptor pharmacology alone explains the risk or that the score is calibrated to clinical outcomes.

The transcriptome-wide association component should also be interpreted as prioritization evidence rather than causal evidence. TWAS associations can reflect linkage disequilibrium between the variants used in expression prediction and variants that influence the phenotype through another mechanism. Colocalization, fine-mapping, replication, and functional experiments are needed before assigning a causal role to a receptor gene or pathway \[10,31,11,32\].

### **4.4 Limitations**

Several limitations should be addressed before the method is used for clinical decision-making. First, the score has not been externally validated against patient-level metabolic outcomes, such as longitudinal weight change, incident diabetes, triglyceride change, or cardiovascular events. The consensus ranking therefore represents computational prioritization, not validated clinical prediction.

Second, the pipeline uses defined daily dose as a standardized dose proxy. Defined daily dose does not necessarily correspond to the prescribed dose, treatment duration, adherence, plasma concentration, or exposure at the receptor. It is useful for harmonization but cannot replace pharmacokinetic or patient-level exposure data.

Third, the default minimum-Ki aggregation selects the strongest observed binding measurement. This may be appropriate for identifying potential pharmacological interactions, but it can overemphasize an isolated measurement and may not represent the central affinity of a drug across assays. Median or hierarchical aggregation should be evaluated in sensitivity analyses.

Fourth, fuzzy matching can introduce uncertainty. A threshold of 80 reduces the likelihood of accepting very poor matches but does not guarantee that the best match is biologically correct. The final workflow should save the accepted ligand, matching score, alternative candidates, and manual adjudication status.

Fifth, receptor mapping contains simplifying assumptions. For example, a generic alpha-1 label is assigned to ADRA1A, a generic alpha-2 label is assigned to ADRA2A, a generic muscarinic label is assigned to CHRM1, and 5-HT3 is mapped to HTR3A. These mappings may be reasonable operational choices but can be incorrect when the source label is not subtype-specific.

Sixth, missing-value imputation can affect both magnitude and rank. Mean or median Ki imputation creates an estimated affinity that was not directly measured. Mean absolute-z and mean-negative-log10-p imputation provide fairness across drugs with incomplete coverage but can also reduce genuine heterogeneity or introduce nonzero TWAS influence where data are absent. The number and location of imputed rows should be reported for every trait and method.

Seventh, the transcriptome-wide association significance indicator uses a nominal p-value threshold of 0.05. Because the significance bonus is applied to metabolic receptor rows, this may increase scores based on findings that would not survive genome-wide or multi-trait correction. A future implementation should support false-discovery-rate and Bonferroni-adjusted thresholds and should separate exploratory from confirmatory analyses.

Eighth, p-value scaling and z-score scaling have different statistical properties. The Logp method ignores direction, while the z-score methods use absolute z scores and therefore also ignore direction. In addition, p values can be sensitive to sample size and LD reference quality. Summary-based TWAS methods can be affected by mismatches between the study population and the reference LD panel \[10,33\].

Finally, the supplied code is a notebook-oriented artifact rather than a fully packaged software release. It includes repeated function definitions, shell commands, hard-coded paths, and several sequential versions of the pipeline. These features should be refactored before publication. A clean implementation should isolate each scoring method, use explicit configuration files, and include tests confirming that the formulas and normalization rules produce the expected values.

### **4.5 Future extensions**

The pipeline can be expanded in several directions. Additional traits could include body mass index, fasting glucose, insulin resistance, hemoglobin A1c, coronary artery disease, blood pressure, and broader metabolic syndrome components. Multi-trait approaches could be used to summarize correlated outcomes while preserving trait-specific results.

The transcriptome-wide association stage could be extended to multi-tissue models. MultiXcan and related methods jointly analyze predicted expression across tissues while accounting for correlations between tissue-specific prediction models \[33\]. This may be useful when a relevant regulatory context is uncertain or when gene expression effects are shared across tissues.

Colocalization and fine-mapping should be integrated before causal interpretation. Methods such as COLOC and eCAVIAR can help distinguish shared regulatory signals from distinct signals in linkage disequilibrium \[12,13\]. A future release could require colocalization support before allowing a transcriptome-wide association result to influence the primary score.

The receptor weights could also be estimated empirically. Candidate approaches include regularized regression using clinical metabolic outcomes, hierarchical Bayesian models, meta-analytic receptor weighting, and machine-learning models that incorporate receptor affinity, dose, treatment duration, and patient-level covariates. Such approaches would require a sufficiently large and well-characterized clinical training dataset and should be evaluated against the current literature-derived weights.

The software could be distributed as a versioned Python package or workflow using a container, continuous-integration tests, and a workflow manager such as Snakemake or Nextflow. A web-based interface could provide interactive rankings, method comparisons, receptor-level decompositions, and downloadable audit files. Any such interface should clearly label the scores as research prioritization indices rather than clinical predictions.

## **5\. Conclusion**

The Ki–defined dose–transcriptome-wide association study pipeline provides a transparent framework for integrating pharmacological and genetically informed evidence into antipsychotic metabolic-risk indices. Across 52 drugs, five lipid traits, and four scaling methods, drug rankings were highly stable. The mean cross-method Spearman correlation was 0.997, only 38 of 260 drug–trait combinations showed rank shifts of at least three positions, clozapine remained the highest-ranked drug under every method, and HRH1 and HTR2A were the most consistent receptor-level contributors. The principal effect of the transcriptome-wide association transformations was to rescale scores rather than substantially reorder drugs. Logp is a reasonable primary method when evidence-weighted scaling is desired, while the baseline method should be retained as a sensitivity analysis. These scores should be used for comparative prioritization and hypothesis generation, not as clinical risk estimates or evidence of causality.

## **6\. Data and Code Availability**

The code and generated outputs used for this manuscript should be deposited in a public repository before submission. The final repository should include the cleaned pipeline, environment specification, parameter files, input-file metadata, output CSV files, figures, and a persistent DOI. 

The supplied output structure includes per-trait and per-method risk-score files, long-format receptor-level files, method-comparison CSV files, consensus rankings, rank-shift files, amplification factors, Friedman-test results, and figures. Input databases should be redistributed only when licensing permits. If the Ki database cannot be redistributed, the repository should provide its name, version, source, access instructions, and a checksum or equivalent version identifier.

## **7\. Ethics Statement**

This analysis used secondary pharmacological, classification, and transcriptome-wide association resources. No new human participants, animals, clinical specimens, or identifiable patient-level data were collected for the present computational analysis.

## **8\. Funding**

Funding information was not included in the supplied materials. The authors should replace this sentence with the verified funding statement before submission.

## **9\. Declaration of Competing Interest**

Not-applicable.

## **10\. Declaration of Generative AI and AI-Assisted Technologies in the Manuscript Preparation Process**

Not-applicable.

# **References**

1. De Hert M, Detraux J, van Winkel R, Yu W, Correll CU. Metabolic and cardiovascular adverse effects associated with antipsychotic drugs. *Nat Rev Endocrinol*. 2012;8:114-126. doi:10.1038/nrendo.2011.156

2. Newcomer JW. Second-generation (atypical) antipsychotics and metabolic effects: a comprehensive literature review. *CNS Drugs*. 2005;19(suppl 1):1-93. doi:10.2165/00023210-200519001-00001

3. Vancampfort D, Stubbs B, Mitchell AJ, et al. Risk of metabolic syndrome and its components in people with schizophrenia and related psychotic disorders, bipolar disorder, and major depressive disorder: a systematic review and meta-analysis. *World Psychiatry*. 2015;14(3):339-347. doi:10.1002/wps.20252

4. Alberti KGM, Eckel RH, Grundy SM, et al. Harmonizing the metabolic syndrome: a joint interim statement of the International Diabetes Federation Task Force on Epidemiology and Prevention; National Heart, Lung, and Blood Institute; American Heart Association; World Heart Federation; International Atherosclerosis Society; and International Association for the Study of Obesity. *Circulation*. 2009;120(16):1640-1645. doi:10.1161/CIRCULATIONAHA.109.192644

5. Allison DB, Mentore JL, Heo M, et al. Antipsychotic-induced weight gain: a comprehensive research synthesis. *Am J Psychiatry*. 1999;156(11):1686-1696. doi:10.1176/ajp.156.11.1686

6. Pillinger T, McCutcheon RA, Vano L, et al. Comparative effects of 18 antipsychotics on metabolic function in patients with schizophrenia, predictors of metabolic dysregulation, and association with psychopathology: a systematic review and network meta-analysis. *Lancet Psychiatry*. 2020;7(1):64-77. doi:10.1016/S2215-0366(19)30416-X

7. Burschinski A, Schneider-Thoma J, Chiocchia V, et al. Metabolic side effects in persons with schizophrenia during mid- to long-term treatment with antipsychotics: a network meta-analysis of randomized controlled trials. *World Psychiatry*. 2023;22(1):116-128. doi:10.1002/wps.21036

8. Gamazon ER, Wheeler HE, Shah KP, et al. A gene-based association method for mapping traits using reference transcriptome data. *Nat Genet*. 2015;47:1091-1098. doi:10.1038/ng.3367

9. Gusev A, Ko A, Shi H, et al. Integrative approaches for large-scale transcriptome-wide association studies. *Nat Genet*. 2016;48:245-252. doi:10.1038/ng.3506

10. Barbeira AN, Dickinson SP, Bonazzola R, et al. Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics. *Nat Commun*. 2018;9:1825. doi:10.1038/s41467-018-03621-1

11. Wainberg M, Sinnott-Armstrong N, Mancuso N, et al. Opportunities and challenges for transcriptome-wide association studies. *Nat Genet*. 2019;51:592-599. doi:10.1038/s41588-019-0385-z

12. Giambartolomei C, Vukcevic D, Schadt EE, et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genet*. 2014;10(5):e1004383. doi:10.1371/journal.pgen.1004383

13. Hormozdiari F, van de Bunt M, Segrè AV, et al. Colocalization of GWAS and eQTL signals detects target genes. *Am J Hum Genet*. 2016;99(6):1245-1260. doi:10.1016/j.ajhg.2016.10.003

14. Emerging Risk Factors Collaboration. Major lipids, apolipoproteins, and risk of vascular disease. *JAMA*. 2009;302(18):1993-2000. doi:10.1001/jama.2009.1619

15. Ference BA, Ginsberg HN, Graham I, et al. Low-density lipoproteins cause atherosclerotic cardiovascular disease. 1\. Evidence from genetic, epidemiologic, and clinical studies. A consensus statement from the European Atherosclerosis Society Consensus Panel. *Eur Heart J*. 2017;38(32):2459-2472. doi:10.1093/eurheartj/ehx144

16. World Health Organization Collaborating Centre for Drug Statistics Methodology. ATC/DDD Index. Accessed July 15, 2026\. [https://atcddd.fhi.no/atc\_ddd\_index/](https://atcddd.fhi.no/atc_ddd_index/)

17. GTEx Consortium. Genetic effects on gene expression across human tissues. *Nature*. 2017;550:204-213. doi:10.1038/nature24277

18. 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature*. 2015;526:68-74. doi:10.1038/nature15393

19. Bland JM, Altman DG. Statistical methods for assessing agreement between two methods of clinical measurement. *Lancet*. 1986;327(8476):307-310. doi:10.1016/S0140-6736(86)90837-8

20. Spearman C. The proof and measurement of association between two things. *Am J Psychol*. 1904;15(1):72-101. doi:10.2307/1412159

21. Kendall MG. A new measure of rank correlation. *Biometrika*. 1938;30(1-2):81-93. doi:10.1093/biomet/30.1-2.81

22. Jaccard P. Étude comparative de la distribution florale dans une portion des Alpes et du Jura. *Bull Soc Vaud Sci Nat*. 1901;37:547-579. doi:10.5169/seals-266450

23. Friedman M. The use of ranks to avoid the assumption of normality implicit in the analysis of variance. *J Am Stat Assoc*. 1937;32(200):675-701. doi:10.1080/01621459.1937.10503522

24. Demšar J. Statistical comparisons of classifiers over multiple data sets. *J Mach Learn Res*. 2006;7:1-30. Accessed July 15, 2026\. [https://www.jmlr.org/papers/v7/demsar06a.html](https://www.jmlr.org/papers/v7/demsar06a.html)

25. Kroeze WK, Hufeisen SJ, Popadak BA, et al. H1-histamine receptor affinity predicts short-term weight gain for typical and atypical antipsychotic drugs. *Neuropsychopharmacology*. 2003;28:519-526. doi:10.1038/sj.npp.1300027

26. Matsui-Sakata A, Ohtani H, Sawada Y. Receptor occupancy-based analysis of the contributions of various receptors to antipsychotics-induced weight gain and diabetes mellitus. *Drug Metab Pharmacokinet*. 2005;20(5):368-378. doi:10.2133/dmpk.20.368

27. Nasrallah HA. Atypical antipsychotic-induced metabolic side effects: insights from receptor-binding profiles. *Mol Psychiatry*. 2008;13:27-35. doi:10.1038/sj.mp.4002066

28. Reynolds GP, Kirk SL. Metabolic side effects of antipsychotic drug treatment—pharmacological mechanisms. *Pharmacol Ther*. 2010;125(1):169-179. doi:10.1016/j.pharmthera.2009.10.010

29. Roth BL, Sheffler DJ, Kroeze WK. Magic shotguns versus magic bullets: selectivity and multisite pharmacology of atypical antipsychotics. *Nat Rev Drug Discov*. 2004;3:353-359. doi:10.1038/nrd1345

30. Yuen JWY, Kim DD, Procyshyn RM, et al. A focused review of the metabolic side-effects of clozapine. *Front Endocrinol*. 2021;12:609240. doi:10.3389/fendo.2021.609240

31. Mancuso N, Freund MK, Johnson R, et al. Probabilistic fine-mapping of transcriptome-wide association studies. *Nat Genet*. 2019;51:675-682. doi:10.1038/s41588-019-0367-1

32. Zhu Z, Zhang F, Hu H, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. *Nat Genet*. 2016;48:481-487. doi:10.1038/ng.3538

33. Barbeira AN, Bonazzola R, Ganguly P, et al. Integrating tissue-specific mechanisms into GWAS summary data analysis. *Nat Genet*. 2021;53:120-130. doi:10.1038/s41588-020-00744-0

