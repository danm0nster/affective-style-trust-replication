[![CC BY 4.0][cc-by-shield]][cc-by]

# Replication data and code for  _The Affective Style of Politics_

This repository contains data and R code to replicate the results in the article
> Nielsen, J. H., & Mønster, D. (2021). The Affective Style of Politics: Evidence from Surveys and Laboratory Experiments. https://doi.org/10.17605/OSF.IO/FJ6DZ

## Data
The following three data files are included:
1. `study1_dk.csv`: Data from Study 1, Denmark.
2. `study1_us.csv`: Data from Study 1, US.
3. `study2.csv`: Data from Study 2.

If you use any of the datasets above, please cite our article, and this repository.

For detailed information about the variables in the data sets, see Nielsen & Mønster (2021).

In addition, there are some data used in a table and plot in the file `iaps_db.csv`. These are taken from Lang et al. (2008).

## Replication scripts
To replicate the results, we provide two scripts&mdash;one for each of the two studies&mdash;that should be run under [R version 4.1](https://cran.r-project.org) or later.

To ensure replicability of the scripts, we use the [groundhog package](https://groundhogr.com/). Using groundhog guarantees that you will use the same versions of all packages that we used when producing the results in the paper. Therefore, you only need to install the `groundhog` package, which will the take care of installing the correct versions of all other packages without interfering with your current R installation. You can install `groundhog` with the following command:
```
install.packages('groundhog')
```

Running the two scripts will generate the tables (in HTML format) and figures (in PDF format) in the article that rely on data. Some tables are split into several parts&mdash;for each of the two countries, Denmark and the US, included in study 1; and valence and arousal or the three affective styles in study 2. In some cases the generated tables are formatted slightly differently from the final tables in the article or online appendix, but the information is the same.

Run the scripts in a new R session, or restart the R session in RStudio, to avoid conflicts with other versions of required packages that may already be loaded.

The scripts use R (R Core Team, 2021) and the following packages
* `groundhog`: Simonsohn & Gruson, 2021.
* `dplyr`: Wickham et al., 2021.
* `table1`: Rich, 2021.
* `markdown`: Allaire et al., 2019.
* `psych`: Revelle 2021.
* `kableExtra`: Zhu, 2021.
* `lme4`: Bates et al., 2015.
* `MuMIn`: Bartón, 2020.
* `lmerTest`: Kuznetsova et al., 2017.
* `texreg`: Leifeld, 2013.
* `effectsize`: Ben-Shachar et al., 2020.
* `tidyr`: Wickham, 2021.
* `ggplot2`: Wickham, 2016.
* `scales`: Wickham & Seidel, 2020.
* `directlabels`: Hocking, 2021.
* `emmeans`: Lenth, 2021.

## References
Allaire, JJ, Jeffrey Horner, Yihui Xie, Vicent Marti and Natacha Porte (2019). markdown: Render Markdown with the C Library 'Sundown'. https://CRAN.R-project.org/package=markdown

Bartoń, Kamil (2020). MuMIn: Multi-Model Inference.  https://CRAN.R-project.org/package=MuMIn

Ben-Shachar M, Lüdecke D, Makowski D (2020). effectsize: Estimation of Effect Size Indices and Standardized Parameters. Journal of Open Source Software, 5(56), 2815. doi: 10.21105/joss.02815

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48. DOI: 10.18637/jss.v067.i01.

Hocking, Toby Dylan (2021). directlabels: Direct Labels for Multicolor Plots. https://CRAN.R-project.org/package=directlabels

Kuznetsova A, Brockhoff PB, Christensen RHB (2017). “lmerTest Package: Tests in Linear Mixed Effects Models.” _Journal of Statistical Software_, 82(13), 1-26. doi: 10.18637/jss.v082.i13 (URL: https://doi.org/10.18637/jss.v082.i13).

Lang, P. J., Bradley, M. M., & Cuthbert, B. N. (1997). International affective picture system (IAPS): Technical manual and affective ratings. NIMH Center for the Study of Emotion and Attention, 1(39-58), 3.

Leifeld, Philip (2013). texreg: Conversion of Statistical Model Output in R to LaTeX and HTML Tables. Journal of Statistical Software, 55(8), 1-24. URL http://dx.doi.org/10.18637/jss.v055.i08.

Lenth, Russell V. (2021). emmeans: Estimated Marginal Means, aka Least-Squares Means. https://CRAN.R-project.org/package=emmeans

Nielsen, J. H., & Mønster, D. (2021). The Affective Style of Politics: Evidence from Surveys and Laboratory Experiments. OSF preprint. [DOI: 10.17605/OSF.IO/FJ6DZ](https://doi.org/10.17605/OSF.IO/FJ6DZ)

R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/.

Revelle, W. (2021) psych: Procedures for Personality and Psychological Research, Northwestern University, Evanston, Illinois, USA, https://CRAN.R-project.org/package=psych

Rich, Benjamin (2021). table1: Tables of Descriptive Statistics in HTML.  https://CRAN.R-project.org/package=table1

Simonsohn, Uri and Hugo Gruson (2021). groundhog: The Simplest Solution to Version-Control for CRAN Packages. https://CRAN.R-project.org/package=groundhog

Wickham, H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.

Wickham, Hadley (2021). tidyr: Tidy Messy Data. https://CRAN.R-project.org/package=tidyr

Wickham, Hadley and Dana Seidel (2020). scales: Scale Functions for Visualization.  https://CRAN.R-project.org/package=scales

Wickham, Hadley, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation.  https://CRAN.R-project.org/package=dplyr

Zhu, Hao (2021). kableExtra: Construct Complex Table with 'kable' and Pipe Syntax.  https://CRAN.R-project.org/package=kableExtra

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
