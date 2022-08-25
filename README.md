panstripe
================

<!-- badges: start -->
[![R-CMD-check](https://github.com/gtonkinhill/panstripe/workflows/R-CMD-check/badge.svg)](https://github.com/gtonkinhill/panstripe/actions)
[![DOI](https://zenodo.org/badge/353592486.svg)](https://zenodo.org/badge/latestdoi/353592486)
<!-- badges: end -->


<p align="center">
<img src="https://github.com/gtonkinhill/panstripe/blob/main/inst/vignette-supp/panstripe.png" alt="alt text" width="500">
</p>

Panstripe improves the post processing of bacterial pangenome analyses.
In particular it aims to replace the dubious but popular pangenome
accumulation curve. The package is currently under development so
frequent changes are to be expected.

-   <a href="#installation" id="toc-installation">Installation</a>
-   <a href="#quick-start" id="toc-quick-start">Quick Start</a>
-   <a href="#citation" id="toc-citation">Citation</a>
-   <a href="#comparing-pangenomes" id="toc-comparing-pangenomes">Comparing
    Pangenomes</a>
-   <a href="#open-vs-closed" id="toc-open-vs-closed">Open vs Closed</a>
-   <a href="#rate-vs-size" id="toc-rate-vs-size">Rate vs Size</a>
-   <a href="#output" id="toc-output">Output</a>
-   <a href="#alternative-models" id="toc-alternative-models">Alternative
    models</a>
-   <a href="#plots" id="toc-plots">Plots</a>
    -   <a href="#pangenome-fit" id="toc-pangenome-fit">Pangenome fit</a>
    -   <a href="#tree-with-presenceabsence"
        id="toc-tree-with-presenceabsence">Tree with presence/absence</a>
    -   <a href="#inferred-ancestral-states"
        id="toc-inferred-ancestral-states">Inferred ancestral states</a>
    -   <a href="#tsne" id="toc-tsne">tSNE</a>
    -   <a href="#accumulation-curves" id="toc-accumulation-curves">Accumulation
        curves</a>
    -   <a href="#etymology" id="toc-etymology">Etymology</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Installation

`panstripe` is currently available on GitHub. It can be installed with
`devtools`

``` r
install.packages("remotes")
remotes::install_github("gtonkinhill/panstripe")
```

If you would like to also build the vignette with your installation run:

``` r
remotes::install_github("gtonkinhill/panstripe", build_vignettes = TRUE)
```

## Quick Start

Panstripe takes as input a phylogeny and gene presence/absence matrix in
Rtab format as is produced by
[Panaroo](https://gtonkinhill.github.io/panaroo/#/),
[Roary](https://github.com/sanger-pathogens/Roary),
[PIRATE](https://github.com/SionBayliss/PIRATE) and other prokaryotic
pangenome tools.

``` r
library(panstripe)
library(ape)
set.seed(1234)

### NOTE: here we load example files from the panstripe package. You should
### replace these variables with the relevant paths to the files you are using.
phylo.file.name <- system.file("extdata", "tree.newick", package = "panstripe")
rtab.file.name <- system.file("extdata", "gene_presence_absence.Rtab", package = "panstripe")
###

# Load files
pa <- read_rtab(rtab.file.name)
tree <- read.tree(phylo.file.name)

# Run panstripe
fit <- panstripe(pa, tree)
fit$summary
#> # A tibble: 6 × 7
#>   term  estimate std.error statistic   p.value `bootstrap CI …` `bootstrap CI …`
#>   <chr>    <dbl>     <dbl>     <dbl>     <dbl>            <dbl>            <dbl>
#> 1 istip   1.49      0.282      5.29   8.03e- 7           0.987             1.98 
#> 2 core    2.97      0.247     12.0    9.69e-21           2.43              3.43 
#> 3 depth   0.0572    0.0619     0.923  3.58e- 1          -0.0605            0.185
#> 4 isti…  -1.93      0.431     -4.48   2.10e- 5          -2.71             -1.06 
#> 5 p       1.25     NA         NA     NA                  1.14              1.44 
#> 6 phi     2.58     NA         NA     NA                  2.06              3.82

# Plot results
plot_pangenome_fits(fit, include_data = TRUE)
```

![](inst/vignette-supp/unnamed-chunk-4-1.png)<!-- -->

A significant p-value for the `tip` term indicates that there is a
different rate of gene exchange at the tips of the phylogeny compared
with the internal branches. This is usually driven by annotation errors
or highly mobile elements that do not persist long enough to be observed
in multiple genomes.

A significant p-value for the `core` term indicates that there is a
significant association between the core genome branch length and the
number of gene exchange events.

The `depth` term is less interesting but indicates when there is a
difference in our ability to detect older gene exchange events.

## Citation

To cite panstripe please use:

## Comparing Pangenomes

As `panstripe` uses a GLM framework it is straightforward to compare the
slope terms between datasets. The `compare_pangenomes` function
considers the interaction term between the data sets and the `core`,
`tip` and `depth` terms.

> **IMPORTANT:** *It is important that the phylogenies for the
> pangenomes being compared are on the same scale. This can be achieved
> most easily by using time scaled phylogenies. Alternatively, it is
> possible to build a single phylogeny and partition it into clades or
> use SNP scaled phylogenies.*

Here, we simulate two pangenomes with different gene exhange rates, fit
the model using `panstripe` and run the `compare_pangenomes` function.

``` r
# Simulate a fast gene gain/loss rate with error
sim_fast <- simulate_pan(rate = 0.001, ngenomes = 200)
sim_slow <- simulate_pan(rate = 5e-04, ngenomes = 200)

# Run panstripe
fit_fast <- panstripe(sim_fast$pa, sim_fast$tree)
fit_slow <- panstripe(sim_slow$pa, sim_slow$tree)

# Compare the fits
result <- compare_pangenomes(fit_fast, fit_slow)
result$summary
#> # A tibble: 4 × 7
#>   term    estimate std.error statistic p.value `bootstrap CI …` `bootstrap CI …`
#>   <chr>      <dbl>     <dbl>     <dbl>   <dbl>            <dbl>            <dbl>
#> 1 depth    -0.0578    0.0268    -2.16  3.14e-2          -0.0990         -0.00988
#> 2 istip    -0.0464    0.173     -0.268 7.88e-1          -0.346           0.199  
#> 3 core     -0.683     0.165     -4.15  3.74e-5          -1.00           -0.380  
#> 4 disper…  NA        NA          0.439 5.08e-1          NA              NA
```

A significant p-value for the `tip` term indicates that the two
pangenomes differ in the rates of gene presence and absence assigned to
the tips of the phylogeny. This is usually driven by either differences
in the annotation error rates between data sets or differences in the
gain and loss of highly mobile elements that do not persist long enough
to be observed in multiple genomes.

A significant p-value for the `core` term indicates that two data sets
have different rates of gene gain and loss.

The `depth` term is less interesting but indicates when there is a
difference in our ability to detect older gene exchange events between
the two pangenomes.

The `dispersion` parameter indicates if there is a significant
difference in the dispersion of the two pangenomes. This suggests that
the relationship between the rate of gene exchange and the size of each
event differs in the two pangenomes. Here, the p-value is obtained using
a Likelihood Ratio Test.

## Open vs Closed

The definition of what constitutes an open or closed pangenome is
somewhat ambiguous. We prefer to consider whether there is evidence for
a temporal signal in the pattern of gene gain and loss. This prevents
annotation errors leading to misleading results.

After fitting a `panstripe` model the significance of the temporal
signal can be assessed by considering the coefficient and p-value of the
‘core’ term in the model. The uncertainty of this estimate can also be
investigated by looking at the bootstrap confidence intervals of the
core term.

Let’s simulate a ‘closed’ pangenome

``` r
sim_closed <- simulate_pan(rate = 0, ngenomes = 100)
fit_closed <- panstripe(sim_closed$pa, sim_closed$tree)

fit_closed$summary
#> # A tibble: 6 × 7
#>   term  estimate std.error statistic   p.value `bootstrap CI …` `bootstrap CI …`
#>   <chr>    <dbl>     <dbl>     <dbl>     <dbl>            <dbl>            <dbl>
#> 1 istip  1.53e+0   1.43e-1   10.7     3.09e-21            1.14            1.93  
#> 2 core  -7.23e+2   3.80e+4   -0.0190  9.85e- 1         -884.           -613.    
#> 3 depth -4.52e-2   3.74e-2   -1.21    2.29e- 1           -0.149           0.0563
#> 4 isti…  7.23e+2   3.80e+4    0.0190  9.85e- 1          613.            884.    
#> 5 p      1.05e+0  NA         NA      NA                   1.00            1.60  
#> 6 phi    3.99e-1  NA         NA      NA                   0.272           1.03
```

The p-value indicates that there is not a significant association
between core branch lengths and gene gain/loss. This is typical of
species that undergo very little to no recombination such as *M.
tuberculosis*.

## Rate vs Size

While comparing the `core` parameter of the model identifies differences
in the association between branch lengths and gene gain and loss it does
not indicate whether this is driven by higher rates of recombination or
simply larger recombination events involving more genes.

The `panstripe` model allows these two scenarios to be investigated by
allowing the dispersion parameter in the model to be different for each
pangenome.

Here, we simulate two data sets with the same recombination rate but
where each recombination event differs in the number of genes involved.

``` r
sim_large <- simulate_pan(rate = 0.001, ngenomes = 100, mean_trans_size = 4)
sim_small <- simulate_pan(rate = 0.001, ngenomes = 100, mean_trans_size = 3)

fit_large <- panstripe(sim_large$pa, sim_large$tree)
fit_small <- panstripe(sim_small$pa, sim_small$tree)

# Compare the fits
result <- compare_pangenomes(fit_large, fit_small)
result$summary
#> # A tibble: 4 × 7
#>   term    estimate std.error statistic p.value `bootstrap CI …` `bootstrap CI …`
#>   <chr>      <dbl>     <dbl>     <dbl>   <dbl>            <dbl>            <dbl>
#> 1 depth     0.0858    0.0490     1.75   0.0807          -0.0209           0.206 
#> 2 istip    -0.181     0.297     -0.609  0.543           -0.803            0.358 
#> 3 core     -0.365     0.226     -1.62   0.106           -0.864            0.0888
#> 4 disper…  NA        NA          6.14   0.0132          NA               NA
```

## Output

The `panstripe` function generates a list with the following attributes

#### summary

A table indicating a subset of the inferred parameters of the GLM. The
p-values and bootstrap confidence intervals can be used to determine
whether each term in the model is significantly associated with gene
gain and loss.

-   **core** indicates whether the branch lengths in the phylogeny are
    associated with gene gain and loss.

-   **tip** indicates associations with genes observed to occur on the
    tips of the phylogeny. These are usually driven by a combination of
    annotation errors and depending upon the temporal sampling density
    also highly mobile elements that are not observed in multiple
    genomes.

-   **depth** indicates whether the rate of gene gain and loss changes
    significantly with the depth of a branch. Typically, our ability to
    detect gene exchange events reduces for older ancestral branches.

-   **p** the inferred index parameter of the underlying Tweedie
    distribution used in the GLM

-   **phi** the inferred dispersion parameter of the GLM

#### model

The output of fitting the GLM model. This object can be used to predict
how many gene gains and losses we would expect to observe given the
parameters of a particular branch.

#### data

A table with the data used to fit the model. The `acc` column indicates
the inferred number of gene gain/loss events for a branch; the `core`
column is the branch length taken from the phylogeny; the `istip` column
indicates whether the branch occurs at the tip of the phylogeny and the
`depth` column indicates the distance from the root node to the branch.

#### ci_samples

A ‘boot’ object, generated by the
[boot](https://cran.r-project.org/web/packages/boot/boot.pdf) package.
Can be used to investigate the uncertainty in the parameter estimates.

#### tree/pa

The original data provided to the `panstripe` function.

## Alternative models

The default Panstripe model assumes a Compound Poisson (Tweedie)
distribution as implemented in the
[Tweedie](https://cran.r-project.org/web/packages/tweedie/) R package.
In some cases there may be insufficient data or the data may not
reliably fit the Tweedie distribution. This is often the case when there
is nearly no gene exchange events inferred to have occurred at the
internal branches of the phylogeny.

To help account for these issues and to add flexibility to the package
it is also possible to fit alternative distributions in place of the
Tweedie model. The most common alternative that we suggest is to use a
Gaussian distribution. This is usually more robust and is less likely to
run into convergence issues in the model fit. In our tests both of these
distributions generally give very similar results.

Panstripe can be run using an alternative distribution as

``` r
fit_gaussian <- panstripe(pa, tree, family = "gaussian")
```

## Plots

Panstripe includes a number of useful plotting functions to help with
interpretation of the output of
[panaroo](https://gtonkinhill.github.io/panaroo/#/).

### Pangenome fit

A simple plot of the fit of the model along with the input data points
can be made by running

``` r
plot_pangenome_fits(fit_fast)
```

![](inst/vignette-supp/unnamed-chunk-9-1.png)<!-- -->

By default this will plot just the model fit. It is also possible to
include the data points of the branches located at the tips of the
phylogeny.

``` r
plot_pangenome_fits(fit_fast, include_data = TRUE)
```

![](inst/vignette-supp/unnamed-chunk-10-1.png)<!-- -->

The function can also take a named list as input allowing for easy
comparisons between data sets. Optionally the plot can be trimmed to
cover the same range when considering multiple pangenome curves.

``` r
plot_pangenome_fits(list(fast = fit_fast, slow = fit_slow), trim = TRUE)
```

![](inst/vignette-supp/unnamed-chunk-11-1.png)<!-- -->

The parameters and accompanying error bars can be plotted as

``` r
plot_pangenome_params(list(fast = fit_fast, slow = fit_slow))
```

![](inst/vignette-supp/unnamed-chunk-12-1.png)<!-- -->

A plot of the residuals of the regression can also be generated using

``` r
plot_residuals(fit_fast)
```

![](inst/vignette-supp/unnamed-chunk-13-1.png)<!-- -->

It is also possible to generate a cumulative plot of the branch lengths
versus the number of gene gain and loss events at each node of the input
phylogeny. As Panstripe models each branch individually and accounts for
the depth of the branch it is not possible to add the Panstripe model
fit to this plot.

``` r
plot_pangenome_cumulative(list(fast = fit_fast, slow = fit_slow))
```

![](inst/vignette-supp/unnamed-chunk-14-1.png)<!-- -->

### Tree with presence/absence

A plot of the phylogeny and the corresponding gene presence/absence
matrix can be made using the `plot_tree_pa` function. This function also
takes an optional vector of gene names to include in the plot.

``` r
# look at only those genes that vary
variable_genes <- colnames(pa)[apply(pa, 2, sd) > 0]

plot_tree_pa(tree = tree, pa = pa, genes = variable_genes, label_genes = FALSE, cols = "black")
```

![](inst/vignette-supp/unnamed-chunk-15-1.png)<!-- -->

### Inferred ancestral states

The `plot_gain_loss` function allows for the visualisation of the fitted
gene gain and loss events on the given phylogeny. The enrichment for
events at the tips of a tree is often driven by a combination of highly
mobile elements and annotation errors.

``` r
plot_gain_loss(fit)
```

![](inst/vignette-supp/unnamed-chunk-16-1.png)<!-- -->

### tSNE

The tSNE dimension reduction technique can be used to investigate
evidence for clusters within the pangenome.

``` r
plot_tsne(pa)
```

![](inst/vignette-supp/unnamed-chunk-17-1.png)<!-- -->

The [Mandrake](https://github.com/johnlees/mandrake) method can also be
used as an alternative to tSNE.

### Accumulation curves

While we do not recommend the use of accumulation curves as they do not
account for population structure, sampling bias or annotation errors we
have included a function to plot them to make it easier for users to
compare methods.

``` r
plot_acc(list(fast = sim_fast$pa, slow = sim_slow$pa))
```

![](inst/vignette-supp/unnamed-chunk-18-1.png)<!-- -->

### Etymology

The name panstripe is an adaptation of ‘pinstripe’, the name of a
[long-nosed potoroo](https://en.wikipedia.org/wiki/Long-nosed_potoroo)
who was a villain in the playstation game [Crash
Bandicoot](https://crashbandicoot.fandom.com/wiki/Pinstripe_Potoroo).
The name was chosen as the program relies on the output of panaroo
(named after the potoroo).
