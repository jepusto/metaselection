# Practice Facilitation Meta-Analysis

Meta-analytic dataset containing results of primary studies examining
the effect of practice facilitation on the uptake of evidence-based
practices (EBPs) in primary care settings.

## Usage

``` r
practice_facilitation
```

## Format

A data frame with 23 rows and 3 variables:

- author:

  first author and publication year of primary study report.

- score:

  score on a scale from 0 to 12, for which higher scores correspond to
  higher quality of the study methods.

- design:

  study design, with CCT = controlled clinical trial, C-RCT = cluster
  randomized controlled trial, RCT = randomized controlled trial.

- allocation_concealed:

  indicator for allocation concealment.

- blinded:

  indicator for whether study was single- or double-blinded.

- intent_to_treat:

  indicator for whether study adhered to intent-to-treat principle.

- outcome:

  description of outcome measure.

- follow_up:

  months of follow-up.

- retention_pct:

  percentage of sample retained at follow-up.

- SMD:

  effect size in form of Hedges' g.

- SE:

  corresponding variance of the effect size.

## Source

Table 1 of Baskerville et al. (2012;
[doi:10.1370/afm.1312](https://doi.org/10.1370/afm.1312) ).

## References

Baskerville NB, Liddy C, Hogg W (2012). “Systematic review and
meta-analysis of practice facilitation within primary care settings.”
*Annals of Family Medicine*, **10**(1), 63–74. ISSN 1544-1717.
[doi:10.1370/afm.1312](https://doi.org/10.1370/afm.1312) .
