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

  First author and publication year of primary study report

- score:

  Scored on a scale from 0 to 12, in which the higher the score, the
  higher the quality of the study methods.

- design:

  Study design, with CCT = controlled clinical trial, C-RCT = cluster
  randomized controlled trial, RCT = randomized controlled trial.

- allocation_concealed:

  Indicator for allocation concealment.

- blinded:

  Indicator for whether study was single- or double-blinded.

- intent_to_treat:

  Indicator for whether study adhered to intent-to-treat principle.

- outcome:

  Description of outcome measure.

- follow_up:

  Months of follow-up.

- retention_pct:

  Percentage of sample retained at follow-up.

- SMD:

  effect size in form of Hedges' g

- SE:

  corresponding variance of the effect size

## Source

[Table 1 of Baskerville et al. (2012)](https://doi.org/10.1370/afm.1312)

## References

Baskerville NB, Liddy C, Hogg W (2012). “Systematic review and
meta-analysis of practice facilitation within primary care settings.”
*Annals of Family Medicine*, **10**(1), 63–74. ISSN 1544-1717.
[doi:10.1370/afm.1312](https://doi.org/10.1370/afm.1312) .
