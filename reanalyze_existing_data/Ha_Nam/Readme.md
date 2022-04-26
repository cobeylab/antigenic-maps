# Background

**Definition of antigenic distance**

The antigenic distance between influenza strains is typically estimated from serological panel data estimated using antisera from naive lab ferrets infected with a single reference strain.
Each ferret generates a primary antibody response specific to the infecting strain, i.
The titer of each antiserum against a panel of heterologous strains indicates how well antibodies specific to strain i cross-react with other strains (j).
If an antiserum's heterologous titer to strain j is much less than the heterologous titer to strain i, then we conclude that the antigenic distance $d_{ij}$ is large, and that i and j weakly cross-react.
Formally, if $s_{ii}$ represents the homologous log$_2$ titer between antiserum i and strain i, and $s_{ij}$ represents the heterologous log$_2$ titer to j, then antigenic distance is defined as

$d_{ij} = s_{ii}-s_{ij}$ (equation 1)

(Neher et al. 2016). Titers are often reported on the log$_2$ scale to reflect 2-fold serial dilutions used in titer measurement. Equation 1 summarized the conceptual foundaiton for many distance calculations, although some studies use slight variations of this definition, e.g. by adjusting for serum potency and strain avidity (Bedford et al. 2014; Smith et al. 2004; Cai et al. 2010).


**Interpretation**

Strictly speaking, antigenic distances should be interpreted as the extent to which antibodies in a primary ferret antiserum rasied against strain i cross-react with strain j.
But in practice, antigenic distance estimates are widely used to model population susceptibilty, to estimate the novelty (and implicitly the fitness) of new variants, and to guide vaccine strain selection.
In these contexts, ferret distance estimates are used as a proxy for the expected level of cross-reactivity between a human antibody reponse to a recent strain i, and a new variant j.
There is a need to systematically assess whether ferret distances accurately and precisely estimate individual cross-reactivity (and heterogeneity in it) against new antigenic variants of influenza, Sars-Cov-2, etc.


# Approach

We calculate **individual distances** using sera collected from individuals who were recently infected by or vaccinated with a known strain. Specifically, this analysis focuses on data from the Ha Nam cohort, available in the supplement of Fonville et al. 2014. This dataset includes sera from:

* 12 individuals from a household cohort study with PCR-confirmed infections (all unvaccinated)
* 106 individuals vaccinated in 1995 with A/Nangchang933/95
* 128 individuals vaccinatd in 1997 with A/Sydney/5/97
* 82 individuals vacinated in 2009 with A/Perth16/2009
* 82 individuals vaccinated in 2010 with A/Vietnam/53/2010

In these analyses, the individual distance is equal to [homologous logtiter]-[heterologous logtiter]. 
The homologous strain is either the vaccine strain, or an isolate that circulated in the same year (and ideally in the same location) as the observed infection.
The heterologous strain is an arbitrary strain in the serological panel.


We also consider 56 individuals from a household cohort study with **no known history of PCR-confirmed infection** (all unvaccinated)

In these analyses, the individual distance is equal to the median of logtiter distances between recently circulating strains, and the test strain.

In all anlyses the **ferret distance** is the Euclidean distance between strain coordinates published in table S3 for the same pair, or pairs of strains used to calculate individual distances.
