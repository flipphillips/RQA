# TODO

![logo](RQA/Documentation/icon.png)

## Utilitron

- $RQAVersion


## Time-series related

- Check for `RegularlySampledQ`
- Even better, allow for non time series and allow for re-sampling option.

## RM

- Change routines to ues sparse matrix, single representation of distance map rather than the binary version.
- Alternately, use distance map but keep distance in rm
- ML of classification typology (typology, can be classified as homogeneous, periodic, drift, and disrupted)

## Measures and Variations of RQA

- Older measures like CORM?
- perpendicular recurrence plot

## Validate

- Variation against Webber numbers?

## Parameter estimation

- Finish implementation of the Estimation of Embedding Dimensions. I based this on an older Kennel reference, but there is a new kid in town:
  
  Kennel, M. B., & Abarbanel, H. D. I. (2002). False neighbors and false strands: a reliable minimum embedding dimension algorithm. *Physical Review. E, Statistical, Nonlinear, and Soft Matter Physics,* 66(2 Pt 2), 026209–026218. https://doi.org/10.1103/PhysRevE.66.026209

- Neighborhood and metric estimation:
    "because the neighbourhood of xi does not have to be the same as that of xj . This property leads to an asymmetric RP."
    "fixed amount of nearest neighbours (FAN)"

    N. Marwan, M.C. Romano, M. Thiel, J. Kurths, Recurrence plots for the analysis of complex systems. Phys. Rep. 438(5–6), 237–329 (2007)
    N. Marwan, How to avoid potential pitfalls in recurrence plot based data analysis. Int. J. Bifurcat. Chaos 21(4), 1003–1017 (2011)

- Tau estimation:
    "A “rule of thumb” for the choice of the threshold " is to select it as a few per cent (not larger than 10 %) of the maximum phase space diameter "
    "Zbilut et al. [7] have suggested to choose " so that the recurrence point density is approximately 1 %. For noisy periodic processes"

    N. Marwan, M.C. Romano, M. Thiel, J. Kurths, Recurrence plots for the analysis of complex systems. Phys. Rep. 438(5–6), 237–329 (2007)
    M. Thiel, M.C. Romano, J. Kurths, R. Meucci, E. Allaria, F.T. Arecchi, Influence of observational noise on the recurrence quantification analysis. Physica D 171(3), 138–152 (2002)
    L. Matassini, H. Kantz, J.A. Hołyst, R. Hegger, Optimizing of recurrence plots for noise reduction. Phys. Rev. E 65(2), 021102 (2002)
    S. Schinkel, O. Dimigen, N. Marwan, Selection of recurrence threshold for signal detection. Eur. Phys. J. Spec. Top. 164(1), 45–53 (2008)
