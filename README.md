## Synopsis

This R code implements the weighted and centered least squares method described in Boruvka, Almirall, Witkiewitz and Murphy (2016). The method is useful for examining proximal or lagged effects of time-varying treatment on a time-varying outcome conditional on baseline or time-varying moderators. The methodology is motivated by mHealth applications (where treatments aimed at behavior change and maintenance are provided in near real time, e.g., via a smartphone) using data arising from a micro-randomized trial (MRT). However, the methodology is applicable to any setting with time-varying treatment, outcome and candidate moderators.

## Motivation for this Code

The method is straightforward to implement, using over-the-counter weighted least squares software (i.e., lm.wt in R).  In the method, the weights are defined as a ratio of two treatment probabilities. The impetus for providing this code is to provide analysts/researchers with the calculations necessary to correct standard errors for small samples and for estimated weights (the numerator of the weights are estimated, or both the numerator and denominator of the weights are estimated).

## Files

Three files are included:

1. weighting-and-centering-estimator-known-randomization-probs.R

This code implements the method for settings where the numerator treatment probabilities are estimated, but the denominator treatment probabilities are known.  This implementation is most likely to be used with data arising from a micro-randomized trial, where the randomizations probabilities are known.

2. weighting-and-centering-estimator-estimated-randomization-probs.R

This code implements the method for settings where both the numerator and denominator treatment probabilities are estimated.  This implementation is most likely to be used with observational study data or any data where the denominator probabilities are not known.

3. fakeMRTdata.txt

This file contains data from a simulated micro-randomized trial. This file can be used to test run the code above.

## Contributors or Testers

Daniel Almirall, Audrey Boruvka, Walter Dempsey, Brook Luers, Timothy NeCamp, Nick Seewald (in alphabetical order).

