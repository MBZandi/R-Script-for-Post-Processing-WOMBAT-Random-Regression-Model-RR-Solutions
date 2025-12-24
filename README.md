# R-Script-for-Post-Processing-WOMBAT-Random-Regression-Model-RR-Solutions
This R script extracts animal breeding values (BVs) from a WOMBAT random regression analysis using Legendre polynomials. It is specifically derived and tested on the results of WOMBAT Example 3 (random regression example in the WOMBAT distribution).

In this example:

The trait analyzed is weight.
The random regression control variable is age (days or similar time points).
A quadratic Legendre polynomial (order 2 → 3 coefficients) is fitted for the additive genetic effect (animal).
Heterogeneous residual variances and other effects are modeled as shown in the original wombat.par file.

The script computes the predicted breeding value trajectory for each animal across all unique ages present in the data and then calculates a total (summed) breeding value — a common summary statistic used in post-processing of random regression models when a single overall BV is desired.
Requirements

R installed.
A file named legendre_Coeff.R containing a function Legdre(p) that returns Legendre polynomial values (at least orders 0–2) for a standardized value p ∈ [-1, 1].
Input files (from WOMBAT Example 3) in your working directory:
mrrtst.dat → data file.
RnSoln_animal.dat → solutions file for the additive genetic (animal) random regression coefficients.


Output

Console output showing progress and previews.
Optional files:
plyleg.txt: Matrix of Legendre basis functions evaluated at standardized ages.
ID_RR_BV_wombat.txt: Table with animal ID and total summed breeding value (BV).
