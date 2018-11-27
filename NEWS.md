# CHANGES IN AnaCoDa 0.1.2

## BUG FIXES
- fixed a bug were the scaling of observed phi values was used inconsitently, causing problems with estimates of Aphi and Sepsilon

## NEW FEATURES
- Added SCUO calculation and improved getCSPEstimates to include reference codons

# CHANGES IN AnaCoDa 0.1.1

## BUG FIXES
- fixed problem with getCSPEstimates where log scaling was falsely enabled

- fixed problem where the grouplist was not stored by writeParameterObject

## NEW FEATURES
- Added functions to calculate the Codon Adaptation Index, Effective Number of Codons and selection coefficients.

- Allow to set initial phi values based on observed phi values stored in genome object.


