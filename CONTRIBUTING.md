# Contribution Guide

## Bugs

For bug reports please submit an issue or pull request.

## Adding new features

One area that `CurvilinearGrids` will benefit the most from is the addition of new discretization schemes. Currently, only the `MontoneExplicitGradientScheme6thOrder` scheme is implemented. New schemes could include central schemes (4th order, 6th order, etc.) or upwind schemes like WENO or WCNS. Not all upwind schemes will work for curvilinear transformations however, since they may not be freestream-preserving. Submit an issue if you have a question. 