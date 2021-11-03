# Potential Issues
This list contains potential bugs & problems detected in the original code while refactoring.

* [check_element](https://github.com/RECETOX/recetox-xMSannotator/blob/master/R/check_element.R) detects `C` if formula contains `Cl` -> No compounds are filtered because they don't contain carbon
* [multilevelannotationstep5](https://github.com/RECETOX/recetox-xMSannotator/blob/c7ad3bb2f4e7cc6aa35d8a0804fe6ef7a8389729/R/multilevelannotationstep5.R#L77-L101): 
   * the `.combine = rbind` parameter of `foreach::foreach()` function causes a twisted output, collapsing a `174x60` 
  matrix to a `61x10` matrix -> 61 are the non-empty rows, but the entries in rows that have more than 10 columns filled get discarded;
   * the issue described above may cause `multilevelannotationstep5` to compute different outputs depending on whether 
  the loop is executed _in parallel_ (`%dopar%`) or _sequentially_ (`%do%`). The results may differ in **chemical scores** or the **number 
  of annotations**. When chemical scores differ, _sequential_ execution will systematically output lower scores.
* in [multilevelannotationstep4](https://github.com/RECETOX/recetox-xMSannotator/blob/master/R/multilevelannotationstep4.R) (and functions called by it) are some suspicious conditions where they check for `curdata$score < 10` and `curdata$score > 10`, but never for actual value `10`
* in [get_chemscorev1.6.71](https://github.com/RECETOX/recetox-xMSannotator/blob/d1f454cd4f4b0f6640965cd17aebc19f3217ef95/R/get_chemscorev1.6.71.R#L136) in the isotope detection, the RT clustering is recomputed with a fixed RT diff of 10 seconds, changing the whole previous RT clustering.
* in [get_chemscorev1.6.71](https://github.com/RECETOX/recetox-xMSannotator/blob/d1f454cd4f4b0f6640965cd17aebc19f3217ef95/R/get_chemscorev1.6.71.R#L312-L327) this part for the calculation is never called under current tests.
* In [get_chemscorev1.6.71](https://github.com/RECETOX/recetox-xMSannotator/blob/d1f454cd4f4b0f6640965cd17aebc19f3217ef95/R/get_chemscorev1.6.71.R#L675) the correlation matrix is apparently only read on the diagonal, which implies that values will always be 1. The function could also just return a matrix of 1s with the size of `mzid_cur` x `mzid_cur`
