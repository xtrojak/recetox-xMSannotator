# Simple Annotation
* The input database is expected to contain only the compounds without adduct information ([example](https://umsa.cerit-sc.cz/library/list#folders/F1c84aa7fc4490e6d/datasets/8ee788c99983ff96)) and not also the adduct related information ([example](https://github.com/kuppal2/xMSannotator/blob/master/data/hmdbAllinf.rda)).
* The code computing the matching was rewritten in C++ for performance reasons and also performs the computation of expected mz value for a compound-adduct-pair (see [here](https://github.com/RECETOX/recetox-xMSannotator/blob/a410b00449db6e4bbc07339a51ef32336f405ce1/xmsannotator/src/match_by_mass.cpp#L22-L47))

# RT Clustering
* The RT-based clustering as the second step of feature grouping has been reworked as described [here](https://github.com/RECETOX/recetox-xMSannotator/pull/17)
* The clustering by this new method results in a finer grouping of peaks into more individual clusters, see the example image below.
![RT Clustering](images/rt_clustering.png)

# Isotopic Pattern Detection
* Instead of checking only the isotopic patterns of [M+H] or [M-H] annotated peaks, this modified version checks all adducts for isotopes.
* Only the second most abundant isotope is checked in the raw data.