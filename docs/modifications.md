# Modifications in recetox-xMSannotator

## Simple Annotation
* The input database is expected to contain only the compounds without adduct information ([example](https://umsa.cerit-sc.cz/library/list#folders/F1c84aa7fc4490e6d/datasets/8ee788c99983ff96)) and not also the adduct related information ([example](https://github.com/kuppal2/xMSannotator/blob/master/data/hmdbAllinf.rda)).
* The code computing the matching was rewritten in C++ for performance reasons and also performs the computation of expected mz value for a compound-adduct-pair (see [here](https://github.com/RECETOX/recetox-xMSannotator/blob/a410b00449db6e4bbc07339a51ef32336f405ce1/xmsannotator/src/match_by_mass.cpp#L22-L47))

## RT Clustering
* The RT-based clustering as the second step of feature grouping has been reworked as described [here](https://github.com/RECETOX/recetox-xMSannotator/pull/17)
* The clustering by this new method results in a finer grouping of peaks into more individual clusters, see the example image below.
![RT Clustering](images/rt_clustering.png)

## Isotopes Detection
* The peaks corresponding to the isotopes of a given compound are detected based on several criteria:
  1. **mass-to-charge ratio**
  2. **mass defect**,
  3. **peak intensity**,
  4. **retention time**.
* For each annotated feature a **peak table** is filtered to find peaks that: 
  1. lie in the same **RT cluster** as the feature, 
  2. have **RT** within `rt_tolerance` from the feature's **RT**,
  3. have **mass defect** within `mass_defect_tolerance` from that of the feature,
* After such peaks are identified, their intensities are compared to theoretical expected intensities of corresponding isotopes.
Peaks that differ in intensities more than `intensity_deviation_tolerance` (relative scale) from expected values are filtered out.
* Lastly, the detected isotopic peaks are appended to the **annotation** table