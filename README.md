# This repository contains:

## (1) R Package to execute the vocal diagnosability test by Isler, Isler & Whitney (1998).
### (2) Other relevant R scripts for Sin et al. (2022) in "SupplementaryMaterials_Sin-et-al._2022" folder. 

#### Installation of (1) in R

	install_github("dbssyck/tanysiptera")

#### Citation

Formulas for (1) were adapted from Isler et al. (1998). The package was prepared for Sin et al. (2022). If you reuse this package please cite:

- Isler, M. L., Isler, P. R., & Whitney, B. M. (1998). Use of vocalizations to establish species limits in antbirds (Passeriformes: Thamnophilidae). The Auk, 115(3), 577-590.

- Sin, Y. C. K., Eaton, J. A., Hutchinson, R. O., & Rheindt, F. E. (2022). Re-assessing species limits in a morphologically cryptic Australasian kingfisher lineage (Coraciiformes: Halcyonidae) using bioacoustic data. Biological Journal of the Linnean Society, 136(4), 489-505.

#### About

The vocal diagnosability test by Isler et al. (1998) (henceforth Isler test) was designed as a conservative method to examine vocal differences between populations of birds. Whlie originally crafted to study differentiation in Antbirds, the Isler test has since been adopted to inspect divergences between many other groups of birds (both Passerines and non-Passerines) across various regions.

The Isler test assigns diagnosability to pairs of populations which continuous vocal characters satisfy the following two criteria:

1) Ranges of measurements for the two populations do not overlap
2) Measurements meets the following criterion


![](https://latex.codecogs.com/svg.image?\overline{x}_{a}&plus;{t}_{a}{SD}_{a}\leq&space;\overline{x}_{b}&plus;{t}_{b}{SD}_{b})


where
![](https://latex.codecogs.com/svg.image?\overline{x}_{i}) = 
mean;
![](https://latex.codecogs.com/svg.image?{SD}_{i}) =
standard deviation; and
![](https://latex.codecogs.com/svg.image?{t}_{i}) =
t-score at 97.5 percentile of the t distribution for n-1 degrees of freedom, with population
![](https://latex.codecogs.com/svg.image?a)
having the smaller set of measurements and population
![](https://latex.codecogs.com/svg.image?b)
with larger set of measurements.

#### Usage

	vocal_diagnose(data, min.sample, prose)

data:	Dataframe containing vocal parameters. Each row should include parameters from one vocal sample (i.e., one individual). Column 1 must contain taxa/population assignment of the sample, and subsequent columns must contain measurements of one vocal parameter each.

min.sample:	Default = 2. Minimum sample size allowed for each taxa/populations for Isler test to be conducted.

prose:	Default = FALSE. If TRUE, the result will be a dataframe with rownames spelled out as "Taxa 1 vs Taxa 2" for each pairwise comparison. If FALSE, the result will be a dataframe with Taxa 1 in first column and Taxa 2 in second column for each pairwise comparison.

#### Outputs

diagnose$diagnosability:	informs whether each pairwise Isler test is diagnosable. NA = no, "1" = yes.

diagnose$inspect:	provides information on the specific criteria that was met in the test
"1": only criterion 1 was met
"2": only criterion 2 was met
"3": both criteria 1 and 2 were met (i.e., the pairs are diagnosable according to the Isler test)

diagnose$param.sum:	provides information on the number of pairs that could were diagnosable by each respective vocal parameter.

#### Example of Data Format

|Taxa|Centre.Frequency..Hz.|Bandwidth.of.first.note..Hz.|Duration.of.first.note..seconds.|Shape.of.second.note|Total.number.of.notes|Position.of.note.with.lowest.centre.frequency
| ------------- |:-------------:| -----:| -----:| -----:| -----:| -----:|
|galatea|1894.933333|648.5333333|0.749|0.10809768|33.66666667|0.038618157
|minor|1900.325|316.9875|0.457208333|0.179217632|24.91666667|0.052062848
|browningi|2257.914286|1425.571429|0.311857143|0|34.28571429|0.082531203
|boanensis|1787.25|175.05|0.029|0.369444445|33.5|0.005280042
|browningi|2402.358333|1445.858333|0.5075|0|25|0.068023823
|doris|1937.533333|1176.1|0.439|0|32.66666667|0.088103372

#### To use the test file

	test <- read.table(system.file ("extdata", "test_data.txt", package = "tanysiptera"), header = TRUE) 
