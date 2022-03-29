## This repository contains:

### 1) R Package to execute the vocal diagnosability test by Isler, Isler & Whitney (1998).
### 2) Other relevant R scripts for Sin et al. (2022) (in "Supplementary Materials" folder). 

#### Installation in R:

	install_github("dbssyck/tanysiptera")

#### Citation:

Formulas for this was adapted from Isler et al. (1998) and this R package was prepared for Sin et al. (2022). If you reuse this package please cite:

- Isler, M. L., Isler, P. R., & Whitney, B. M. (1998). Use of vocalizations to establish species limits in antbirds (Passeriformes: Thamnophilidae). The Auk, 115(3), 577-590.

- Sin, Y.C.K., Eaton, J. A., Hutchinson, R. O. & Rheindt, F. E. (2022). Re-assessing species limits in a morphologically cryptic australasian kingfisher lineage (Coraciiformes: Halcyonidae) using bioacoustic data. Biological Journal of the Linnean Society.

#### About

The vocal diagnosability test by Isler et al. (1998) (henceforth Isler criterion) was designed as a conservative method to examine vocal differences between populations of birds. Whlie originally crafted to study differentiation in Antbirds, the Isler criterion has since been adopted to inspect divergences between many other groups of birds (both Passerines and non-Passerines) across various regions.

The Isler criterion assigns diagnosability to pairs of populations which continuous vocal characters satisfy the following two criteria:

1) Ranges of measurements for the two populations do not overlap
2) Measurements meets the following criterion
x_a + t_a * SD_a <= x_b + t_b * SD_b
where subscripts refer to the two populations: (a) with smaller set of measurements and (b) with larger set of measurements. x = mean; SD = standard deviation; and t = t-score at 97.5 percentile of the t distribution for n-1 degrees of freedom.

#### Usage

	vocal_diagnose(dataframe, min.sample, verbose)

Your dataframe should contain X rows x Y columns.

Each row should include parameters from one vocal sample (i.e., one individual). The first column must contain taxa/population assignment, and subsequent columns should each contain vocal parameters.

Row and column names need not be specified, but setting column names corresponding to each vocal parameter will allow results to be read easily.








