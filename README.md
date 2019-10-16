# DAWG-Fall-2019

DAWG will be working with a 16S rRNA amplicon dataset provided by Dr. Judie Howrylak, a pulmonologist with Penn State Hershey. 
In total there are ~ 207 samples taken from 84 patients presenting with the same illness when entering the hospital. Sequencing was performed with the 515/806 primer set. 

60 patients had 2 samples taken from them, stool and trachael at day 1
3 Patient had one sample taken from them, either stool or trachael at day 1

18 patients had 4 samples taken from them, stool and trachael at day 1 and again at day 5
4 patients had 3 samples taken from them, stool and trachael at day 1 and either stool or trachael at day 5

The goal of this project is to determine if there are biomarkers that can be used to predict patient outcome (survival or mortality) once they are admitted to see the researcher. It  is unclear which set of samples, stool or tracheal, that will provide us with better predictions but Dr. Howrylak hypothesizes that the tracheal swabs will be better predictors. 

Metadata for this dataset includes:
gender: male = 1, female = 0
Race: white = 0, black = 1, hispanic = 2
Age = numerical
ARDS: no ARDS = 0, diagnosis with ARDS = 1 (acute respiratory distress syndrome)
SIRS: no SIRS = 0, SIRS = 1 (systemic inflammatory response syndrome)
Sepsis: no sepsis = 0, sepsis = 1; 30 samples from patients who did not present with sepsis
Sever sepsis: no severe sepsis = 0, sepsis = 1
Shock: no sepsis shock = 0, sepsis shock = 1
APACHE: numerical value for APACHE II disease severity score (Acute Physiologic Assessment and Chronic Health Evaluation)
inhosp: no in-hospital mortality = 0, in-hospital mortality = 1
Month: no mortality in 1 month = 0, mortality in 1 month = 1; 109 samples from patients who had month = 0, 99 samples from patients who had month = 1
Year: no mortality in 1 year = 0, mortality in 1 year = 1
Death: number of death free days; 69 samples from patients who died within 14 days of being admitted; 
ICU: number of ICU-Free days
Vent: number of ventilator-free days
Location: location of the patient prior to admission to the ICU, home = 0, nursing home = 1, rehab facility = 2, outside hospital = 3


Packages that you will need installed in R include:
install.packages("BiocManager")
vegan
metacoder
selbal
phyloseq
randomForest
rfUtilities
caret
plyr
tibble
Hmisc
devtools
ecodist
RVAideMemoire

Some packages will need to installed using the following:
BiocManager::install("package")

To install selbal use this code:
library(devtools) 

install_github(repo = "UVic-omics/selbal")

