Magnetic Resonance Spectroscopy analysis

---
1H spectroscopy in the SHANK3 mouse model for Autism spectrum disorder 

---

<p align="center">
<img src= "https://www.smartnets-etn.eu/wp-content/uploads/2020/08/dondersru-1.png" 
alt="donderslogo" width=300>

<p align="center">
<img src="https://upload.wikimedia.org/wikipedia/commons/d/dd/Radboud_university_medical_center_logo.png?20170808070641" 
alt="radboudumclogo" width=300>
 
<p align="center">
MSc Cognitive Neuroscience Thesis
<p align="center" >Channelle Tham </p>
</p>

## Abstract
<!--- > [!IMPORTANT]\
> We ask that all users read our [legal disclaimer](https://github.com/simple-icons/simple-icons/blob/develop/DISCLAIMER.md) before using icons from Simple Icons. ---> 

Research on human phenotyping has broadened our knowledge of the genetic foundations of neurodevelopmental disorders, notably autism spectrum disorder. However, moving beyond correlations in human cohorts presents difficulties. In this study, we examine neurometabolic differences linked to autism using the SHANK3 mouse model. We use proton (1H) magnetic resonance spectroscopy to track neurochemical alterations in the transgenic mouse models. This non-invasive method measures the concentrations of in vivo neurometabolites.

We acquired single voxel PRESS scans in homologous brain areas—the cingulate cortex and thalamus—across three genotypes: SHANK3+/+ (WT), SHANK3-/+ (HET), and SHANK3-/- (KO). We imaged mixed-sex mice during adolescence (30 days post-natal) or early adulthood (70 days post-natal) using a Bruker BioSpec 11.7T scanner. We processed all spectral data using Spectroscopy Analysis Tools (SPANT), with rigorous quality control measures and visual inspection by two analysts. We compared metabolite concentrations using effect size analysis with unpaired Hedge’s g, examining how genotype affects metabolite levels in wild-type, heterozygous, and homozygous mice.

## Section 1: Preprocessing and data conversion using spec2nii

### 1A: Using spec2nii for Bruker (FID) scans 
<source ~/.bashrc>
module load anaconda3 <br>
conda activate spec2nii <br>
cd /project/4180000.24/test <br>

spec2nii bruker -m 2DSEQ 2DSEQ_FILE_or_DIR 
<p align="left" > OR </p>
spec2nii bruker -m FID FID_FILE_or_DIR
<p
</p>

```html
example: spec2nii bruker -m FID -o /project/4180000.24/test ./20221114_134901_aRi001_1_1_1/18/fid
```

### 1B: Using spec2nii for Siemens twix (.dat) scans
No conversion required, change format of to "twix"
```html
example: mrs_data <- read_mrs('~/filepath.dat', format = "twix")
```


## Section 2: Visual analysis and spectra production using SPANT
### 2A: Installing SPANT for R/Rstudio on Donders High Performance Computer Cluster (HPC)
```html
cd ~/R/x86_64-pc-linux-gnu-library/4.1`
rm -rf 00LOCK*
```

R version: 4.1.0 <br>
RStudio version: 1.4.1717 <br>
Memory requirement:  8GB-64GB <br>
Select "load preinstalled packages" <br>

```html
install.packages(“spant”)
library(spant)
```
> [!Note]\
> If there is a non-zero exit code error when installing SPANT package, delete "00LOCK-spant" folder from /R/x86_64-pc-linux-gnu-library/4.1 

### 2B: How to manually plot spectra  

 ```html
mrs_data <- read_mrs(file_name, format = "nifti")
mrs_proc <- hsvd_filt(mrs_data, xlim = c(8, 6), scale = "ppm") #|> shift(-1.90)
plot(mrs_proc, xlim = c(4, 0.5))
spectrum_file <- paste0(file_name, "_spectrum.png")

basis <- sim_basis_1h_brain_press(mrs_data)
print(basis)

stackplot(basis, xlim = c(5.5, 0.5), labels = basis$names, y_offset = 10)
fit_res <- fit_mrs(mrs_proc, basis, opts = abfit_opts(noise_region = c(6, 8)))
plot(fit_res)
stdev_FID<-fit_res$res_tab
t_stdev_FID <-t(stdev_FID)
print(t_stdev_FID)

file_path <- "stdev_output.txt"
write.table(t_stdev_FID, file = file_path, col.names = TRUE, row.names = TRUE) 

amps <- fit_amps(fit_res)
print(amps)
result <- amps
t_result <- t(result)
print(t_result) #transposes and prints concentration values as a txt file

input_file_path <- "t_result.txt"
output_file_path <- "amps_tCR_correction.txt"
data <- read.table(input_file_path, header = TRUE, sep = "\t")
reference_metabolite <- "tCR"
reference_index <- which(colnames(data) == reference_metabolite)
data[, -reference_index] <- data[, -reference_index] / data[, reference_index]
write.table(data, file = output_file_path, sep = "\t", row.names = FALSE)
results_file <- paste0(file_name, "tCR_ampsoutput.txt")

```
### 2B Example: 
```r
install.packages("spant")
library(spant)
```
```r
spant 2.18.0

Attaching package: 'spant'

The following object is masked from 'package: stats':

    sd
```

```r
mrs_data <- read_mrs('test/FID_001_18.nii.gz',format='nifti')
plot(mrs_data)
mrs_proc<- hsvd_filt(mrs_data,xlim = c(7,6),scale = 'ppm') |> shift(-1.90)
plot(mrs_proc,xlim=c(4.5,0.5)) 
[image of blue spec here]
![](https:/github.com/grandjeanlab/shank3mrs/blob/main/raw00118.png)


```

```r
basis <- sim_basis_1h_brain_press(mrs_proc)
print(basis)
[image of basis set parameters here]
```

```r
stackplot(basis, xlim = c(4, 0.5), labels = basis$names, y_offset = 5)
[basis stackplot here]
```

```html
fit_res <- fit_mrs(mrs_proc, basis, opts = abfit_opts(noise_region = c(6, 8)))
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |======================================================================| 100%
```

```r
plot(fit_res)
[put fitted/observed spectrum here]
fit_res$res_tab
[insert output here]
```

```r
amps <- fit_amps(fit_res)
amps

[amps]
```


### 2C:  Automatically importing data from the working directory to plot spectra 

 ```html
file_name <- list.files("~/project/test/codetest", pattern = "\\.nii\\.gz$", full.names = TRUE)

a <- function(file_name)
{
  mrs_data <- read_mrs(file_name, format = 'nifti')
  mrs_proc <- hsvd_filt(mrs_data, xlim = c(8, 6), scale = 'ppm') #|> shift(-1.90)
  plot(mrs_proc, xlim = c(4, 0.5))

  basis <- sim_basis_1h_brain_press(mrs_data)
  stackplot(basis, xlim = c(5.5, 0.5), labels = basis$names, y_offset = 10)
  fit_res <- fit_mrs(mrs_proc, basis, opts = abfit_opts(noise_region = c(6, 8)))
  plot(fit_res)

  fit_res$res_tab
  stdev_FID<-fit_res$res_tab
  t_stdev_FID <-t(stdev_FID)
  print(t_stdev_FID)
  file_path <- (file_name, "stdev_output.txt")
  write.table(t_stdev_FID, file = file_path, col.names = TRUE, row.names = TRUE) 

  amps <- fit_amps(fit_res)
  result <- amps
  t_result <- t(result)
  print(t_result)
  
  spectrum_file <- paste0(file_name, "_spectrum.png")
  results_file <- paste0(file_name, "amps_output.txt")

  input_file_path <- (file_name, "amps_output.txt")
  output_file_path <- (filename, "amps_tCR_correction.txt")
  data <- read.table(input_file_path, header = TRUE, sep = "\t")
  reference_metabolite <- "tCR"
  reference_index <- which(colnames(data) == reference_metabolite)
  data[, -reference_index] <- data[, -reference_index] / data[, reference_index]
  write.table(data, file = output_file_path, sep = "\t", row.names = FALSE)
  
  dev.copy(png, file = spectrum_file)
  dev.off()
  
  write.table(t_result, file = results_file, sep = "\t", row.names = FALSE)
}

for (file_name in file_list)
{
  a(file_name)
}
```

### 2D: Automatically importing data from a CSV file to plot spectra 
```html
# Read the CSV file
csv_file <- "xlsxcoding.csv" 
patient_data <- read.csv(csv_file)

# Set the working directory where the patient files are saved
setwd("~/project/test/codetest") 

# Iterate through each row in the dataframe
for (i in 1:nrow(patient_data)) 
{
  patient_id <- patient_data$file_name[i]
  phase_shift <- as.numeric(patient_data$PhaseShift[i])
  
  if (patient_data$mrsthal_exclude == "Y" | patient_data$mrspfc_exclude == “Y”)  #checks if the subject should be excluded
{
    print(paste(“skipping patient", patient_id...))
    next  # Skip to the next row
 }
  
file_name <- paste0(patient_id, ".nii.gz")
if (file.exists(file_name)) {
    mrs_data <- read_mrs(file_name, format = "nifti")
    mrs_proc <- hsvd_filt(mrs_data, xlim = c(8, 6), scale = "ppm") |> shift(-patient_data$phaseshift[i])
    plot(mrs_proc, xlim = c(4, 0.5))
    
    png(filename = paste0(patient_id, "_spectrum.png"))
    plot(mrs_proc, xlim = c(4, 0.5))

    basis <- sim_basis_1h_brain_press(mrs_data)
    print(basis)
    
    stackplot(basis, xlim = c(5.5, 0.5), labels = basis$names, y_offset = 10)
    fit_res <- fit_mrs(mrs_proc, basis, opts = abfit_opts(noise_region = c(6, 8)))
    plot(fit_res)

    fit_res$res_tab
    stdev_FID<-fit_res$res_tab
    t_stdev_FID <-t(stdev_FID)
    print(t_stdev_FID)
    file_path <- (file_name, "stdev_output.txt")
    write.table(t_stdev_FID, file = file_path, col.names = TRUE, row.names = TRUE)

    amps <- fit_amps(fit_res)
    print(amps)
    
    result <- amps
    t_result <- t(result)
    print(t_result)

write.table(t_result, file = paste0(patient_id, "_t_result.txt"), sep = "\t", quote = FALSE)

input_file_path <- (patient_id, "_t_result.txt")

data <- read.table(input_file_path, header = TRUE, sep = "\t")
reference_metabolite <- "tCR"
reference_index <- which(colnames(data) == reference_metabolite)
data[, -reference_index] <- data[, -reference_index] / data[, reference_index]
write.table(data, file = output_file_path, sep = "\t", row.names = FALSE)

write.table(t_result_corrected, file = paste0(patient_id, “_tCRcorrected.txt"), sep = "\t", quote = FALSE)
    dev.off()
  } 
else 
{
    print(paste("patient file", file_name, "not found. Skipping..."))
  }
}
```
## Section 3: Data analysis of metabolite concentration values
### 3A: Computing effect size with Hedge's g
```html
We normalised the relative concentration ratios of metabolites to the signal of total creatine (phosphocreatine + creatine = tCr), as seen in Section 2D.
```
```html
Using the website estimationstats.com, we compared the normalised values from transgenic groups with shared control, providing standardized mean differences, 95% confidence intervals, effect sizes, and standard deviations
```
### 3B: Using the modelbased package to run a contrast analysis 
```html
install.packages ("tidyverse"); install.packages ("glue"); install.packages ("knitr"); install.packages("data.table")
install.packages("lme4"); install.packages("multcomp"); install.packages("parameters"); install.packages("effectsize"); install.packages("performance"); install.packages("emmeans"); install.packages("modelbased"); install.packages("effsize"); install.packages("ggpubr"); install.packages("ggplot2"); install.packages("forestplot")

library("tidyverse"); library("glue"); library("knitr"); library("data.table"); library("lme4"); library("multcomp"); library("parameters"); library("effectsize"); library("performance"); library("emmeans"); library("modelbased"); library("ggpubr"); library("ggplot2"); library("forestplot"); library(effsize)

pfc_file <- "linmod_pfc.csv" 
pfc_data <- read.csv(pfc_file)

thal_file <- "linmod_thal.csv"
thal_data <- read.csv(thal_file)

metab_tCR <- c("CrCH2_tCR", "Ala_tCR", "Asp_tCR", "Cr_tCR", "GABA_tCR", "Glc_tCR", "Gln_tCR", "GSH_tCR", "Glu_tCR", "GPC_tCR", "Ins_tCR", "Lac_tCR", "Lip09_tCR", "Lip13a_tCR", "Lip13b_tCR", "Lip20_tCR", "MM09_tCR", "MM12_tCR", "MM14_tCR", "MM17_tCR", "MM20_tCR", "tLM20_tCR")
metab_stdev <- c("CrCH2_stdev", "Ala_stdev", "Asp_stdev", "Cr_stdev", "GABA_stdev", "Glc_stdev", "Gln_stdev", "GSH_stdev", "Glu_stdev", "GPC_stdev", "Ins_stdev", "Lac_stdev", "Lip09_stdev", "Lip13a_stdev", "Lip13b_stdev", "Lip20_stdev", "MM09_stdev", "MM12_stdev", "MM14_stdev", "MM17_stdev", "MM20_stdev", "tLM20_stdev")

#default
mod_[metab] <- lm([metab_tCR] ~ genotype + sex + age, data = df, weights = 1/[metab_stdev])
contrasts <- estimate_contrasts(mod_[metab], at= "genotype")

#example
mod_CrCH2 <- lm(CrCH2_tCR ~ genotype + sex + age, data = thal_data, weights = 1/CrCH2_stdev)
contrasts_CrCH2 <- estimate_contrasts(mod_CrCH2, at= "genotype")
```


## Contributors
<ul style=“list-style-type:circle”>
<li>  Channelle Tham (Animal data analysis)  </li>
<li>  Alejandro Rivera-Olvera (Animal data acquisition)  </li>
<li> Sabrina van Heukelum (Animal data acquisition)  </li>
<li> Andor Veltien (Lab technician/hardware)  </li>
<li> Nicolaas Puts (Human data acquisition)  </li>
<li> Viola Hollestein (Human data acquisition)  </li>
<li> Judith Homberg (Supervision)  </li>
<li> Joanes Grandjean (Daily supervision)  </li>
</ul>

## Sources 
````html
Clarke WT, Bell TK, Emir UE, Mikkelsen M, Oeltzschner G, Shamaei A, Soher BJ, Wilson M. NIfTI-MRS: A standard data format for magnetic resonance spectroscopy. Magn Reson Med. 2022. doi: 10.1002/mrm.29418.
````
