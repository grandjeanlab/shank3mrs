
## Section 1: Preprocessing and data conversion using spec2nii

### 1A: Using spec2nii for Bruker (FID) scans 
#### Activating the Anaconda environment to install spec2nii
````html
source ~/.bashrc
module load anaconda3 
conda activate spec2nii 
cd /project/4180000.24/test 
````

#### Converting Bruker scans to  NIfTI-MRS
- Format	= Bruker
- File extension = fid OR 2dseq
- SVS = Yes
- MRSI = Yes
- Automatic orientation = Yes 
  
```html
spec2nii bruker -m 2DSEQ 2DSEQ_FILE_or_DIR 
spec2nii bruker -m FID FID_FILE_or_DIR
```


### 1B: Using spec2nii for Siemens twix (.dat) scans
No conversion needed for SPANT, you can use .dat file for processing </p>
  EXAMPLE: mrs_data <- read_mrs('~/filepath.dat', format = "twix")

- Format	= Siemens
- File extension = .dat
- SVS = Yes
- MRSI = Partial handling (see Note below) 
- Automatic orientation = Yes


> [!NOTE]
> As spec2nii is not a reconstruction program, it cannot convert MRSI data. Far too little information is held in the twix headers to reconstruct arbitrary k,t-space data. However, if passed a file containing MRSI data spec2nii will attempt to create an empty NIfTI-MRS file with the correct orientation information, data shape, and header information. This empty file can then have data inserted from an offline reconstruction routine.
> Source: https://github.com/wtclarke/spec2nii

## Section 2: Spectra production using SPANT
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
> [!TIP]
> If there is a non-zero exit code error when installing SPANT package, delete "00LOCK-spant" folder from /R/x86_64-pc-linux-gnu-library/4.1 

### 2B: How to manually plot spectra  

 ```html
mrs_data <- read_mrs(file_name, format = "nifti")
mrs_proc <- hsvd_filt(mrs_data, xlim = c(8, 6), scale = "ppm") #|> shift(-1.90)
plot(mrs_proc, xlim = c(4, 0.5))
spectrum_file <- paste0(file_name, "_spectrum.png")

basis <- sim_basis_1h_brain_press(mrs_proc)
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
mrs_proc<- hsvd_filt(mrs_data,xlim = c(7,6),scale = 'ppm') |> shift(-1.90)
plot(mrs_proc,xlim=c(4,0.5))
spectrum_file <- paste0(FID_001_18, "_spectrum.png")
```
![spec00118](https://github.com/grandjeanlab/shank3mrs/blob/633d6ea44bc26910416707b100b0b4120177e577/figure/spec00118.png)


```r
basis <- sim_basis_1h_brain_press(mrs_proc)
print(basis)
```
![setpara00118](https://github.com/grandjeanlab/shank3mrs/blob/774bc6e61120d08f3f1636870d5f26772a04ab8a/setpara00118.png)


```r
stackplot(basis, xlim = c(5.5, 0.5), labels = basis$names, y_offset = 5)
```
![stackplot00118](https://github.com/grandjeanlab/shank3mrs/blob/9c1c2039ed6f01ae2600355a71edcb2efa4b7d5c/stackplot00118.png)

```r
fit_res <- fit_mrs(mrs_proc, basis, opts = abfit_opts(noise_region = c(6, 8)))
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |======================================================================| 100%

plot(fit_res)
```
![fitted00118](https://github.com/grandjeanlab/shank3mrs/blob/774bc6e61120d08f3f1636870d5f26772a04ab8a/fitted00118.png)

```r
stdev_FID<-fit_res$res_tab
t_stdev_FID <-t(stdev_FID)
print(t_stdev_FID)
file_path <- "stdev_output.txt"
write.table(t_stdev_FID, file = file_path, col.names = TRUE, row.names = TRUE) 
```

```r
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
![amps00118.png](https://github.com/grandjeanlab/shank3mrs/blob/774bc6e61120d08f3f1636870d5f26772a04ab8a/amps00118.png) 

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

