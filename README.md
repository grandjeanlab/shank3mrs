# Magnetic Resonance Spectroscopy analysis

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
Human phenotyping research enriches our understanding of the biological correlates of neurodevelopmental disorders. However, progressing beyond correlations in human cohorts is challenging. Thus, a comparative analysis was performed using the SHANK3 mouse model of Autism spectrum disorder (ASD). Neurochemical profiling via magnetic resonance spectroscopy (MRS) fingerprinted the metabolite impacts in the SHANK3 mouse model across the cingulate cortex and thalamus. MRS was chosen since it facilitates trans-species comparison since identical metabolites are recorded in homologous brain regions. 

The aim was to uncover shared metabolic alterations in mice, offering a comprehensive understanding of neurodevelopmental metabolic changes. Single voxel PRESS scans were acquired in homologous brain areas—cingulate cortex and thalamus—across three genotypes: SHANK3+/+ (WT), SHANK3-/+ (HET), and SHANK3-/- (KO). Mice, of mixed sexes (Male/Female 1:1), were imaged during adolescence (30 days) or early adulthood (70 days) using a Bruker BioSpec 11.7T. All spectra underwent processing via Spectroscopy Analysis Tools (SPANT) and visual inspection by two analysts, with rigorous quality control measures. 

Effect size with unpaired Hedge’s g compared metabolite concentration of WT relative to HET and KO. A significant effect size for Inositol in the Thalamus (-0.93 [95%CI -1.63, -0.0751]) suggests a plausible difference between WT and KO. Minimal evidence was found for Glutamate's involvement in the thalamus and cingulate cortex. Unexpectedly, Inositol levels were lower in the thalamus, potentially linked to reduced thalamic inflammation in SHANK3 KO mice. Given shared metabolites, alignment between SHANK3 mouse data and human findings was hypothesized. Expanding the study to include group comparisons and a translational perspective could enhance understanding of neurometabolic alterations underlying autism spectrum disorder.


## Contributions
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

## Section 1: Preprocessing

### 1A: Using spec2nii for Bruker (FID) scans 

```html
source ~/.bashrc
module load anaconda3
conda activate spec2nii
cd /project/4180000.24/test
spec2nii bruker -m FID -o /project/4180000.24/test ./20221114_134901_aRi001_1_1_1/18/fid
```

### 1B: Using spec2nii for Siemens twix (.dat) scans
```html
No conversion is needed, change format of to "twix"
example:
mrs_data <- read_mrs('~/project/test/human/141793096301/sub-141793096301_PRESS_ACC_unsuppressed.dat', format = "twix")
```

## Section 2: Using SPANT to producing spectra 
### 2A: Installing SPANT for R/Rstudio on Donders Computer Cluster (HPC)
```html
cd ~/R/x86_64-pc-linux-gnu-library/4.1`
rm -rf 00LOCK*
rstudio # R version 4.1.0, RStudio version 1.4.1717, 8GB-64GB, load preinstalled packages 
install.packages(“spant”)
library(spant)
```
> [!Note]\
> If there is a non-zero exit code error when installing SPANT package, delete "00LOCK-spant" folder from /R/x86_64-pc-linux-gnu-library/4.1 

### 2B: How to manually plot fitted/observed spectrum with basis plot information 
 ```html
mrs_data <- read_mrs(file_name, format = "nifti")
#example --> mrs_data <- read_mrs('test/FID_001_18.nii.gz')
mrs_proc <- hsvd_filt(mrs_data, xlim = c(8, 6), scale = "ppm") #|> shift(-1.90)
plot(mrs_proc, xlim = c(4, 0.5))

basis <- sim_basis_1h_brain_press(mrs_data)
print(basis)

stackplot(basis, xlim = c(5.5, 0.5), labels = basis$names, y_offset = 10)
fit_res <- fit_mrs(mrs_proc, basis, opts = abfit_opts(noise_region = c(6, 8)))
plot(fit_res)

amps <- fit_amps(fit_res)
print(amps)

result <- amps
t_result <- t(result)
print(t_result) #transposes and prints concentration values as a txt file

spectrum_file <- paste0(file_name, "_spectrum.png")
results_file <- paste0(file_name, "_results.txt")
```

### 2C: How to automatically read from the working directory and plot fitted/observed spectrum 
 ```html
file_list <- list.files("~/project/test/codetest", pattern = "\\.nii\\.gz$", full.names = TRUE)

a <- function(file_name) {
  mrs_data <- read_mrs(file_name, format = 'nifti')
  mrs_proc <- hsvd_filt(mrs_data, xlim = c(8, 6), scale = 'ppm') #|> shift(-1.90)
  plot(mrs_proc, xlim = c(4, 0.5))
  basis <- sim_basis_1h_brain_press(mrs_data)
  stackplot(basis, xlim = c(5.5, 0.5), labels = basis$names, y_offset = 10)
  fit_res <- fit_mrs(mrs_proc, basis, opts = abfit_opts(noise_region = c(6, 8)))
  plot(fit_res)
  amps <- fit_amps(fit_res)
  result <- amps
  t_result <- t(result)
  print(t_result)
  
  spectrum_file <- paste0(file_name, "_spectrum.png")
  results_file <- paste0(file_name, "_results.txt")
  
  dev.copy(png, file = spectrum_file)
  dev.off()
  
  write.table(t_result, file = results_file, sep = "\t", row.names = FALSE)
}

for (file_name in file_list) {
  a(file_name)
}
```

### 2D: How to automatically read from a CSV file and plot fitted/observed spectrum 
 ```html
WORK IN PROGRESS
```

## Sources 
````html
Clarke WT, Bell TK, Emir UE, Mikkelsen M, Oeltzschner G, Shamaei A, Soher BJ, Wilson M. NIfTI-MRS: A standard data format for magnetic resonance spectroscopy. Magn Reson Med. 2022. doi: 10.1002/mrm.29418.
````
