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


## Section 1: Preprocessing

### 1A: Using spec2nii for Bruker (FID) scans 

<source ~/.bashrc><br>
module load anaconda3 <br>
conda activate spec2nii <br>
cd /project/4180000.24/test <br>

spec2nii bruker -m 2DSEQ 2DSEQ_FILE_or_DIR 
<p align="left" > OR </p>
spec2nii bruker -m FID FID_FILE_or_DIR
<p
</p>

```html
spec2nii bruker -m FID -o /project/4180000.24/test ./20221114_134901_aRi001_1_1_1/18/fid
```

### 1B: Using spec2nii for Siemens twix (.dat) scans
```html
No conversion is needed, change format of to "twix"
example:
mrs_data <- read_mrs('~/filepath.dat', format = "twix")
```


## Section 2: Using SPANT to producing spectra 
### 2A: Installing SPANT for R/Rstudio on Donders Computer Cluster (HPC)
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

### 2B: How to manually plot fitted/observed spectrum with basis plot information 
 ```html
mrs_data <- read_mrs(file_name, format = "nifti")
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
### 2B Example: 
```html
install.packages("spant")
library(spant)
```
```html
spant 2.18.0

Attaching package: 'spant'

The following object is masked from 'package: stats':

    sd
```

```html
mrs_data <- read_mrs('test/FID_001_18.nii.gz',format='nifti')
plot(mrs_data)
mrs_proc<- hsvd_filt(mrs_data,xlim = c(7,6),scale = 'ppm') |> shift(-1.90)
plot(mrs_proc,xlim=c(4.5,0.5)) 
[image of blue spec here]
```

```html
basis <- sim_basis_1h_brain_press(mrs_proc)
print(basis)
[image of basis set parameters here]
```

```html
stackplot(basis, xlim = c(4, 0.5), labels = basis$names, y_offset = 5)
[stackplot here]
```

```html
fit_res <- fit_mrs(mrs_proc, basis, opts = abfit_opts(noise_region = c(6, 8)))
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |======================================================================| 100%
```

```html
plot(fit_res)
[put fitted/observed spectrum here]
```

```html
amps <- fit_amps(fit_res)
amps

X.CrCH2   3.049333
Ala       0.000000
Asp       1.944614
Cr        5.870648
GABA     14.704498
Glc       0.000000
Gln       0.000000
GSH       2.952613
Glu      34.348136
GPC       2.887518
Ins      10.837776
Lac       1.507300
Lip09     0.000000
Lip13a    4.604497
Lip13b    0.000000
Lip20    14.135484
MM09     14.331642
MM12      3.456480
MM14     30.517761
MM17     15.648993
MM20     89.169940
NAA       4.593830
NAAG      2.316995
PCh       2.194539
PCr       4.701732
sIns      0.000000
Tau      33.200250
tNAA      6.910825
tCr      10.572380
tCho      5.082057
Glx      34.348136
tLM09    14.331642
tLM13    38.578739
tLM20   103.305424
```

### 2C: How to automatically read from the working directory and plot fitted/observed spectrum 
 ```html
file_list <- list.files("~/project/test/codetest", pattern = "\\.nii\\.gz$", full.names = TRUE)

a <- function(file_name)
{
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

### 2D: How to automatically read from a CSV file and plot fitted/observed spectrum 
 ```html
WORK IN PROGRESS
```
### 2E: How to manually standardize metabolite values by tCr
```html
input_file_path <- "amps_output.txt"
output_file_path <- "amps_tCR_correction.txt"
data <- read.table(input_file_path, header = TRUE, sep = "\t")
reference_metabolite <- "tCR"
reference_index <- which(colnames(data) == reference_metabolite)
data[, -reference_index] <- data[, -reference_index] / data[, reference_index]
write.table(data, file = output_file_path, sep = "\t", row.names = FALSE)
```

## Quality Control Criteria 
insert green, red and orange spectrums and explanations on how to QC
<table border="1">
  <tr>
    <th></th>
    <th>Spectrum Example</th>
    <th>Notes</th>
  </tr>
  <tr>
    <td>Green</td>
    <td>...</td>
  </tr>
  <tr>
    <td>Orange</td>
    <td>...</td>
  </tr>
  <tr>
    <td>Red</td>
    <td>...</td>
  </tr>
</table>

</body>

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
