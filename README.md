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

Human phenotyping research has historically been used to enrich our understanding of the biological correlates of neurodevelopmental disorders such as autism spectrum disorder. However, advancing beyond correlations in human cohorts presents challenges. The current study presents a comparative analysis using the SHANK3 mouse model to investigate neurometabolic changes associated with autism spectrum disorder. 
We conducted neurochemical profiling via magnetic resonance spectroscopy, which fingerprinted the metabolite changes in the SHANK3 mouse model across the cingulate cortex and thalamus. 
By selecting magnetic resonance spectroscopy, we enabled cross-species comparison and examined identical metabolites in homologous brain regions across different species.
We aimed to provide a comprehensive understanding of the neurodevelopmental alterations associated with autism spectrum disorder. 

To achieve this, our primary objective is to discover shared neurometabolic changes and loss of function in SHANK3 mice. We acquired single voxel PRESS scans in homologous brain areas—cingulate cortex and thalamus—across three genotypes: SHANK3+/+ (WT), SHANK3-/+ (HET), and SHANK3-/- (KO). We imaged mixed-sex mice (Male/Female 1:1) during adolescence (30 days) or early adulthood (70 days) using a Bruker BioSpec 11.7T with CryoProbe. 
All spectral data underwent processing via Spectroscopy Analysis Tools (SPANT) with rigorous quality control measures and visual inspection by two analysts. Effect size with unpaired Hedge’s g compared metabolite concentration of WT relative to HET and KO mice. 


## Section 1: Preprocessing

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
No conversion is needed, change format of to "twix"
```html
example: mrs_data <- read_mrs('~/filepath.dat', format = "twix")
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
spectrum_file <- paste0(file_name, "_spectrum.png")

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

input_file_path <- "t_result.txt"
output_file_path <- "amps_tCR_correction.txt"
data <- read.table(input_file_path, header = TRUE, sep = "\t")
reference_metabolite <- "tCR"
reference_index <- which(colnames(data) == reference_metabolite)
data[, -reference_index] <- data[, -reference_index] / data[, reference_index]
write.table(data, file = output_file_path, sep = "\t", row.names = FALSE)
results_file <- paste0(file_name, "amps_output.txt")

```
### Example: 
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
!(https://github.com/grandjeanlab/shank3mrs/blob/main/raw00118.png)


```

```html
basis <- sim_basis_1h_brain_press(mrs_proc)
print(basis)
[image of basis set parameters here]
```

```html
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

```html
plot(fit_res)
[put fitted/observed spectrum here]
```

```html
amps <- fit_amps(fit_res)
amps

[amps]
```

### 2C: How to automatically read from the working directory and plot fitted/observed spectrum 
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
    # Read MRS data
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

<!--- > ## Quality Control Criteria 
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

</body> ----->

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
