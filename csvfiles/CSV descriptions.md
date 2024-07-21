## Description of contents in CSV files 

### Name	
- FID + Subject ID + Thalamus/CC scan number 

### Genotype	
- WT = Wild-Type (SHANK3+/+)
- HET = Heterozygous (SHANK3-/+)
- KO = Homozygous Knockout (SHANK3-/-)

### Age	
- PND30 = Postnatal days 30, adolescence
- PND70 = Postnatal days 70, adulthood 

### Sex	
- Male or female
  
### Weight	
- Mice were weighed in grams prior to the data acquisition session
  
### Rotation	
- Used to adjust the phase orientation of the chemical shift axis or the spectrum itself to correct for distortions or misalignments

### Shift	
- Represents the amount by which the spectral data is shifted along the chemical shift axis, typically in parts per million (ppm).
- A negative value means the spectrum will be shifted to a lower ppm value (to the left).

### Thal/CC file
- This corresponds the the number of the thalamus/cc file in the database 

### Mrsthal/cc exclude
Unrecoverable scans were excluded from the analysis using a binary exclusion scheme
- Y = yes
- N = no

### Exclusion reason 
Explains the reason why the scan was excluded from the study
- dead = died during scan preparation or during data acquisition
- inc = animal was included
- ppr = poor peak resolution
- lac = abnormally high lactate (Lac) peaks
- naa = abnormally low NAA peaks

> [!Note]\
> The "E" written in the name column indicates exclusion, and the exclusion reason column specifies the exact reason for exclusion

### Metabolite abbreviation  
- CrCH2 = Creatine methylene
- Ala = Alanine
- Asp = Aspartate
- Cr = Creatine
- GABA = Gamma-Aminobutyric Acid
- Glc = Glucose
- Gln = Glutamine
- GSH = Glutathione
- Glu = Glutamate
- GPC = Glycerophosphocholine
- Ins = Inositol
- Lac = Lactate
- Lip09 = Lipid 09
- Lip13a = Lipid 13a
- Lip13b = Lipid 13b
- Lip20 = Lipid 20
- MM09 = Metabolite Marker 09
- MM12 = Metabolite Marker 12
- MM14 = Metabolite Marker 14
- MM17 = Metabolite Marker 17
- MM20 = Metabolite Marker 20
- tLM20 = Total Lipid Marker 20

### Metabolites_stdev
- Refers to the absolute metabolite concentration divided by the standard deviation

### Metabolites_tCr
- Refers to the absolute metabolite concentration corrected by the corresponding tCr value


