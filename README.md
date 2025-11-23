# **Welcome to Multi-omic Oncometabolism GPS**

**Authors:**  
[Higor Almeida Cordeiro Nogueira](https://www.researchgate.net/profile/Higor-Cordeiro-Nogueira)   
[Emanuell de Souza Rodrigues](https://www.researchgate.net/profile/Emanuell-Rodrigues-De-Souza)   
[Victor dos Santos Lopes](https://www.linkedin.com/in/victor-lopes-880604377)   
[Enrique Medina-Acosta](https://www.researchgate.net/profile/Enrique-Medina-Acosta)

---

## 🌐 **About the Project**

**Multi-omic OncoMetabolismGPS** is an interactive Shiny application developed as part of the research associated with the pre-print **A Multi-Omic Atlas of Convergent and Divergent Metabolic Regulatory Circuitries in Cancer**, available at **[BioRxiv](https://www.biorxiv.org/content/10.1101/2025.11.15.688631v1).**

The application enables exploration of multi-omic metabolic signatures and regulatory circuitries across 33 tumor types, integrating:
- Genomic, transcriptomic, epigenomic, proteomic, mutational, and phenotypic data  
- Immune contexture profiles and clinical outcomes  
- Directionality of regulatory relationships (convergent vs. divergent)

This platform supports hypothesis generation, biomarker discovery, and investigation of metabolic therapeutic vulnerabilities in cancer.

---

## 🧠 **Conceptual Framework**

<p align="center">
  <img src="https://github.com/BioCancerInformatics/Multi-omic-Oncometabolism-GPS/blob/main/Multi-omic-oncometabolismGPS/www/Figure_1.png" width="1000">
</p>

---

## 🔗 **Access the Online Application**
🚀 **Run in Browser**  
https://oncometabolismgps.shinyapps.io/Multi-omicOncometabolismGPSShiny/

---

## 🎧 **Podcast Overview**
Short discussion about the scientific motivations behind **OncoMetabolismGPS**.

<p align="center">
  <audio controls>
    <source src="https://raw.githubusercontent.com/BioCancerInformatics/Multi-omic-Oncometabolism-GPS/main/Multi-omic-oncometabolismGPS/www/Podcast.m4a" type="audio/mp4">
    Your browser does not support the audio element.
  </audio>
</p>

<p align="center">
  <a href="https://raw.githubusercontent.com/BioCancerInformatics/Multi-omic-Oncometabolism-GPS/main/Multi-omic-oncometabolismGPS/www/Podcast.m4a"
     download="Podcast.m4a"
     style="
       font-size:18px; 
       background:#2c3e50; 
       color:white; 
       padding:10px 20px; 
       border-radius:6px; 
       text-decoration:none;">
     ⬇️ Download Podcast
  </a>
</p>

---

## 🎥 **Project Presentation Video**

<p align="center">
  <video width="800" controls>
    <source src="https://raw.githubusercontent.com/BioCancerInformatics/Multi-omic-Oncometabolism-GPS/main/Multi-omic-oncometabolismGPS/www/Presentation.mp4" type="video/mp4">
    Your browser does not support the video tag.
  </video>
</p>

---

## 💻 **Run Locally in R**

To launch this tool locally in R, download **Multi-omic-oncometabolismGPS paste**, modify the path to the parent directory of the source directory, and run the code.

```r
library(shiny)
setwd("/path/to/parent/dir/of/source/")
runApp()
```

---

## 📌 Citation

If you use this repository, its datasets, analytical pipelines, figures, methods, or conceptual framework in your research, please cite:

**Nogueira, H. A. C.; Souza, E. R.; Lopes, V. S.; Medina-Acosta, E.** (2025) 
**A Multi-Omic Atlas of Convergent and Divergent Metabolic Regulatory Circuitries in Cancer.** 
Preprint. https://doi.org/10.1101/2025.11.15.688631

Citing this work supports the continued development of this multi-omic atlas.

---

## 🐞 Bug Reports

Please open an **issue** on GitHub or contact:  
📧 **[Higor Almeida Cordeiro Nogueira](higoralmeida1995@gmail.com)**  

---

## ⚙️ Tested Environment

```
R version 4.3.1 (2023-06-16)
Platform: x86_64-w64-mingw32 (64-bit)

Core Shiny Framework

shiny – 1.11.1
shinycssloaders – 1.1.0

Data Handling & Utilities

dplyr – 1.1.4
tidyr – 1.3.1
stringr – 1.5.1
tibble – 3.2.1
purrr – 1.1.0
readr – 2.1.5
readxl – 1.4.5
glue – 1.8.0
memoise – 2.0.1
qs – 0.27.3

Plotting and Visualization

ggplot2 – 3.5.2
ggsci – 3.2.0
cowplot – 1.2.0
grid – 4.3.1
gridtext – 0.1.5
igraph – 2.1.4
tidygraph – 1.3.1
ggraph – 2.2.1

survival – 3.8.3
survminer – 0.5.0
fmsb – 0.7.6
ggradar – 0.2

Interactive Tables

DT – 0.33

UCSC Xena Tools

UCSCXenaShiny – 2.2.0
UCSCXenaTools – 1.6.1

Deployment

rsconnect – 1.5.0

Development Tools

devtools – 2.4.5
```

