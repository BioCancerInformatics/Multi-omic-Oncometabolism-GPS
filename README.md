# Welcome to Multi-omic Oncometabolism GPS Shiny

**[Higor Almeida Cordeiro Nogueira](https://www.researchgate.net/profile/Higor-Cordeiro-Nogueira), [Emanuell de Souza Rodrigues](https://www.researchgate.net/profile/Emanuell-Rodrigues-De-Souza), [Victor dos Santos Lopes](), [Enrique Medina-Acosta](https://www.researchgate.net/profile/Enrique-Medina-Acosta)**

<p align="center">
  <img src="https://github.com/CancerRCD/CancerRCDShiny/blob/bbb0be5495097c0bedfc711565c8279061b50306/www/Figure%203_HRF.png" width="1000">
</p>

**Multi-omic-Oncometabolism-GPS** is an R Shiny application designed for researchers to explore results from Omic layer- and metabolic pathway-specific signatures for pan-cancer biomarker and therapeutic target discovery.   

---

## 🔗 Useful Links
- 🔥 [Online App](https://oncometabolismgps.shinyapps.io/Multi-omic-oncometabolism-GPS/)  

---

### ▶️ Run Locally
To launch this tool locally in R, download this repository, modify the path to the parent directory of the source directory, and run the code.

```r
library(shiny)
setwd("/path/to/parent/dir/of/source/")
runApp()
```

---

## 🐞 Bug Reports

Please open an **issue** on GitHub or contact:  
📧 **[Higor Almeida Cordeiro Nogueira](higoralmeida1995@gmail.com)**  

---

## ⚙️ Tested Environment

```
R version 4.3.1 (2023-06-16)
Platform: x86_64-w64-mingw32 (64-bit)

Core Packages

UCSCXenaShiny – version 2.2.0
UCSCXenaTools – version 1.6.1
DT – built under R version 4.3.3
rsconnect – version not specified (attached, with function masking noted)
dplyr – built under R version 4.3.3
fmsb – built under R version 4.3.3
survminer – built under R version 4.3.3
survival – built under R version 4.3.3
tidyr – built under R version 4.3.3
devtools – built under R version 4.3.3
usethis – built under R version 4.3.3
memoise – built under R version 4.3.3
ggsci – built under R version 4.3.3
gridtext – built under R version 4.3.3
cowplot – version not specified (attached, with masking reported)
ggpubr – version not specified (loaded as a required package by survminer)
```

