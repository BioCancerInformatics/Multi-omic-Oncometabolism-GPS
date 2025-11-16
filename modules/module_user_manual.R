# --- Novo UI ---
mod_user_manual_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    tags$style(HTML("      
      .manual-text {
        font-family: 'Segoe UI', 'Roboto', sans-serif;
        font-size: 16px;
        line-height: 1.7;
        color: #333;
        padding: 25px;
        background-color: #ffffff;
        border: 1px solid #E0E0E0;
        border-radius: 12px;
        box-shadow: 0 3px 8px rgba(0,0,0,0.05);
        max-height: 750px;
        overflow-y: auto;
      }
    ")),
    
    h4("📘 User Manual – Cancer Metabolism GPS"),
    
    div(class = "manual-text", 
        pre(
          "User Manual: General Navigation in CancerMetabolismGPS

App Overview
CancerMetabolismGPS is an interactive platform for exploring omic-layer signatures and metabolic pathways across 33 cancer types. The app contains multiple tabs that support:
- Simple gene queries
- Integrative signature explorations
- Identification of biomarkers and therapeutic targets

1. Home
Access:
- Automatically opens on the “Home” tab
What you see:
- Welcome message: “Welcome to CancerMetabolismGPS”
- Brief description of platform purpose
- Image summarizing the app workflow

2. Search Your Gene
Access:
- Click “Search Your Gene” tab
Features:
- Search for any gene or mature miRNA of interest
- Check if it is in database signatures
- View associated nomenclatures
- Obtain automatic interpretation
- Download information

How to use:

Step 1: Data Entry
- Field: Enter Gene → type gene name (use official symbol recommended, e.g., RRM2)
- Button: Search Gene → click after typing
- Status Message:
  - Gene found: present in database
  - Gene not found: absent in database

Step 2: Viewing Information

Left Column: Main Results
- Corresponding Signature:
  - Displays multiomic signature and signature ID (standardized code)
- Button: Download Data → download full data in .csv

Central Column: Corresponding Nomenclatures
- Corresponding Nomenclature:
  - Lists nomenclatures containing searched gene
  - Each code summarizes omic, clinical, biological features
- Interaction:
  - Click a nomenclature to view detailed interpretation

Right Column: Signature Interpretation
- Interpretation includes:
  - Cancer type (abbreviation and full name)
  - Unique signature ID
  - Associated metabolism
  - Metabolic pathway
  - Clinical meaningful interaction
  - Type of metabolic cell death
  - Omic layers and phenotypes traits
  - Tumor vs non-tumor status
  - Prognostic association (Hazard Ratio)
  - Kaplan-Meier prognostic classifications
  - Tumor microenvironment classification
  - Tumor immune profile
- Button: Download Interpretation → download as .txt

Download Features:
- Raw data (.csv)
- Signature interpretation (.txt)

FAQs:

Q: Do I need the signature code?
A: No, just the gene name.

Q: Can I click multiple nomenclatures?
A: Yes, each click updates the interpretation.

Q: Can I download data?
A: Yes, both raw data and interpretation.

Q: What does “No signatures contain this gene” mean?
A: Gene exists but is not part of any multiomic signature.

3. Omic Layer- and Metabolic Pathway-specific Signatures
Access:
- Hover over “Omic Layer- and Metabolic Pathway-specific Signatures” in top menu
- Dropdown with options appears

Sub-tabs:
1. CNV-Specific Signatures → Copy number alterations
2. mRNA-Specific Signatures → Protein-coding gene expression
3. Methylation-Specific Signatures → DNA methylation
4. Mutation-Specific Signatures → Somatic mutations
5. miRNA-Specific Signatures → MicroRNA expression
6. Protein-Specific Signatures → Protein expression
7. Transcript-Specific Signatures → Transcripts

How to use filtering:

1. Filter 1: Omic Layer
- Auto-filled based on dataset
- Fixed but dynamic → updates if dataset changes

2. Filter 2: Cancer Type
- Displays available tumor types (CTAB)
- Selection restricts next filters

3. Filter 3: Metabolism
- Shows metabolic processes linked to signatures (e.g., carbohydrates)
- Active only after selecting cancer type

4. Filter 4: Pathway
- Lists metabolic pathways (e.g., glycolysis)
- Active only after metabolism selection

5. Filter 5: Metabolic Cell Death
- Filters by types (e.g., ferroptosis)
- Active only after pathway selection

6. Filter 6: Signature
- Lists signatures matching previous filters
- Select one → system displays integrative summary

Reset Filters:
- Button: Reset Filters → clears all selections

Download:
- Button: Download Summary Results → .tsv with summary column

Results after filtering:

- Integrative Summary:
  - Signature nomenclature
  - Associated cancer type
  - Molecular class
  - Metabolic processes and pathways

- Molecular Interaction:
  - Highlights clinically significant interactions

- Molecular Correlation:
  - Spearman correlation (ρ), p-value
  - Differential expression (tumor vs non-tumor)

- Survival Analysis:
  - Cox Regression → associations with OS, DSS, DFI, PFI
  - Kaplan-Meier → worst prognosis group

- Microenvironment & Immunophenotype:
  - Tumor microenvironment (e.g., “Hot”, “Cold”)
  - Immune profile

- Metabolic Cell Death Association:
  - Displays if relevant

- Table with Complete Data:
  - Direct review of variables and metrics

FAQs:

Q: Can I apply filters in any order?
A: No, filtering is sequential.

Q: Do I need to apply all filters?
A: Yes, to define the signature and generate summary.

Q: What does “No data available…” mean?
A: No records match selected criteria.

Q: What happens if I click “Reset Filters”?
A: All selections cleared.

4. Top Signatures
Access:
- Hover over “Top Signatures” → dropdown with options

Options:
1. Clinical Meaningful Signatures → highest clinical relevance
2. Top Ranked Signatures → based on quantitative ranking

Features:
- Visualizations and interpretations
- Quick browsing of relevant signatures
- Downloadable results

How to use:
- Select option
- Explore list and interpret data

5. Analysis and Plotting
Access:
- Hover over “Analysis and Plotting” → submenu with analyses

Sub-tabs:
1. Correlation Analysis → omic layers vs phenotypes
2. Tumor vs Normal Analysis → expression comparisons
3. Cox Analysis → survival associations
4. Survival Analysis → Kaplan-Meier by prognosis
5. Immune Infiltrates Analysis → immune context exploration

How to use:
- Select analysis
- Input signature
- View charts and results
- Export as needed

6. About Us
Access:
- Click on “About Us” tab

What you find:
- Platform developer information
- Possible contact details and acknowledgments"
        )
    )
  )
}