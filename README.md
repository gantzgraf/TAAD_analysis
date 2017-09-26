# TAAD_analysis
Code for [some-paper] to define the frequency of causative genetic variants and phenotypic risk factors associated with a genetic aetiology of thoracic aortic aneurysm/dissection (TAAD).

## Getting Started
```bash
git clone https://github.com/superDross/TAAD_analysis
PYTHONPATH=$PYTHONPATH:/full/path/to/TAAD_analysis
# ensure input_files and output directories exist before execution
python3 TAAD_analysis
```

## Input Files
An ```input_files``` directory should exist within ```TAAD_analysis``` containing:
```
input_files/
│   
├── UK_Depth/
│   ├── depth_vs_taadx/
│   └── depth_vs_taadz/
├── Yale_Depth/
│   ├── depth_vs_taadx/
│   └── depth_vs_taadz/
│   
├── UK_All_Variants_Data.csv
├── UK_Most_Damaging_Data.csv
├── UK_Phenotype_Data.csv
├── Yale_All_Variants_Data.csv
├── Yale_Most_Damaging_Data.csv
├── Yale_Phenotype_Data.csv
└── Yale_Survival_Data_Clean.csv
```

## Data Cleaning
![](docs/data_cleaning.png?raw=true)
