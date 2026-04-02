# Description

-   This repo is for my BCB430Y1 project
-   This repo is forked from [scrattch.hicat](https://github.com/AllenInstitute/scrattch.hicat) and modified with one commit [d04b812](https://github.com/DDMMMAA/scrattch.hicat-BCB430-/commit/d04b8121259ad7519733000dcc226d6143c9c10a).
-   Analysis script is in` ~/analysis_script/`
-   Result is in `~/data/Result/`
-   Report and presentation slides is in `~/report`

## Installation

`scrattch.hicat` has several dependencies, including two from BioConductor and one from Github:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")

devtools::install_github("JinmiaoChenLab/Rphenograph")
```

Once these dependencies are installed, this modified `scrattch.hicat` can be installed with:
```
devtools::install_github("DDMMMAA/scrattch.hicat")
```
