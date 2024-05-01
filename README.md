# MATLAB code to replicate main results from "Functional ultrasound imaging of freely moving adult human neural activity through an acoustically transparent cranial prosthetic"
==============================


## Project Organization
------------
    
    ├── LICENSE
    ├── LICENSE.md                      <- license
    ├── README.md                       <- The top-level README for developers using this project
    ├── docs                            <- documentation 
    ├── +w2b                            <- window 2 the brain (w2b) package for commonly used functionality
    ├── acquisition                     <- sample acquisition code for fUSI
    ├── experiments                     <- source code for experiment specific scripts 
    ├── references                      <- data dictionaries, manuals, other explanatory materials
    └── third-party                     <- third party submodules, repositories, functions, etc.

------------
## Associated dataset
Available from CaltechData: [Click here to download](https://doi.org/10.22002/f3y3k-em558). 

Packaged as .zip file in data hierarchy used by MATLAB functions. Once downloaded, extract contents from the .zip file.
   
    ├── human                                   Human data
    |   └──S{*A}R{*B}.mat                       {*A} is the session number and {*B} is the run number.
    ├── in vitro                                In vitro data
    |   └──In_vitro_Doppler_Data_STM.mat        Combined in vitro data
    ├── rodent                                  Rodent data
    |   └──YYYYMMDD                             YYYY=Year, MM=Month, DD=Day
    |       └──{Run name}                       
    |           └──Dop.mat                      Doppler data
    |           └──UF.mat                       Metadata


------------
## Other information
Tested on MATLAB R2021a and R2023a on Windows 10, and MATLAB R2021b and R2022b on Mac Sonoma 14.1.2
Please send feedback and suggestions to: [sumner.norman@gmail.com](mailto:sumner.norman@gmail.com), [crabut@caltech.edu](mailto:crabut@caltech.edu), or [wsgriggs@gmail.com](mailto:wsgriggs@gmail.com)

Zenodo archive - 
[![DOI](https://zenodo.org/badge/755770495.svg)](https://zenodo.org/doi/10.5281/zenodo.10645590)

==============================

### In publications, please reference:
Rabut, C., Norman, S.L., Griggs, W.S., Russin, J.R., Jann, K., Christopoulos, V., Liu, C., Andersen, R.A., and Shapiro, M.G. Functional ultrasound imaging of freely moving adult human neural activity through an acoustically transparent cranial prosthetic. Science Translational Medicine. _Accepted._ (See BioRxiv version [here](https://www.biorxiv.org/content/10.1101/2023.06.14.544094v1))

------------
## Getting set up

1. Clone the repository 

    ```bash
    $ git clone https://github.com/wsgriggs2/window-to-the-brain
    ```
    *note: If you are interested in contributing, please get in touch with the authors (Claire, Sumner, or Whitney). In the meantime, feel free to submit issues and/or feature requests.*

2. Download paired dataset from [CaltechData](https://doi.org/10.22002/f3y3k-em558) and unzip in a location that is convenient.

3. Run `run_me_first.m` within `.\experiments`. This will set up any necessary pathing.

4. Try out any of the scripts within the `.\experiments` folder!

## Dependencies

This repository makes use of multiple MATLAB toolboxes including Bioinformatics Toolbox, Curve Fitting Toolbox, Image Processing Toolbox, Parallel Processing Toolbox, Signal Processing Toolbox, and Statistics and Machine Learning Toolbox.


