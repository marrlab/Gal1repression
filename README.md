# Gal1repression

contains code and data accompanying 

__Modeling of Gal1 expression in single yeast cells reveals repression memory__

Lea Schuh<sup>1,2,3,4</sup>, Poonam Bheda<sup>3</sup>, Igor Kukhtevich<sup>3</sup>, Robert Schneider<sup>3,5,6,\*</sup>, Carsten Marr<sup>1,2,\*</sup>

<sub><sup>
<sup>1</sup>Institute of Computational Biology, Helmholtz Zentrum München - German Research Center for Environmental Health, 85764 Neuherberg, Germany <br>
<sup>2</sup>Institute of AI for Health, Helmholtz Zentrum München - German Research Center for Environmental Health, 85764 Neuherberg, Germany <br>
<sup>3</sup>Institute of Functional Epigenetics, Helmholtz Zentrum München - German Research Center for Environmental Health, 85764 Neuherberg, Germany <br>
<sup>4</sup>Department of Mathematics, Technical University of Munich, Garching, 85748, Germany <br>
<sup>5</sup>German Center for Diabetes Research (DZD), 85764 Neuherberg, Germany <br>
<sup>6</sup>Faculty of Biology, Ludwig-Maximilians University of Munich, 82152 Planegg-Martinsried, Germany <br>
*correspondence: robert.schneider@helmholtz-muenchen.de, carsten.marr@helmholtz-muenchen.de <br>
</sup></sub>

## Summary

Changing environments force cells to adapt. A faster gene expression response to a repeated stimulus, called reinduction memory, facilitates such adaptation and has been widely studied in different organisms. However, it remains unknown how memory in gene repression is mediated at the single-cell level. Here, we monitored galactokinase 1 (GAL1) expression in single Saccharomyces cerevisiae cells over 16 hours to analyze gene repression kinetics introduced by carbon source switches from galactose to glucose. We describe single-cell repression kinetics with ordinary differential equations and employ model selection to quantitatively compare different repression hypotheses. Repression delay, i.e. the time between carbon source switch and termination of Gal1 expression, is decreased upon a second repression, suggesting the existence of repression memory. Moreover, we found that the previously identified gain-of-reinduction memory mutant elp6Δ also shows gain-of-repression memory, exhibiting a stronger decrease in the repression delay between first and second repression. Together, our work provides the first mechanistic characterization of repression memory at the single-cell level. <br>

Required software and toolboxes:

- MATLAB (R2017b)
- [PESTO](https://github.com/ICB-DCM/PESTO version 1.1.9 for parameter estimation and MCMC sampling

The analysis was performed on a macOS Big Sur version 11.4. <br>

## Folder Structure

Here you can find the list of all folders contained in this repository:

- [Data](Data)
- [Figures](Figures)
- [Functions](Functions)
- [Segmentation](Segmentation)

## How to run the code

The important scripts are (i) [main_script.m](main_script.m) for running the workflow containing the data preprocessing, parameter estimation and statistical analysis (see Workflow) and (ii) plotFigureX.m to plot the single figure panels of figure X (see Figures).

## Segmentation

The raw MATLAB structures containing the information extracted from segmenting single budding yeast cells can be found in the folder [Segmentation](Segmentation), where the file name ExpX_posY_Z_A summarizes: X - experiment name, Y - position, Z - strain, A - strain name. This folder will be made available upon publication.

## Data

Single cell information important for our analysis (namely cell ID, mother cell ID, detection frame, relative GFP intensities per time and cell area per time) was extracted and summarized in MATLAB structures named SExpX_posY_Z_A according to their raw MATLAB files in [Segmentation](Segmentation). The MATLAB structure and excel sheet containing the total GFP of all computed non-dividing cells, can also be found in the folder [Data](Data). This folder will be made available upon publication.

## Workflow

The workflow to the whole manuscript can be found in the script [main_script.m](main_script.m).

## Figures
plotFigureX.m scripts recreate all figure panels of figure X. The created figure files as pdfs are saved to and can be found in the folder [Figures](Figures). You can also find the full figures in this folder. 

## Functions
All functions used by the [main_script.m](main_script.m) can be found here. More information on the functions can be found in the scripts themselves eg. input, output. 

## Results 
All outputs generated by the [main_script.m](main_script.m) apart from the data structures are saved here eg. estimated parameters for each total GFP trace for repressions 1 and 2 and wildtype and elp6 as well as the outputs of the statistical analysis.


