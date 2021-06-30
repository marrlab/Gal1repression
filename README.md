# Gal1repression

contains code and data accompanying 

__XXX__

Lea Schuh<sup>1,2</sup>, Poonam Bheda<sup>3</sup>, Igor Kukhtevich<sup>3</sup>, Robert Schneider<sup>3,4,5,\*</sup>, Carsten Marr<sup>1,\*</sup>

<sub><sup>
<sup>1</sup>Institute of Computational Biology, Helmholtz Zentrum München - German Research Center for Environmental Health, 85764 Neuherberg, Germany <br>
<sup>2</sup>Department of Mathematics, Technical University of Munich, Garching, 85748, Germany <br>
<sup>3</sup>Institute of Functional Epigenetics, Helmholtz Zentrum München - German Research Center for Environmental Health, 85764 Neuherberg, Germany <br>
<sup>4</sup>German Center for Diabetes Research (DZD), 85764 Neuherberg, Germany <br>
<sup>5</sup>Faculty of Biology, Ludwig-Maximilians University of Munich, 82152 Planegg-Martinsried, Germany <br>
*correspondence: robert.schneider@helmholtz-muenchen.de, carsten.marr@helmholtz-muenchen.de <br>
</sup></sub>

## Summary

Changing environments require cells to adapt their gene expression. Transcriptional memory facilitates such an adapted cellular response to the repeated exposure of a stimulus. While re-induction memory, a faster response to a repeated stimulus, has been widely studied in different organisms, it remains unknown whether there also exists memory in gene repression. 
Here, we monitor galactokinase 1 (Gal1) expression in single Saccharomyces cerevisiae cells over time to identify single-cell repression kinetics. Using ordinary differential equations, multi-start maximum likelihood optimization, quantitative model selection and statistical analysis, we are able to quantitatively describe single-cell repression kinetics. We found that the repression delay is decreased upon a repeated second repression, suggesting the existence of re-repression memory. Moreover, we found that the previously identified gain-of-re-induction memory mutant Δelp6 also shows gain-of-re-repression memory by exhibiting a stronger decrease in the repression delay between first and second repression compared to wildtype cells. 
Together, our work provides the first quantitative characterization of repression kinetics and re-repression memory at the single-cell level. <br>

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

The raw MATLAB structures containing the information extracted from segmenting single budding yeast cells can be found in the folder [Segmentation](Segmentation), where the file name ExpX_posY_Z_A summarizes: X - experiment name, Y - position, Z - strain, A - strain name.

## Data

Single cell information important for our analysis (namely cell ID, mother cell ID, detection frame, relative GFP intensities per time and cell area per time) was extracted and summarized in MATLAB structures named SExpX_posY_Z_A according to their raw MATLAB files in [Segmentation](Segmentation). The MATLAB structure containing the total GFP of all computed non-dividing cells, can also be found in the folder [Data](Data). All other outputs generated by the [main_script.m](main_script.m) are saved here eg. estimated parameters for each total GFP tra for repressions 1 and 2 and wildtype and elp6 as well as the outputs of the statistical analysis. 

## Workflow

The workflow to the whole manuscript can be found in the script [main_script.m](main_script.m).

## Figures
plotFigureX.m scripts recreate all figure panels of figure X. The created figure files as pdfs are saved to and can be found in the folder [Figures](Figures). You can also find the full figures in this folder. 

## Functions
All functions used by the [main_script.m](main_script.m) can be found here. More information on the functions can be found in the scripts themselves eg. input, output. 


