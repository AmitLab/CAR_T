# source code for CAR-T study 2024
this repo contains source code used to analyze the data for the study:

Single Cell analysis of pre-treatment peripheral immune composition can predict CART
outcomes in patients with Large B-Cell Lymphoma

it contains R scripts to produce visualizations, to run the scripts yourself, please contact the authors of the paper to access the raw data
it also contains source code to create the clinical prediction presented in figure 4, this can be run from intermediate data representation (not raw data), from the data directory. 

## clinical prediction
instruction:
1. clone the repo
2. create python env
   python version used: python=3.9
   use requirements file with PyPi
   pip install -r requirements
3. run the fig4/prediction_notebooks/CAR_T_experiment.ipynb
   1. overview:
      1. this notebook train models with different configuration and produces .pkl with results file
      2. run this 3 times with different configurations, the configurations are just after the imports
      3. results files (.pkl) will be saved to data dir, don't delete them
   2. instruction
      1. edit general configuration (repo directory)
      2. edit experiment configuration, to reproduce paper follow instruction in notebook
4. run the visualization notebook fig4/prediction_notebooks/CAR_T_vis.ipynb
   1. edit the in the general configuration the repo directory
   2. after running this notebook a figure dir with execution date will be created with all figures
