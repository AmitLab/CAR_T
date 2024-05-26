# clinical prediction for CAR-T study 2024

instruction:
1. clone the repo
2. create python env
   use requirments file with pipy (it is recommended to use conda)
   pip install -r requirements
3. run the CAR_T_with_ltr.ipynb
   1. overview:
      1. this notebook train models with different configuration and produces .pkl with results file
      2. run this 8 times with different configurations, the configurations are just after the imports
      3. results files (.pkl) will be saved to data dir, don't delete them
   2. instruction
      1. edit general configuration
      2. edit experiment configuration, to reproduce paper follow instruction in notebook
4. run the visualization notebook
   1. after running this notebook a figure dir with execution date will be created with all figures