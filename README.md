# clinical prediction for CAR-T study 2024

instruction:
1. clone the repo
2. create python env
   python version used: python=3.9
   use requirements file with PyPi
   pip install -r requirements
3. run the CAR_T_experiment.ipynb
   1. overview:
      1. this notebook train models with different configuration and produces .pkl with results file
      2. run this 3 times with different configurations, the configurations are just after the imports
      3. results files (.pkl) will be saved to data dir, don't delete them
   2. instruction
      1. edit general configuration (repo directory)
      2. edit experiment configuration, to reproduce paper follow instruction in notebook
4. run the visualization notebook
   1. edit the in the general configuration the repo directory
   2. after running this notebook a figure dir with execution date will be created with all figures