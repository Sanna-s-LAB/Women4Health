This folder contains the scripts necessary to reproduce the statistical analyses, as well as two additional folders:

1) data: This folder stores all the data required for the analyses. It is also the location where some of the data generated by the pipeline are saved. You will find some datasets with just the rownames and colnames but no individual data for the starting mandatory 
data and empty directories where some data will be saved by the pipeline. It's important to keep the structure of the folder as it is.

2) results: This folder will contain the results of all analyses, including both images and tables. The results are organized into subfolders named after the respective analyses. They are empty but the results of the pipeline will be saved there. It's important to keep the structure of the folder as it is.


Keep in mind that:

- The paths in each script are written according to "complete_pipeline" folder being your working directory. The script "8. tsne.R" and "3. adonis_analyses.R" 
  do not save all the plots automatically so you will be able to generate all plots but some of them ought to 
  be saved manually.

- The script independence_tests.R must not be run directly but must be run
  using the bash version (2. independence_tests.sh) so they must be in the same folder

- Likewise, the script 4. bacterial_residual.py must be run by terminal by the command:
  
  python "4. bacterial_residual.py" True
  
  or
  
  python "4. bacterial_residual.py" False
  
  depending on whether you want to filter the samples by prevalence and mean abundance or not
  (if True the filter will be applied). Check the comments in the script -- you should pay 
  attention to how the script is linked to independence test results.
  Moreover, the script needs the script utilities.py in the same folder.
