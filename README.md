# mandy
## Model-based ANalysis of DYnamic chromatin interactions from PLAC-Seq and HiChIP data

Running this code is straight-forward, but you will need to have the output directories from running [MAPS](https://github.com/ijuric/MAPS) on your data for all replicates and also for the combined data.  
But first install the requirements. These are all R packages that can be downloaded either from CRAN or from Bioconductor:  
CRAN packages: dplyr, data.table, VGAM, timereg, yaml, argparse  
Bioconductor packages: preprocessCore, limma  
  


Generate a file in similar to the **mandy.yaml** file in this repository. In this file, you will need to specify:  
- output_filename: where the output should be written and the name for it  
- chip_filenames: comma-separated list of ChIP peaks for each of the samples you want to compare  
- bedpe_dirs: comma-separated path to the MAPS output of all replicates (the directory where these path points to, should include the **reg_raw** files output from MAPS  
- maps_peaks_filenames: Another comma-separated list! Each of these entries point to one of the bedpe files generated by MAPS (in the same directory as reg_raw files).  
- genome: This should be one of the values *hg* or *mm*. hg for human and mm for mouse (well, that's a lie. if this is set to anything but hg, Mandy will assume it's a mouse. Pun intended)  
- groups: a string like "0,0,1,1,1" (don't forget the double-quotes). This means the first two directories in the **bedpe_dirs** above are biological replicates of one sample and the other three are biological replicates of the second sample. You can use anything else instead of 0 and 1, as long as they show the relationships between the replicates and samples.  
mandy_utils_path: This is the path to the *mandy_utils.R* file that you can download from this repository. It includes all the functions.  

Once everything is set up (i.e. all requirements are installed, you have the MAPS outputs, and the YAML file is created), you can run the script by this command:  
```
Rscript -c yaml_filename
```

If you find anything interesting let us know, or we can wait until you publish it. And hopefully one day you will be able to cite us too when we are published. Good luck!


