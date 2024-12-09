# Field quality control scripts

George Lu, July 2022

## Purpose

This repository contains processing scripts by George Lu used in 2022 Greenland field work. Instructions for their use are below. This folder also contains links to some data that one can test the scripts with (`data_links.md`). It also contains a spreadsheet used to calculate ApRES storage space and power usage (`ApRESCalculator.xlsx`). It is currently calculating the numbers based on the May 2022 deployment of ApRES.  


## Instructions for Installation and Setup of August 2022 check scripts

### Accessing code:
These instructions come with the files saved at: https://github.com/autonomous-phase-sensitive-radar/apres-recipes

To access the latest version of this code, go to the link above, click the button that says `Code`, and click the option to download ZIP. Unzip the folder, and then you will have the code accessible.

### Code requirements:
These scripts require the latest version of MATLAB installed, and in order to use the scripts in the field, MATLAB needs to work offline. This means making sure that the license is properly renewed (for more information regarding these licenses, go to: https://www.cuit.columbia.edu/content/matlab). Additionally, MATLAB needs to have the Signal Processing Toolbox along with the Statistics Toolbox installed. Furthermore, these scripts call upon additional utility scripts from Keith Nicholls. They are included in the Github folder with the processing scripts as of July 18, 2022. When opening MATLAB, the folder containing these utility scripts (named `nicholls_utils`) should be added to the path. 

The directory structure should look like: 
> … whatever directory the processing folder is stored in
>> /processing
>>>…field processing scripts
>>>>/nicholls_utils
>>>>>… utility scripts from Keith in here

All the relevant files outside of the main field processing scripts (of which this instruction document is a part of) should already be on the toughbook used in the May 2022 field trip to Ilulissat. Consequently, MATLAB and its relevant toolboxes should not need to be reinstalled, though it is important to check that it still works offline and that the license is valid. 

### Running the scripts in the field:
- After plugging in the SD card, go to the file `run_script.m`. This is the main script that runs all the other checkscripts. 
- It can be run right away as is, but there are some modifications that can be made to change your results. These specific modifications are detailed in the sections about their respective individual script instructions
- Pressing the run button will start the script, which will then prompt you to select a directory in a popup. Click the SD card location and press open to select the SD card as your folder. The name of the SD card should be the last item in the path. This path can also be hardcoded into the script by setting the variable `SD_card_directory` to the SD card path.
- The script takes in two optional inputs: `site_name` and `mode`
    - `site_name` allows the user to set the logfile name to refer to a specific site
    - `mode` sets whether the script looks at data taken in unattended mode or attended mode. 
- Then, in sequence, the run script will check for data gaps with `check_dates_vs_time.m`, check for clipping/too much attenuation with `check_clipping_attenuation.m`, and then check some example vertical velocity calculations with `check_vertical_velocity.m`. 
- The script will output a log file titled "log-_sitename_-_mm-dd-yy-HH-MM_.txt", and the end of the log file will contain a summary of all the tests. Here the sitename is specified by the user and will default to 'unspecified-site' otherwise. The date and time is automatically generated when the script is run. 

### Running the scripts individually:
There are 3 functions scripts to run: `check_dates_vs_time.m`, `check_clipping_attenuation.m`, and `check_vertical_velocity.m`. They are listed in order of importance/order to run. Their individual instructions are as follows.

`check_dates_vs_time.m`: This script checks to make sure there are no glaring data gaps. It has two settings that allow for either a quick coarse check (in field) or a longer fine check (in town). 
- The function requires a path to the folder, `MyFolder` (make sure that the path is terminated by a '/') and a specification for the `resolution` variable (in the field it would be set at 1 for the quick check; in town, 0). In the main script, the resolution is set to coarse (1) and the path is determined through a prompt to the user.
- The script will generate a plot marking the dates (and times if it's the thorough check) of each file created (or burst for thorough check). One can use this plot to identify if there's any lapses in measurements. 

`check_clipping_attenuation.m`: This script checks to make sure there is no clipping or to make sure that the signal is not too highly attenuated. The main script runs it twice in the field, once for clipping and once for high attenuation. 
- Many inputs are required: `MyFolder` (same as with the other scripts), `max_plots`, `file_spacing`, `burst_spacing`, `clipping`, and `amplitude`. 
- `max_plots` specifies the max number of plots to generate. This is to prevent potentially hundreds of plots popping up in the case that many measurements are flagged. It is preset to 10.
- `file_spacing` and `burst_spacing` specify how many files (roughly 1 file a day) and how many bursts (roughly 96 bursts in a file), respectively, to skip over when iterating through all the bursts. This is to speed up the script runtime in the field. These specifications are preset at 2 (file) and 20 (burst) in the main script. 
- `clipping`  specifies whether we want to detect clipping or high attenuation (prioritizing clipping, since we want conservative attenuation settings for night and winter). Setting it to 1 means we're flagging clipped data; 0 means flagging highly attenuated data. 
- If we are flagging highly attenuated data, specify the amplitude that serves as a cutoff with the variable `amplitude`. Note that clipping happens when the amplitude is greater than 1.25. The preset for this is 0.25. 
- The code will display a warning for which file and burst contains flagged data, and also generate a plot showing the raw voltage, a histogram of voltages, and a range plot generated by that specific measurement. We can use those plots to verify if the measurements are indeed clipped or too attenuated. 

`check_vertical_velocity.m`: This script does a rough check to make sure that vertical velocity measurements make sense. It will be run by the main script in the field and returns the vertical velocities of specified files and bursts.
- The code will compare the first and last file recorded on the SD card. It generates 3 plots: the first one being a vertical velocity between the first and second burst of the first file, the second one being between the first and second burst of the last file, and the third one being between the first burst from the first file and the first burst from the last file. 
- The main things to look at are the first two plots - we want to make sure that those vertical velocities are essentially zero throughout the range, as the measurements are done very close to one another. The third plot of vertical velocities should be nonzero, but there is not really an expected trend. 
- Some additional comparisons/fine tuning can be done (though probably off the ice) by using commands of similar structure to lines (46-50) specifying which bursts to compare and whether or not to average the files that contain those bursts. 
- Furthermore, the file `fmcw_process_config_vsr.m` contains many tuning settings (explained in the file itself) that can modify the results. I have left it as is for the time being. 

