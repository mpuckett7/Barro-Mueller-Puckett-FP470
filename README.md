# CS 470 Final Project

Authors: Gavin Barro, Mason Puckett, Beau Mueller


Running Experiments:
 -> There are testing scripts set up for multiple parallel technologies based on the launch, run, and view shell script examples on the clustor reference
 -> DO NOT execute the run scripts. Use the launch scripts to run a full suite of experiments for both the original implementation and an implementation we adjusted.
 -> You can use the view scripts to view the output of each experiment in the terminal or look through the generated text files

Config.mk Files/Setups:
 -> In CONFIG_FILE_OPTIONS we have premade config make files that change what if any parallel technology is being used when executables are ran.
 -> To select which set up you would like to use go to CONFIG-FILE-OPTIONS and copy the contents from one of the config files
 -> Then paste the contents into BOTH OpenLB/config.mk and openLB_original/config.mk
 -> the make file structure within the examples and library use config.mk at the top level so by changing the contents of those specific files the
    technology will automatically be updated and you don't need to do anything else

