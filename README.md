# heavy_ion_data
Heavy Ion Data Version 1.1
Creator: Ivan Chernyshev, June 20 2017

## General Info
- Heavy Ion Data is a basic data analysis tool designed to analyze collision data from heavy ion collisions to make it easier to draw conclusions from the data 
- Depends heavily on CERN's data analysis framework ROOT
- Uses LaTeX to generate report files and to put data on an easy-to-read surface
- Outputs data as graphs in ATLAS style
- Currently, the program is incomplete, so errors may occur and all users are assumed to understand this fact 

## Techical Specifications
- my_code: the main generator of the data graphs produced by this software
- pt_combiner: a comparator which compares graohs created by my_code
- old: a storage compartment for data created by older versions of the code, as well as old veersions of the code
- atlasstyle-00-03-05: the folder with the libraries for ATLAS's style
