# Similarity hhsearch

Script to construct a graph from .hhr file from a hhsearch analysis

Dependencies :
--------------

- Python 3.5
   - Numpy 1.11.2
   - Matplotlib 1.5.3
   - Pandas 0.19.0
   - Networkx 1.11  

Usage:
------

```
usage: parse_hhmsearch.py [-h] -hhr <path> [-a <file>] [-o <path>] [-c <COVERAGE_MIN>] [-add <COLUMN_NAME> [<COLUMN_NAME> ...]]

-----------------------------------------

  _   _ _   _ ____     ____                 _
 | | | | | | |  _ \   / ___|_ __ __ _ _ __ | |__
 | |_| | |_| | |_) | | |  _| '__/ _` | '_ \| '_ \
 |  _  |  _  |  _ <  | |_| | | | (_| | |_) | | | |
 |_| |_|_| |_|_| \_\  \____|_|  \__,_| .__/|_| |_|
                                     |_|

-----------------------------------------

optional arguments:
  -h, --help            show this help message and exit

General input dataset options:
  -hhr <path>, --hhrfolder <path>
                        Path to the HHR result folder
  -a <file>, --annotationFile <file>
                        File with the information about the annotations
  -o <path>, --output <path>
                        Using <path> for output files (default: HHRFolder directory)
  -c <COVERAGE_MIN>, --coverage_min <COVERAGE_MIN>
                        The minimum coverage for the alignmenent in the results (between 0 and 1)
  -add <COLUMN_NAME> [<COLUMN_NAME> ...], --add_data <COLUMN_NAME> [<COLUMN_NAME> ...]
                        columns added in the definition file that you want to be add in the network information
```

Format annotation file
----------------------
The annotation file need to be a comma separate file with these columns:

name_of_the_hhr_file_without_extension **[TAB]** new_name_you_want

and comment begin with **#**
