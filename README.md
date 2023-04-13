# YClon
Python3 script used to group B cells into clones (clonotypes) to be used in Mac OS, Linux and Windows 10. By default it takes MiAIRR formated .tsv files.

Installation: Make sure you have a Python 3.7.3 or later version and pip3 installed.
              It uses the following python libraries that are not on python standard installation: alive_progress,pandas, numpy, sklearn, hdbscan and sparse_dot_topn.
              To install those libraries you can use pip command:
              
              pip3 install alive_progress pandas numpy scikit-learn
              
Copy YClon.py into a folder. 
[click here with the right button of the mouse and save as... to download YClon.py](https://github.com/jao321/YClon/raw/main/YClon.py)
You can [click here to download the repository](https://github.com/jao321/YClon/archive/refs/heads/main.zip)

              or if you have git installed you can clone the repository with the following command line
              git clone https://github.com/jao321/YClon.git

In a terminal window you can execute Python3 with YClon script. 
If you would like to use the default settings you should follow the name of the script with the filename. If the file is in the same folder as YClon.py there is no need for the path:
            
              python3 path/YClon.py YClon_input_test_airr_only_essential_info.tsv

---
# Arguments
If you would like to change the settings you can use the following arguments

| Command        | Possible inputs                                                                                                                                                                                                        | What does it do                                                                                                                                                                                                                                                                                                |
|----------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| --input        | path to the file: $PATH/YClon_input_test_airr_only_essential_info.tsv                                                                                                                                                  | Gets the path of the file with information  for clonotyping                                                                                                                                                                                                                                                    |
| --thr          | A number between 0 and 1  (default 0.09)                                                                                                                                                                                | Maximum dissimilarity between two sequences to determine if they are part of the same clonotype                                                                                                                                                                                                                   |
| --sequence     | A string with no spaces and special character  (default cdr3)                                                                                                                                                          | Gets the name of the column that contains  the sequence that will be compared                                                                                                                                                                                                                                  |
| --v_gene       | A string with no spaces and special character  (default v_call)                                                                                                                                                        | Gets the name of the column that contains  the sequence V gene annotation for that sequence                                                                                                                                                                                                                    |
| --j_gene       | A string with no spaces and special character  (default j_call)                                                                                                                                                        | Gets the name of the column that contains  the sequence J gene annotation for that sequence                                                                                                                                                                                                                    |
| --seq_id       | A string with no spaces and special character  (default sequence_id)                                                                                                                                                   | Gets the name of the column that contains  the unique identifier of for that sequence                                                                                                                                                                                                                          |
| --sep          | A character or escape sequence  (default '\t')                                                                                                                                                                         | Determines what should be used to read the input file                                                                                                                                                                                                                                                          |
| --kmer_length  | A number between 1 and 9 (default 3)                                                                                                                                                                                   | Determines the length to which the sequences should be interpreted into kmers                                                                                                                                                                                                                                  |
| --out          | Name or Path for the output file. path to the file: $PATH/output_file_name name of the file: output_file_name  By default the output file has the same name of the input file with the addition of "_YClon_clonotyped" | Gets the name that the output file should have. you should not include any extension to the name  If no path is provided, the output file will be  saved in the same directory as the input file.  It will be added "_YClon_clonotyped" to the name provided as well as the same extension of the  input file  |
| --dir_out      | Path for the directory where the output file should be saved $PATH/directory_name  By default the output file is saved in the same directory of the input file                                                         | Gets the path to where the output file  should be saved.                                                                                                                                                                                                                                                       |
| --low_memory   |                                                                                                                                                                                                                        | Runs YClon using less memory by collapsing  identical cdr3 sequences into one.  The output file has every sequence as the input as the collapsed sequences are expanded after clonotyping                                                                                                                      |
| --short_output |                                                                                                                                                                                                                        | Saves output file with only essential information: Sequence identifier V gene annotation J gene annotation Sequence that were compared Clone identifier                                                                                                                                                        |
| --folder       | Path for the folder with multiple files that should be  clonotyped $PATH/directory_name                                                                                                                                | Runs YClon in every file of a specific folder                                                                                                                                                                                                                                                                  |
| --version      |                                                                                                                                                                                                                        | Returns which version of YClon the user have                                                                                                                                                                                                                                                                   |
              
Examples on how to run YClon with different arguments

Standard, but changing the threshold of dissimilarity from 0.09 to 0.15
                
                python3 YClon.py --input YClon_input_test_airr_only_essential_info.tsv --thr 0.15


Saving only the essential columns
                
                python YClon.py --input YClon_input_test_airr_only_essential_info.tsv --short_output


The output file will be written in the same folder as the input file, but with "_YClon_clonotyped" add to its name.

# Plug and play
---

You can also run YClon via the Colab notebook:

https://colab.research.google.com/drive/1vQFF8e2hxMAeDSJNooVx6zaJ2hzPVbuH?usp=sharing
