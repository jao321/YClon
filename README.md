# YClon
Python3 script used to group B cells into clones (clonotypes) to be used in Mac OS X, Linux and Windows 10. By default it takes MiAIRR formated .tsv files.

Installation: Make sure you have a Python 3.7.3 or later version and pip3 installed.
              It uses the following python libraries that are not standard: alive_progress,pandas, numpy, sklearn, hdbscan and sparse_dot_topn. 
              In order to install those libraries you can use pip command:
              
              pip3 install alive_progress pandas numpy sklearn hdbscan sparse_dot_topn
              
Copy YClon.py into a folder. 

In a terminal you can execute Python3 with YClon script. 
If you would like to use the default settings you should follow the name of the script with the filename. If the file is in the same folder as YClon.py there is no need for the path:
            
              Python3 path/YClon.py YClon_input_test_airr_only_essential_info.tsv

If you would like to change the settings you can use the following arguments

usage: YClon.py [--input INPUT_FILE] 
                [--method {AHAM (default), HDBSCAN}]
                [--thr THRESHOLD_VALUE (0.09 as default)]
                [--sequence SEQUENCE_COLUMN_NAME (cdr3 as default)]
                [--v_gene V GENE COLUMN NAME]
                [--j_gene J GENE COLUMN NAME]
                [--seq_id SEQUENCE ID COLUMN NAME]
                [--sep SEPARATOR ON THE FILE]
                
                Python3 path/YClon.py --input YClon_input_test_airr_only_essential_info.tsv --thr 0.15 --method HDBSCAN

The output file will be written in the same folder as the input file, but with "_YClon_clonotyped" add to its name.
