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


The file to be clonotyped should have at least 4 columns with information of a unique sequence identifier, V gene annotation, J gene annotation and CDR3 or other column that maybe compared. Like the example bellow

| sequence_id | v_call      | j_call   | cdr3                                        |
|-------------|-------------|----------|---------------------------------------------|
| 6048431_1   | IGHV1-69*19 | IGHJ3*02 | gcgagccatattttaggatattgtaatagtaccagctgctat  |
| 26048436_1  | IGHV4-34*01 | IGHJ3*02 | gcgagagtgggaacagcatattgtggtggtgactgctatgatg |



In a terminal window you can execute Python3 with YClon script. 
If you would like to use the default settings you should follow the name of the script with the filename. If the file is in the same folder as YClon.py there is no need for the path:
            
              python3 path/YClon.py YClon_input_test_airr_only_essential_info.tsv

---
# Arguments
If you would like to change the settings you can use the following arguments

| Command        | Possible inputs                                                                                                                                                                                                                                          | What does it do                                                                                                                                                                                                                                                                                                                             |
|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| --input        | path to the file: <br> $PATH/input_file.tsv                                                                                                                                                                                                              | Gets the path of the file with information <br> for clonotyping                                                                                                                                                                                                                                                                             |
| --thr          | A number between 0 and 1 <br> (default 0.9)                                                                                                                                                                                                              | Maximum similarity between two sequences to <br> determine if they are part of the same clonotype                                                                                                                                                                                                                                           |
| --sequence     | A string with no spaces and special character <br> (default cdr3)                                                                                                                                                                                        | Gets the name of the column that contains <br> the sequence that will be compared                                                                                                                                                                                                                                                           |
| --v_gene       | A string with no spaces and special character <br> (default v_call)                                                                                                                                                                                      | Gets the name of the column that contains <br>  the sequence V gene annotation for that sequence                                                                                                                                                                                                                                            |
| --j_gene       | A string with no spaces and special character <br> (default j_call)                                                                                                                                                                                      | Gets the name of the column that contains <br> the sequence J gene annotation for that sequence                                                                                                                                                                                                                                             |
| --seq_id       | A string with no spaces and special character <br> (default sequence_id)                                                                                                                                                                                 | Gets the name of the column that contains <br> the unique identifier of for that sequence                                                                                                                                                                                                                                                   |
| --sep          | A character or escape sequence <br> (default '\t')                                                                                                                                                                                                       | Determines what should be used to read the <br>input file                                                                                                                                                                                                                                                                                   |
| --kmer_length  | A number between 1 and 9<br> (default 3)                                                                                                                                                                                                                 | Determines the length to which the sequences <br> should be interpreted into kmers                                                                                                                                                                                                                                                          |
| --out          | Name or Path for the output file.<br> path to the file:<br> $PATH/output_file_name<br> name of the file:<br> output_file_name <br><br>  By default the output file has the<br> same name of the input file<br> with the addition  of "_YClon_clonotyped" | Gets the name that the output file should have.<br> You should not include any extension to the name. <br>  If no path is provided, the output file will be <br> saved in the same directory as the input file. <br><br> It will be added "_YClon_clonotyped" to the name<br> provided as well as the same extension of the <br>input file  |
| --dir_out      | Path for the directory where the output<br> file should be saved<br> $PATH/directory_name<br> <br>  By default the output file is <br>saved in the samedirectory of the input file                                                                       | Gets the path to where the output file <br> should be saved.                                                                                                                                                                                                                                                                                |
| --low_memory   |                                                                                                                                                                                                                                                          | Runs YClon using less memory by collapsing <br> identical cdr3 sequences into one.<br> <br>  The output file has every sequence as the input<br> as the collapsed sequences are expanded after<br> clonotyping                                                                                                                              |
| --short_output |                                                                                                                                                                                                                                                          | Saves output file with only essential information:<br> Sequence identifier<br> V gene annotation<br> J gene annotation<br> Sequence that were compared<br> Clone identifier                                                                                                                                                                 |
| --folder       | Path for the folder with multiple files<br> that should be clonotyped<br> $PATH/directory_name                                                                                                                                                           | Runs YClon in every file of a specific folder                                                                                                                                                                                                                                                                                               |
| --version      |                                                                                                                                                                                                                                                          | Returns which version of YClon the user have                                                                                                                                                                                                                                                                                                |
              
Examples on how to run YClon with different arguments

Standard, but changing the threshold of dissimilarity from 0.09 to 0.15
                
                python3 YClon.py --input YClon_input_test_airr_only_essential_info.tsv --thr 0.15


Saving only the essential columns
                
                python YClon.py --input YClon_input_test_airr_only_essential_info.tsv --short_output


The output file will be written in the same folder as the input file, but with "_YClon_clonotyped" add to its name.

# Output file and what we can do with it
---

A standard YClon output file will be a .tsv and it will look something like this. Where you have two new columns attached to the columns present in the original input file: 

- clone_id: This is a unique identifier for the clone of which that particular sequence is part
- clones_seq_count : This is the number of sequences present in that clone

| sequence_id | v_call      | d_call      | j_call   | cdr3                                                                     | clone_id | clone_seq_count |
|-------------|-------------|-------------|----------|--------------------------------------------------------------------------|----------|-----------------|
| 6048431_1   | IGHV1-69*19 | IGHD2-2*02  | IGHJ3*02 | gcgagccatattttaggatattgtaatagtaccagctgctatacgccaacccttagggcggcttttgatatc | 2        | 1               |
| 6048436_1   | IGHV4-34*01 | IGHD2-21*02 | IGHJ3*02 | gcgagagtgggaacagcatattgtggtggtgactgctatgatgatgcttttgatatc                | 11       | 5               |
| 6048436_2   | IGHV4-34*01 | IGHD2-21*02 | IGHJ3*02 | gcgagagtgggaacagcatattgtggtggtgactgctatgatgatgcttttgatatc                | 11       | 5               |
| 6048436_3   | IGHV4-34*01 | IGHD2-21*02 | IGHJ3*02 | gcgagagtgggaacagcatattgtggtggtgactgctatgatgatgcttttgatatc                | 11       | 5               |
| 6048437_1   | IGHV1-2*02  | IGHD6-19*01 | IGHJ4*02 | gcgagagaaagggggagcagtgaatggtacggcgttgactac                               | 20       | 9               |
| 6048437_2   | IGHV1-2*02  | IGHD6-19*01 | IGHJ4*02 | gcgagagaaagggggagcagtgaatggtacggcgttgactac                               | 20       | 9               |
| 6048437_3   | IGHV1-2*02  | IGHD6-19*01 | IGHJ4*02 | gcgagagaaagggggagcagtgaatggtacggcgttgactac                               | 20       | 9               |
| 6048437_4   | IGHV1-2*02  | IGHD6-19*01 | IGHJ4*02 | gcgagagaaagggggagcagtgaatggtacggcgttgactac                               | 20       | 9               |
| 6048437_5   | IGHV1-2*02  | IGHD6-19*01 | IGHJ4*02 | gcgagagaaagggggagcagtgaatggtacggcgttgactac                               | 20       | 9               |
| 6048437_6   | IGHV1-2*02  | IGHD6-19*01 | IGHJ4*02 | gcgagagaaagggggagcagtgaatggtacggcgttgactac                               | 20       | 9               |
| 6048437_7   | IGHV1-2*02  | IGHD6-19*01 | IGHJ4*02 | gcgagagaaagggggagcagtgaatggtacggcgttgactac                               | 20       | 9               |
| 6048437_8   | IGHV1-2*02  | IGHD6-19*01 | IGHJ4*02 | gcgagagaaagggggagcagtgaatggtacggcgttgactac                               | 20       | 9               |
| 6048437_9   | IGHV1-2*02  | IGHD6-19*01 | IGHJ4*02 | gcgagagaaagggggagcagtgaatggtacggcgttgactac                               | 20       | 9               |
| 6048439_1   | IGHV1-8*01  | IGHD3-10*01 | IGHJ5*02 | gcgagagctccccccggacgaggtcggggaataatgcaagccggggtctggttcgacccc             | 78       | 2               |
| 6048441_1   | IGHV3-23*01 | IGHD4-17*01 | IGHJ6*03 | gcgaaagatcacggttccaacccctactactactactacatggacgtc                         | 84       | 8               |


As an example showed below we can see the the first sequence is on clone named "2" and it is a clone with only one sequence.
Meanwhile the second sequence is part of the clone "11" and this clone has 5 sequences

| sequence_id | v_call      | d_call      | j_call   | cdr3                                                                     | clone_id | clone_seq_count |
|-------------|-------------|-------------|----------|--------------------------------------------------------------------------|----------|-----------------|
| 6048431_1   | IGHV1-69*19 | IGHD2-2*02  | IGHJ3*02 | gcgagccatattttaggatattgtaatagtaccagctgctatacgccaacccttagggcggcttttgatatc | 2        | 1               |
| 6048436_1   | IGHV4-34*01 | IGHD2-21*02 | IGHJ3*02 | gcgagagtgggaacagcatattgtggtggtgactgctatgatgatgcttttgatatc                | 11       | 5               |


What kind of information you can gather from it?

We can attest the clonal diversity via the command line:
                python3 YClon.py diversity YClon_input_test_airr_only_essential_info_YClon_clonotyped.tsv
                
This will output the Simpson diversity, Shannon diversity indexes for that particular repertoire and the Shannon evenesss of it.

    ---------------------------------------------
    DIVERSITY REPORT


    Simpson diversity:0.9998301997772924


    Shannon diverity: 9.600194781822921


    Shannon eveness: 0.00040002478360860543
    ---------------------------------------------

You can also use an accessory R script that is also available in this repository
To do so, you can use the following code:

                Rscript rarefaction_coverage.R YClon_input_test_airr_only_essential_info_YClon_clonotyped.tsv
![alt text](https://raw.githubusercontent.com/jao321/YClon/main/YClon_input_test_airr_only_essential_info_YClon_clonotyped_RAREFACTION.png)

# Plug and play
---

You can also run YClon via the Colab notebook:

https://colab.research.google.com/drive/1vQFF8e2hxMAeDSJNooVx6zaJ2hzPVbuH?usp=sharing
