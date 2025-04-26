


![alt text](https://raw.githubusercontent.com/jao321/YClon/main/yclon_logo.png)

YClon is a fast computational method for grouping B cell receptor (BCR) clones, which requires the user to provide a BCR annotated repertoire in the Minimal Information Adaptive Immune Receptor Repertoire (MiAIRR) format. This tool is particularly useful for analyzing antibody repertoires that contain a large number of sequences, such as those with more than 1,000,000 sequences, and for analyzing multiple repertoires simultaneously.

The first step of YClon involves subsetting the data based on shared V and J gene segments and CDR3 length. The next step involves decomposing the sequences into n-grams of size 3, which are then counted and used to form a vector. To determine the similarity of sequences, the cosine similarity between the vectors is calculated. In order to further enhance the accuracy of the grouping process, the n-gram frequency is normalized using the TF-IDF weighting method. 

The Python3 script implemented for that can be used in Mac OS, Linux and Windows 10. To install and run this software, please ensure that you have Python 3.8 or a later version. 

Furthermore, to use this tool, you will need to install several Python libraries that are not included in the standard Python distribution. These libraries include alive_progress, pandas, numpy and scikit-learn(1.4 or higher). To install these libraries, you can use the pip command in your terminal or command prompt on your operating system:

              
              pip3 install alive_progress pandas numpy scikit-learn
              
Copy YClon.py into a folder. 

You can download the [script file(click here with the right button of the mouse and save as...)](https://raw.githubusercontent.com/jao321/YClon/main/YClon.py) or a [folder](https://github.com/jao321/YClon/archive/refs/heads/main.zip) that contains the script, the README file and an example MiAIRR file (YClon_input_test_airr_only_essential_info.tsv). Alternatively, if you have git installed you can clone the repository with the following command line:

              
              git clone https://github.com/jao321/YClon.git


Once you have installed the required libraries, you can run YClon by providing an annotated BCR file in the MiAIRR format as input. It is important to note that the input file should have a minimum of four columns: a unique sequence identifier, V gene annotation, J gene annotation, and CDR3 sequence or another comparable column. For further guidance, please refer to the example provided below:

| sequence_id | v_call      | j_call   | cdr3                                        |
|-------------|-------------|----------|---------------------------------------------|
| 6048431_1   | IGHV1-69*19 | IGHJ3*02 | gcgagccatattttaggatattgtaatagtaccagctgctat  |
| 26048436_1  | IGHV4-34*01 | IGHJ3*02 | gcgagagtgggaacagcatattgtggtggtgactgctatgatg |


To execute YClon script, you can open a terminal window and run Python3 with the script file. 
            
              python3 path/to/YClon.py YClon_input_test_airr_only_essential_info.tsv

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
| --metric   |            "hamming"(default), "kmer"                                                                                                                                                                                                                                    | Choose if you want to run YClon with kmer comparison or hamming distance                                                                                                                             |
| --short_output |                                                                                                                                                                                                                                                          | Saves output file with only essential information:<br> Sequence identifier<br> V gene annotation<br> J gene annotation<br> Sequence that were compared<br> Clone identifier                                                                                                                                                                 |
| --folder       | Path for the folder with multiple files<br> that should be clonotyped<br> $PATH/directory_name                                                                                                                                                           | Runs YClon in every file of a specific folder                                                                                                                                                                                                                                                                                               |
| --version      |                                                                                                                                                                                                                                                          | Returns which version of YClon the user have                                                                                                                                                                                                                                                                                                |
              
Examples on how to run YClon with different arguments

Standard, but changing the threshold of dissimilarity from 0.09 to 0.15
                
                python3 YClon.py --input YClon_input_test_airr_only_essential_info.tsv --thr 0.15


Saving only the essential columns
                
                python YClon.py --input YClon_input_test_airr_only_essential_info.tsv --short_output


The output file will be written in the same folder as the input file, but with "_YClon_clonotyped" add to its name.

    Clonotyping |████████████████████████████████████████| 3059/3059 [100%] in 8.6s (356.57/s) 
    Writing output file |████████████████████████████████████████| 82928/82928 [100%] in 0.6s (129002.16/s) 
    The work was completed in: 33.000

The prompt showed at the terminal has three information:
3059 is the number of subsets that has the same V gene, J gene and CDR3 length.
82928 is the number of sequences in the repertoire.
33.000 is how many seconds it took to complete the job.

# Output file and its uses
---

A standard output file generated by YClon will be in .tsv format and will contain  two new additional columns in addition to the columns present in the original input file:
-	clone_id: This column provides a unique identifier for the clone to which a particular sequence belongs.
-	clones_seq_count : This column indicates the number of sequences present in the respective clone. 



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


For instance, in the example below, we observe that the first sequence belongs to a clone labeled as "2," which consists of a single sequence. Conversely, the second sequence is associated with the clone "11," which comprises a total of five sequences.

| sequence_id | v_call      | d_call      | j_call   | cdr3                                                                     | clone_id | clone_seq_count |
|-------------|-------------|-------------|----------|--------------------------------------------------------------------------|----------|-----------------|
| 6048431_1   | IGHV1-69*19 | IGHD2-2*02  | IGHJ3*02 | gcgagccatattttaggatattgtaatagtaccagctgctatacgccaacccttagggcggcttttgatatc | 2        | 1               |
| 6048436_1   | IGHV4-34*01 | IGHD2-21*02 | IGHJ3*02 | gcgagagtgggaacagcatattgtggtggtgactgctatgatgatgcttttgatatc                | 11       | 5               |


What kind of information you can gather from it?

The clonal diversity can be assessed by utilizing the command line tool:

                python3 YClon.py diversity YClon_input_test_airr_only_essential_info_YClon_clonotyped.tsv
                
This will output the Simpson diversity, Shannon diversity indexes for that particular repertoire and the Shannon evenesss of it.

    ---------------------------------------------
    DIVERSITY REPORT


    Simpson diversity:0.9998301997772924


    Shannon diversity: 9.600194781822921


    Shannon evenness: 0.00040002478360860543
    ---------------------------------------------


Additionally, you have the option to utilize an accompanying R script, which is also provided in this repository. However, please ensure that you have R installed on your computer beforehand. To run the R script, you can execute the following code snippet:

                Rscript rarefaction_coverage.R YClon_input_test_airr_only_essential_info_YClon_clonotyped.tsv
![alt text](https://raw.githubusercontent.com/jao321/YClon/main/YClon_input_test_airr_only_essential_info_YClon_clonotyped_RAREFACTION.png)

# Plug and play
---

You can also run YClon via the Colab notebook:

https://colab.research.google.com/drive/1vQFF8e2hxMAeDSJNooVx6zaJ2hzPVbuH?usp=sharing



# Cite this work
---
    @article{GERVASIO2023113576,
    title = {Yclon: Ultrafast clustering of B cell clones from high-throughput immunoglobulin repertoire sequencing data},
    journal = {Journal of Immunological Methods},
    volume = {523},
    pages = {113576},
    year = {2023},
    issn = {0022-1759},
    doi = {https://doi.org/10.1016/j.jim.2023.113576},
    url = {https://www.sciencedirect.com/science/article/pii/S0022175923001588},
    author = {João Gervásio and Alice Ferreira and Liza F. Felicori},
    keywords = {B cell clonotyping method, Antibody repertoire sequencing, Agglomerative clustering},
    abstract = {Motivation
    The next-generation sequencing technologies have transformed our understanding of immunoglobulin (Ig) profiles in various immune states. Clonotyping, which groups Ig sequences into B cell clones, is crucial in investigating the diversity of repertoires and changes in antigen exposure. Despite its importance, there is no widely accepted method for clonotyping, and existing methods are computationally intensive for large sequencing datasets.
    Results
    To address this challenge, we introduce YClon, a fast and efficient approach for clonotyping Ig repertoire data. YClon uses a hierarchical clustering approach, similar to other methods, to group Ig sequences into B cell clones in a highly sensitive and specific manner. Notably, our approach outperforms other methods by being more than 30 to 5000 times faster in processing the repertoires analyzed. Astonishingly, YClon can effortlessly handle up to 2 million Ig sequences on a standard laptop computer. This enables in-depth analysis of large and numerous antibody repertoires.
    Availability and implementation
    YClon was implemented in Python3 and is freely available on GitHub.}
    }

