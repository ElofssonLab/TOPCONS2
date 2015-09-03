## TOPCONS2

This is the standalone version of web-server http://topcons.net.
This software package is supposed to be run on Ubuntu x64 system.
It might also work on other Linux boxes but have not been tested.

If you are interested in running TOPCONS2 on other systems, please contact
Arne Elofsson (arne@bioinfo.se)

## Description

TOPCONS2 is an updated version of the widely used TOPCONS for predicting
membrane protein topologies using consensus prediction.  It is faster yet
more accurate than the old TOPCONS according to our solid benchmarking.
Moreover, it predicts not only the trans-membrane helices, but also the
location of signal peptide

TOPCONS2 may be used for the purpose of academic research only.

## Reference
Tsirigos, K.D., Peters, C., Shu, N., Kall, L., Elofsson, A., 2015. The TOPCONS
web server for consensus prediction of membrane protein topology and signal
peptides. Nucleic Acids Res. 43, W401-W407

## Installation and usage:

1. Check out the software from the github by

    `$ git clone https://github.com/ElofssonLab/TOPCONS2`

2. Download the database for TOPCONS2 from
    http://topcons.net/static/download/topcons2_database.zip
   and unzip it by 

    `$ unzip topcons2_database.zip`

3. Change to the folder 'topcons2_webserver' and create a soft link to the
   downloaded database

    `$ ln -s /path/to/the/downloaded/database database`

4. Install dependencies if not installed

    *    perl-Moose (e.g. sudo apt-get install perl-Moose)
    *    bioperl    (e.g. cpan > install  CJFIELDS/BioPerl-1.6.924.tar.gz )
    *    IPC        (e.g. cpan > install IPC::Run)
    *    kalign     (e.g. sudo apt-get install kalign2)
    *    hmmer3.0   (note that hmmscan should be compatible with the pfam database
                     otherwise, you may encounter format incompatible problem)

5. Test the topcons2 workflow

   change to the folder 'topcons2_webserver/test'
   and run the following commands 

    `$ ../workflow/pfam_workflow.py one_seq.fasta rst1 ../tools/blast-2.2.26/ ../database/blast/uniref90.fasta`

    `$ ../workflow/pfam_workflow.py multiple_seqs.fasta rst2 ../tools/blast-2.2.26/ ../database/blast/uniref90.fasta`

   The example results can be found in the folder 'rst_one_seq' and
   'rst_multiple_seqs' for the example fasta file one_seq.fasta and
   multiple_seqs.fasta respectively.

* Description of the output results
    If the input contains one sequence and the out_path is "one_seq"
    The tree view of all output files will be

    <pre>
    one_seq
    |-- seq_0              # if there are more than one sequence in the input
                           # file, the subfolders will be named as seq_1, seq_2 ...
    |   |-- DG1.txt
    |   |-- OCTOPUS
    |   |   `-- query.top
    |   |-- PolyPhobius
    |   |   `-- query.top
    |   |-- SCAMPI_MSA
    |   |   `-- query.top
    |   |-- SPOCTOPUS
    |   |   `-- query.top
    |   |-- Topcons
    |   |   |-- reliability.final
    |   |   |-- reliability.txt
    |   |   |-- topcons.gnu
    |   |   |-- topcons.large.png
    |   |   |-- topcons.png
    |   |   |-- topcons.top          # the predicted topology for TOPCONS
    |   |   |-- total_image.gnu
    |   |   |-- total_image.large.png
    |   |   `-- total_image.png
    |   |-- dg.txt
    |   |-- nicetop.html
    |   |-- philius
    |   |   `-- query.top
    |   `-- seq.fa
    `-- time.txt
    </pre>

