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

The software is open source and licensed under the GPL license.

## Reference
Tsirigos, K.D., Peters, C., Shu, N., Kall, L., Elofsson, A., 2015. The TOPCONS
web server for consensus prediction of membrane protein topology and signal
peptides. Nucleic Acids Res. 43, W401-W407

## Installation and usage:

1. Check out the software from the github by

        git clone https://github.com/ElofssonLab/TOPCONS2

2. Download the database for TOPCONS2 from
    http://topcons.net/static/download/topcons2_database.zip
   and unzip it by 

        unzip topcons2_database.zip

3. Change to the folder 'topcons2_webserver' and create a soft link to the
   downloaded database

        ln -s /path/to/the/downloaded/database database

    You may need to update the BLAST database for uniref90 by the script [install_blastdb.sh](topcons2_webserver/install_blastdb.sh) with the following command:


        ./install_blastdb.sh database

4. Install dependencies if not installed

    *    Cmake      (for installation of modhmm, e.g. sudo apt-get install cmake)
    *    perl-Moose (e.g. sudo apt-get install perl-Moose)
    *    bioperl-1.6.924   (e.g. cpan > install  CJFIELDS/BioPerl-1.6.924.tar.gz )
    *    biopython  (e.g. sudo pip install biopython )
    *    IPC        (e.g. cpan > install IPC::Run)
    *    kalign     (e.g. sudo apt-get install kalign)
    *    hmmer3.0   (note that hmmscan should be compatible with the pfam database
                     otherwise, you may encounter format incompatible problem
                     http://hmmer.org/download.html)
    *    Gnuplot    (e.g. sudo apt-get install gnuplot)
    *    Java       (make sure that the command java is in the PATH, 
                     e.g. sudo apt-get install default-jre)
    *    convert from ImageMagick (e.g. sudo apt-get install imagemagick)
    *    xsltproc   (e.g. sudo apt-get install xsltproc)
    *    gengetopt  (e.g. sudo apt-get install gengetopt)
    *    awk, sort, head

    Note that the commands hmmscan, kalign, gnuplot, convert, sort, awk and
    head should be in the PATH

5. Install modhmm

   change to the folder `topcons2_webserver/predictors/source/modhmm`

        bash fresh_install.sh /path/to/topcons2_webserver/predictors

6. Test the topcons2 workflow

   change to the folder `topcons2_webserver/test`
   and run the following commands 

        ../workflow/pfam_workflow.py one_seq.fasta rst1 ../tools/blast-2.2.26/ ../database/blast/uniref90.fasta

        ../workflow/pfam_workflow.py multiple_seqs.fasta rst2 ../tools/blast-2.2.26/ ../database/blast/uniref90.fasta

   The example results can be found in the folder `rst_one_seq` and
   `rst_multiple_seqs` for the example fasta file `one_seq.fasta` and
   `multiple_seqs.fasta` respectively.

* Description of the output results
    If the input contains one sequence and the out_path is "one_seq"
    The tree view of all output files will be

    <pre>
    └── one_seq
    ├── seq_0
    │   ├── DG1.txt
    │   ├── dg.txt
    │   ├── Homology
    │   │   ├── query.fa.total_aligns
    │   │   └── query.top
    │   ├── nicetop.html
    │   ├── OCTOPUS
    │   │   ├── NN_PRF_FILES
    │   │   │   └── query.prf
    │   │   └── query.top
    │   ├── philius
    │   │   └── query.top
    │   ├── PolyPhobius
    │   │   └── query.top
    │   ├── SCAMPI_MSA
    │   │   └── query.top
    │   ├── seq.fa
    │   ├── SPOCTOPUS
    │   │   ├── NN_PRF_FILES
    │   │   │   └── query.nnprf
    │   │   └── query.top
    │   └── Topcons
    │       ├── reliability.final
    │       ├── reliability.txt
    │       ├── topcons.gnu
    │       ├── topcons.large.png
    │       ├── topcons.png
    │       ├── topcons.top
    │       ├── total_image.gnu
    │       ├── total_image.large.png
    │       └── total_image.png
    └── time.txt
    </pre>


## Run TOPCONS2 with Docker

1. Download the database for TOPCONS2 from
    http://topcons.net/static/download/topcons2_database.zip
   and unzip it by 

        unzip topcons2_database.zip

    and saved as `/data/topcons2_database`

    Update the BLAST database for uniref90 by the script [install_blastdb.sh](topcons2_webserver/install_blastdb.sh) with the following command:

        ./install_blastdb.sh /data/topcons2_database/blast

2. Pull the Docker image by 

        docker pull nanjiang/topcons2

    or you can also build the Docker image locally by 

        docker build -t topcons2 .

    within the cloned folder, i.e. `TOPCONS2/`

3. Run Docker container by 

        docker run -v /data/topcons2_database:/data/topcons2_database -it nanjiang/topcons2

4. Test run

        cd /app/topcons2/test

        ../run_topcons2.sh one_seq.fasta -outpath rst1

    The result will be available at `rst1/one_seq/seq_0`

    If you want to run TOPCONS2 Docker from outside of the container, suppose
    you want to output the result to the folder `/scratch`, then start the
    container by 

        docker run -e USER_ID=$(id -u $USER) -v /data/topcons2_database:/data/topcons2_database -v /scratch:/scratch -it  --name topcons2 -d nanjiang/topcons2

    Now you can test run the TOPCONS2 using the example sequence provided in
    the package by

        docker exec --user user topcons2 script /dev/null -c "/app/topcons2/run_topcons2.sh /app/topcons2/test/one_seq.fasta -outpath  /scratch/rst1"

    To use you own sequence, you can copy your query sequence, e.g. `yourseq.fasta` to `/scratch` and then run
    the following command in the shell terminal

        docker exec --user user topcons2 script /dev/null -c "/app/topcons2/run_topcons2.sh /scratch/yourseq.fasta -outpath  /scratch/rst_yourseq"

    The prediction result will be output to `/scratch/rst_yourseq` given successful
    run.

    **Note** also that the output files are owned by your current user, i.e. `$USER`

## Run TOPCONS2 with Singularity

You can also run TOPCONS2 with Singularity 

1. Database preparation: do the same as step 1 for when running with Docker.

2. Pull the Docker image by 

        singularity pull docker://nanjiang/topcons2

3. Suppose you have saved the database at `/data/topcons2_database` and you
   have a user writable folder `/scratch`, then you can run TOPCONS2 with the
   following command

        singularity exec  -B /data:/data -B /scratch:/scratch topcons2.img /app/topcons2/run_topcons2.sh /app/topcons2/test/one_seq.fasta -outpath /scratch/rst1

   The prediction result will be output to `/scratch/rst1` given successful run.


## Run only the sub-predictors OCTOPUS and SPOCTOPUS
If you only need to run OCTOPUS and SPOCTOPUS within the TOPCONS2 package, you
need to install the whole package using the procedure described above for
TOPCONS2. Then, use the script `pfam_workflow_octopus.py`.

Examples:
change to the folder 'topcons2_webserver/test' and run the following commands 

        ../run_octopus.sh -outpath rst2 multiple_seqs.fasta

The result of predicted topologies in Fasta format can be found in
`rst1/multiple_seqs.OCTOPUS.topfa` and  `rst1/multiple_seqs.SPOCTOPUS.topfa`

If you do not need the individual output files nor the ANN output, you can
run the commands with the "-remove-individual-files" flag, that is

        ../run_octopus.sh -outpath rst2 -RM multiple_seqs.fasta

## Run only the sub-predictors OCTOPUS and SPOCTOPUS using Docker

Run OCTOPUS and SPOCTOPUS using Docker is similar to run the whole TOPCONS2
workflow.

start the container by (assume you have the write permission to the folder
`/scratch`)

        docker run -e USER_ID=$(id -u $USER) -v /data/topcons2_database:/data/topcons2_database -v /scratch:/scratch -it  --name topcons2 -d nanjiang/topcons2

And make your test run by

        docker exec --user user topcons2 script /dev/null -c "/app/topcons2/run_octopus.sh /app/topcons2/test/one_seq.fasta -outpath  /scratch/rst2"

