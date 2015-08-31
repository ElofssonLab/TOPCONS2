TOPCONS2
========

Github for the TOPCONS2

This is the standalone version of web-server http://topcons.net
This software package is supposed to be run on Ubuntu x64 system.
It might also work on other Linux boxes but have not been tested.

If you are interested in running TOPCONS2 on other systems, please contact
Arne Elofsson (arne@bioinfo.se)

Installation and usage:

1. Check out the software from the github by

    $ git clone https://github.com/ElofssonLab/TOPCONS2

2. Download the database for TOPCONS2 from
    http://topcons.net/static/download/topcons2.0_database.zip
   and unzip it by 
    $ unzip topcons2.0_database.zip

3. Change to the folder topcons2_webserver and create a soft link to the
   downloaded database
    $ ln -s /path/to/the/downloaded/database database

4. Install dependencies if not installed
    a)    perl-Moose (e.g. sudo apt-get install perl-Moose)
    b)    bioperl    (e.g. cpan > install  CJFIELDS/BioPerl-1.6.924.tar.gz )
    c)    IPC        (e.g. cpan > install IPC::Run)
    d)    kalign     (e.g. sudo apt-get install kalign2)
    e)    hmmer3.0   (note that hmmscan should be compatible with the pfam database
                      otherwise, you may encounter format incompatible problem)

5. Test topcons2 workflow
   change to the folder test
   run the command

    $ ../workflow/pfam_workflow.py t1.fa myout1 ../tools/blast-2.2.26/ ../database/blast/uniref90.fasta

   You will find the results in the folder myout1
