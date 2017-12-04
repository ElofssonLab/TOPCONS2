FROM nanjiang/common-ubuntu
LABEL maintainer "Nanjiang Shu (nanjiang.shu@nbis.se)"
LABEL version "1.1"


#================================
# Install basics
#===============================
RUN apt-get  update -y
RUN apt-get install -y apt-utils  \
                       curl wget bc \
                       python python-dev python-pip \
                       build-essential  \
                       libxml2-dev libxslt1-dev libsqlite3-dev zlib1g-dev   \
                       cmake  \
                       gnuplot \
                       kalign \
                       imagemagick \
                       gengetopt \
                       xsltproc \
                       perlbrew \
                       locales-all \
                       default-jre

RUN apt-get install -y blast2 ncbi-blast+
#================================
# Install python package 
#===============================
RUN pip install --upgrade pip && pip install biopython==1.70

#================================
#  Add topcons2 source code
#===============================
WORKDIR /app
# add the source code to WORKDIR /app
ADD . ./topcons2

RUN mkdir -p /app/download /data/  /scratch/

# link data
RUN cd /app/topcons2 &&\
        ln -s /data/topcons2_database database

# building modhmm
RUN cd /app/topcons2/predictors/source/modhmm && \
    bash fresh_install.sh /app/topcons2/predictors/

#================================
# Install HMMER 
#===============================
RUN cd /app/download && \
        wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz -O  hmmer-3.1b2-linux-intel-x86_64.tar.gz && \
        tar -xvzf  hmmer-3.1b2-linux-intel-x86_64.tar.gz  && \
        cd hmmer-3.1b2-linux-intel-x86_64 && \
        ./configure && \
        make && \
        make install
#===============================
# install perl moudles
#===============================
RUN perlbrew install-cpanm
RUN /root/perl5/perlbrew/bin/cpanm Moose &&\
    /root/perl5/perlbrew/bin/cpanm IPC::Run &&\
    /root/perl5/perlbrew/bin/cpanm --force CJFIELDS/BioPerl-1.6.924.tar.gz

#================================
# Setting ENVs
#===============================
ENV PERL5LIB "/app/topcons2/database/pfam_seq/PfamScan"
ENV USER_DIRS "/app"

RUN rm -rf /app/download

CMD ["/bin/bash" ]
