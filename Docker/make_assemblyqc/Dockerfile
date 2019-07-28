FROM ubuntu
RUN apt-get update --fix-missing

RUN apt-get install -y build-essential wget git autoconf cmake libjsoncpp-dev
#RUN apt-get install -y python3
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
/bin/bash Miniconda3-latest-Linux-x86_64.sh -f -b -p /opt/conda && \
rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/conda/bin:$PATH
RUN apt-get install -y libpango1.0-dev
RUN apt-get install -y libcairo2-dev

RUN conda create -n LTR_retriever
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install -c conda-forge perl perl-text-soundex
RUN conda install -c bioconda cd-hit
RUN conda install -c bioconda/label/cf201901 repeatmasker

RUN git clone https://github.com/oushujun/LTR_retriever.git
RUN echo "source activate LTR_retriever" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH


RUN python3 -m pip install biopython
RUN python3 -m pip install statistics
RUN python3 -m pip install pandas
RUN python3 -m pip install plotly
RUN python3 -m pip install docopt
RUN python3 -m pip install numpy
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y libgsl-dev
RUN apt-get install -y libboost-graph-dev
RUN apt-get install -y libmysql++-dev
#RUN apt-get install -y libbamtools-dev
#RUN apt-get install -y libsuitesparse-dev liblpsolve55-dev
RUN apt-get install -y libbz2-dev liblzma-dev
RUN apt-get install -y libncurses5-dev
RUN apt-get install -y libssl-dev libcurl4-openssl-dev
RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.3.1.tar.gz
RUN mkdir /opt/augustus
RUN mv augustus-3.3.1.tar.gz /opt/augustus
WORKDIR /opt/augustus
RUN tar zxvf augustus-3.3.1.tar.gz

RUN git clone git://github.com/pezmaster31/bamtools.git /opt/bamtools
RUN mkdir /opt/bamtools/build
WORKDIR /opt/bamtools/build
RUN cmake -DCMAKE_INSTALL_PREFIX=/opt/bamtools/ ..
RUN make
RUN make install
WORKDIR /opt/augustus/augustus-3.3.1/auxprogs/bam2hints
RUN sed -i 's|INCLUDES = /usr/include/bamtools|INCLUDES = /opt/bamtools/include/bamtools|' Makefile
RUN sed -i 's|LIBS = -lbamtools -lz|LIBS = /opt/bamtools/lib/libbamtools.a -lz|' Makefile
WORKDIR /opt/augustus/augustus-3.3.1/auxprogs/filterBam/src
RUN sed -i 's|BAMTOOLS = /usr/include/bamtools|BAMTOOLS = /opt/bamtools/include/bamtools|' Makefile
RUN sed -i 's|INCLUDES = -I$(BAMTOOLS) -Iheaders -I./bamtools|INCLUDES = -I/opt/bamtools/include/bamtools -Iheaders -I./bamtools|' Makefile
RUN sed -i 's|LIBS = -lbamtools -lz|LIBS = /opt/bamtools/lib/libbamtools.a -lz|' Makefile


RUN git clone https://github.com/samtools/htslib.git /opt/htslib
WORKDIR /opt/htslib
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install

RUN git clone https://github.com/samtools/samtools.git /opt/samtools
WORKDIR /opt/samtools
RUN autoheader
RUN autoconf -Wno-syntax
RUN ./configure
RUN make
RUN make install

RUN git clone https://github.com/samtools/bcftools.git /opt/bcftools
WORKDIR /opt/bcftools
RUN autoheader
RUN autoconf -Wno-syntax
RUN ./configure
RUN make
RUN make install

RUN git clone https://github.com/samtools/tabix.git /opt/tabix
WORKDIR /opt/tabix
RUN make
WORKDIR /opt/augustus/augustus-3.3.1/auxprogs/bam2wig
RUN sed -i 's|SAMTOOLS=$(TOOLDIR)/samtools/|SAMTOOLS=/opt/samtools/|' Makefile
RUN sed -i 's|HTSLIB=$(TOOLDIR)/htslib/|HTSLIB=/opt/htslib/|' Makefile
RUN sed -i 's|BCFTOOLS=$(TOOLDIR)/bcftools/|BCFTOOLS=/opt/bcftools/|' Makefile
RUN sed -i 's|TABIX=$(TOOLDIR)/tabix/|TABIX=/opt/tabix/|' Makefile
WORKDIR /opt/augustus/augustus-3.3.1/auxprogs
RUN sed -i 's|cd bam2wig; make;|#cd bam2wig; make;|' Makefile
RUN sed -i 's|cd bam2wig; make clean;|#cd bam2wig; make clean;|' Makefile
WORKDIR /opt/augustus/augustus-3.3.1/
RUN make
RUN make install

RUN apt-get install -y ncbi-blast+

RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz
RUN tar zxf hmmer.tar.gz
WORKDIR hmmer-3.2.1
RUN ./configure --prefix /opt/hmmer
RUN make
RUN make install

ENV PATH="/opt/augustus/augustus-3.3.1/bin:${PATH}"
ENV AUGUSTUS_CONFIG_PATH="/opt/augustus/augustus-3.3.1/config"
RUN git clone https://gitlab.com/ezlab/busco.git /opt/busco
WORKDIR /opt/busco
RUN python setup.py install
RUN mv /opt/busco/config/config.ini.default /opt/busco/config/config.ini
WORKDIR /opt/busco/config
RUN sed -i 's|path = /home/osboxes/BUSCOVM/hmmer/hmmer-3.1b2-linux-intel-ia32/binaries/|path = /opt/hmmer/bin/|' config.ini
RUN sed -i 's|path = /home/osboxes/BUSCOVM/augustus/augustus-3.2.2/bin/|path = /opt/augustus/augustus-3.3.1/bin/|' config.ini
RUN sed -i 's|path = /home/osboxes/BUSCOVM/augustus/augustus-3.2.2/scripts/|path = /opt/augustus/augustus-3.3.1/scripts/|' config.ini

RUN ["chmod", "+x", "/opt/busco/scripts/run_BUSCO.py"]

RUN mkdir /opt/contamination_check
ADD contamination.py /opt/contamination_check
ADD BtLog.py /opt/contamination_check
ADD BtIO.py /opt/contamination_check
ADD blobtools.py /opt/contamination_check
ADD Univec_modified.fasta /opt/contamination_check
ADD Univec_modified.fasta.nhr /opt/contamination_check
ADD Univec_modified.fasta.nin /opt/contamination_check
ADD Univec_modified.fasta.nsq /opt/contamination_check
ADD taxon_n /opt/contamination_check
ADD NG.py /opt
ADD assembly_stats.py /opt
ADD pipeline.sh /opt
RUN wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/ltr-finder/LTR_FINDER.x86_64-1.0.5.tar.gz
RUN mv LTR_FINDER.x86_64-1.0.5.tar.gz /opt
WORKDIR /opt
RUN tar zxvf LTR_FINDER.x86_64-1.0.5.tar.gz
COPY run_ltrfinder_in_pieces.pl /opt/LTR_FINDER.x86_64-1.0.5
COPY cut.pl /opt/LTR_FINDER.x86_64-1.0.5
RUN ["chmod", "+x", "/opt/LTR_FINDER.x86_64-1.0.5/run_ltrfinder_in_pieces.pl"]
RUN wget http://genometools.org/pub/binary_distributions/gt-1.5.10-Linux_x86_64-64bit-complete.tar.gz
RUN tar zxvf gt-1.5.10-Linux_x86_64-64bit-complete.tar.gz

ENTRYPOINT ["bash", "/opt/pipeline.sh"]
