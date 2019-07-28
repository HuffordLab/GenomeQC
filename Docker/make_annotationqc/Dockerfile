FROM ubuntu
RUN apt-get update --fix-missing

RUN apt-get install -y build-essential wget git autoconf cmake libjsoncpp-dev
RUN apt-get install -y python3
RUN apt-get install -y python3-pip
RUN python3 -m pip install statistics
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

ADD gff3_stats.py /opt
ADD pipeline.sh /opt
ENTRYPOINT ["bash", "/opt/pipeline.sh"]
