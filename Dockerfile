FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
    curl \
    wget \
    ca-certificates \
    vim \
    sudo \
    git \
    bzip2 \
 && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.3.31-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
ENV PATH /opt/conda/bin:$PATH

RUN conda install -c anaconda numpy && conda clean -ya
RUN conda install -c bioconda/label/cf201901 bwa && conda clean -ya
RUN conda install -c bioconda samtools==1.7 && conda clean -ya
RUN conda install -c bioconda bedtools==2.17.0 && conda clean -ya
RUN conda install -c bioconda seqtk && conda clean -ya
RUN conda install -c conda-forge pathlib && conda clean -ya

RUN wget --quiet https://drive5.com/cgi-bin/upload3.py?license=2019110501245926689 -O /usr/bin/usearch \
 && chmod +x /usr/bin/usearch
#COPY ./usearch /usr/bin

# Create a working directory
RUN mkdir /STRsearch
WORKDIR /STRsearch

COPY ./STRsearch /STRsearch
RUN chmod +x /STRsearch/pipeline.py
RUN sed -i 's$YOUR PATH/bwa$/opt/conda/bin/bwa$g' /STRsearch/conf.py
RUN sed -i 's$YOUR PATH/samtools$/opt/conda/bin/samtools$g' /STRsearch/conf.py
RUN sed -i 's$YOUR PATH/bamToFastq$/opt/conda/bin/bamToFastq$g' /STRsearch/conf.py
RUN sed -i 's$YOUR PATH/seqtk$/opt/conda/bin/seqtk$g' /STRsearch/conf.py
RUN sed -i 's$YOUR PATH/usearch$/usr/bin/usearch$g' /STRsearch/conf.py

#RUN git clone https://github.com/AnJingwd/STRsearch.git

ENTRYPOINT ["python3","/STRsearch/pipeline.py"]
CMD ["-h"]
