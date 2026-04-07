# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install R and all system dependencies first
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    gnupg \
    apt-transport-https \
    ca-certificates \
    fastqc \
    python3 \
    python3-pip \
    curl \
    wget \
    git \
    build-essential \
    unzip \
    tar \
    openjdk-17-jre-headless \
    cutadapt \
    libbz2-dev \
    zlib1g-dev \
    libboost-all-dev \
    perl \
    libperlio-gzip-perl \
    cpanminus \
    # Add these for TransDecoder Perl dependencies
    libdb-dev \
    liburi-escape-xs-perl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Add R repository and install R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install required Python packages
RUN pip3 install --no-cache-dir \
    pandas requests multiqc biopython \
    matplotlib plotly goatools numpy scipy scikit-learn statsmodels seaborn

# Create data directories for databases
RUN mkdir -p /data/go /data/uniprot /data/contaminants

# Download GO database
RUN wget -O /data/go/go-basic.obo https://purl.obolibrary.org/obo/go/go-basic.obo && \
    chmod 644 /data/go/go-basic.obo
    
# Set environment variable for GO database
ENV GO_OBO_PATH="/data/go/go-basic.obo"

# Copy pre-built DIAMOND databases for out-of-the-box functionality
COPY data/uniprot/uniprot.dmnd /data/uniprot/uniprot.dmnd
COPY data/contaminants/contaminants_uniprot.dmnd /data/contaminants/contaminants_uniprot.dmnd

# Set database environment variables
ENV UNIPROT_DB="/data/uniprot/uniprot.dmnd"
ENV CONTAMINANTS_DB="/data/contaminants/contaminants_uniprot.dmnd"

# Set file permissions for databases
RUN chmod 644 /data/uniprot/uniprot.dmnd && \
    chmod 644 /data/contaminants/contaminants_uniprot.dmnd

# Install Diamond
RUN wget https://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz && \
    tar xzf diamond-linux64.tar.gz && \
    mv diamond /usr/local/bin/ && \
    chmod +x /usr/local/bin/diamond && \
    rm -rf diamond-linux64.tar.gz

# Install required Perl modules for TransDecoder
RUN cpanm --notest URI::Escape

# Install TransDecoder
RUN wget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.5.0.tar.gz && \
    tar xzf TransDecoder-v5.5.0.tar.gz && \
    mv TransDecoder-TransDecoder-v5.5.0 /opt/TransDecoder && \
    cd /opt/TransDecoder && \
    make && \
    ln -s /opt/TransDecoder/TransDecoder.LongOrfs /usr/local/bin/ && \
    ln -s /opt/TransDecoder/TransDecoder.Predict /usr/local/bin/ && \
    rm -rf TransDecoder-v5.5.0.tar.gz

# Set TransDecoder environment variables
ENV PERL5LIB="/opt/TransDecoder/PerlLib:$PERL5LIB"
ENV PATH="/opt/TransDecoder:/opt/TransDecoder/util:$PATH"

# Install Salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz && \
    tar xzf salmon-1.10.0_linux_x86_64.tar.gz && \
    mkdir -p /opt/salmon && \
    cp -r salmon-latest_linux_x86_64/* /opt/salmon/ && \
    ln -s /opt/salmon/bin/salmon /usr/local/bin/salmon && \
    chmod +x /opt/salmon/bin/salmon && \
    rm -rf salmon-1.10.0_linux_x86_64.tar.gz salmon-latest_linux_x86_64

# Configure Salmon wrapper
RUN mv /usr/local/bin/salmon /usr/local/bin/salmon-original && \
    echo '#!/bin/bash' > /usr/local/bin/salmon && \
    echo 'LD_LIBRARY_PATH="/opt/salmon/lib:$LD_LIBRARY_PATH" /usr/local/bin/salmon-original "$@"' >> /usr/local/bin/salmon && \
    chmod +x /usr/local/bin/salmon

# Install required R packages
RUN Rscript -e 'install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "RColorBrewer", "reshape2", "optparse"), repos="https://cloud.r-project.org/", dependencies=TRUE)'

# Install BiocManager and topGO
RUN Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org/")' && \
    Rscript -e 'BiocManager::install("topGO", dependencies=TRUE)'

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# Create app directory for scripts (not /workspace - that's for user data)
RUN mkdir -p /app/bin/annotation \
             /app/bin/contamination \
             /app/bin/enrichment \
             /app/modules/explot/cli \
             /app/modules/explot/categories \
             /app/modules/landscape

# Copy scripts and modules to app directory
COPY bin/ /app/bin/
COPY sceptr/ /app/sceptr/
COPY modules/ /app/modules/

# Make scripts executable
RUN find /app -type f -name "*.py" -exec chmod +x {} \; && \
    find /app -type f -name "*.R" -exec chmod +x {} \;

# Add app paths to environment (not workspace paths)
ENV PATH="/app/bin:/app/bin/annotation:/app/bin/contamination:/app/bin/enrichment:/app/modules/explot/cli:$PATH"

# Simple entrypoint
RUN echo '#!/bin/bash' > /usr/local/bin/entrypoint.sh && \
    echo 'exec "$@"' >> /usr/local/bin/entrypoint.sh && \
    chmod +x /usr/local/bin/entrypoint.sh

# Verification
RUN echo "Verifying installations:" && \
    which salmon && salmon --version && \
    which TransDecoder.LongOrfs && \
    which diamond && diamond version && \
    which fastqc && fastqc --version && \
    echo "Database verification:" && \
    ls -lh /data/uniprot/uniprot.dmnd && \
    ls -lh /data/contaminants/contaminants_uniprot.dmnd && \
    ls -lh /data/go/go-basic.obo && \
    echo "All tools and databases verified"

# Working directory for user data
WORKDIR /workspace

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]

LABEL version="1.0.0" \
      description="SCEPTR - Statistical Characterisation of Expression Profiles in Transcriptomes" \
      maintainer="James McCabe" \
      databases.included="true" \
      databases.uniprot="Swiss-Prot (curated)" \
      databases.contaminants="Custom curated" \
      databases.go="go-basic.obo"
