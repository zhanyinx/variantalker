FROM broadinstitute/gatk:4.5.0.0

LABEL maintainer="variantalker"
LABEL description="GATK 4.5.0.0 with PostgreSQL support for variant annotation caching"
LABEL version="1.0.0"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Install system dependencies including locales
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-dev \
    postgresql-client \
    libpq-dev \
    gcc \
    locales \
    && rm -rf /var/lib/apt/lists/*

# Generate UTF-8 locale
RUN locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# Update environment variables to use the generated locale
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# Upgrade pip
RUN pip3 install --no-cache-dir --upgrade pip

# Install Python dependencies for PostgreSQL and data processing
RUN pip3 install --no-cache-dir \
    psycopg[binary]==3.0.18

# Verify installations
RUN python3 -c "import psycopg; print(f'psycopg version: {psycopg.__version__}')" && \
    python3 -c "import pandas; print(f'pandas version: {pandas.__version__}')" && \
    python3 -c "import numpy; print(f'numpy version: {numpy.__version__}')" && \
    gatk --version

# Set working directory
WORKDIR /workspace

# Default command
CMD ["/bin/bash"]
