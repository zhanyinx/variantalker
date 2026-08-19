FROM griffithlab/civicpy

RUN apt-get update -y
RUN apt-get install -y tabix
RUN pip install bgzip
RUN civicpy update
RUN cp -r /root/.civicpy /home/
