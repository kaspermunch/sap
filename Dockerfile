FROM python:2.7.9

RUN apt-get update && apt-get install -y clustalw

ADD MANIFEST.in ez_setup.py setup.cfg setup.py README /code/sap/
ADD SAP/ /code/sap/SAP/
ADD icons/ /code/sap/icons/
ADD ext/ /code/sap/ext/
RUN cd code/sap/ && python -u setup.py install

ADD server /code/server

RUN ln -s /usr/bin/clustalw /usr/bin/clustalw2

RUN cd /code/server && pip install -r requirements.txt

RUN mkdir /opt/blast \
      && curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz \
      | tar -zxC /opt/blast --strip-components=1

ENV PATH /opt/blast/bin:$PATH
