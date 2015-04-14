FROM python:2.7.9

RUN apt-get update && apt-get install -y ncbi-blast+ clustalw

ADD MANIFEST.in ez_setup.py setup.cfg setup.py README /code/sap/
ADD SAP/ /code/sap/SAP/
ADD icons/ /code/sap/icons/
ADD ext/ /code/sap/ext/
RUN cd code/sap/ && python -u setup.py install

ADD server /code/server
RUN cd /code/server && pip install -r requirements.txt
