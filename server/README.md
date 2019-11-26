

# Docker

Run the miniconda docker container and run bash

    docker run -it continuumio/miniconda3:4.7.12 /bin/bash

In the docker container: Create and export a sap environment in the miniconda container

    conda create --name sap -c conda-forge -c bioconda -c biobuilds -c anaconda python=2.7 blast clustalw flask celery redis redis-py gunicorn flask-mail gcc_linux-64
    conda activate sap
    pip install biopython==1.75
    conda env export > environment.yml

Outside the docker container (read id off prompt in docker conatainer):

    docker cp b98e8292bb04:/environment.yml .

Close continuumio/miniconda3 container.

Delete old so files built on OSX (becuase the docker container links to these files and thinks they are up to date.)

    python setup.py clean

Build the container using develop.yml so that I can edit files on my own file system even though sap runs in the container

    docker-compose build

Take container down

    docker-compose down

Run bash in the container

    docker-compose exec web /bin/bash
    cd ..
    conda activate sap
    sap --project mytest --email kaspermunch@birc.au.dk tests/query.fasta
    sap --compile '(COI[Gene Name]) AND barcode[Keyword] AND Aves[Orgn]' --database local_database.fasta --email kaspermunch@birc.au.dk

Start the web service on port 7000

    docker-compose up


