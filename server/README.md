

# Docker server

Run the miniconda docker container and run bash

    docker run -it continuumio/miniconda3:4.7.12 /bin/bash

In the docker container: Create and export a sap environment in the miniconda container

    conda create --name sap -c conda-forge -c bioconda -c biobuilds -c anaconda python=2.7 blast clustalw flask celery redis redis-py gunicorn flask-mail gcc_linux-64 biopython=1.75 vine=1.3
    conda activate sap
    conda env export > environment.yml

Outside the docker container (read id off prompt in docker conatainer):

    docker cp b98e8292bb04:/environment.yml .

Close continuumio/miniconda3 container.

Delete old so files built on OSX (becuase the docker container links to these files and thinks they are up to date.)

    python2.7 setup.py clean

Build the container using develop.yml so that I can edit files on my own file system even though sap runs in the container

    docker-compose build

Tag images (replace VERSION for the version e.g. 1-9-9):

    docker image tag sap_web kaspermunch/sap_web:version-VERSION
    docker image tag sap_web kaspermunch/sap_web:latest

    docker image tag sap_worker kaspermunch/sap_worker:version-VERSION
    docker image tag sap_worker kaspermunch/sap_worker:latest

Login to DockerHub:

    docker login

Push docker images (replace VERSION for the version e.g. 1-9-9):

    docker push kaspermunch/sap_web:version-VERSION
    docker push kaspermunch/sap_web:latest

    docker push kaspermunch/sap_worker:version-VERSION
    docker push kaspermunch/sap_worker:latest

Testing the server:

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

# Docker for distribution

**NB:** The following needs to be run in the root dir of the project.

Run the miniconda docker container and run bash

    docker run -it continuumio/miniconda3:4.7.12 /bin/bash

In the docker container: Create and export a sap environment in the miniconda container

    conda create --name sap -c conda-forge -c bioconda -c biobuilds -c anaconda python=2.7 blast clustalw flask celery redis redis-py gunicorn flask-mail gcc_linux-64 biopython=1.75 vine=1.3
    conda activate sap
    conda env export > environment.yml

Outside the docker container (read id off prompt in docker conatainer):

    docker cp b98e8292bb04:/environment.yml .

Close continuumio/miniconda3 container.

Delete old so files built on OSX (becuase the docker container links to these files and thinks they are up to date.)

    python2.7 setup.py clean

Build Docker image (replace VERSION for the version e.g. 1-9-9):

    docker build -t kaspermunch/sap:latest -t kaspermunch/sap:version-VERSION  .

Login to DockerHub:

    docker login

Push docker image (replace VERSION for the version e.g. 1-9-9):

    docker push kaspermunch/sap:version-VERSION
    docker push kaspermunch/sap:latest

Run bash in image interactively:

    docker run -it --entrypoint /bin/bash kaspermunch/sap

Run docker image (with Docker):

    docker run --rm -v $PWD:/code/sap kaspermunch/sap:latest [ARGUMENTS]

On a machine without access to Docker, the image can be pulled using Singularity (you only need to do this once):

    singularity pull sap docker://kaspermunch/sap:latest

Then just run the image as any other executable:

    ./sap [ARGUMENTS]


# Developing in container

In VScode, run "Open workspace in container"

Open terminal and run

    python setup.py develop

Then

    activate sap

Then you can run sap on terminal to test changes to code.