
# Statistical Assignment Package (SAP)

SAP makes statistical assingment of an unknown DNA sequence to taxa represented by
sequences in GenBank. Each assignment is associated with a posterior probility that serve
as a confidence in the assignment.

### Web service

SAP is available as a web service too at: [services.birc.au.dk/sap](https://services.birc.au.dk/sap).

### Installation

SAP is now distributed as a Docker image allowing it to run on both Mac, Windows, and Linux. 
To use SAP in this way, you need to download and install <a href="https://www.docker.com/get-started">Docker Desktop</a>.
Having done that, you run sap from the command line using this command on Mac/Linux:</p>

```
docker run --rm -v $PWD:/code/sap kaspermunch/sap:latest
```

and this one on Windows:

```
docker run --rm -v %CD%:/code/sap kaspermunch/sap:latest
```

The first time you do this it will pull the image from DockerHub and run it. On subsequent runs it will run the cached image.

See the [documentation](https://services.birc.au.dk/sap) for how to use sap on the command line.</p>

Where it is not possible to install Docker, as is often the case on shared computing facilities, the image can be pulled using [Singularity](https://sylabs.io/singularity) (you only need to do this once):

```
singularity pull sap docker://kaspermunch/sap:latest
```

Then just run the image as any other executable:

```
./sap [ARGUMENTS]
```

You can find more documentation and detail on how to install and run SAP on the [web service](https://services.birc.au.dk/sap) under Downloads.

### Issues and feature requests

You can report a bug or suggest a feature by [opening an issue](https://github.com/kaspermunch/sap/issues).
