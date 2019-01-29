# NaRLE
Nanopore Run Length Encoding: analysis pipeline

## Description

NaRLE is a pipeline for the analysis of ONT "Nanopore" data using Run Length Encoding

## Installation


### Dependencies ###

```
python2
pip
docker
virtualenv
```

### Quick Start ###

```
git clone git@github.com:UCSC-nanopore-cgl/NaRLE.git
cd NaRLE
virtualenv -p python2 venv
. venv/bin/activate
pip install -e .
toil-narle --help
```


### Docker Initialization ###

All tools used in the workflow must be run from docker images.  In an ideal world, these are up-to-date on dockerhub, but it is advisable that you build them manually.

The images are tagged with the current git version of the NaRLE repository (so the state of the docker file can be determined through git history).  As such, it is recommended that you only publish docker images which were created from an unmodified (or fully committed) repository.

There are custom entrypoints which store the executable's output in a log file as well as resource usage debugging information.  These require that the current working directory is mounted onto '/data'

#### Quick Start ####

```
cd docker
ls | xargs -n1 -I{} bash -c "cd {} &amp;&amp; make"
```
