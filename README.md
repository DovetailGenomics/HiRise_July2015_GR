# HiRise

This release of HiRise is intended to allow readers of our paper
“Chromosome-scale shotgun assembly using an in vitro method for
long-range linkage” to reproduce the assemblies described therein.
For information on how to access the latest version of HiRise, please
visit http://www.dovetailgenomics.com/services.

## Setup

**NOTE:    
Where applicable, `/path/to/hirise` will refer to the path of the
`hirise` directory that contains this README.**

### Using Vagrant

A `VagrantFile` has been included for convenience. This will set up a
virtual machine with all the dependencies necessary to run HiRise. To
use this file, VirtualBox and Vagrant will need to be installed. To
set up a *Ubuntu 14.04* virtual machine and install all the necessary
dependencies: 
``` 
$ cd /path/to/hirise && vagrant up 
```

Once the virtual machine is ready:
```
$ vagrant ssh
vagrant@vagrant-ubuntu-trusty-64:~$ cd /vagrant && bash hirise_commands.bash
```

### Using Docker

A `Dockerfile` has also been included for proponents of Docker. To set
up a *Ubuntu 14.04* docker image and install all the necessary
dependencies:

```
$ docker build -t hirise .
```

Once the docker image has been built:
```
$ docker run --rm -ti -v `pwd`:/opt/hirise hirise /bin/bash
root@d540664db6a2:/opt/release-2.0.5# cd /opt/hirise && bash hirise_commands.bash
```

### Manual Setup

The basic dependencies will first need to be installed. On
Ubuntu/Debian, `apt-get` can be used, as follows:
```
$ sudo apt-get update
$ sudo apt-get -y upgrade
$ sudo apt-get install -y build-essential cmake libboost-all-dev liblapack-dev gfortran samtools python3-dev python3-pip
```

The required Python packages may be installed as such:
```
$ sudo pip3 install -r /path/to/hirise/requirements.txt
```

HiRise uses `merauder`, the gapclosing module of the Meraculous assembler for gap closing. This may be downloaded and installed as follows:
```
$ wget http://downloads.sourceforge.net/project/meraculous20/release-2.0.5.tgz
$ tar -zxvf release-2.0.5.tgz
$ mkdir $HOME/meraculous && cd release-2.0.5 && ./install.sh $HOME/meraculous
$ export PATH=$HOME/meraculous/bin:$PATH
```

The HiRise scripts will need to be added to the `PATH` as well (be sure to use the correct path to the `hirise` directory):
```
$ export PATH=/path/to/hirise/scripts:$PATH
```

HiRise can then be run as follows:
```
$ cd /path/to/hirise
$ wget https://s3-us-west-2.amazonaws.com/dovetail-public-data1/examples.tgz && tar -zxvf examples.tgz
$ bash hirise_commands.bash
```

## Preparing bam files for HiRise

HiRise expects sorted, indexed BAM files containing alignments of
Chicago reads to a starting assembly to be scaffolded. Read mapping
quality scores (MAPQ scores) should be computed based on the alignment
of each read, independent of its paired end read, and the MQ tag
should be set to provide the MAPQ score of the paired read. (The
samblaster --addMateTags command provides one way to add these tags.)


