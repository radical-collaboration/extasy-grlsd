# ExTASY (Gromacs, LSDMap)
 

## Requirements

* Gromacs and lsdmap to be installed on Blue Waters (currently installed
and made public from Vivek's account)

## EnTK, Radical Pilot Installation

```
pip install radical.pilot
git clone https://github.com/radical-collaboration/extasy-grlsd.git
git clone https://github.com/radical-cybertools/radical.entk.git
cd radical.entk
git checkout feature/gpu
pip install .
apt-get install rabbitmq-server

```
* Note: For the current version, you will have to install RabbitMQ. 
[This link](http://radicalentk-06.readthedocs.io/en/arch-v0.6/install.html) provides two methods in which
you can install RabbitMQ.

You will need to use docker to run rabbitMQ for this project.
For setting up RabbitMQ with Docker use [This Link](http://radicalentk-06.readthedocs.io/en/arch-v0.6/install.html)

Example docker commands to run rabbitmq command 

```
sudo docker run -d --name rabbit-1  -P rabbitmq:3
```
* Note: that if you will run rabbitmq you will need to make a different name for each rabbitmq-server 

## Setting up access to HPCs

Currently, all packages and permissions are setup for Blue Waters.

[Blue Waters](https://bluewaters.ncsa.illinois.edu/user-guide)
requires GSISSH access. Instructions to setup gsissh access for Ubuntu can be 
found [here](https://github.com/vivek-bala/docs/blob/master/misc/gsissh_setup_stampede_ubuntu_xenial.sh/).
Please note that this has been tested only for xenial and trusty (for trusty, 
simple replace 'xenial' with 'trusty' in all commands). Even then, there might 
be some additional steps to setup gsissh correctly for your system. Happy to 
help!


## Setup environment

Firstly, clone the current repository

```
git clone git@github.com:radical-collaboration/extasy-grlsd.git
cd extasy-grlsd
```

Next, you need to set a few environment variables:
```
export RADICAL_ENTK_VERBOSE=info
export RADICAL_PILOT_DBURL="mongo db url"
export RP_ENABLE_OLD_DEFINES=True
export GLOBUS_LOCATION=/usr
```
For profiling:
```
export RADICAL_ENTK_PROFILE=True
export RADICAL_PILOT_PROFILE=True
```

## Executing the script

Setup the walltime, allocation and cores you require in resource_config.rcfg.


The ```extasy_grlsd.py``` script contains the information about the application 
execution workflow and the associated data movement. Please take a look at all 
the comments to understand the various sections.

Execution command: 

```
rabbitmq-server &
python extasy_grlsd.py --RPconfig resource_config.rcfg --Kconfig gromacslsdmap.wcfg
```

The behavior of the RabbitMQ server is visible under http://localhost:15672/#/ with login guest and password guest

## Note

* Hopefully not, but there might be lingering permission issues which will get 
detected once other users start running the code.
