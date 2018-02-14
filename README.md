# ExTASY (Gromacs, LSDMap)

## Installation

* no need to install anything on bluewaters, uses a python installation by us

## Local python installation

In an conda python2.7 installation:
```
conda install -c conda-forge rabbitmq-server tmux pip git
pip install git+https://github.com/radical-cybertools/radical.utils.git@devel
pip install git+https://github.com/radical-cybertools/saga-python.git@feature-gpu
pip install git+https://github.com/radical-cybertools/radical.pilot.git@feature-gpu
pip install git+https://github.com/radical-cybertools/radical.entk.git@feature-gpu
pip install git+https://github.com/radical-cybertools/radical.analytics@devel
```
### Alternative Docker installation

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
This should work without typing any password:
```
gsissh username@bw.ncsa.illinois.edu
```


## Setup environment

Firstly, clone the current repository

```
git clone git@github.com:radical-collaboration/extasy-grlsd.git
cd extasy-grlsd
```

Next, you need to set a few environment variables, you can replace the RADICAL_PILOT_DBURL with your own mongoDB on mlab:
```
export RADICAL_ENTK_VERBOSE=info
export RP_ENABLE_OLD_DEFINES=True
export GLOBUS_LOCATION='/usr/' #assuming gsissh is at /usr/bin/gsissh
export RADICAL_ENTK_PROFILE=True
export RADICAL_PILOT_PROFILE=True
export SAGA_PTY_SSH_TIMEOUT=300
export RADICAL_PILOT_DBURL='mongodb://eh22:a3Qv*zs0@ds141209.mlab.com:41209/clementigroup'
```

Start the rabbitmq server

```
rabbitmq-server &
```

The behavior of the RabbitMQ server is visible under http://localhost:15672/#/ with login guest and password guest. If you need to restart the rabbitmq server type:
```
rabbitmqctl stop
rabbitmq-server &
```

## Executing the script

Setup the walltime, allocation and cores you require in resource_config.rcfg and all settings in the used settings_*.wcfg.
If you want to start a new adaptive sampling set start_iter to 0, if you want to extend the last adaptive sampling with more iterations set start_iter to the next iteration to run. 

Execution command for Ala2 "Alanine dipeptide", for longer simulations best to run inside tmux on an machine which can run undisturbed for long times:

```
python extasy_grlsd.py --Kconfig settings_ala2.wcfg
```

Execution command for Ala12 "Alanine12": 

```
python extasy_grlsd.py --Kconfig settings_ala12.wcfg
```


## Your own system
The MD simulation is in openmm, you have to inp_files:
* the gromacs structure
* a copy of run-openmm-ala12.py with any changes for you system (number of steps, forcefield,...)
* copy settings_ala12.wcfg and change the filenames



## Results
* output directory  will have for each iteration some output
* full output on bluewaters on at "remote_output_directory" as set in settings_*.wcfg


## Profiling
run ```python analytics_timing.py```, this gives information how much time which steps took.

## Notes 
The ```extasy_grlsd.py``` script contains the information about the application
execution workflow and the associated data movement. Please take a look at all
the comments to understand the various sections. 

