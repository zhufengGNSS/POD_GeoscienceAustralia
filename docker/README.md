# Docker Instructions

## Follow these instructions to install docker:
https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04

## Build with just packages
To build the image with just the required packages run the docker build command

```
docker build -t pod:0.1 .
```

## Build with packages and compiled pod code

```
cd ..
docker build -t pod:0.1 -f docker/dockerfile-fullbuild .
```


## Run 

To run and enter the container

```
docker run it  pod:0.1 /bin/bash
```

To run just a command:
```
docker run -d --name day1 pod:0.1 /root/pod/bin/pod -h
```

To view the logs
```
docker logs day1
```

