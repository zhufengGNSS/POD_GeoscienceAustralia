FROM ubuntu:disco

RUN apt-get update
#==============================================================================
# install the build infrastructure to build an ami
#==============================================================================
RUN apt-get install -y packer 
RUN apt-get install -y unzip
RUN apt-get install -y jq 
RUN apt-get install -y wget  
RUN apt-get install -y curl  
RUN apt-get install -y python  
RUN wget --quiet https://releases.hashicorp.com/terraform/0.11.3/terraform_0.11.3_linux_amd64.zip \
  && unzip terraform_0.11.3_linux_amd64.zip \
  && mv terraform /usr/bin \
  && rm terraform_0.11.3_linux_amd64.zip
RUN curl "https://s3.amazonaws.com/aws-cli/awscli-bundle.zip" -o "awscli-bundle.zip" \
  && unzip awscli-bundle.zip \
  && ./awscli-bundle/install -i /usr/local/aws -b /usr/local/bin/aws
#==============================================================================
# install the dependencies to compile the pea
#==============================================================================
RUN apt-get install -y gcc 
RUN apt-get install -y make
RUN apt-get install -y libboost-dev
RUN apt-get install -y liblapack3
RUN ln -s /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/liblapack.so
RUN apt-get install -y libblas3
RUN ln -s /usr/lib/x86_64-linux-gnu/libblas.so.3 /usr/lib/x86_64-linux-gnu/libblas.so
#==============================================================================
# Make the directory structure, copy the code and compile
#==============================================================================
RUN mkdir -p /data/acs/pod
COPY ./ /data/acs/pod/
WORKDIR /data/acs/pod/
WORKDIR /data/acs/pod/src/
RUN make
#WORKDIR /data/acs/pde/
#RUN ./TST01-download.sh
#WORKDIR /data/acs/pde/src/test/
#RUN make
#RUN ./test_antenna
