FROM ubuntu:22.04 AS spython-base
LABEL Author PhilipSmith
LABEL Version v0.1.0
RUN TZ=Etc/UTC && \
ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
echo $TZ > /etc/timezone
RUN apt-get -y update
RUN apt-get install -y git-all
RUN apt-get clean
RUN git clone https://github.com/Phil9S/gel-sig-workflow.git
