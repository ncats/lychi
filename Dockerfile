FROM openjdk:8-jdk as buildContainer

RUN apt-get update
RUN apt-get install -y maven

WORKDIR /opt/app
COPY . /opt/app

RUN chmod +x ./make.sh
RUN bash ./make.sh

CMD bash