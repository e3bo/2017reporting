FROM r-base
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN install2.r --error packrat
RUN mkdir -p /home/docker/
COPY ./packrat /home/docker/packrat
COPY ./.Rprofile /home/docker/
RUN /usr/bin/Rscript -e "setwd(\"/home/docker\"); source(\".Rprofile\"); packrat::restore()"
