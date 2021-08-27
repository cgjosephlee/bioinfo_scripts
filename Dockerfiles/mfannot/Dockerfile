FROM nbeck/mfannot:latest

ENV LC_ALL=C

WORKDIR /usr/local/bin
RUN sed -i 's/getlogin/getlogin || "docker"/g' HMMsearchCombiner HMMsearchWrapper grab-fasta grab-seq mf2sqn umac /mfannot/mfannot

VOLUME /project
WORKDIR /project

CMD ["mfannot", "-h"]
