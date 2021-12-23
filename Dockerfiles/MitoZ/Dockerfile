FROM guanliangmeng/mitoz:2.4-alpha

RUN mv /root/.etetoolkit /

RUN apt-get update
RUN apt-get install -y patch

COPY mitoZ.patch /app/release_MitoZ_v2.4-alpha
RUN cd /app/release_MitoZ_v2.4-alpha && \
    patch -p1 < mitoZ.patch
RUN ln -s /app/release_MitoZ_v2.4-alpha /app/MitoZ
