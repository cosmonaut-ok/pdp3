FROM debian:stretch-slim

WORKDIR /tmp

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

RUN echo 'deb-src http://deb.debian.org/debian stretch main' >> /etc/apt/sources.list
RUN apt-get update -y -qq

RUN apt-get -y -qq install apt-utils build-essential git imagemagick pandoc pandoc-citeproc libtool libtool-bin wget

# because mainstream openjdk does not want to be installed w/o some (pseudo)graphics
RUN apt-get -y -qq install default-jdk-headless || true
RUN apt-get -qq -y build-dep hdf5 || true

RUN wget -qO- "https://www.hdfgroup.org/package/source-bzip/?wpdmdl=12594&refresh=5bbc7778635b21539078008" | tar xjf -

WORKDIR /tmp/hdf5-1.10.3

RUN ./autogen.sh
RUN ./configure --enable-cxx --enable-build-mode=production --prefix=/usr
RUN make
RUN make install
