FROM ubuntu

ARG USERNAME=foo

RUN apt-get update \
 && apt-get -y upgrade \
 && apt-get -y install git make libopenmpi-dev \
 && apt-get -y install sudo pkgconf \
 && apt-get -y install python3 automake cmake \
 && apt-get -y install fuse libfuse-dev pandoc

RUN \
  # sshd
  apt-get -y install --no-install-recommends \
    openssh-server \
  # sshd_config
  && printf '%s\n' \
    'PermitRootLogin yes' \
    'PasswordAuthentication yes' \
    'PermitEmptyPasswords yes' \
    'UsePAM no' \
    > /etc/ssh/sshd_config.d/auth.conf \
  # ssh_config
  && printf '%s\n' \
    'Host *' \
    '    StrictHostKeyChecking no' \
    > /etc/ssh/ssh_config.d/ignore-host-key.conf

RUN useradd -m -s /bin/bash $USERNAME \
 && echo "$USERNAME ALL=(ALL:ALL) NOPASSWD: ALL" > /etc/sudoers.d/$USERNAME \
 # delete passwd
 && passwd -d $USERNAME

USER $USERNAME
RUN cd \
 && git clone -c feature.manyFiles=true https://github.com/spack/spack.git \
 && . spack/share/spack/setup-env.sh \
 && spack external find automake autoconf libtool cmake m4 pkgconf \
 && spack install mochi-margo ^mercury~boostsys ^libfabric fabrics=rxm,sockets,tcp,udp \
 && echo ". spack/share/spack/setup-env.sh" >> .bashrc \
 && echo "spack load mochi-margo" >> .bashrc

#COPY --chown=$USERNAME:$USERNAME chfs /home/$USERNAME/chfs

RUN cd \
 && . spack/share/spack/setup-env.sh \
 && spack load mochi-margo \
 && git clone https://github.com/otatebe/chfs.git \
 && cd chfs \
 && autoreconf -i \
 && ./configure \
 && make \
 && sudo make install \
 && sudo ldconfig \
 && make clean

COPY --chown=$USERNAME:$USERNAME . /home/$USERNAME/argot-io

USER $USERNAME
RUN cd \
 && . spack/share/spack/setup-env.sh \
 && spack load mochi-margo \
# && git clone https://github.com/tsukuba-ccs/argot-io.git \
 && cd argot-io \
 && autoreconf -i \
 && ./configure \
 && make \
 && sudo make install \
 && make clean

COPY --chown=$USERNAME:$USERNAME dev/hosts /home/$USERNAME/hosts
COPY --chown=$USERNAME:$USERNAME dev/test.sh /home/$USERNAME/test.sh
