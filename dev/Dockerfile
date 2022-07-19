FROM ubuntu

ARG USERNAME=foo

RUN apt-get update \
 && apt-get -y upgrade \
 && apt-get -y install git make libopenmpi-dev \
 && apt-get -y install sudo pkg-config

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

COPY --chown=$USERNAME:$USERNAME . /home/$USERNAME/argot-io

USER $USERNAME
RUN cd \
# && git clone https://github.com/tsukuba-ccs/argot-io.git \
 && cd argot-io \
 && autoreconf -i \
 && ./configure \
 && make \
 && sudo make install \
 && make clean

COPY --chown=$USERNAME:$USERNAME dev/hosts /home/$USERNAME/hosts
COPY --chown=$USERNAME:$USERNAME dev/test.sh /home/$USERNAME/test.sh
