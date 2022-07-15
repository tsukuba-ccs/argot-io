# Docker containers to test ARGOT-IO

## How to use

- build image

      % docker compose build

- execute containers

      % docker compose up -d

  This executes four contaiers using the docker image.  You can login to all containers by ssh.  For details, see docker-compose.yml.

- login to a container and execute ARGOT-IO

      % ssh 172.30.3.2
      % sh test.sh

  output files are created in work directory.  Performance result is in work/io_bench-out.

- shutdown containers

      % docker compose down
