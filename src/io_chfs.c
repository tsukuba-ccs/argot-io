#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <margo.h>
#include "io_ops.h"
#include <chfs.h>

void
ch_init()
{
	chfs_init(NULL);
}

void
ch_term()
{
	chfs_term();
}

int
ch_write_open(char *fn, void **r)
{
	int fd;

	fd = chfs_create(fn, O_WRONLY, 0644);
	if (fd < 0) {
		perror(fn);
		exit(EXIT_FAILURE);
	}
	*r = (void *)(intptr_t)fd;
	return (0);
}

int
ch_read_open(char *fn, void **r)
{
	int fd;

	fd = chfs_open(fn, O_RDONLY);
	if (fd < 0) {
		perror(fn);
		exit(EXIT_FAILURE);
	}
	*r = (void *)(intptr_t)fd;
	return (0);
}

int
ch_write(void *b, int s, void *r)
{
	int f = (intptr_t)r;
	int n;

	n = chfs_write(f, b, s);
	if (n < 0) {
		perror("write");
		exit(EXIT_FAILURE);
	}
	return (n);
}

int
ch_read(void *b, int s, void *r)
{
	int f = (intptr_t)r;
	int n;

	n = chfs_read(f, b, s);
	if (n < 0) {
		perror("read");
		exit(EXIT_FAILURE);
	}
	return (n);
}

int
ch_close(void *r)
{
	int f = (intptr_t)r;

	return (chfs_close(f));
}

int
ch_mkdir(char *fn, int mode)
{
	return (chfs_mkdir(fn, mode));
}

int
ch_unlink(char *fn)
{
	return (chfs_unlink(fn));
}

struct io_ops io_ops_chfs = {
	.init = ch_init,
	.term = ch_term,
	.write_open = ch_write_open,
	.read_open = ch_read_open,
	.write = ch_write,
	.read = ch_read,
	.close = ch_close,
	.mkdir = ch_mkdir,
	.unlink = ch_unlink
};

struct io_ops *
io_chfs(void)
{
	return (&io_ops_chfs);
}
