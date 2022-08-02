#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include "io_ops.h"
#include "config.h"

void
posix_init(void)
{
}

void
posix_term(void)
{
}

int
posix_write_open(char *fn, void **r)
{
	FILE *f;

	f = fopen(fn, "w");
	if (f == NULL)
		perror(fn), exit(EXIT_FAILURE);
	*r = f;
	return (0);
}

int
posix_read_open(char *fn, void **r)
{
	FILE *f;

	f = fopen(fn, "r");
	if (f == NULL)
		perror(fn), exit(EXIT_FAILURE);
	*r = f;
	return (0);
}

int
posix_write(void *b, int s, void *r)
{
	FILE *f = r;

	return (fwrite(b, s, 1, f));
}

int
posix_read(void *b, int s, void *r)
{
	FILE *f = r;

	return (fread(b, s, 1, f));
}

int
posix_close(void *r)
{
	FILE *f = r;
	int rv;

	rv = fclose(f);
	if (rv == EOF)
		return (-1);
	return (0);
}

int
posix_mkdir(char *fn, int mode)
{
	return (mkdir(fn, mode));
}

int
posix_unlink(char *fn)
{
	return (unlink(fn));
}

struct io_ops io_ops_posix = {
	.init = posix_init,
	.term = posix_term,
	.write_open = posix_write_open,
	.read_open = posix_read_open,
	.write = posix_write,
	.read = posix_read,
	.close = posix_close,
	.mkdir = posix_mkdir,
	.unlink = posix_unlink
};

struct io_ops *
io_posix(void)
{
	return (&io_ops_posix);
}

struct io_ops *io_ops;
#ifdef HAVE_LIBGFARM
struct io_ops *io_gfarm(void);
#endif
#ifdef HAVE_CHFS
struct io_ops *io_chfs(void);
#endif

void
init_io_ops(char *api)
{
	if (strcmp(api, "posix") == 0)
		io_ops = io_posix();
#ifdef HAVE_LIBGFARM
	else if (strcmp(api, "gfarm") == 0)
		io_ops = io_gfarm();
#endif
#ifdef HAVE_CHFS
	else if (strcmp(api, "chfs") == 0)
		io_ops = io_chfs();
#endif
	if (io_ops == NULL)
		fprintf(stderr, "%s: not supported\n", api), exit(EXIT_FAILURE);
}
