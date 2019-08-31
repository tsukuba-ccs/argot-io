#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "io_ops.h"
#include <gfarm/gfarm.h>

void
gfarm_init()
{
	gfarm_initialize(NULL, NULL);
}

void
gfarm_term()
{
	gfarm_terminate();
}

/* XXX - internal function */
gfarm_error_t gfarm_realpath_by_gfarm2fs(const char *, char **);

int
gfarm_write_open(char *fn, void **r)
{
	GFS_File f;
	char *pi = NULL;
	gfarm_error_t e;

	e = gfarm_realpath_by_gfarm2fs(fn, &pi);
	if (e == GFARM_ERR_NO_ERROR)
		fn = pi;
	e = gfs_pio_create(fn, GFARM_FILE_WRONLY, 0644, &f);
	if (e != GFARM_ERR_NO_ERROR) {
		fprintf(stderr, "%s: %s\n", fn, gfarm_error_string(e));
		exit(EXIT_FAILURE);
	}
	free(pi);
	*r = f;
	return (0);
}

int
gfarm_read_open(char *fn, void **r)
{
	GFS_File f;
	char *pi = NULL;
	gfarm_error_t e;

	e = gfarm_realpath_by_gfarm2fs(fn, &pi);
	if (e == GFARM_ERR_NO_ERROR)
		fn = pi;
	e = gfs_pio_open(fn, GFARM_FILE_RDONLY, &f);
	if (e != GFARM_ERR_NO_ERROR) {
		fprintf(stderr, "%s: %s\n", fn, gfarm_error_string(e));
		exit(EXIT_FAILURE);
	}
	free(pi);
	*r = f;
	return (0);
}

int
gfarm_write(void *b, int s, void *r)
{
	GFS_File f = r;
	int n;
	gfarm_error_t e;

	e = gfs_pio_write(f, b, s, &n);
	if (e != GFARM_ERR_NO_ERROR) {
		fprintf(stderr, "write: %s\n", gfarm_error_string(e));
		exit(EXIT_FAILURE);
	}
	return (n);
}

int
gfarm_read(void *b, int s, void *r)
{
	GFS_File f = r;
	int n;
	gfarm_error_t e;

	e = gfs_pio_read(f, b, s, &n);
	if (e != GFARM_ERR_NO_ERROR) {
		fprintf(stderr, "read: %s\n", gfarm_error_string(e));
		exit(EXIT_FAILURE);
	}
	return (n);
}

int
gfarm_close(void *r)
{
	GFS_File f = r;
	gfarm_error_t e;

	e = gfs_pio_close(f);
	if (e == GFARM_ERR_NO_ERROR)
		return (0);
	errno = gfarm_error_to_errno(e);
	return (-1);
}

int
gfarm_mkdir(char *fn, int mode)
{
	gfarm_error_t e;

	e = gfs_mkdir(fn, mode);
	if (e == GFARM_ERR_NO_ERROR)
		return (0);
	errno = gfarm_error_to_errno(e);
	return (-1);
}

int
gfarm_unlink(char *fn)
{
	gfarm_error_t e;

	e = gfs_unlink(fn);
	if (e == GFARM_ERR_NO_ERROR)
		return (0);
	errno = gfarm_error_to_errno(e);
	return (-1);
}

struct io_ops io_ops_gfarm =
{
	.init = gfarm_init,
	.term = gfarm_term,
	.write_open = gfarm_write_open,
	.read_open = gfarm_read_open,
	.write = gfarm_write,
	.read = gfarm_read,
	.close = gfarm_close,
	.mkdir = gfarm_mkdir,
	.unlink = gfarm_unlink
};

struct io_ops *
io_gfarm(void)
{
	return(&io_ops_gfarm);
}
