struct io_ops {
	void (*init)(void);
	void (*term)(void);
	int (*write_open)(char *, void **);
	int (*read_open)(char *, void **);
	int (*write)(void *, int, void *);
	int (*read)(void *, int, void *);
	int (*close)(void *);
	int (*mkdir)(char *, int);
};

extern struct io_ops *io_ops;

void init_io_ops(char *);
