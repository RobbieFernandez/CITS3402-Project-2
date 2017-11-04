#ifndef SITE_H
#define SITE_H

typedef struct site_struct SITE;

struct site_struct {
	int row;
	int col;
	int n;
	int e;
	int s;
	int w;
	int populated; 
};

enum route {
	ROWS,
	COLS,
	BOTH
};

typedef enum route ROUTE;

typedef struct results RESULTS;

struct results {
	int size;
	int percolates;
};

#endif
