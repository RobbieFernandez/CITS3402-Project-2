#include "site.h"

#ifndef SITE_STACK_H
#define SITE_STACK_H

typedef struct site_stack SITE_STACK;

struct site_stack {
	SITE* site;
	SITE_STACK* next;
	int size;
};

extern void  push(SITE_STACK* stack, SITE* site);
extern SITE* pop(SITE_STACK* stack);
extern void  init_stack(SITE_STACK* stack_ptr);
extern bool  is_empty(SITE_STACK* stack_ptr);
extern void free_stack(SITE_STACK* stack_ptr);

#endif
