#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "stack.h"
#include "site.h"

// This stack is NOT thread safe.
// Only access a stack from a single thread at a time.

// Push a site onto the stack.
void push(SITE_STACK* stack, SITE* site) {
	// Copy the current top stack item, then edit the previous one.
	// This way we can edit the stack pointer without returning anything.
	SITE_STACK* stack_cpy = (SITE_STACK*) malloc(sizeof(SITE_STACK));
	memcpy(stack_cpy, stack, sizeof(SITE_STACK));
	stack->next = stack_cpy;
	stack->site = site;
	stack->size = stack->next->size + 1;
}

// Pop a site off the stack and return it.
SITE* pop(SITE_STACK* stack) {
	// Copy "next" into the stack pointer.
	// Free "next"
	SITE* site = NULL;
	if (!is_empty(stack)) {
		site = stack->site;
		SITE_STACK* next = stack->next;
		memcpy(stack, next, sizeof(SITE_STACK));
		if (next != NULL)
			free(next);
	}
	return site;
}

// Initialise a stack's pointers to NULL and size to 0.
void init_stack(SITE_STACK* stack_ptr) {
	#pragma omp critical
	{	
		stack_ptr->next = NULL;
		stack_ptr->site = NULL;
		stack_ptr->size = 0;
	}
}

inline bool is_empty(SITE_STACK* stack_ptr) {
	return stack_ptr->size == 0;
}

// Free the memory allocated to the stack.
// Doesn't free the SITE pointers in the stack.
void free_stack(SITE_STACK* stack_ptr) {
	#pragma omp critical
	{
		while(stack_ptr != NULL) {
			SITE_STACK* next_ptr = stack_ptr->next;
			free(stack_ptr);
			stack_ptr = next_ptr;
		}
	}
}

// Print the pointers for each site in the stack.
// Doesn't alter the stack.
void print_stack(SITE_STACK* stack) {
	int n = stack->size;
	for(int i=0; i<n; i++) {
		printf("%p\n", (void*)stack->site);
		stack = stack->next;
	}
}
