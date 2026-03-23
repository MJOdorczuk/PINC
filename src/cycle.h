/**
 * @file		cycle.h
 * @brief		PINC cycle.
 * @author		Michal Jan Odorczuk <michaljo@uio.no>
 *
 * PINC cycle implementation.
 */

#ifndef CYCLE_H
#define CYCLE_H

#include "core.h"

typedef struct {
    void (*acc)();
    void (*distr)();
    void (*solverInterface)();
    void (*extractEmigrants)();
    void (*collide)();
    void (*solve)();
    void (*solverAlloc)();
    void (*solverFree)();
} Methods;

// TODO: Add the description
void cycle(dictionary *ini);

#endif