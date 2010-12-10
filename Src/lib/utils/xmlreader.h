/*
 * This source code is part of the PyQuante Quantum Chemistry suite.
 *  
 * Written by Gabriele Lanaro, 2009-2010
 * Copyright (c) 2009-2010, Gabriele Lanaro
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 */
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include "../contracted-gto.h"
#include "../primitive-gto.h"

ContractedGTO **parseFile(char *fn, int *dim);
ContractedGTO **parseBasisSet(xmlDocPtr doc, xmlNodePtr cur, int *ndim);
ContractedGTO *parseCgto(xmlDocPtr doc, xmlNodePtr cur,
			 double x, double y, double z,
			 int l, int m, int n);

