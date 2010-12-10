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
#include <stdio.h>
#include <stdlib.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include "../contracted-gto.h"
#include "../primitive-gto.h"
#include "xmlreader.h"


/* Reading from the file */
ContractedGTO ** parseFile(char *fn, int *dim) {

    xmlDocPtr doc;
    xmlNodePtr cur;

    doc = xmlReadFile(fn, NULL, 0);
    cur = xmlDocGetRootElement(doc);

    if (cur == NULL) {
        fprintf(stderr,"empty document\n");
        xmlFreeDoc(doc);
        return NULL;
    }

    return parseBasisSet(doc,cur, dim);
}

/* Parse cgtos */
/*
  Returns an array of ContractedGTOs
*/
ContractedGTO **parseBasisSet(xmlDocPtr doc, xmlNodePtr cur, int *ndim) {

    int l,m,n;
    float x,y,z;
    char *origin,*powers;
    ContractedGTO ** cgtos = NULL;
    int dim=0;

    cur = cur->children;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (const xmlChar *) "cgbf")) {
            origin = (char *) xmlGetProp(cur, (xmlChar *)"origin");
            powers = (char *) xmlGetProp(cur, (xmlChar *)"powers");

            sscanf(origin, "(%f,%f,%f)", &x,&y,&z);
            sscanf(powers, "(%d,%d,%d)", &l,&m,&n);

            cgtos = realloc(cgtos, sizeof(ContractedGTO *)*(dim+1));

            cgtos[dim] = parseCgto(doc,cur,
                                   x,y,z,// Additional info needed
                                   l,m,n);

            dim++;
        }
        cur = cur->next;
    }
    *ndim = dim;
    return cgtos;
}


/* Parse cgto */
ContractedGTO *parseCgto(xmlDocPtr doc, xmlNodePtr cur, double x, double y, double z,
                         int l, int m, int n) {
    PrimitiveGTO **pgtos = NULL;
    double *coefs = NULL;
    int dim = 0;

    char *alpha_str,*coef_str;
    float alpha, coef;
    ContractedGTO *retcgto;

    cur = cur->children;

    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (const xmlChar *) "prim")) {
            /* Fetching */
            alpha_str =(char *) xmlGetProp(cur, (xmlChar *)"exp");
            coef_str = (char *) xmlGetProp(cur, (xmlChar *)"coeff");
            sscanf(alpha_str,"%f",&alpha);
            sscanf(coef_str,"%f",&coef);

            /* Creating objects */
            pgtos = realloc(pgtos, sizeof(PrimitiveGTO)*(dim+1));

            pgtos[dim] = primitive_gto_new(alpha,x,y,z,l,m,n,coef) ;
            dim++;
        }
        cur = cur->next;
    }
    retcgto =  contracted_gto_new(x,y,z,l,m,n,0);// We don't need atid in these tests, maybe
    contracted_gto_from_primitives(retcgto,pgtos, dim);

    return retcgto;
}

