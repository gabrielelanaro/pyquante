#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include "../contracted-gto.h"
#include "../primitive-gto.h"

ContractedGTO **parseFile(char *fn, int *dim);
ContractedGTO **parseBasisSet(xmlDocPtr doc, xmlNodePtr cur, int *ndim);
ContractedGTO *parseCgto(xmlDocPtr doc, xmlNodePtr cur,
			 double x, double y, double z,
			 int l, int m, int n);

