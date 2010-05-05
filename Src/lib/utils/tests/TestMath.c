#include <assert.h>
#include <setjmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <utils/CuTest.h>
#include <utils/math.h>

#define ACC_EPS 0.00000001
/*-------------------------------------------------------------------------*
 *  Test
 *-------------------------------------------------------------------------*/

/* Tests the boys function */
void TestFm(CuTest* tc)
{
  CuAssertDblEquals(tc, Fm(0,0), 1, ACC_EPS);
  CuAssertDblEquals(tc, Fgamma(0,1), 0.7468241328124271, ACC_EPS);
  CuAssertDblEquals(tc, Fgamma(0.5,6.8), 0.073447516533246701, ACC_EPS);
}

/*-------------------------------------------------------------------------*
 * main
 *-------------------------------------------------------------------------*/

CuSuite* GetSuite(void)
{
  CuSuite* suite = CuSuiteNew();
  
  SUITE_ADD_TEST(suite, TestFm);
  
  return suite;
}
