#include <assert.h>
#include <setjmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <utils/CuTest.h>

/*-------------------------------------------------------------------------*
 *  Test
 *-------------------------------------------------------------------------*/

void Test(CuTest* tc)
{
	CuAssert(tc, "test should pass", 1 == 0 + 1);
}

/*-------------------------------------------------------------------------*
 * main
 *-------------------------------------------------------------------------*/

CuSuite* GetSuite(void)
{
	CuSuite* suite = CuSuiteNew();

	SUITE_ADD_TEST(suite, Test);

	return suite;
}
