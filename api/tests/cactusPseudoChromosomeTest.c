/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusReferenceTestShared.h"

static void testTeardown() {
	if(!nestedTest) {
		cactusReferenceTestSharedTeardown();
	}
}

static void testSetup() {
	if(!nestedTest) {
		cactusReferenceTestSharedSetup();
	}
}

void testPseudoChromosome_construct(CuTest* testCase) {
	testSetup();
	assert(testCase != NULL);
	//tested by setup code.
	testTeardown();
}

void testPseudoChromosome_getName(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, pseudoChromosome_getName(pseudoChromosome1) != NULL_NAME);
	CuAssertTrue(testCase, pseudoChromosome_getName(pseudoChromosome2) != NULL_NAME);
	CuAssertTrue(testCase, pseudoChromosome_getName(pseudoChromosome1) != pseudoChromosome_getName(pseudoChromosome2));
	testTeardown();
}

void testPseudoChromosome_get5End(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, pseudoChromosome_get5End(pseudoChromosome1) == end1);
	CuAssertTrue(testCase, pseudoChromosome_get5End(pseudoChromosome2) == end7);
	testTeardown();
}

void testPseudoChromosome_get3End(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, pseudoChromosome_get3End(pseudoChromosome1) == end6);
	CuAssertTrue(testCase, pseudoChromosome_get3End(pseudoChromosome2) == end8);
	testTeardown();
}

void testPseudoChromosome_getReference(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, pseudoChromosome_getReference(pseudoChromosome1) == reference);
	CuAssertTrue(testCase, pseudoChromosome_getReference(pseudoChromosome2) == reference);
	testTeardown();
}

void testPseudoChromosome_getPseudoAdjacencyNumber(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, pseudoChromosome_getPseudoAdjacencyNumber(pseudoChromosome1) == 3);
	CuAssertTrue(testCase, pseudoChromosome_getPseudoAdjacencyNumber(pseudoChromosome2) == 1);
	testTeardown();
}

void testPseudoChromosome_getPseudoAdjacencyByIndex(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, pseudoChromosome_getPseudoAdjacencyByIndex(pseudoChromosome1, 0) == pseudoAdjacency1);
	CuAssertTrue(testCase, pseudoChromosome_getPseudoAdjacencyByIndex(pseudoChromosome1, 1) == pseudoAdjacency2);
	CuAssertTrue(testCase, pseudoChromosome_getPseudoAdjacencyByIndex(pseudoChromosome1, 2) == pseudoAdjacency3);
	testTeardown();
}

void testPseudoChromosome_pseudoAdjacencyIterator(CuTest* testCase) {
	testSetup();
	PseudoChromsome_PseudoAdjacencyIterator *iterator;
	iterator = pseudoChromosome_getPseudoAdjacencyIterator(pseudoChromosome1);

	CuAssertTrue(testCase, pseudoChromosome_getNextPseudoAdjacency(iterator) == pseudoAdjacency1);
	CuAssertTrue(testCase, pseudoChromosome_getNextPseudoAdjacency(iterator) == pseudoAdjacency2);
	CuAssertTrue(testCase, pseudoChromosome_getNextPseudoAdjacency(iterator) == pseudoAdjacency3);
	CuAssertTrue(testCase, pseudoChromosome_getNextPseudoAdjacency(iterator) == NULL);

	PseudoChromsome_PseudoAdjacencyIterator *iterator2;
	iterator2 = pseudoChromosome_copyPseudoChromosomeIterator(iterator);

	CuAssertTrue(testCase, pseudoChromosome_getPreviousPseudoAdjacency(iterator) == pseudoAdjacency3);
	CuAssertTrue(testCase, pseudoChromosome_getPreviousPseudoAdjacency(iterator) == pseudoAdjacency2);
	CuAssertTrue(testCase, pseudoChromosome_getPreviousPseudoAdjacency(iterator) == pseudoAdjacency1);
	CuAssertTrue(testCase, pseudoChromosome_getPreviousPseudoAdjacency(iterator) == NULL);

	CuAssertTrue(testCase, pseudoChromosome_getPreviousPseudoAdjacency(iterator2) == pseudoAdjacency3);
	CuAssertTrue(testCase, pseudoChromosome_getPreviousPseudoAdjacency(iterator2) == pseudoAdjacency2);
	CuAssertTrue(testCase, pseudoChromosome_getPreviousPseudoAdjacency(iterator2) == pseudoAdjacency1);
	CuAssertTrue(testCase, pseudoChromosome_getPreviousPseudoAdjacency(iterator2) == NULL);

	pseudoChromosome_destructPseudoAdjacencyIterator(iterator);
	pseudoChromosome_destructPseudoAdjacencyIterator(iterator2);

	testTeardown();
}

void testPseudoChromosome_serialisation(CuTest* testCase) {
	testSetup();
	int32_t i;
	void *vA = binaryRepresentation_makeBinaryRepresentation(pseudoChromosome1,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))pseudoChromosome_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	int32_t index1 = pseudoAdjacency_getIndex(pseudoAdjacency1);
	int32_t index2 = pseudoAdjacency_getIndex(pseudoAdjacency2);
	int32_t index3 = pseudoAdjacency_getIndex(pseudoAdjacency3);

	pseudoChromosome_destruct(pseudoChromosome1);
	void *vA2 = vA;
	pseudoChromosome1 = pseudoChromosome_loadFromBinaryRepresentation(&vA2, reference);
	pseudoAdjacency1 = pseudoChromosome_getPseudoAdjacencyByIndex(pseudoChromosome1, index1);
	pseudoAdjacency2 = pseudoChromosome_getPseudoAdjacencyByIndex(pseudoChromosome1, index2);
	pseudoAdjacency3 = pseudoChromosome_getPseudoAdjacencyByIndex(pseudoChromosome1, index3);
	free(vA);

	nestedTest = 1;

	testPseudoChromosome_getName(testCase);
	testPseudoChromosome_get5End(testCase);
	testPseudoChromosome_get3End(testCase);
	testPseudoChromosome_getReference(testCase);
	testPseudoChromosome_getPseudoAdjacencyNumber(testCase);
	testPseudoChromosome_getPseudoAdjacencyByIndex(testCase);
	testPseudoChromosome_pseudoAdjacencyIterator(testCase);
	testPseudoChromosome_construct(testCase);

	nestedTest = 0;
	testTeardown();
}

CuSuite* cactusPseudoChromosomeTestSuite(void) {
	CuSuite* suite = CuSuiteNew();

	SUITE_ADD_TEST(suite, testPseudoChromosome_getName);
	SUITE_ADD_TEST(suite, testPseudoChromosome_get5End);
	SUITE_ADD_TEST(suite, testPseudoChromosome_get3End);
	SUITE_ADD_TEST(suite, testPseudoChromosome_getReference);
	SUITE_ADD_TEST(suite, testPseudoChromosome_getPseudoAdjacencyNumber);
	SUITE_ADD_TEST(suite, testPseudoChromosome_getPseudoAdjacencyByIndex);
	SUITE_ADD_TEST(suite, testPseudoChromosome_pseudoAdjacencyIterator);
	SUITE_ADD_TEST(suite, testPseudoChromosome_serialisation);
	SUITE_ADD_TEST(suite, testPseudoChromosome_construct);

	return suite;
}
