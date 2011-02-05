/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_REFERENCE_H_
#define CACTUS_REFERENCE_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Cactus reference ordering functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * A reference ordering structure for the non-stub ends in the flower.
 *
 * There are three types of end in any flower:
 * 'attached stub-ends', 'free stub-ends' and 'block ends'.
 *
 * Each block has two 'block ends', that represent the two ends of the block.
 *
 * 'Stub-ends' represent telomeres, or the end of blocks that are not observed in the considered
 * flower (but may be observed in a parent flower). Stub ends tell us about the boundaries
 * of the sequences in a flower. Stub-ends are further broken down into 'free stub-ends' and 'attached stub-ends'.
 *
 * Free stub-ends are treated as missing information, i.e. the ends of blocks or telomeres
 * which we did not observe in any flower. They are not part of the reference structure.,
 *
 * Attached stub ends must be present in any parent flower, and if they represent a block end, the
 * block must be present in a parent flower.
 *
 * A reference ordering orders all the attached stub ends and block ends in a flower
 * into a set of linear "pseudo chromosomes".
 * Let (a, b) represent a "pseudo-adjacency" between two ends a and b, so called because
 * there may not necessarily exist an actual adjacency between caps representing
 * the two ends.
 *
 * A reference ordering can be represented as a 2-d list of lists of pseudo adjacencies.
 *
 * ( ( (L_{0, 0}, R_{0, 0}), (L_{0, 1}, R_{0, 1}) ... (L_{0, N}, R_{0, N}) ),
 *   ( (L_{1, 0}, R_{1, 0}), (L_{1, 1}, R_{1, 1}) ... (L_{1, O}, R_{1, O}) ),
 *   ...
 *   ( (L_{M, 0}, R_{M, 0}), (L_{M, 1}, R_{M, 1}) ... (L_{M, P}, R_{M, P}) ) )
 *
 * Where pair (L_{i, j}, R_{i, j}) represents the jth pseudo-adjacency, 0 <= j <= N, in the i_th  pseudo-chromosome, 0 <= i <= M.
 *
 * A reference ordering REF is valid for a flower only if :
 * (1) the first, L_{i, 0}, and last, R_{i, N}, ends
 * of the first and last pseudo-adjacencies of each pseudo-chromosome
 * in REF are attached stub ends.
 * Call this property of REF "properly terminated"
 * (2) for each pseudo adjacency (a, b) both a and b are in the same group,
 * call this property of REF "group respecting".
 * (3) for each pair of contiguous psuedo-adjacencies
 *  ((L_{i, j}, R_{i, j}), (L_{i, j+1}, R_{i, j+1}))
 *  in a pseudo-chromosome in REF, R_{i, j} L_{i, j+1} are opposite ends of the same block.
 *  Call this property of REF "block respecting".
 *
 * A reference ordering REF is valid for flower NET if and only if REF contains only
 * the complete set of attached stub end and block ends for NET and is properly terminated,
 * group respecting and block respecting.
 */
Reference *reference_construct(Flower *flower);

/*
 * Destructor for reference.
 */
void reference_destruct(Reference *reference);

/*
 * Gets flower associated with reference.
 */
Flower *reference_getFlower(Reference *reference);

/*
 * Returns the number of pseudo-chromosomes in the reference.
 */
int32_t reference_getPseudoChromosomeNumber(Reference *reference);

/*
 * Get pseudo-chromosome with the given name.
 */
PseudoChromosome *reference_getPseudoChromosome(Reference *reference, Name name);

/*
 * Gets the first pseudo chromosome in the list of pseudo-chromosomes.
 */
PseudoChromosome *reference_getFirst(Reference *reference);

/*
 * Gets an iterator to iterate through the pseudo-chromosomes in the reference.
 */
Reference_PseudoChromosomeIterator *reference_getPseudoChromosomeIterator(Reference *reference);

/*
 * Gets the next pseudo-chromosome from the iterator.
 */
PseudoChromosome *reference_getNextPseudoChromosome(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator);

/*
 * Gets the previous pseudo-chromosome from the iterator.
 */
PseudoChromosome *reference_getPreviousPseudoChromosome(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator);

/*
 * Duplicates the iterator.
 */
Reference_PseudoChromosomeIterator *reference_copyPseudoChromosomeIterator(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator);

/*
 * Destructs the iterator.
 */
void reference_destructPseudoChromosomeIterator(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator);

/*
 * The following is checked by this function.. which creates an assert error if it fails.
 *
 * The every attached stub and block end is in exactly one pseudo-adjacency in the reference.
 * That free stubs are not in the reference.
 * That every pseudo-chromosome contains one or more pseudo-adjacencies.
 * That order of ends in the pseudo chromosome is valid (i.e. that the 5 and 3 prime ends
 * match and that internal pseudo-adjacencies span blocks).
 * That the two ends in every pseudo-adjacency are in the same descendant group.
 * That all ends in the reference are in the positive orientation.
 * That for all pseudo adjacencies whose end are in a link, the link contains both ends of the
 * pseudo adjacency.
 */
void reference_check(Reference *reference);

#endif
