
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007-2008, J. Craig Venter Institute
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static const char *rcsid = "$Id: AS_ALN_bruteforcedp.C 4371 2013-08-01 17:19:47Z brianwalenz $";

#include "AS_global.H"
#include "AS_ALN_bruteforcedp_parallel.H"
#include "AS_UTL_reverseComplement.H"
#include <omp.h>

#undef DEBUG

#define MATCH            0
#define GAPA             1
#define GAPB             2
#define STOP             3

//  3, -2, -2 original
//  5, -4, -3 suggested
//  3, -6, -4 got around one problem, but made many more
//
#define MATCHSCORE       3
#define GAPSCORE        -6
#define MISMATCHSCORE   -4

#define SLOP             10

//  Scores in the dynamic programming matrix are unsigned ints, currently 30 bits wide.
//
//  The smallest score is 0.
//  The zero score is 2^29 (512 million)
//  The largest score is 2^30-1 (1 billion).
//
//  Edges that we don't want to end up at are initialized to DP_NEGT (2^28, 256 million) which will
//  take an awful lot of alignment to make it decent, or to underflow.
//
//  Edges that we want to end up at are set to DP_ZERO (2^29, 512 million).

#define DP_NEGT          (1 << 28)  //  A very negative value
#define DP_ZERO          (1 << 29)  //  Zero


void
alignLinker(char           *alignA,
            char           *alignB,
            char           *stringA,
            char           *stringB,
            dpCell *M2,
            alignLinker_s  *a) {

  int32 lenA = strlen(stringA);
  int32 cols = lenA;
  int32 lenB = strlen(stringB);
  int32 rows = lenB;
  int32 maxDiagLen = (rows < cols) ? rows : cols;
  int32 topMatrixCellIndex = rows * cols - 1;

//  rows = 9;
//  cols = 10;

  memset(a, 0, sizeof(alignLinker_s));

  if ((lenA > AS_READ_MAX_NORMAL_LEN) || (lenB > AS_READ_MAX_NORMAL_LEN)) {
    fprintf(stderr, "alignLinker()-- Reads too long.  %d or %d > %d\n", lenA, lenB, AS_READ_MAX_NORMAL_LEN);
    return;
  }

  int32 maxCell = topMatrixCellIndex;
  int32 scoreMax  = 0;
  int32 maxCellJ = rows;
  int32 maxCellI = cols;

  // http://stackoverflow.com/a/2112951/412582
  #pragma omp parallel default(none) shared(maxCell, scoreMax, maxCellI, maxCellJ) firstprivate(stringA, stringB, lenA, lenB, rows, cols, maxDiagLen, topMatrixCellIndex, M2)
  {

    int32 numTimesteps = rows + cols - 1;

    int32 oldOldDiagStart = 0;
    int32 oldDiagStart = 0;
    int32 diagStart = 0;

    for (int32 timestep=0; timestep < numTimesteps; timestep++) {

      int nthreads = omp_get_num_threads();
      int tid = omp_get_thread_num();

      int32 rowStart;
      int32 rowEnd;

      if (timestep < rows) {
        rowStart = timestep;
      } else {
        rowStart = rows - 1;
      }

      if (timestep < cols) {
        rowEnd = 0;
      } else {
        rowEnd = timestep - cols + 1;
      }

      int32 diagLen = rowStart - rowEnd + 1;
      int32 row = rowStart;
      int32 diagPos = 0;

      int32 col = (timestep-row);

      while (row >= rowEnd && col % nthreads != tid) {
        row--;
        col = (timestep-row);
        diagPos++;
      }

      while (row >= rowEnd) {

        int32 up = topMatrixCellIndex;
        int32 left = topMatrixCellIndex;
        int32 upLeft = topMatrixCellIndex;

        if (timestep < maxDiagLen) {
          if (diagLen - 1 == 0) {
            // this is 0,0
          } else if (diagPos == 0) {
            // no left or up left
            up = oldDiagStart + diagPos;
          } else if (diagPos == diagLen - 1) {
            // no up or up left
            left = oldDiagStart + diagPos - 1;
          } else {
            up = oldDiagStart + diagPos;
            left = up - 1;
            upLeft = oldOldDiagStart + diagPos - 1;
          }

        } else if (diagLen == maxDiagLen) {

          // if we're not in the upper triangle and this is true, it means the matrix is not square
          // so the else part of the next conditional means cols > rows

          if (rows > cols) {

            up = oldDiagStart + diagPos;

            // diagPos == 0 has no left or up left

            if (diagPos != 0) {

              left = oldDiagStart + diagPos - 1;
              upLeft = oldOldDiagStart + diagPos -1;

            }

          } else {

            left = oldDiagStart + diagPos;

            // diagPos == diagLen - 1 has no up or up left

            if (diagPos != diagLen - 1) {

              up = oldDiagStart + diagPos + 1;

              if (timestep == maxDiagLen) {

                // this is the first diag after the upper triangle
                upLeft = oldOldDiagStart + diagPos;

              } else {

                upLeft = oldOldDiagStart + diagPos + 1;

              }

            }
          }

        } else {

          // lower triangle

          up = oldDiagStart + diagPos + 1;
          left = oldDiagStart + diagPos;
          if (rows >= cols && diagLen == maxDiagLen - 1) {
            upLeft = oldOldDiagStart + diagPos;
          } else {
            upLeft = oldOldDiagStart + diagPos + 1;
          }

        }

        int32 i = col + 1;
        int32 j = row + 1;

        int upLeftScore = DP_ZERO;
        int leftScore = DP_ZERO;
        int upScore = DP_ZERO;

        if (upLeft != topMatrixCellIndex) {
          upLeftScore = M2[upLeft].score;
        }

        if (left != topMatrixCellIndex) {
          leftScore = M2[left].score;
        }

        if (up != topMatrixCellIndex) {
          upScore = M2[up].score;
        }

        //  Pick the max of these

        int ulMatrixScore = upLeftScore + ((stringA[i-1] == stringB[j-1]) ? MATCHSCORE : MISMATCHSCORE);
        int lfMatrixScore = leftScore + GAPSCORE;  //  Gap in B
        int upMatrixScore = upScore + GAPSCORE;  //  Gap in A

        //  For unitig consensus, if the consensus sequence (stringA) is lowercase, count it as a
        //  match, otherwise ignore the gap it induces in stringB.

        if (stringA[i-1] == 'a')
          if (stringB[j-1] == 'A')
            ulMatrixScore = upLeftScore + MATCHSCORE;
          else {
            ulMatrixScore = upLeftScore;
            lfMatrixScore = leftScore;
          }

        if (stringA[i-1] == 'c')
          if (stringB[j-1] == 'C')
            ulMatrixScore = upLeftScore + MATCHSCORE;
          else {
            ulMatrixScore = upLeftScore;
            lfMatrixScore = leftScore;
          }

        if (stringA[i-1] == 'g')
          if (stringB[j-1] == 'G')
            ulMatrixScore = upLeftScore + MATCHSCORE;
          else {
            ulMatrixScore = upLeftScore;
            lfMatrixScore = leftScore;
          }

        if (stringA[i-1] == 't')
          if (stringB[j-1] == 'T')
            ulMatrixScore = upLeftScore + MATCHSCORE;
          else {
            ulMatrixScore = upLeftScore;
            lfMatrixScore = leftScore;
          }

        int32 diagArrayPos = diagStart + diagPos;

        //  Set score to the smallest value possible; we will then ALWAYS pick an action below.

        dpCell newCell;

        newCell.score = 0;
        newCell.action = MATCH;
        newCell.pred = topMatrixCellIndex;

        if (newCell.score < ulMatrixScore) {
          newCell.score  = ulMatrixScore;
          newCell.action = MATCH;
          newCell.pred = upLeft;
        }

        if (newCell.score < lfMatrixScore) {
          newCell.score  = lfMatrixScore;
          newCell.action = GAPB;
          newCell.pred = left;
        }

        if (newCell.score < upMatrixScore) {
          newCell.score  = upMatrixScore;
          newCell.action = GAPA;
          newCell.pred = up;
        }

        M2[diagArrayPos] = newCell;


        row -= nthreads;
        col = (timestep-row);
        diagPos += nthreads;

      }

      #pragma omp barrier

      #pragma omp single
      {

        int32 diagStartScore = M2[diagStart].score;
        int32 diagEndScore = M2[diagStart + diagLen - 1].score;

        if (timestep < maxDiagLen) {

          if (diagLen == maxDiagLen) {

            if (rows == cols) {

              // both left and right

              if (diagStartScore > scoreMax) {
                maxCell = diagStart;
                scoreMax = diagStartScore;
                maxCellJ = rowStart + 1;
                maxCellI = timestep - rowStart + 1;
              }

              if (diagEndScore > scoreMax) {
                maxCell = diagStart + diagLen - 1;
                scoreMax = diagEndScore;
                maxCellJ = rowEnd + 1;
                maxCellI = timestep - rowEnd + 1;
              }

            } else if (rows > cols) {

              // right

              if (diagEndScore > scoreMax) {
                maxCell = diagStart + diagLen - 1;
                scoreMax = diagEndScore;
                maxCellJ = rowEnd + 1;
                maxCellI = timestep - rowEnd + 1;
              }

            } else {

              // left

              if (diagStartScore > scoreMax) {
                maxCell = diagStart;
                scoreMax = diagStartScore;
                maxCellJ = rowStart + 1;
                maxCellI = timestep - rowStart + 1;
              }

            }

          }

        } else if (diagLen == maxDiagLen) {

          // if we're not in the upper triangle and this is true, it means the matrix is not square
          // so the else part of the next conditional means cols > rows

          if (rows > cols) {

            if (rowStart == rows - 1) {

              // left
              if (diagStartScore > scoreMax) {
                maxCell = diagStart;
                scoreMax = diagStartScore;
                maxCellJ = rowStart + 1;
                maxCellI = timestep - rowStart + 1;
              }

            }

            // right

            if (diagEndScore > scoreMax) {
              maxCell = diagStart + diagLen - 1;
              scoreMax = diagEndScore;
              maxCellJ = rowEnd + 1;
              maxCellI = timestep - rowEnd + 1;
            }

          } else {

            if (rowEnd == 0) {

              // right

              if (diagEndScore > scoreMax) {
                maxCell = diagStart + diagLen - 1;
                scoreMax = diagEndScore;
                maxCellJ = rowEnd + 1;
                maxCellI = timestep - rowEnd + 1;
              }
            }

            // left

            if (diagStartScore > scoreMax) {
              maxCell = diagStart;
              scoreMax = diagStartScore;
              maxCellJ = rowStart + 1;
              maxCellI = timestep - rowStart + 1;
            }

          }

        } else {

          if (diagStart == topMatrixCellIndex) {

            // this is the last cell
            // could use either, using left here
            if (diagStartScore > scoreMax) {
              maxCell = diagStart;
              scoreMax = diagStartScore;
              maxCellJ = rowStart + 1;
              maxCellI = timestep - rowStart + 1;
            }

          } else {

            // left and right

            if (diagStartScore > scoreMax) {
              maxCell = diagStart;
              scoreMax = diagStartScore;
              maxCellJ = rowStart + 1;
              maxCellI = timestep - rowStart + 1;
            }

            if (diagEndScore > scoreMax) {
              maxCell = diagStart + diagLen - 1;
              scoreMax = diagEndScore;
              maxCellJ = rowEnd + 1;
              maxCellI = timestep - rowEnd + 1;
            }

          }

        }

      }

      oldOldDiagStart = oldDiagStart;
      oldDiagStart = diagStart;
      diagStart += diagLen;

    }

  }

  int32 endI=maxCellI, curI=maxCellI;
  int32 endJ=maxCellJ, curJ=maxCellJ;

  int32  alignLen  = 0;

  int32  nGapA     = 0;
  int32  nGapB     = 0;
  int32  nMatch    = 0;
  int32  nMismatch = 0;

  int32 curCellPos = maxCell;

  do {

    dpCell curCell = M2[curCellPos];

    switch (curCell.action) {

      case MATCH:
        alignA[alignLen] = stringA[curI-1];
        alignB[alignLen] = stringB[curJ-1];

        if (alignA[alignLen] == alignB[alignLen]) {
          alignA[alignLen] = tolower(alignA[alignLen]);
          alignB[alignLen] = tolower(alignB[alignLen]);
          nMatch++;
        } else {
          alignA[alignLen] = toupper(alignA[alignLen]);
          alignB[alignLen] = toupper(alignB[alignLen]);
          nMismatch++;
        }

        curI--;
        curJ--;
        alignLen++;

        break;
      case GAPA:
        alignA[alignLen] = '-';
        alignB[alignLen] = stringB[curJ-1];
        alignLen++;
        nGapA++;
        curJ--;
        break;
      case GAPB:
        alignA[alignLen] = stringA[curI-1];
        alignB[alignLen] = '-';
        alignLen++;
        nGapB++;
        curI--;
        break;

    }

    curCellPos = curCell.pred;

  } while (curCellPos != topMatrixCellIndex);

  alignA[alignLen] = 0;
  alignB[alignLen] = 0;

  reverse(alignA, alignB, alignLen);

  a->matches    = nMatch;
  a->alignLen   = alignLen;
  a->begI       = curI;
  a->begJ       = curJ;
  a->endI       = endI;
  a->endJ       = endJ;
  a->lenA       = lenA;
  a->lenB       = lenB;
  a->pIdentity  = (double)(nMatch) / (double)(nGapA + nGapB + nMatch + nMismatch);
  a->pCoverageA = (double)(a->endI - a->begI) / (double)(lenA);
  a->pCoverageB = (double)(a->endJ - a->begJ) / (double)(lenB);

}
