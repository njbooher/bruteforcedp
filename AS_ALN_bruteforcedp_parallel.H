
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

#ifndef AS_ALN_BRUTEFORCEDP_PARALLEL
#define AS_ALN_BRUTEFORCEDP_PARALLEL

#include "AS_global.H"

static const char *rcsid_AS_ALN_BRUTEFORCEDP = "$Id: AS_ALN_bruteforcedp.H 4371 2013-08-01 17:19:47Z brianwalenz $";


/*
 * Only need the score for the last 2 diagonals, and in the end following action
 * back through the matrix is all that matters
 *
 * could save space by making pred an offset in the matrix (store as unsigned but treat as negative always)
 * could determine if up-left action based soley on how big the offset is
 * if wanted to add more stupid logic, could also determine left/up by calculating index
 * might be able to get away with 16bit int then
 *
 * could also store diagonal start indicies as iterate, and diagPos of maxCell
 * then with stored action and current diagPos can calculate index of stored action cell in that direction
 *
 */
typedef struct {
  uint64  score  : 30;
  uint64  action : 2;
  uint64  pred   : 32;
} dpCell;


typedef struct {
  int32   matches;
  int32   alignLen;
  int32   begI, begJ;
  int32   endI, endJ;
  int32   lenA, lenB;

  double  pIdentity;
  double  pCoverageA;
  double  pCoverageB;
} alignLinker_s;


void
alignLinker(char           *alignA,
            char           *alignB,
            char           *stringA,
            char           *stringB,
            dpCell *M2,
            alignLinker_s  *a);

#endif  //  AS_ALN_BRUTEFORCEDP_PARALLEL
