/***************************************************************************
 *   Copyright (C) 2022 by                                                 *
 *   BUI Quang Minh <m.bui@anu.edu.au>                                     *
 *   Nhan Ly-Trong <trongnhan.uit@gmail.com>                               *
 *   Chris Bielow <chris.bielow@fu-berlin.de>                              *
 *                                                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#pragma once

#include <utils/tools.h>

 //
// Some math & matrix functions
//
// To be converted to SSE/AVX eventually
//


template <int length>
RealNumType dotProduct(const RealNumType* const vec1, const RealNumType* const vec2)
{
  RealNumType result{0};
  for (int j = 0; j < length; ++j)
  {
    result += vec1[j] * vec2[j];
  }
  return result;
}


template <int length>
RealNumType matrixEvolve(const RealNumType* const vec1,
                         const RealNumType* const vec2,
                         const RealNumType* mutation_mat_row,
                         const RealNumType total_blength)
{
  RealNumType result{ 0 };
  for (int i = 0; i < length; ++i, mutation_mat_row += length)
  {
    // NHANLT NOTE:
    // tot2: likelihood of i evolves to j
    // tot2 = (1 + mut[i,i] * total_blength) * lh(seq2,i) + mut[i,j] * total_blength * lh(seq2,j)
    RealNumType tot2 = dotProduct<length>(mutation_mat_row, vec2);

    // NHANLT NOTE:
    // tot = the likelihood of observing i * the likelihood of i evolves to j
    result += vec1[i] * (vec2[i] + total_blength * tot2);
  }
  return result;
}

template <int length>
RealNumType matrixEvolveRoot(const RealNumType* const vec1,
                             const RealNumType* const vec2,
                             const int seq1_state,
                             const RealNumType* model_root_freqs, 
                             const RealNumType* transposed_mut_mat_row, 
                             const RealNumType* mutation_mat_row, 
                             const RealNumType total_blength, 
                             const RealNumType seq1_region_plength_observation2node)
{
  RealNumType result{ 0 };
  for (int i = 0; i < length; ++i, mutation_mat_row += length)
  {
    // NHANLT NOTE: UNSURE
    // tot2: likelihood that we can observe seq1_state elvoving from i (from root) (account for the fact that the observation might have occurred on the other side of the phylogeny with respect to the root)
    // tot2 = root_freqs[seq1_state] * (1 + mut[seq1_state,seq1_state] * plength_observation2node) + root_freqs[i] * mut[i,seq1_state] * plength_observation2node
    RealNumType tot2;

    if (seq1_state == i)
      tot2 = model_root_freqs[i] * (1.0 + transposed_mut_mat_row[i] * seq1_region_plength_observation2node);
    else
      tot2 = model_root_freqs[i] * (transposed_mut_mat_row[i] * seq1_region_plength_observation2node);

    // NHANLT NOTE:
    // tot3: likelihood of i evolves to j
    // tot3 = (1 + mut[i,i] * total_blength) * lh(seq2,i) + mut[i,j] * total_blength * lh(seq2,j)
    RealNumType tot3 = dotProduct<length>(mutation_mat_row, vec2);
    result += tot2 * (vec2[i] + total_blength * tot3);
  }
  return result;
}
