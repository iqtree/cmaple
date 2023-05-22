/****************************************************************************
 *   Copyright (C) 2022 by
 *   Nhan Ly-Trong <trongnhan.uit@gmail.com>
 *   Chris Bielow <chris.bielow@fu-berlin.de>
 *   Nicola De Maio <demaio@ebi.ac.uk>
 *   BUI Quang Minh <m.bui@anu.edu.au>
 *
 *
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the
 *   Free Software Foundation, Inc.,
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 ***************************************************************************/

#pragma once

#include "tools.h"

 //
// Some math & matrix functions
//
// To be converted to SSE/AVX eventually
//


template <cmaple::StateType length>
cmaple::RealNumType dotProduct(const cmaple::RealNumType* const vec1, const cmaple::RealNumType* const vec2)
{
    cmaple::RealNumType result{0};
    for (cmaple::StateType j = 0; j < length; ++j)
  {
    result += vec1[j] * vec2[j];
  }
  return result;
}

template <cmaple::StateType length>
cmaple::RealNumType sumMutationByLh(const cmaple::RealNumType* const vec1, const cmaple::RealNumType* const vec2)
{
    cmaple::RealNumType result{0};
    for (cmaple::StateType j = 0; j < length; ++j)
        result += (vec1[j] > 0.1 ? vec2[j] : 0);
    
  return result;
}


template <cmaple::StateType length>
cmaple::RealNumType matrixEvolve(const cmaple::RealNumType* const vec1,
                                 const cmaple::RealNumType* const vec2,
                                 const cmaple::RealNumType* mutation_mat_row,
                                 const cmaple::RealNumType total_blength)
{
    cmaple::RealNumType result{ 0 };
    for (cmaple::StateType i = 0; i < length; ++i, mutation_mat_row += length)
  {
    // NHANLT NOTE:
    // tot2: likelihood of i evolves to j
    // tot2 = (1 + mut[i,i] * total_blength) * lh(seq2,i) + mut[i,j] * total_blength * lh(seq2,j)
      cmaple::RealNumType tot2 = dotProduct<length>(mutation_mat_row, vec2);

    // NHANLT NOTE:
    // tot = the likelihood of observing i * the likelihood of i evolves to j
    result += vec1[i] * (vec2[i] + total_blength * tot2);
  }
  return result;
}

template <cmaple::StateType length>
cmaple::RealNumType matrixEvolveRoot(const cmaple::RealNumType* const vec2,
                                     const cmaple::StateType seq1_state,
                                     const cmaple::RealNumType* model_root_freqs,
                                     const cmaple::RealNumType* transposed_mut_mat_row,
                                     const cmaple::RealNumType* mutation_mat_row,
                                     const cmaple::RealNumType total_blength,
                                     const cmaple::RealNumType seq1_region_plength_observation2node)
{
    cmaple::RealNumType result{ 0 };
    for (cmaple::StateType i = 0; i < length; ++i, mutation_mat_row += length)
  {
    // NHANLT NOTE: UNSURE
    // tot2: likelihood that we can observe seq1_state elvoving from i (from root) (account for the fact that the observation might have occurred on the other side of the phylogeny with respect to the root)
    // tot2 = root_freqs[seq1_state] * (1 + mut[seq1_state,seq1_state] * plength_observation2node) + root_freqs[i] * mut[i,seq1_state] * plength_observation2node
      cmaple::RealNumType tot2;

    if (seq1_state == i)
      tot2 = model_root_freqs[i] * (1.0 + transposed_mut_mat_row[i] * seq1_region_plength_observation2node);
    else
      tot2 = model_root_freqs[i] * (transposed_mut_mat_row[i] * seq1_region_plength_observation2node);

    // NHANLT NOTE:
    // tot3: likelihood of i evolves to j
    // tot3 = (1 + mut[i,i] * total_blength) * lh(seq2,i) + mut[i,j] * total_blength * lh(seq2,j)
      cmaple::RealNumType tot3 = dotProduct<length>(mutation_mat_row, vec2);
    result += tot2 * (vec2[i] + total_blength * tot3);
  }
  return result;
}

template <cmaple::StateType length>
cmaple::RealNumType updateVecWithState(cmaple::RealNumType* const update_vec, const cmaple::StateType seq1_state,
                               const cmaple::RealNumType* const vec,
                               const cmaple::RealNumType factor)
{
    cmaple::RealNumType result{0};
    for (cmaple::StateType i = 0; i < length; ++i)
  {
    if (i == seq1_state)
      update_vec[i] *= (1.0 + vec[i] * factor);
    else
      update_vec[i] *= vec[i] * factor;
    result += update_vec[i];
  }
  return result;
}

template <cmaple::StateType length>
void setVecWithState(cmaple::RealNumType* const set_vec, const cmaple::StateType seq1_state,
  const cmaple::RealNumType* const vec,
  const cmaple::RealNumType factor)
{
  for (cmaple::StateType i = 0; i < length; ++i)
    set_vec[i] = vec[i] * factor;

    set_vec[seq1_state] += 1.0;
}

template <cmaple::StateType length>
void updateCoeffs(cmaple::RealNumType* const root_freqs, cmaple::RealNumType* const transposed_mut_mat_row, cmaple::RealNumType* const likelihood, cmaple::RealNumType* const mutation_mat_row, const cmaple::RealNumType factor, cmaple::RealNumType& coeff0, cmaple::RealNumType& coeff1)
{
    for (cmaple::StateType i = 0; i < length; ++i)
    {
        coeff0 += root_freqs[i] * transposed_mut_mat_row[i] * factor * likelihood[i];
        coeff1 += mutation_mat_row[i] * likelihood[i];
    }
}

template <cmaple::StateType length>
void setVecByProduct(cmaple::RealNumType* const set_vec, const cmaple::RealNumType* const vec1, const cmaple::RealNumType* const vec2)
{
    for (cmaple::StateType j = 0; j < length; ++j)
        set_vec[j] = vec1[j] * vec2[j];
}

/* NHANLT: I'm not sure if there is an AVX instruction to reset all entries of a vector to zero */
template <cmaple::StateType length>
void resetVec(cmaple::RealNumType* const set_vec)
{
    for (cmaple::StateType i = 0; i < length; ++i)
        set_vec[i] = 0;
}

template <cmaple::StateType length>
cmaple::RealNumType resetLhVecExceptState(cmaple::RealNumType* const set_vec, const cmaple::StateType state, const cmaple::RealNumType state_lh)
{
    resetVec<length>(set_vec);
    
    set_vec[state] = state_lh;
    
    return state_lh;
}
