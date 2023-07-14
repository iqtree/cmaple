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
 //
// Some math & matrix functions  (some using SSE/AVX/NEON for significantly better speed)
//

// let's include some x64 instrinsics from SIMDe; SIMDe will translate to Neon for ARM automatically
#include <simde/x86/sse2.h>
#include <simde/x86/sse4.1.h>
#include <simde/x86/avx.h>

#include <utils/tools.h>

#include <assert.h>
#include <vector>

// Multiply 4 floats by another 4 floats.
// inspired by https://stackoverflow.com/a/59495197
template<int offsetRegs>
inline simde__m128 mul4(const float* p1, const float* p2)
{
  constexpr int lanes = offsetRegs * 4;
  const simde__m128 a = simde_mm_loadu_ps(p1 + lanes);
  const simde__m128 b = simde_mm_loadu_ps(p2 + lanes);
  return simde_mm_mul_ps(a, b);
}

// Multiply 4 doubles by another 4 doubles.
template<int offsetRegs>
inline simde__m256d mul4(const double* p1, const double* p2)
{
  constexpr int lanes = offsetRegs * 4;
  const simde__m256d a = simde_mm256_loadu_pd(p1 + lanes);
  const simde__m256d b = simde_mm256_loadu_pd(p2 + lanes);
  return simde_mm256_mul_pd(a, b);
}

// sum up 4 floats (SSE)
// see https://stackoverflow.com/a/59495197
inline float horiz_sum(simde__m128 v) {
  // Add 4 values into 2
  const simde__m128 r2 = simde_mm_add_ps(v, simde_mm_movehl_ps(v, v));
  // Add 2 lower values into the final result
  const simde__m128 r1 = simde_mm_add_ss(r2, simde_mm_movehdup_ps(r2));
  // Return the lowest lane of the result vector.
  // The intrinsic below compiles into noop, modern compilers return floats in the lowest lane of xmm0 register.
  return simde_mm_cvtss_f32(r1);
}

// sum up 4 doubles (AVX)
// see https://stackoverflow.com/a/49943540
inline double horiz_sum(simde__m256d v) {
  simde__m128d vlow = simde_mm256_castpd256_pd128(v);
  simde__m128d vhigh = simde_mm256_extractf128_pd(v, 1); // high 128
  vlow = simde_mm_add_pd(vlow, vhigh);     // reduce down to 128

  simde__m128d high64 = simde_mm_unpackhi_pd(vlow, vlow);
  return  simde_mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}


// Compute dot product of float vectors
template <int length, typename RealType>
inline RealType dotProduct(const RealType* p1, const RealType* p2)
{ // this is the default and 'fine' for DNA etc, when the number of operations is small. For AA's (length==20) there are specialized versions available.
  RealType result{ 0 };
  for (int j = 0; j < length; ++j)
  {
    result += p1[j] * p2[j];
  }
  return result;
}

template <>
inline float dotProduct<20>(const float* p1, const float* p2)
{
  const int count = 20;
  const float* const p1End = p1 + count;

  // Process all 20 values. Nothing to add yet, just multiplying.
  auto dot0 = mul4<0>(p1, p2);
  auto dot1 = mul4<1>(p1, p2);
  auto dot2 = mul4<2>(p1, p2);
  auto dot3 = mul4<3>(p1, p2);
  auto dot4 = mul4<4>(p1, p2);

  // 20 to 4
  const auto dot01 = simde_mm_add_ps(dot0, dot1);
  const auto dot23 = simde_mm_add_ps(dot2, dot3);
  const auto dot401 = simde_mm_add_ps(dot4, dot01);
  const auto dot23401 = simde_mm_add_ps(dot23, dot401);

  return horiz_sum(dot23401);
}

template <>
inline double dotProduct<20>(const double* p1, const double* p2)
{
  const int count = 20;
  const double* const p1End = p1 + count;

  // Process all 20 values. Nothing to add yet, just multiplying.
  auto dot0 = mul4<0>(p1, p2);
  auto dot1 = mul4<1>(p1, p2);
  auto dot2 = mul4<2>(p1, p2);
  auto dot3 = mul4<3>(p1, p2);
  auto dot4 = mul4<4>(p1, p2);

  // 20 to 4
  const auto dot01 = simde_mm256_add_pd(dot0, dot1);
  const auto dot23 = simde_mm256_add_pd(dot2, dot3);
  const auto dot401 = simde_mm256_add_pd(dot4, dot01);
  const auto dot23401 = simde_mm256_add_pd(dot23, dot401);

  return horiz_sum(dot23401);
}


template <StateType length>
RealNumType sumMutationByLh(const RealNumType* const vec1, const RealNumType* const vec2)
{
    RealNumType result{0};
    for (StateType j = 0; j < length; ++j)
        result += (vec1[j] > 0.1 ? vec2[j] : 0);
    
  return result;
}


template <StateType length>
RealNumType matrixEvolve(const RealNumType* const vec1,
                         const RealNumType* const vec2,
                         const RealNumType* mutation_mat_row,
                         const RealNumType total_blength)
{
  RealNumType result{ 0 };
  for (StateType i = 0; i < length; ++i, mutation_mat_row += length)
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

template <StateType length>
RealNumType matrixEvolveRoot(const RealNumType* const vec2,
                             const StateType seq1_state,
                             const RealNumType* model_root_freqs, 
                             const RealNumType* transposed_mut_mat_row, 
                             const RealNumType* mutation_mat_row, 
                             const RealNumType total_blength, 
                             const RealNumType seq1_region_plength_observation2node)
{
  RealNumType result{ 0 };
  for (StateType i = 0; i < length; ++i, mutation_mat_row += length)
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

template <StateType length>
RealNumType updateVecWithState(RealNumType* const update_vec, const StateType seq1_state,
  const RealNumType* const vec,
  const RealNumType factor)
{
  RealNumType result{0};
  for (StateType i = 0; i < length; ++i)
  {
    if (i == seq1_state)
      update_vec[i] *= (1.0 + vec[i] * factor);
    else
      update_vec[i] *= vec[i] * factor;
    result += update_vec[i];
  }
  return result;
}

template <StateType length>
void setVecWithState(RealNumType* const set_vec, const StateType seq1_state,
  const RealNumType* const vec,
  const RealNumType factor)
{
  for (StateType i = 0; i < length; ++i)
    set_vec[i] = vec[i] * factor;

    set_vec[seq1_state] += 1.0;
}

template <StateType length>
void updateCoeffs(RealNumType* const root_freqs, RealNumType* const transposed_mut_mat_row, RealNumType* const likelihood, RealNumType* const mutation_mat_row, const RealNumType factor, RealNumType& coeff0, RealNumType& coeff1)
{
    for (StateType i = 0; i < length; ++i)
    {
        coeff0 += root_freqs[i] * transposed_mut_mat_row[i] * factor * likelihood[i];
        coeff1 += mutation_mat_row[i] * likelihood[i];
    }
}

template <StateType length>
void setVecByProduct(RealNumType* const set_vec, const RealNumType* const vec1, const RealNumType* const vec2)
{
    for (StateType j = 0; j < length; ++j)
        set_vec[j] = vec1[j] * vec2[j];
}

/* NHANLT: I'm not sure if there is an AVX instruction to reset all entries of a vector to zero */
template <StateType length>
void resetVec(RealNumType* const set_vec)
{
    for (StateType i = 0; i < length; ++i)
        set_vec[i] = 0;
}

template <StateType length>
RealNumType resetLhVecExceptState(RealNumType* const set_vec, const StateType state, const RealNumType state_lh)
{
    resetVec<length>(set_vec);
    
    set_vec[state] = state_lh;
    
    return state_lh;
}
