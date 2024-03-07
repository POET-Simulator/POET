// //  Time-stamp: "Last modified 2023-03-31 14:59:49 mluebke"

/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2022 Marco De Lucia, Max Luebke (GFZ Potsdam)
**
** POET is free software; you can redistribute it and/or modify it under the
** terms of the GNU General Public License as published by the Free Software
** Foundation; either version 2 of the License, or (at your option) any later
** version.
**
** POET is distributed in the hope that it will be useful, but WITHOUT ANY
** WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
** A PARTICULAR PURPOSE. See the GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License along with
** this program; if not, write to the Free Software Foundation, Inc., 51
** Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef HASHFUNCTIONS_H_
#define HASHFUNCTIONS_H_

#include <cstdint>

namespace poet {

// Sum of POET interpreted as ASCII
constexpr uint32_t HASH_SEED = 80 + 79 + 69 + 84;

uint64_t Murmur2_64A(int len, const void *key);
uint64_t md5_sum(int len, const void *key);

} // namespace poet

#endif // HASHFUNCTIONS_H_
