//  Time-stamp: "Last modified 2023-03-27 14:50:53 mluebke"
/*
**-----------------------------------------------------------------------------
** MurmurHash2 was written by Austin Appleby, and is placed in the public
** domain. The author hereby disclaims copyright to this source code.
**-----------------------------------------------------------------------------
**
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

#include "poet/HashFunctions.hpp"

#include <cstddef>
#include <openssl/evp.h>
#include <openssl/md5.h>
#include <stdexcept>

#if defined(_MSC_VER)

#define BIG_CONSTANT(x) (x)

// Other compilers

#else // defined(_MSC_VER)

#define BIG_CONSTANT(x) (x##LLU)

#endif // !defined(_MSC_VER)

// HACK: I know this is not a good practice, but this will do it for now!
EVP_MD_CTX *ctx = NULL;

void poet::initHashCtx(const EVP_MD *md) {
  if (ctx == NULL) {
    ctx = EVP_MD_CTX_new();
    EVP_DigestInit_ex(ctx, md, NULL);
  }
}

void poet::freeHashCtx() {
  EVP_MD_CTX_free(ctx);
  ctx = NULL;
}

uint64_t poet::hashDHT(int key_size, const void *key) {
  unsigned char sum[MD5_DIGEST_LENGTH];
  uint32_t md_len;
  uint64_t retval, *v1, *v2;

  // calculate md5 using MD5 functions
  EVP_DigestUpdate(ctx, key, key_size);
  EVP_DigestFinal_ex(ctx, sum, &md_len);

  if (md_len != MD5_DIGEST_LENGTH) {
    throw std::runtime_error("Something went wrong during MD5 hashing!");
  }

  // divide hash in 2 64 bit parts and XOR them
  v1 = (uint64_t *)&sum[0];
  v2 = (uint64_t *)&sum[8];
  retval = *v1 ^ *v2;

  return retval;
}

//-----------------------------------------------------------------------------
// MurmurHash2, 64-bit versions, by Austin Appleby

// The same caveats as 32-bit MurmurHash2 apply here - beware of alignment
// and endian-ness issues if used across multiple platforms.

// 64-bit hash for 64-bit platforms
// objsize: 0x170-0x321: 433

uint64_t poet::Murmur2_64A(int len, const void *key) {
  const uint64_t m = BIG_CONSTANT(0xc6a4a7935bd1e995);
  const int r = 47;

  uint64_t h = HASH_SEED ^ (len * m);

  const uint64_t *data = (const uint64_t *)key;
  const uint64_t *end = data + (len / 8);

  while (data != end) {
    uint64_t k = *data++;

    k *= m;
    k ^= k >> r;
    k *= m;

    h ^= k;
    h *= m;
  }

  const unsigned char *data2 = (const unsigned char *)data;

  switch (len & 7) {
  case 7:
    h ^= uint64_t(data2[6]) << 48;
  case 6:
    h ^= uint64_t(data2[5]) << 40;
  case 5:
    h ^= uint64_t(data2[4]) << 32;
  case 4:
    h ^= uint64_t(data2[3]) << 24;
  case 3:
    h ^= uint64_t(data2[2]) << 16;
  case 2:
    h ^= uint64_t(data2[1]) << 8;
  case 1:
    h ^= uint64_t(data2[0]);
    h *= m;
  };

  h ^= h >> r;
  h *= m;
  h ^= h >> r;

  return h;
}
