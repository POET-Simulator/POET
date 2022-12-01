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

#include "poet/HashFunctions.hpp"

#include <cstddef>
#include <openssl/evp.h>
#include <openssl/md5.h>
#include <stdexcept>

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

uint64_t poet::hashDHT(int key_size, void *key) {
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
