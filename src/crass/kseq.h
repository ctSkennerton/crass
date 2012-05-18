/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

/* Last Modified: 12APR2009 */

/* De-macro'd by the ACE-team 18MAY2012 wattup! */

#ifndef AC_KSEQ_H
  #define AC_KSEQ_H

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

typedef struct //__kstream_t
{
	char *buf;
	int begin, end, is_eof;
	gzFile f;
} kstream_t;

typedef struct //__kstring_t
{
        size_t l, m;
        char *s;
} kstring_t;

typedef struct
{
	kstring_t name, comment, seq, qual;
	int last_char;
	kstream_t *f;
} kseq_t;

kstream_t *ks_init(gzFile f);
void ks_destroy(kstream_t *ks);
int ks_getc(kstream_t *ks);
int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret);
kseq_t *kseq_init(gzFile fd);
void kseq_rewind(kseq_t *ks);
void kseq_destroy(kseq_t *ks);
int kseq_read(kseq_t *seq);

#endif