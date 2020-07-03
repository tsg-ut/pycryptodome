#include "endianess.h"
#include <assert.h>

void convert_to_radix255(uint32_t out[9], const uint64_t in[4]);
void convert_from_radix255(uint64_t out[4], const uint32_t in[9]);
void reduce_25519(uint64_t x[4]);
void cswap(uint32_t a[10], uint32_t b[10], uint32_t c[10], uint32_t d[10], unsigned cond);
void invert(uint32_t out[10], const uint32_t x[10]);
void ladder(uint8_t shared_secret[32], const uint8_t *k, size_t len, const uint8_t pubkey[32]);

static const uint64_t modulus[4] =  { 0xffffffffffffffedULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL, 0x7fffffffffffffffULL };
static const uint64_t modulus2[4] = { 0xffffffffffffffdaULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL };
static const uint64_t hundhund[4] = { 0xe08063f1e8753fb4ULL, 0x29e492f797f6605cULL, 0x1f6de7b30d1327efULL, 0x534a930de945ebf3ULL };

static const uint32_t modulus_32[10] = { 0x3ffffed, 0x1ffffff, 0x3ffffff, 0x1ffffff, 0x3ffffff, 0x1ffffff, 0x3ffffff, 0x1ffffff, 0x3ffffff, 0x1ffffff };
static const uint32_t hundhund_32[10] = { 0x753fb4, 0x18fc7a, 0xb9c10, 0x1bcbfb3, 0xa7924b, 0x11327ef, 0x2f3d986, 0x17e63ed, 0xde945e, 0x14d2a4c };

void test_to(void)
{
    uint32_t out[10];
    uint64_t in[4];

    memset(out, 0xAA, sizeof out);
    memcpy(in, modulus, sizeof modulus);
    convert_to_radix255(out, in);

    assert(out[0] == 0x3ffffed);
    assert(out[1] == 0x1ffffff);
    assert(out[2] == 0x3ffffff);
    assert(out[3] == 0x1ffffff);
    assert(out[4] == 0x3ffffff);
    assert(out[5] == 0x1ffffff);
    assert(out[6] == 0x3ffffff);
    assert(out[7] == 0x1ffffff);
    assert(out[8] == 0x3ffffff);
    assert(out[9] == 0x1ffffff);

    memset(out, 0xAA, sizeof out);
    memcpy(in, hundhund, sizeof hundhund);
    convert_to_radix255(out, in);

    assert(out[0] == 0x753fb4);
    assert(out[1] == 0x18fc7a);
    assert(out[2] == 0xb9c10);
    assert(out[3] == 0x1bcbfb3);
    assert(out[4] == 0xa7924b);
    assert(out[5] == 0x11327ef);
    assert(out[6] == 0x2f3d986);
    assert(out[7] == 0x17e63ed);
    assert(out[8] == 0xde945e);
    assert(out[9] == 0x14d2a4c);

    in[0] = 0xAAAAAAAAAAAAAAAA;
    in[1] = 0xBBBBBBBBBBBBBBBB;
    in[2] = 0xCCCCCCCCCCCCCCCC;
    in[3] = 0xDDDDDDDDDDDDDDDD;
    convert_to_radix255(out, in);
    assert(out[0] == 0x2aaaaaa);
    assert(out[1] == 0xaaaaaa);
    assert(out[2] == 0x3777555);
    assert(out[3] == 0x1dddddd);
    assert(out[4] == 0x2eeeeee);
    assert(out[5] == 0xcccccc);
    assert(out[6] == 0x2666666);
    assert(out[7] == 0x1bbb999);
    assert(out[8] == 0x1dddddd);
    assert(out[9] == 0x3777777);
}

void test_from(void)
{
    uint64_t out[4];

    memset(out, 0xAA, sizeof out);
    convert_from_radix255(out, modulus_32);

    assert(out[0] == modulus[0]);
    assert(out[1] == modulus[1]);
    assert(out[2] == modulus[2]);
    assert(out[3] == modulus[3]);

    memset(out, 0xAA, sizeof out);
    convert_from_radix255(out, hundhund_32);

    assert(out[0] == hundhund[0]);
    assert(out[1] == hundhund[1]);
    assert(out[2] == hundhund[2]);
}

void test_reduce(void)
{
    uint64_t x[4];

    /** 0 **/

    x[0] = x[1] = x[2] = x[3] = 0;
    reduce_25519(x);
    assert(x[0] == 0);
    assert(x[1] == 0);
    assert(x[2] == 0);
    assert(x[3] == 0);

    memcpy(x, modulus, sizeof x);
    reduce_25519(x);
    assert(x[0] == 0);
    assert(x[1] == 0);
    assert(x[2] == 0);
    assert(x[3] == 0);

    memcpy(x, modulus2, sizeof x);
    reduce_25519(x);
    assert(x[0] == 0);
    assert(x[1] == 0);
    assert(x[2] == 0);
    assert(x[3] == 0);

    /** 1 **/

    x[0] = 1;
    x[1] = x[2] = x[3] = 0;
    reduce_25519(x);
    assert(x[0] == 1);
    assert(x[1] == 0);
    assert(x[2] == 0);
    assert(x[3] == 0);

    memcpy(x, modulus, sizeof x);
    x[0]++;
    reduce_25519(x);
    assert(x[0] == 1);
    assert(x[1] == 0);
    assert(x[2] == 0);
    assert(x[3] == 0);

    memcpy(x, modulus2, sizeof x);
    x[0]++;
    reduce_25519(x);
    assert(x[0] == 1);
    assert(x[1] == 0);
    assert(x[2] == 0);
    assert(x[3] == 0);

    /** 38 **/

    x[0] = 38;
    x[1] = x[2] = x[3] = 0;
    reduce_25519(x);
    assert(x[0] == 38);
    assert(x[1] == 0);
    assert(x[2] == 0);
    assert(x[3] == 0);

    x[0] = 0x13;
    x[1] = x[2] = 0;
    x[3] = 0x8000000000000000ULL;
    reduce_25519(x);
    assert(x[0] == 38);
    assert(x[1] == 0);
    assert(x[2] == 0);
    assert(x[3] == 0);
}

int check(uint32_t x[10], uint32_t s)
{
    unsigned i;
    for (i=0; i<10; i++)
        if (x[i] != s) {
            return 0;
        }
    return 1;
}

void test_cswap(void)
{
    uint32_t a[10];
    uint32_t b[10];
    uint32_t c[10];
    uint32_t d[10];

    memset(a, 0xAA, sizeof a);
    memset(b, 0x55, sizeof b);
    memset(c, 0x77, sizeof c);
    memset(d, 0x11, sizeof d);
    cswap(a, b, c, d, 0);

    assert(check(a, 0xAAAAAAAA));
    assert(check(b, 0x55555555));
    assert(check(c, 0x77777777));
    assert(check(d, 0x11111111));

    cswap(a, b, c, d, 1);
    assert(check(c, 0xAAAAAAAA));
    assert(check(d, 0x55555555));
    assert(check(a, 0x77777777));
    assert(check(b, 0x11111111));
}

void test_invert(void)
{
    uint64_t in[4];
    uint32_t x[10];
    uint64_t out[4];

    in[0] = 1;
    in[1] = in[2] = in[3] = 0;
    convert_to_radix255(x, in);
    invert(x, x);
    convert_from_radix255(out, x);
    reduce_25519(out);
    assert(out[0] == 1);
    assert(out[1] == 0);
    assert(out[2] == 0);
    assert(out[3] == 0);

    in[0] = 2;
    in[1] = in[2] = in[3] = 0;
    convert_to_radix255(x, in);
    invert(x, x);
    convert_from_radix255(out, x);
    reduce_25519(out);
    assert(out[0] == 0xfffffffffffffff7);
    assert(out[1] == 0xffffffffffffffff);
    assert(out[2] == 0xffffffffffffffff);
    assert(out[3] == 0x3fffffffffffffff);

    in[0] = 0xAAAAAAAAAAAAAAAA;
    in[1] = 0xBBBBBBBBBBBBBBBB;
    in[2] = 0xCCCCCCCCCCCCCCCC;
    in[3] = 0xDDDDDDDDDDDDDDDD;
    convert_to_radix255(x, in);
    invert(x, x);
    convert_from_radix255(out, x);
    reduce_25519(out);
    assert(out[0] == 0x6cf8847ba332c4c7);
    assert(out[1] == 0x984028b39c0b5e92);
    assert(out[2] == 0x2404af2276fdd005);
    assert(out[3] == 0x22336ebc77628108);
}

void test_ladder(void)
{
    uint8_t scalar[32] = "\xa5\x46\xe3\x6b\xf0\x52\x7c\x9d\x3b\x16\x15\x4b\x82\x46\x5e\xdd\x62\x14\x4c\x0a\xc1\xfc\x5a\x18\x50\x6a\x22\x44\xba\x44\x9a\xc4";
    uint8_t pubkey[32] = "\xe6\xdb\x68\x67\x58\x30\x30\xdb\x35\x94\xc1\xa4\x24\xb1\x5f\x7c\x72\x66\x24\xec\x26\xb3\x35\x3b\x10\xa9\x03\xa6\xd0\xab\x1c\x4c";
    uint8_t expout[32] = "\xc3\xda\x55\x37\x9d\xe9\xc6\x90\x8e\x94\xea\x4d\xf2\x8d\x08\x4f\x32\xec\xcf\x03\x49\x1c\x71\xf7\x54\xb4\x07\x55\x77\xa2\x85\x52";
    uint8_t out[32];

    scalar[0] &= 248;
    scalar[31] &= 127;
    scalar[31] |= 64;
    ladder(out, scalar, 32, pubkey);
    assert(0 == memcmp(out, expout, 32));
}

int main(void)
{
    test_to();
    test_from();
    test_reduce();
    test_cswap();
    test_invert();
    test_ladder();
    return 0;
}
