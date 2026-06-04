/* Measures the miniz uncompress() step CoolProp pays at first fluid load
   (FluidLibrary.cpp:49). Compress all_fluids.json with miniz, then time the
   inflate into a 7x buffer, exactly as the runtime does. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "miniz.h"

static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1.0e6;
}

int main(int argc, char** argv) {
    const char* path = (argc > 1) ? argv[1] : "dev/all_fluids.json";
    FILE* f = fopen(path, "rb");
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char* raw = (unsigned char*)malloc(n);
    if (fread(raw, 1, n, f) != (size_t)n) {
        fprintf(stderr, "read fail\n");
        return 1;
    }
    fclose(f);

    /* Compress with miniz (zlib) to produce the embedded .z blob. */
    mz_ulong clen = compressBound(n);
    unsigned char* comp = (unsigned char*)malloc(clen);
    if (compress(comp, &clen, raw, n) != MZ_OK) {
        fprintf(stderr, "compress fail\n");
        return 1;
    }
    printf("json %.2f MB  ->  zlib blob %.2f MB\n\n", n / (1024.0 * 1024.0), clen / (1024.0 * 1024.0));

    /* Time uncompress() into a 7x buffer, as FluidLibrary.cpp does. */
    const int iters = 20;
    double best = 1e30, sum = 0, worst = 0;
    unsigned char* out = (unsigned char*)malloc((size_t)clen * 7);
    for (int i = 0; i < iters; ++i) {
        mz_ulong outlen = (mz_ulong)clen * 7;
        double t0 = now_ms();
        int code = uncompress(out, &outlen, comp, clen);
        double dt = now_ms() - t0;
        if (code != MZ_OK || outlen != (mz_ulong)n) {
            fprintf(stderr, "uncompress fail %d\n", code);
            return 1;
        }
        if (dt < best) best = dt;
        if (dt > worst) worst = dt;
        sum += dt;
    }
    printf("miniz uncompress   min %6.2f ms   mean %6.2f ms   max %6.2f ms   (%d runs)\n", best, sum / iters, worst, iters);
    return 0;
}
