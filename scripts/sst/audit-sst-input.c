#include <eccodes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

static void checked(int code, const char *operation) {
    if (code) {
        fprintf(stderr, "FAIL: %s: %s\n", operation, codes_get_error_message(code));
        exit(2);
    }
}
static long integer(codes_handle *h, const char *key) {
    long value;
    checked(codes_get_long(h, key, &value), key);
    return value;
}
static double *values(codes_handle *h, size_t *count) {
    checked(codes_get_size(h, "values", count), "values count");
    double *result = malloc(*count * sizeof(*result));
    if (!result) { perror("malloc"); exit(2); }
    size_t decoded = *count;
    checked(codes_get_double_array(h, "values", result, &decoded), "decode values");
    if (decoded != *count) { fprintf(stderr, "FAIL: incomplete decode\n"); exit(2); }
    return result;
}
static void grid(codes_handle *h, char result[128]) {
    size_t length = 128;
    char type[128];
    checked(codes_get_string(h, "gridType", type, &length), "grid type");
    if (strcmp(type,"reduced_gg")) { fprintf(stderr,"FAIL: expected reduced Gaussian grid\n"); exit(2); }
    /* GRIB-1 rounds the global east bound differently (359.929 vs 359.930).
       Actual global reduced-grid point locations are determined by N and pl.
       Check the global domain, all scan/latitude keys and every row length. */
    long east = integer(h,"longitudeOfLastGridPoint");
    if (east < 359800 || east >= 360000) { fprintf(stderr,"FAIL: not a global grid\n"); exit(2); }
    const char *keys[] = {"N","Nj","numberOfDataPoints","scanningMode", "latitudeOfFirstGridPoint",
                         "latitudeOfLastGridPoint","longitudeOfFirstGridPoint"};
    uint64_t digest = UINT64_C(14695981039346656037);
    for (size_t k=0; k<sizeof(keys)/sizeof(keys[0]); k++) {
        uint64_t item=(uint64_t)integer(h,keys[k]);
        for (int byte=0; byte<8; byte++) { digest=(digest^(item&255))*UINT64_C(1099511628211); item>>=8; }
    }
    size_t rows;
    checked(codes_get_size(h,"pl",&rows),"row count");
    long *pl=malloc(rows*sizeof(*pl));
    if (!pl) { perror("malloc row list"); exit(2); }
    checked(codes_get_long_array(h,"pl",pl,&rows),"row lengths");
    for (size_t k=0; k<rows; k++) {
        uint64_t item=(uint64_t)pl[k];
        for (int byte=0; byte<8; byte++) { digest=(digest^(item&255))*UINT64_C(1099511628211); item>>=8; }
    }
    free(pl);
    snprintf(result,128,"global_reduced_gg_%016" PRIx64,digest);
}
int main(int argc, char **argv) {
    if (argc != 3 && argc != 4) {
        fprintf(stderr, "Usage: %s SURFACE_INIT CLIMATE_INIT [MASK_DATE]\n", argv[0]);
        return 2;
    }
    long requested_mask_date = argc == 4 ? strtol(argv[3], NULL, 10) : 0;
    FILE *initial = fopen(argv[1], "rb");
    if (!initial) { perror(argv[1]); return 2; }
    double *mask = NULL, *sst = NULL;
    size_t mask_count = 0, sst_count = 0;
    char mask_grid[128], sst_grid[128];
    int err = 0;
    codes_handle *h;
    long initial_date = 0;
    while ((h = codes_handle_new_from_file(NULL, initial, PRODUCT_GRIB, &err))) {
        long parameter = integer(h, "paramId");
        if (parameter == 172 && (!requested_mask_date || integer(h,"dataDate") == requested_mask_date)) {
            if (mask) {
                size_t duplicate_count;
                char duplicate_grid[128];
                double *duplicate = values(h, &duplicate_count);
                grid(h, duplicate_grid);
                if (duplicate_count != mask_count || strcmp(duplicate_grid, mask_grid) ||
                    memcmp(mask, duplicate, mask_count * sizeof(*mask))) {
                    fprintf(stderr, "FAIL: multiple different masks require model-read selection\n");
                    return 2;
                }
                printf("INFO: duplicate mask at date=%ld is numerically identical\n", integer(h,"dataDate"));
                free(duplicate);
            } else {
                mask = values(h, &mask_count);
                grid(h, mask_grid);
            }
        } else if (parameter == 34) {
            if (sst) { fprintf(stderr, "FAIL: duplicate initial SST\n"); return 2; }
            sst = values(h, &sst_count);
            grid(h, sst_grid);
            initial_date = integer(h, "dataDate");
        }
        codes_handle_delete(h);
    }
    checked(err, "read initial GRIB");
    fclose(initial);
    if (!mask || !sst || mask_count != sst_count || strcmp(mask_grid, sst_grid)) {
        fprintf(stderr, "FAIL: initial SST and mask are missing or on different grids\n");
        return 2;
    }
    size_t sea_count = 0, invalid_initial = 0;
    double initial_min = INFINITY, initial_max = -INFINITY;
    for (size_t i = 0; i < mask_count; i++) {
        if (!isfinite(mask[i]) || mask[i] < 0 || mask[i] > 1) {
            fprintf(stderr, "FAIL: invalid mask at index %zu\n", i); return 2;
        }
        /* The current UPDCLIE implementation updates SST where LSM <= 0.5. */
        if (mask[i] <= 0.5) {
            sea_count++;
            if (!isfinite(sst[i]) || sst[i] < 250 || sst[i] > 330) invalid_initial++;
            if (sst[i] < initial_min) initial_min = sst[i];
            if (sst[i] > initial_max) initial_max = sst[i];
        }
    }
    printf("Initial date=%ld mask_date=%ld points=%zu ocean_points=%zu ocean_SST_min=%.10g max=%.10g invalid=%zu grid=%s\n",
           initial_date, requested_mask_date, sst_count, sea_count, initial_min, initial_max, invalid_initial, sst_grid);
    FILE *climate = fopen(argv[2], "rb");
    if (!climate) { perror(argv[2]); return 2; }
    size_t frames = 0, invalid_total = 0;
    while ((h = codes_handle_new_from_file(NULL, climate, PRODUCT_GRIB, &err))) {
        if (integer(h, "paramId") == 139) {
            size_t count;
            double *temperature = values(h, &count);
            char climate_grid[128];
            grid(h, climate_grid);
            if (count != sst_count || strcmp(climate_grid, sst_grid)) {
                fprintf(stderr, "FAIL: climate temperature is on a different grid\n"); return 2;
            }
            size_t invalid = 0, nonfinite = 0;
            double minimum = INFINITY, maximum = -INFINITY, difference = 0;
            for (size_t i = 0; i < count; i++) {
                if (mask[i] > 0.5) continue;
                if (!isfinite(temperature[i])) nonfinite++;
                if (!isfinite(temperature[i]) || temperature[i] < 250 || temperature[i] > 330) invalid++;
                if (temperature[i] < minimum) minimum = temperature[i];
                if (temperature[i] > maximum) maximum = temperature[i];
                if (fabs(temperature[i] - sst[i]) > difference) difference = fabs(temperature[i] - sst[i]);
            }
            printf("Climate date=%ld time=%ld ocean_SST_min=%.10g max=%.10g outside_250_330K=%zu nonfinite=%zu max_difference_from_initial=%.10g\n",
                   integer(h, "dataDate"), integer(h, "dataTime"), minimum, maximum, invalid, nonfinite, difference);
            frames++;
            invalid_total += invalid;
            free(temperature);
        }
        codes_handle_delete(h);
    }
    checked(err, "read climate GRIB");
    fclose(climate);
    free(mask);
    free(sst);
    printf("SUMMARY climate_frames=%zu invalid_ocean_values=%zu input_unchanged=yes\n", frames, invalid_total);
    if (!frames || invalid_initial) return 2;
    return invalid_total ? 1 : 0;
}
