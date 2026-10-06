#ifndef SST_GRIB_UTILS_H
#define SST_GRIB_UTILS_H
#include <eccodes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

static inline void sst_check(int err, const char *operation) {
    if (err) { fprintf(stderr,"FAIL: %s: %s\n",operation,codes_get_error_message(err)); exit(2); }
}
static inline long sst_long(codes_handle *h, const char *key) {
    long v; sst_check(codes_get_long(h,key,&v),key); return v;
}
static inline double *sst_values(codes_handle *h, size_t *count) {
    sst_check(codes_get_size(h,"values",count),"values size");
    double *v=malloc(*count*sizeof(*v));
    if (!v) { perror("malloc"); exit(2); }
    size_t n=*count;
    sst_check(codes_get_double_array(h,"values",v,&n),"decode values");
    if (n != *count) { fprintf(stderr,"FAIL: incomplete decode\n"); exit(2); }
    return v;
}
static inline void sst_grid(codes_handle *h, char result[128]) {
    char type[128]; size_t length=sizeof(type);
    sst_check(codes_get_string(h,"gridType",type,&length),"grid type");
    if (strcmp(type,"reduced_gg")) { fprintf(stderr,"FAIL: expected reduced Gaussian grid\n"); exit(2); }
    long east=sst_long(h,"longitudeOfLastGridPoint");
    if (east < 359800 || east >= 360000) { fprintf(stderr,"FAIL: expected global O1280 domain\n"); exit(2); }
    const char *keys[]={"N","Nj","numberOfDataPoints","scanningMode","latitudeOfFirstGridPoint",
                        "latitudeOfLastGridPoint","longitudeOfFirstGridPoint"};
    uint64_t hash=UINT64_C(14695981039346656037);
    for (size_t k=0;k<sizeof(keys)/sizeof(keys[0]);k++) {
        uint64_t v=(uint64_t)sst_long(h,keys[k]);
        for (int b=0;b<8;b++) { hash=(hash^(v&255))*UINT64_C(1099511628211); v>>=8; }
    }
    size_t rows; sst_check(codes_get_size(h,"pl",&rows),"row count");
    long *pl=malloc(rows*sizeof(*pl));
    if (!pl) { perror("malloc"); exit(2); }
    sst_check(codes_get_long_array(h,"pl",pl,&rows),"row lengths");
    for (size_t k=0;k<rows;k++) {
        uint64_t v=(uint64_t)pl[k];
        for (int b=0;b<8;b++) { hash=(hash^(v&255))*UINT64_C(1099511628211); v>>=8; }
    }
    free(pl);
    snprintf(result,128,"global_reduced_gg_%016" PRIx64,hash);
}
typedef struct {
    double *mask, *sst;
    size_t count;
    long date, time, mask_date;
    char geometry[128];
} sst_source;
static inline sst_source sst_load(const char *path, long expected_mask_date) {
    sst_source s={0};
    FILE *f=fopen(path,"rb");
    if (!f) { perror(path); exit(2); }
    size_t nm=0,ns=0; char mg[128]; int err=0; codes_handle *h;
    while ((h=codes_handle_new_from_file(NULL,f,PRODUCT_GRIB,&err))) {
        long param=sst_long(h,"paramId");
        if (param==172 && !s.mask) {
            /* SUGRIDG / GRID_IN uses the first matching field (LREADALL=false). */
            s.mask_date=sst_long(h,"dataDate");
            if (s.mask_date != expected_mask_date) { fprintf(stderr,"FAIL: first model mask date differs\n"); exit(2); }
            s.mask=sst_values(h,&nm); sst_grid(h,mg);
        } else if (param==34) {
            if (s.sst) { fprintf(stderr,"FAIL: duplicate initial SST\n"); exit(2); }
            s.sst=sst_values(h,&ns); sst_grid(h,s.geometry);
            s.date=sst_long(h,"dataDate"); s.time=sst_long(h,"dataTime");
        }
        codes_handle_delete(h);
    }
    sst_check(err,"read initial file"); fclose(f);
    if (!s.mask || !s.sst || nm!=ns || strcmp(mg,s.geometry)) {
        fprintf(stderr,"FAIL: source SST/mask geometry mismatch\n"); exit(2);
    }
    s.count=ns;
    size_t sea=0;
    for (size_t i=0;i<ns;i++) {
        if (!isfinite(s.mask[i]) || s.mask[i]<0 || s.mask[i]>1) { fprintf(stderr,"FAIL: invalid mask\n"); exit(2); }
        if (s.mask[i]<=0.5) {
            sea++;
            if (!isfinite(s.sst[i]) || s.sst[i]<250 || s.sst[i]>330) { fprintf(stderr,"FAIL: invalid initial ocean SST\n"); exit(2); }
        }
    }
    if (!sea) { fprintf(stderr,"FAIL: no ocean points\n"); exit(2); }
    printf("Source date=%ld time=%ld first_model_mask_date=%ld ocean_points=%zu geometry=%s\n",
           s.date,s.time,s.mask_date,sea,s.geometry);
    return s;
}
#endif
