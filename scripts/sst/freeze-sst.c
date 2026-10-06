#define _POSIX_C_SOURCE 200809L
#include "sst-grib-utils.h"
#include <fcntl.h>
#include <unistd.h>

int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr,"Usage: %s SURFACE_INIT ORIGINAL_CLIMATE NEW_CLIMATE FIRST_MASK_DATE EXPECTED_FRAMES\n",argv[0]);
        return 2;
    }
    sst_source source=sst_load(argv[1],strtol(argv[4],NULL,10));
    long expected=strtol(argv[5],NULL,10);
    if (expected<1) return 2;
    FILE *in=fopen(argv[2],"rb");
    if (!in) { perror(argv[2]); return 2; }
    int fd=open(argv[3],O_WRONLY|O_CREAT|O_EXCL,0600);
    if (fd<0) { perror(argv[3]); return 2; }
    FILE *out=fdopen(fd,"wb");
    if (!out) { perror("fdopen"); return 2; }
    int err=0; codes_handle *h; size_t frames=0,records=0;
    while ((h=codes_handle_new_from_file(NULL,in,PRODUCT_GRIB,&err))) {
        if (sst_long(h,"paramId")==139) {
            char geometry[128]; sst_grid(h,geometry);
            size_t n; double *v=sst_values(h,&n);
            if (n!=source.count || strcmp(geometry,source.geometry)) { fprintf(stderr,"FAIL: climate grid mismatch\n"); return 2; }
            double missing; sst_check(codes_get_double(h,"missingValue",&missing),"missing value");
            for (size_t i=0;i<n;i++) {
                /* This package's climate bitmap is exactly the static ocean mask. */
                if ((source.mask[i]<=0.5) == (v[i]==missing)) {
                    fprintf(stderr,"FAIL: original SST bitmap/model mask mismatch at %zu\n",i); return 2;
                }
                if (source.mask[i]<=0.5) v[i]=source.sst[i];
            }
            sst_check(codes_set_double_array(h,"values",v,n),"encode fixed initial SST");
            size_t decoded; double *check=sst_values(h,&decoded);
            if (decoded!=n) return 2;
            double max_error=0;
            for (size_t i=0;i<n;i++) {
                double delta=fabs(check[i]-v[i]);
                if (delta>max_error) max_error=delta;
                if (check[i]!=v[i]) {
                    fprintf(stderr,"FAIL: repacking changed an intended value at %zu by %.10g K\n",i,delta); return 2;
                }
            }
            printf("Fixed forcing date=%ld time=%ld points=%zu exact_initial_SST=yes land_unchanged=yes max_repack_error=%.10g\n",
                   sst_long(h,"dataDate"),sst_long(h,"dataTime"),n,max_error);
            free(v); free(check); frames++;
        }
        const void *message; size_t length;
        sst_check(codes_get_message(h,&message,&length),"GRIB message");
        if (fwrite(message,1,length,out)!=length) { perror("write climate"); return 2; }
        codes_handle_delete(h); records++;
    }
    sst_check(err,"read original climate"); fclose(in);
    if (frames!=(size_t)expected) { fprintf(stderr,"FAIL: expected %ld SST frames, got %zu\n",expected,frames); return 2; }
    if (fflush(out) || fsync(fd) || fclose(out)) { perror("flush climate"); return 2; }
    free(source.mask); free(source.sst);
    printf("PASS: records=%zu SST_frames=%zu; original input untouched; staged output requires independent validation\n",records,frames);
    return 0;
}
