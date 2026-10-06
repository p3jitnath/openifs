#include "sst-grib-utils.h"

static void same_section(codes_handle *a,codes_handle *b,int section,const unsigned char *am,
                         const unsigned char *bm,size_t an,size_t bn) {
    char offset_key[40],length_key[40];
    snprintf(offset_key,sizeof(offset_key),"offsetSection%d",section);
    snprintf(length_key,sizeof(length_key),"section%dLength",section);
    long ao=sst_long(a,offset_key),bo=sst_long(b,offset_key);
    long al=sst_long(a,length_key),bl=sst_long(b,length_key);
    if (ao<0 || bo<0 || al<0 || bl<0 || al!=bl || (size_t)ao+(size_t)al>an ||
        (size_t)bo+(size_t)bl>bn || memcmp(am+ao,bm+bo,(size_t)al)) {
        fprintf(stderr,"FAIL: GRIB metadata/bitmap section %d changed\n",section); exit(2);
    }
}
int main(int argc,char **argv) {
    if (argc!=6) {
        fprintf(stderr,"Usage: %s SURFACE_INIT ORIGINAL_CLIMATE FROZEN_CLIMATE FIRST_MASK_DATE EXPECTED_FRAMES\n",argv[0]); return 2;
    }
    sst_source source=sst_load(argv[1],strtol(argv[4],NULL,10));
    long expected=strtol(argv[5],NULL,10);
    FILE *a=fopen(argv[2],"rb"),*b=fopen(argv[3],"rb");
    if (!a || !b) { perror("open climate"); return 2; }
    int ae=0,be=0; size_t records=0,frames=0,untouched=0; codes_handle *ah,*bh;
    long previous_date=0;
    while ((ah=codes_handle_new_from_file(NULL,a,PRODUCT_GRIB,&ae))) {
        bh=codes_handle_new_from_file(NULL,b,PRODUCT_GRIB,&be);
        if (!bh) { fprintf(stderr,"FAIL: truncated replacement\n"); return 2; }
        const void *am,*bm; size_t an,bn;
        sst_check(codes_get_message(ah,&am,&an),"original message");
        sst_check(codes_get_message(bh,&bm,&bn),"replacement message");
        if (sst_long(ah,"paramId")!=139) {
            if (an!=bn || memcmp(am,bm,an)) { fprintf(stderr,"FAIL: non-SST record %zu changed\n",records+1); return 2; }
            untouched++;
        } else {
            if (sst_long(ah,"edition")!=1 || sst_long(bh,"edition")!=1) return 2;
            /* PDS=date/parameter/units/levels; GDS=all grid geometry; bitmap=land mask. */
            for (int s=1;s<=3;s++) same_section(ah,bh,s,am,bm,an,bn);
            char geometry[128]; sst_grid(bh,geometry);
            if (strcmp(geometry,source.geometry)) return 2;
            size_t oldn,newn; double *old=sst_values(ah,&oldn),*v=sst_values(bh,&newn);
            if (oldn!=source.count || newn!=oldn) return 2;
            size_t ocean=0,land=0;
            for (size_t i=0;i<newn;i++) {
                if (source.mask[i]<=0.5) {
                    if (!isfinite(v[i]) || v[i]!=source.sst[i]) {
                        fprintf(stderr,"FAIL: SST is not exactly initial at point %zu\n",i); return 2;
                    }
                    ocean++;
                } else {
                    if (v[i]!=old[i]) { fprintf(stderr,"FAIL: land value changed at %zu\n",i); return 2; }
                    land++;
                }
            }
            long date=sst_long(bh,"dataDate");
            if (date<=previous_date || date<source.date || sst_long(bh,"dataTime")!=1200) {
                fprintf(stderr,"FAIL: unexpected forcing dates or times\n"); return 2;
            }
            previous_date=date;
            printf("Verified date=%ld time=1200 exact_initial_ocean_values=%zu unchanged_land_values=%zu metadata_and_bitmap_identical=yes\n",date,ocean,land);
            free(old); free(v); frames++;
        }
        codes_handle_delete(ah); codes_handle_delete(bh); records++;
    }
    sst_check(ae,"original EOF");
    bh=codes_handle_new_from_file(NULL,b,PRODUCT_GRIB,&be);
    sst_check(be,"replacement EOF");
    if (bh) { fprintf(stderr,"FAIL: extra replacement record\n"); return 2; }
    fclose(a); fclose(b); free(source.mask); free(source.sst);
    if (expected<1 || frames!=(size_t)expected) { fprintf(stderr,"FAIL: wrong SST frame count\n"); return 2; }
    printf("PASS: records=%zu SST_frames=%zu non_SST_records_byte_identical=%zu; SST constant in time, spatial pattern preserved\n",records,frames,untouched);
    return 0;
}
