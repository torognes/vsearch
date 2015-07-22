TESTS += \
./tests/check_vsearch.o \
./tests/test_align_simd_nuc.o \
./tests/test_align_simd_aa.o \
./tests/test_db.o \
./tests/test_maps.o \
./tests/test_dprofile_fill_nuc.o

TEST_DEPS += \
./tests/tests.h