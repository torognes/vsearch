TESTS += \
./tests/check_vsearch.o \
./tests/test_align_simd_nuc.o \
./tests/test_align_simd_aa.o \
./tests/test_db.o \
./tests/test_maps.o \
./tests/test_dprofile_fill_nuc.o \
./tests/test_matrices.o \
./tests/helper_functions.o

TEST_DEPS += \
./tests/tests.h \
./tests/helper_functions.h