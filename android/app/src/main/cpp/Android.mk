LOCAL_PATH:= $(call my-dir)
export MAINDIR:= $(LOCAL_PATH)

include $(CLEAR_VARS)

#include $(MAINDIR)/clapack/Android.mk
include $(MAINDIR)/lrslib/Android.mk
#include $(MAINDIR)/jfftw/Android.mk

LOCAL_PATH := $(MAINDIR)

# include $(CLEAR_VARS)
# LOCAL_MODULE:= lapack
# # LOCAL_STATIC_LIBRARIES := tmglib clapack blas f2c
# LOCAL_STATIC_LIBRARIES := clapack blas f2c
# LOCAL_EXPORT_C_INCLUDES := $(LOCAL_C_INCLUDES)
# LOCAL_EXPORT_LDLIBS := $(LOCAL_LDLIBS)
# include $(BUILD_STATIC_LIBRARY)
#
# include $(CLEAR_VARS)
# LOCAL_MODULE:= testlapack
# LOCAL_SRC_FILES:= testclapack.cpp
# LOCAL_STATIC_LIBRARIES := lapack
# include $(BUILD_SHARED_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE := roomshape
LOCAL_SRC_FILES:= roomshape.cpp ext_funcs.cpp prune_echoes.cpp trilaterate_back.cpp
LOCAL_STATIC_LIBRARIES := lrs
LOCAL_C_INCLUDES := $(MAINDIR)/Eigen
include $(BUILD_SHARED_LIBRARY)
