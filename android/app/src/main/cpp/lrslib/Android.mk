LOCAL_PATH:= $(call my-dir)

include $(CLEAR_VARS)
include $(LOCAL_PATH)/Config.mk

LOCAL_MODULE:= lrs
LOCAL_CFLAGS := -O3 -DLRSLONG
ifeq ($(TARGET_ARCH),arm)
	LOCAL_CFLAGS += -DB32
endif
ifeq ($(TARGET_ARCH),x86)
	LOCAL_CFLAGS += -DB32
endif
LOCAL_SRC_FILES:= $(ALLOBJ)  

include $(BUILD_STATIC_LIBRARY)
