LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE := Hough
LOCAL_SRC_FILES = hough.cpp
LOCAL_CPPFLAGS    := -std=c++11 -D_REENTRANT -fexceptions -frtti -w -Wall

LOCAL_C_INCLUDES += $(LOCAL_PATH) \
                    $(LOCAL_PATH)/../eigen/eigen3/build/include/

LOCAL_EXPORT_C_INCLUDES += $(LOCAL_PATH)
LOCAL_EXPORT_LDLIBS := -lz

LOCAL_STATIC_LIBRARIES := Eigen

include $(BUILD_STATIC_LIBRARY)