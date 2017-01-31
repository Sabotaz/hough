// **********************************************************************************
//
// BSD License.
// This file is part of a Hough Transformation tutorial
// see: http://www.keymolen.com/2013/05/hough-transformation-c-implementation.html
//
// Copyright (c) 2013, Bruno Keymolen, email: bruno.keymolen@gmail.com
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// Redistributions in binary form must reproduce the above copyright notice, this
// list of conditions and the following disclaimer in the documentation and/or other
// materials provided with the distribution.
// Neither the name of "Bruno Keymolen" nor the names of its contributors may be
// used to endorse or promote products derived from this software without specific
// prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// **********************************************************************************

#ifndef HOUGH_H_
#define HOUGH_H_

#include <vector>
#include <Eigen/Core>

#include <android/log.h>

#define LOG_TAG "keymolen"
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO, LOG_TAG, __VA_ARGS__)
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, __VA_ARGS__)
#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, __VA_ARGS__)

#define HOUGH_UNIT_ 100

namespace keymolen {

    struct box {
        float minx;
        float maxx;
        float miny;
        float maxy;

	    inline Eigen::Vector2f project(Eigen::Vector2f pt) {
	        Eigen::Vector2f p = {minx, miny};
	        return (pt - p)*HOUGH_UNIT_;
        }

	    inline Eigen::Vector2f unproject(Eigen::Vector2f pt) {
	        Eigen::Vector2f p = {minx, miny};
	        return pt/HOUGH_UNIT_ + p;
        }

        inline float getW() { return (maxx - minx)* HOUGH_UNIT_; }
        inline float getH() { return (maxy - miny)* HOUGH_UNIT_; }

    } typedef box;

	class Hough {
	public:
		Hough();
		virtual ~Hough();
	public:
		void Transform(std::vector<Eigen::Vector2f> points, box bornes);
		std::vector< std::pair< std::pair< Eigen::Vector2f, Eigen::Vector2f >, std::vector<Eigen::Vector2f> > > GetLines(int threshold);
		std::vector<Eigen::Vector2f>** GetAccu(int *w, int *h);
	private:
		std::vector<Eigen::Vector2f>** _accu;
		int _accu_w;
		int _accu_h;
		int _img_w;
		int _img_h;
		box _bornes;
	};

}

#endif /* HOUGH_H_ */
