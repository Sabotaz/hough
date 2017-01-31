// **********************************************************************************
//
// BSD License.
// This file is part of a Hough Transformation tutorial,
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

#include "hough.h"
#include <cmath>
#include <iostream>
#include <string.h>
#include <stdlib.h>

#define DEG2RAD 0.017453293f

namespace keymolen {

	Hough::Hough():_accu(0), _accu_w(0), _accu_h(0), _img_w(0), _img_h(0) {

	}

	Hough::~Hough() {
		if(_accu) {
             for(int r=0;r<_accu_h;r++) {
                 delete[] _accu[r];
            }
             delete[] _accu;
         }
	}

	void Hough::Transform(std::vector<Eigen::Vector2f> points, box b) {

		_bornes = b;

		_img_w = _bornes.getW();
		_img_h = _bornes.getH();

		//Create the accu
		double hough_h = sqrt(2.0) * std::max(_img_w, _img_h) / 2.0;
		_accu_h = hough_h * 2.0 + 1; // -r -> +r
		_accu_w = 180;

		if(_accu_h == 0) return;

		if(_accu) {
            for(int r=0;r<_accu_h;r++) {
                delete[] _accu[r];
			}
            delete[] _accu;
        }
		//_accu = (std::vector<Eigen::Vector2f>*)std::malloc(_accu_h * _accu_w * sizeof(std::vector<Eigen::Vector2f>));
		_accu = new std::vector<Eigen::Vector2f>*[_accu_h];
		for(int r=0;r<_accu_h;r++) {
		    _accu[r] = new std::vector<Eigen::Vector2f>[_accu_w];
			for(int t=0;t<_accu_w;t++) {
                std::vector<Eigen::Vector2f> vect;
                _accu[r][t] = vect;
            }
        }

		double center_x = _img_w/2;
		double center_y = _img_h/2;

        for (auto pt : points) {
            auto proj = _bornes.project(pt);
            for(int t=0;t<180;t++) {
                double r = (proj(0) - center_x) * std::cos(t * DEG2RAD) + (proj(1) - center_y) * std::sin(t * DEG2RAD);
                _accu[(int)(round(r + hough_h))][t].push_back(pt);
            }
        }
	}

	std::vector< std::pair< std::pair< Eigen::Vector2f, Eigen::Vector2f >, std::vector<Eigen::Vector2f> > > Hough::GetLines(int threshold) {
	    std::vector< std::pair< std::pair< Eigen::Vector2f, Eigen::Vector2f >, std::vector<Eigen::Vector2f> > > lines;

		if(_accu == 0)
			return lines;

		for(int r=0;r<_accu_h;r++) {
			for(int t=0;t<_accu_w;t++) {
				if(_accu[r][t].size() >= threshold) {
					//Is this point a local maxima (9x9)
					int loc = 15 / 2;
					int max = _accu[r][t].size();
					for(int ly=-loc;ly<=loc;ly++) {
						for(int lx=-loc;lx<=loc;lx++) {
							if( (ly+r>=0 && ly+r<_accu_h) && (lx+t>=0 && lx+t<_accu_w)  ) {
								if( _accu[r+ly][t+lx].size() > max ) {
									max = _accu[r+ly][t+lx].size();
									ly = lx = 5;
								}
							}
						}
					}
					if(max > _accu[r][t].size())
						continue;

					float x1, y1, x2, y2;
					x1 = y1 = x2 = y2 = 0;

					if(t >= 45 && t <= 135) {
						//y = (r - x cos(t)) / sin(t)
						x1 = 0;
						y1 = ((r-(_accu_h/2)) - ((x1 - (_img_w/2) ) * std::cos(t * DEG2RAD))) / std::sin(t * DEG2RAD) + (_img_h / 2);
						x2 = _img_w - 0;
						y2 = ((r-(_accu_h/2)) - ((x2 - (_img_w/2) ) * std::cos(t * DEG2RAD))) / std::sin(t * DEG2RAD) + (_img_h / 2);
					} else {
						//x = (r - y sin(t)) / cos(t);
						y1 = 0;
						x1 = ((r-(_accu_h/2)) - ((y1 - (_img_h/2) ) * std::sin(t * DEG2RAD))) / std::cos(t * DEG2RAD) + (_img_w / 2);
						y2 = _img_h - 0;
						x2 = ((r-(_accu_h/2)) - ((y2 - (_img_h/2) ) * std::sin(t * DEG2RAD))) / std::cos(t * DEG2RAD) + (_img_w / 2);
					}
                    Eigen::Vector2f p1 = _bornes.unproject({x1, y1});
                    Eigen::Vector2f p2 = _bornes.unproject({x2, y2});
                    std::pair< Eigen::Vector2f, Eigen::Vector2f > line = std::make_pair(p1, p2);
					lines.push_back(std::make_pair(line, _accu[r][t]));

				}
			}
		}

		LOGI("lines: %d, threshold: %d", lines.size(), threshold);
		return lines;
	}

	std::vector<Eigen::Vector2f>** Hough::GetAccu(int *w, int *h) {
		*w = _accu_w;
		*h = _accu_h;

		return _accu;
	}
}
