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
#define CLAMP(x, min, max) x<min?min:(x>max?max:x)

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

    void Hough::recomputeBox(pcl::PointCloud<Eigen::Vector2f>::Ptr cloud) {

        box b = _bornes;

        for (auto pt : cloud->points) {
            b.minx = std::min(b.minx, pt.x());
            b.maxx = std::max(b.maxx, pt.x());
            b.miny = std::min(b.miny, pt.y());
            b.maxy = std::max(b.maxy, pt.y());
        }

        if (b != _bornes)
            updateBox(b);
    }

    void Hough::updateBox(box newbox) {
		int new_img_w = newbox.getW();
		int new_img_h = newbox.getH();

		//Create the accu
		int new_hough_h = sqrt(2.0) * std::max(new_img_w, new_img_h) / 2.0 + 1;
		int new_accu_h = new_hough_h * 2.0 + 1; // -r -> +r
		int new_accu_w = 180;

		std::vector<Eigen::Vector2f>** new_accu = new std::vector<Eigen::Vector2f>*[new_accu_h];
		for(int r=0;r<new_accu_h;r++) {
		    new_accu[r] = new std::vector<Eigen::Vector2f>[new_accu_w];
			for(int t=0;t<new_accu_w;t++) {
                std::vector<Eigen::Vector2f> vect;
                new_accu[r][t] = vect;
            }
        }

		if(_accu) {
            // copy old _accu
            int dw = (new_accu_w - _accu_w)/2;
            int dh = (new_accu_h - _accu_h)/2;

            for(int r=0;r<_accu_h;r++) {
                for(int t=0;t<_accu_w;t++) {
                    auto cell = GetAccuCell(r, t);
                    new_accu[r+dh][t+dw].insert(new_accu[r+dh][t+dw].end(), cell->begin(), cell->end());

                }
            }

            for(int r=0;r<_accu_h;r++) {
                delete[] _accu[r];
			}
            delete[] _accu;
        }
        _accu = new_accu;
        _bornes = newbox;

		_img_w = new_img_w;
		_img_h = new_img_h;

		_hough_h = new_hough_h;
		_accu_h = new_accu_h;
		_accu_w = new_accu_w;
    }

    void Hough::AddPointCloud(pcl::PointCloud<Eigen::Vector2f>::Ptr cloud) {

        if (cloud->points.size() == 0) return;

        recomputeBox(cloud);

		if(_accu_h == 0) return;

		double center_x = _img_w/2;
		double center_y = _img_h/2;

        for (auto pt : cloud->points) {
            auto proj = _bornes.project(pt);
            for(int t=0;t<180;t++) {
                double r = (proj(0) - center_x) * std::cos(t * DEG2RAD) + (proj(1) - center_y) * std::sin(t * DEG2RAD);
                auto cell = GetAccuCell(CLAMP((int)(round(r + _hough_h)), 0, _accu_h-1), t);
                cell->push_back(pt);
            }
        }
	}

    std::vector<Eigen::Vector2f>* Hough::GetAccuCell(int h, int w) {
        return &_accu[h][w];
    }

    void Hough::RemovePointCloud(pcl::PointCloud<Eigen::Vector2f>::Ptr cloud) {

        if (cloud->points.size() == 0) return;

		if(_accu_h == 0) return;

		double center_x = _img_w/2;
		double center_y = _img_h/2;

        for (auto pt : cloud->points) {
            auto proj = _bornes.project(pt);
            for(int t=0;t<180;t++) {
                double r = (proj(0) - center_x) * std::cos(t * DEG2RAD) + (proj(1) - center_y) * std::sin(t * DEG2RAD);
                auto cell = GetAccuCell(CLAMP((int)(round(r + _hough_h)), 0, _accu_h-1), t);
                auto it = std::find(cell->begin(), cell->end(), pt);
                if (it != cell->end()) {
                  std::swap(*it, cell->back());
                  cell->pop_back();
                }
            }
        }
	}

    std::vector<std::pair<Eigen::Vector2f, Eigen::Vector2f>> Hough::GetSegments(int threshold, float dist_between_parts) {
	    std::vector<std::pair<Eigen::Vector2f, Eigen::Vector2f>> segments;
	    for (auto eye : GetAccuEyes(threshold)) {
	        auto cell = GetAccuCell(eye.first, eye.second);

            if (cell->size()== 0) continue;

            auto line = ToParametrizedLine(eye.first, eye.second);

            auto sort = [&line](const Eigen::Vector2f &a, const Eigen::Vector2f &b) {
                Eigen::Vector2f p1 = line.projection(a);
                Eigen::Vector2f p2 = line.projection(b);
                return (p1 - line.origin()).norm() < (p2 - line.origin()).norm();
            };

            if (not std::is_sorted(cell->begin(), cell->end(), sort))
                std::sort(cell->begin(), cell->end(), sort);

            Eigen::Vector2f start = (*cell)[0];
            for (int i = 0; i < cell->size()-1; ++i ) {
                if (((*cell)[i+1]-(*cell)[i]).norm() > dist_between_parts) {
                    segments.push_back(std::make_pair(start,(*cell)[i]));
                    start = (*cell)[i+1];
                }
            }
            segments.push_back(std::make_pair(start,(*cell)[cell->size()-1])); // add last segment
	    }
	    return segments;
    }

	std::vector<Eigen::ParametrizedLine<float, 2>> Hough::GetLines(int threshold) {
	    std::vector<Eigen::ParametrizedLine<float, 2>> lines;
	    for (auto eye : GetAccuEyes(threshold)) {
	        lines.push_back(ToParametrizedLine(eye.first, eye.second));
	    }
	    return lines;
	}

	Eigen::ParametrizedLine<float, 2> Hough::ToParametrizedLine(int r, int t) {

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
        Eigen::ParametrizedLine<float, 2> line = Eigen::ParametrizedLine<float, 2>::Through(p1,p2);
        return line;

	}

	std::vector<std::pair<int,int>> Hough::GetAccuEyes(int threshold) {
	    std::vector<std::pair<int,int>> eyes;

		if(_accu == 0)
			return eyes;

		for(int r=0;r<_accu_h;r++) {
			for(int t=0;t<_accu_w;t++) {
				if(GetAccuCell(r, t)->size() >= threshold) {
					//Is this point a local maxima (9x9)
					int loc = 15 / 2;
					int max = GetAccuCell(r, t)->size();
					for(int ly=-loc;ly<=loc;ly++) {
						for(int lx=-loc;lx<=loc;lx++) {
							if( (ly+r>=0 && ly+r<_accu_h) && (lx+t>=0 && lx+t<_accu_w)  ) {
								if( GetAccuCell(r+ly, t+lx)->size() > max ) {
									max = GetAccuCell(r+ly, t+lx)->size();
									ly = lx = 5;
								}
							}
						}
					}
					if(max > GetAccuCell(r, t)->size())
						continue;

					eyes.push_back(std::make_pair(r, t));

				}
			}
		}
		return eyes;
	}

	std::vector<Eigen::Vector2f>** Hough::GetAccu(int *w, int *h) {
		*w = _accu_w;
		*h = _accu_h;

		return _accu;
	}
}
