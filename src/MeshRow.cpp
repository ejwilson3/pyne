#include "discretize.h"
#include <stdlib.h>
#include <math.h>
#include <stdexcept>

namespace pyne {


MeshRow::MeshRow(int num_rays, bool grid) {
  if (grid && sqrt(num_rays*num_rays) != num_rays) {
    throw std::runtime_error("Must have a perfect square number of rays for"
                              "a grid.");
  }
  this->num_rays = num_rays;
  this->grid = grid;
  sums.at(0).vol_frac = 0;
  sums.at(0).rel_error = 0;
}

//MeshRow::~MeshRow() {}

void MeshRow::setDimension1(double div1, double div2) {
  d1div1 = div1;
  d1div2 = div2;
}

void MeshRow::setDimension2(double div1, double div2) {
  d2div1 = div1;
  d2div2 = div2;
}

void MeshRow::fireRays(int d3, std::vector<double> divs) {
  d3divs = divs;
  //What is the vol?
  //WATCH THIS
  //vol vol;

  std::vector<double> width;
  for (int i = 0; i < d3divs.size() - 1; i++) {
    width.at(i) = d3divs.at(i+1) - d3divs.at(i);
  }

  for (int i = 0; i < num_rays; i++) {
    int x, y;
    double a, b, c, d;
    startPoints(i);
    //TODO
    //if (!point_in_volume(vol, point, direction) {
    //  vol = find_volume(?);
    //}
    double value;

    //??????????????????????
    //WATCH THIS
    //Placeholder; distance should come from some dagmc thing.
    double distance = (rand() % 300)/100;

    for (int j = 0; j < width.size(); j++) {
      do {
        if (distance >= width.at(j)) {
          value = 1;
          distance -= width.at(j);
        }
        else {
          value = distance/width.at(j);
        }
        //why did I add this?
        //double[2] sample;
        //WATCH THIS
        sums.at(0).vol_frac += value;
        sums.at(0).rel_error += value*value;
      //figure out how to iterate to the next distance...
      } while (false);
    }
  }
}

//WATCH THIS
std::vector<disc_result> MeshRow::getSums() {
  return sums;
}

void MeshRow::startPoints(int iter) {
  double dist_x = d1div2 - d1div1;
  double dist_y = d2div2 - d2div1;
  if (grid) {
    int points = sqrt(num_rays);
    //Iterate from 1 to points, repeat points times.
    int iter_x = iter%points + 1;
    //Iterate slowly from 1 to points once.
    int iter_y = iter/points + 1;

    start_point_x = d1div1 + dist_x/(points + 1)*iter_x;
    start_point_y = d1div1 + dist_y/(points + 1)*iter_y;
  }
  else {
    start_point_x = d1div1 + dist_x*rand();
    start_point_y = d2div1 + dist_y*rand();
  }
}

} //namespace pyne
