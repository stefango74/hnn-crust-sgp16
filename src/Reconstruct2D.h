/*
 * Reconstruct2D.h
 *
 *  Created on: Jun 25, 2015
 *      Author: Stefan Ohrhallinger
 */

#ifndef RECONSTRUCT2D_H_
#define RECONSTRUCT2D_H_

#include <set>
//#include <list>
#include <stack>
#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>
//#include <CGAL/Orthogonal_incremental_neighbor_search.h>
//#include <CGAL/Search_traits.h>
#include <ANN/ANN.h>
#include "Point.h"
#include "circle.h"

#define PI 3.1415926

#ifdef USE_CGAL
struct IPoint
{
  Point p;
  int i = -1;

  IPoint() { p[0] = p[1] = 0.0; }
  IPoint (Point p_P, int index) { p = p_P; i = index; }

  Point point() const { return p; }
  int index() const { return i; }

  int& index() { return i; }

  bool operator==(const IPoint& ip) const
  {
    return (ip.p == p) && (ip.index() == index());
  }

  bool  operator!=(const IPoint& p) const { return ! (*this == p); }
}; //end of class

namespace CGAL
{
  template <>
  struct Kernel_traits<IPoint> {
    struct Kernel {
      typedef double FT;
      typedef double RT;
    };
  };
}

struct IDistance {
  typedef IPoint Query_item;

  double transformed_distance(const IPoint& p1, const IPoint& p2) const {
    double distx= p1.p[0]-p2.p[0];
    double disty= p1.p[1]-p2.p[1];
    return distx*distx+disty*disty;
  }

  template <class TreeTraits>
  double min_distance_to_rectangle(const IPoint& p,
				   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const {
    double distance(0.0), h = p.p[0];
    if (h < b.min_coord(0)) distance += (b.min_coord(0)-h)*(b.min_coord(0)-h);
    if (h > b.max_coord(0)) distance += (h-b.max_coord(0))*(h-b.max_coord(0));
    h=p.p[1];
    if (h < b.min_coord(1)) distance += (b.min_coord(1)-h)*(b.min_coord(1)-h);
    if (h > b.max_coord(1)) distance += (h-b.max_coord(1))*(h-b.min_coord(1));
    return distance;
  }

  template <class TreeTraits>
  double max_distance_to_rectangle(const IPoint& p,
				   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const {
    double h = p.p[0];

    double d0 = (h >= (b.min_coord(0)+b.max_coord(0))/2.0) ?
                (h-b.min_coord(0))*(h-b.min_coord(0)) : (b.max_coord(0)-h)*(b.max_coord(0)-h);

    h=p.p[1];
    double d1 = (h >= (b.min_coord(1)+b.max_coord(1))/2.0) ?
                (h-b.min_coord(1))*(h-b.min_coord(1)) : (b.max_coord(1)-h)*(b.max_coord(1)-h);
    return d0 + d1;
  }

  double new_distance(double& dist, double old_off, double new_off,
		      int /* cutting_dimension */)  const {
    return dist + new_off*new_off - old_off*old_off;
  }

  double transformed_distance(double d) const { return d*d; }

  double inverse_of_transformed_distance(double d) { return std::sqrt(d); }

}; // end of struct Distance

struct Construct_coord_iterator {
  typedef  const double* result_type;

  const Point operator()(const IPoint& p) const
  { return static_cast<const Point>(p.p); }

//  const Point operator()(const IPoint& p, int)  const
//  { return static_cast<const Point>(p.p+3); }
};
#endif

using namespace std;

/*
typedef CGAL::Search_traits<double, IPoint, const Point, Construct_coord_iterator> Traits;
typedef CGAL::Orthogonal_incremental_neighbor_search<Traits, IDistance> NN_incremental_search;
typedef NN_incremental_search::iterator NN_incremental_iterator;
typedef NN_incremental_search::Tree Tree;
*/

typedef enum { CONFORMING, NONCONFORMING, NOISY, OUTLIER } PointsEnum;
typedef enum { UNIJECTIVE, BIJECTIVE } EdgeEnum;

class Reconstruct2D {
private:
	bool isClosed;
	vector<Point> points, projPoints, normals;
	vector<PointsEnum> pClasses;
	map<pair<int, int>, EdgeEnum> visEdgeMap;
	vector<Circle> circles;
	vector<pair<int, int> > arcs;
	void determineNeighbors(ANNkd_tree *kdTree, ANNpointArray ann_points, int kMax,
			vector<pair<int, int> > &manifoldNeighbors);

public:
	Reconstruct2D(const vector<Point> &p_points, bool isClosed);
	void reconstruct();
	void reconstructNoisy();
	map<pair<int, int>, EdgeEnum> getEdgeMap();
	vector<Point> getProjectedPoints();
	vector<Circle> getCircles();
	vector<pair<int, int> > getArcs();
	vector<PointsEnum> getPointClassification();
	vector<Point> getNormals();
};

typedef struct Node *NodePtr;

struct Node
{
	int index;
	vector<pair<NodePtr, int> > neighbors;
};

#endif /* RECONSTRUCT2D_H_ */
