#include "cavc/polyline.hpp"

using namespace cavc;

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  // closed polyline representing a circle going counter clockwise
  double radius = 10.0;
  Polyline<double> ccwCircle;
  ccwCircle.addVertex(0, 0, 1);
  ccwCircle.addVertex(2.0 * radius, 0, 1);
  ccwCircle.isClosed() = true;

  // closed polyline representing a circle going clockwise
  Polyline<double> cwCircle;
  cwCircle.addVertex(0, 0, -1);
  cwCircle.addVertex(2.0 * radius, 0, -1);
  cwCircle.isClosed() = true;

  // path length of polyline
  double length = getPathLength(ccwCircle);
  assert(utils::fuzzy::equal(length, 2.0 * utils::pi<double>() * radius));

  // signed area of closed polyline (area is positive if counter clockwise, negative if clockwise)
  double area = getArea(ccwCircle);
  assert(utils::fuzzy::equal(area, utils::pi<double>() * radius * radius));
  assert(utils::fuzzy::equal(getArea(cwCircle), -area));

  // polyline extents
  AABB<double> extents = getExtents(ccwCircle);
  assert(utils::fuzzy::equal(extents.xMin, 0.0));
  assert(utils::fuzzy::equal(extents.yMin, -radius));
  assert(utils::fuzzy::equal(extents.xMax, 2.0 * radius));
  assert(utils::fuzzy::equal(extents.yMax, radius));

  // Closest point on polyline to a point given
  ClosestPoint<double> closestPoint(ccwCircle, Vector2<double>(radius, 10.0 * radius));
  // index is the starting vertex index of the closest segment (going clockwise so arc starting at
  // the second vertex is closest)
  assert(closestPoint.index() == 1);
  assert(fuzzy::equal(closestPoint.point(), Vector2<double>(radius, radius)));
  assert(utils::fuzzy::equal(closestPoint.distance(), 9.0 * radius));
}
