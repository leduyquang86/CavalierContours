#ifndef CAVC_POLYLINEOFFSETISLANDS_HPP
#define CAVC_POLYLINEOFFSETISLANDS_HPP
#include "polyline.hpp"
#include "polylinecombine.hpp"
#include "polylineoffset.hpp"
#include <unordered_map>
#include <vector>

namespace cavc {
template <typename Real> struct OffsetLoop {
  std::size_t parentLoopIndex;
  Polyline<Real> polyline;
  StaticSpatialIndex<Real> spatialIndex;
};

template <typename Real> struct OffsetLoopSet {
  std::vector<OffsetLoop<Real>> ccwLoops;
  std::vector<OffsetLoop<Real>> cwLoops;
};

template <typename Real> class ParallelOffsetIslands {
public:
  ParallelOffsetIslands() {}
  OffsetLoopSet<Real> compute(OffsetLoopSet<Real> const &input, Real offsetDelta);

private:
  // type to represent a slice point (intersect) on an OffsetLoop
  struct LoopSlicePoint {
    // intersect between loops
    PlineIntersect<Real> intr;
    bool noSliceAfterForIndex1;
  };

  // type to represent a set of slice points (intersects) between two loops
  struct SlicePointSet {
    // index of first loop involved
    std::size_t loopIndex1;
    // index of second loop involved
    std::size_t loopIndex2;
    // all of the slice points (intersects) between the two loops
    std::vector<LoopSlicePoint> slicePoints;
  };

  // get an offset loop by index i, maps to ccw then cw offset loops
  OffsetLoop<Real> &getOffsetLoop(std::size_t i) {
    return i < _ccwOffsetLoops.size() ? _ccwOffsetLoops[i]
                                       : _cwOffsetLoops[i - _ccwOffsetLoops.size()];
  }

  OffsetLoop<Real> const &getParentLoop(std::size_t i) {
    return i < _inputSet->ccwLoops.size() ? _inputSet->ccwLoops[i]
                                           : _inputSet->cwLoops[i - _inputSet->ccwLoops.size()];
  }

  void createOffsetLoops(const OffsetLoopSet<Real> &input, Real absDelta);
  void createOffsetLoopsIndex();
  void createSlicePoints();

  struct DissectionPoint {
    std::size_t otherLoopIndex;
    Vector2<Real> pos;
  };

  struct DissectedSlice {
    // open polyline representing the slice
    Polyline<Real> pline;
    // index of the loop the slice is from
    std::size_t sliceParentIndex;
    // index of the loop that intersected the parent loop to form the start of the slice
    std::size_t startLoopIndex;
    // index of the loop that intersected the parent loop to form the end of the slice
    std::size_t endLoopIndex;
  };

  bool pointOnOffsetValid(std::size_t skipIndex, Vector2<Real> const &pt, Real absDelta);
  void createSlicesFromLoop(std::size_t loopIndex, Real absDelta,
                            std::vector<DissectedSlice> &result);

  OffsetLoopSet<Real> const *_inputSet;

  // counter clockwise offset loops, these surround the clockwise offset loops
  std::vector<OffsetLoop<Real>> _ccwOffsetLoops;
  // clockwise (island) offset loops, these are surrounded by the counter clockwise loops
  std::vector<OffsetLoop<Real>> _cwOffsetLoops;
  std::size_t totalOffsetLoopsCount() { return _ccwOffsetLoops.size() + _cwOffsetLoops.size(); }
  // spatial index of all the offset loops
  std::unique_ptr<StaticSpatialIndex<Real>> _offsetLoopsIndex;
  using IndexPair = std::pair<std::size_t, std::size_t>;
  // set to keep track of already visited pairs of loops when finding intersects
  std::unordered_set<IndexPair, internal::IndexPairHash> _visitedLoopPairs;
  // buffers to use for querying spatial indexes
  std::vector<std::size_t> _queryStack;
  std::vector<std::size_t> _queryResults;
  // slice point sets from intersects between loops
  std::vector<SlicePointSet> _slicePointSets;
  // lookup used to get slice points for a particular loop index (holds indexes to sets in
  // _slicePointSets)
  std::vector<std::vector<std::size_t>> _slicePointsLookup;
  // dissection points used to form slices for a particular loop in createSlicesFromLoop
  std::unordered_map<std::size_t, std::vector<DissectionPoint>> _loopDissectionPoints;
};

template <typename Real>
void ParallelOffsetIslands<Real>::createOffsetLoops(const OffsetLoopSet<Real> &input,
                                                    Real absDelta) {
  // create counter clockwise offset loops
  _ccwOffsetLoops.clear();
  std::size_t parentIndex = 0;
  for (auto const &loop : input.ccwLoops) {
    auto offsets = parallelOffset(loop.polyline, absDelta);
    for (auto &offset : offsets) {
      // must check if orientation inverted (due to collapse of very narrow or small input)
      if (getArea(offset) < Real(0)) {
        continue;
      }
      auto index = createApproxSpatialIndex(offset);
      _ccwOffsetLoops.push_back({parentIndex, std::move(offset), std::move(index)});
    }
    parentIndex += 1;
  }

  // create clockwise offset loops (note counter clockwise loops may result from outward offset)
  _cwOffsetLoops.clear();
  for (auto const &loop : input.cwLoops) {
    auto offsets = parallelOffset(loop.polyline, absDelta);
    for (auto &offset : offsets) {
      auto index = createApproxSpatialIndex(offset);
      if (getArea(offset) < Real(0)) {
        _cwOffsetLoops.push_back({parentIndex, std::move(offset), std::move(index)});
      } else {
        _ccwOffsetLoops.push_back({parentIndex, std::move(offset), std::move(index)});
      }
    }
    parentIndex += 1;
  }
}

template <typename Real> void ParallelOffsetIslands<Real>::createOffsetLoopsIndex() {
  // create spatial index for all offset loop bounding boxes
  _offsetLoopsIndex = std::make_unique<StaticSpatialIndex<Real>>(totalOffsetLoopsCount());
  for (auto const &posC : _ccwOffsetLoops) {
    auto const &i = posC.spatialIndex;
    _offsetLoopsIndex->add(i.minX(), i.minY(), i.maxX(), i.maxY());
  }

  for (auto const &negC : _cwOffsetLoops) {
    auto const &i = negC.spatialIndex;
    _offsetLoopsIndex->add(i.minX(), i.minY(), i.maxX(), i.maxY());
  }
  _offsetLoopsIndex->finish();
}

template <typename Real> void ParallelOffsetIslands<Real>::createSlicePoints() {
  _visitedLoopPairs.clear();
  _slicePointSets.clear();
  _slicePointsLookup.clear();

  // find all intersects between all offsets
  std::size_t totalOffsetCount = totalOffsetLoopsCount();
  _slicePointsLookup.resize(totalOffsetCount);
  PlineIntersectsResult<Real> intrsResults;
  for (std::size_t i = 0; i < totalOffsetCount; ++i) {
    auto const &loop1 = getOffsetLoop(i);
    auto const &index1 = loop1.spatialIndex;
    _queryResults.clear();
    _offsetLoopsIndex->query(index1.minX(), index1.minY(), index1.maxX(), index1.maxY(),
                              _queryResults, _queryStack);

    for (std::size_t j : _queryResults) {
      // skip same index (no self intersects among the offset loops)
      if (i == j) {
        continue;
      }
      // skip reversed index order (would end up comparing the same loops)
      if (_visitedLoopPairs.find({j, i}) != _visitedLoopPairs.end()) {
        continue;
      }
      _visitedLoopPairs.emplace(i, j);

      auto const &loop2 = getOffsetLoop(j);
      intrsResults.intersects.clear();
      intrsResults.coincidentIntersects.clear();
      // finding intersects
      findIntersects(loop1.polyline, loop2.polyline, index1, intrsResults);
      if (intrsResults.hasIntersects()) {
        _slicePointSets.emplace_back();
        auto &slicePointSet = _slicePointSets.back();
        slicePointSet.loopIndex1 = i;
        slicePointSet.loopIndex2 = j;
        for (auto &intr : intrsResults.intersects) {
          slicePointSet.slicePoints.push_back({std::move(intr), false});
        }

        // add coincident start and end points
        if (intrsResults.coincidentIntersects.size() != 0) {
          auto coinSliceResult = sortAndjoinCoincidentSlices(intrsResults.coincidentIntersects,
                                                             loop1.polyline, loop2.polyline);
          for (auto &sp : coinSliceResult.sliceStartPoints) {
            slicePointSet.slicePoints.push_back({std::move(sp), false});
          }
          for (auto &ep : coinSliceResult.sliceEndPoints) {
            slicePointSet.slicePoints.push_back({std::move(ep), true});
          }
        }

        _slicePointsLookup[i].push_back(_slicePointSets.size() - 1);
        _slicePointsLookup[j].push_back(_slicePointSets.size() - 1);
      }
    }
  }
}

template <typename Real>
bool ParallelOffsetIslands<Real>::pointOnOffsetValid(std::size_t skipIndex, const Vector2<Real> &pt,
                                                     Real absDelta) {
  // test distance against input polylines
  std::size_t const inputTotalCount = _inputSet->ccwLoops.size() + _inputSet->cwLoops.size();
  for (std::size_t i = 0; i < inputTotalCount; ++i) {
    if (i == skipIndex) {
      continue;
    }
    auto const &parentLoop = getParentLoop(i);
    if (!internal::pointValidForOffset(parentLoop.polyline, absDelta, parentLoop.spatialIndex, pt,
                                       _queryStack)) {
      return false;
    }
  }

  return true;
}

template <typename Real>
void ParallelOffsetIslands<Real>::createSlicesFromLoop(std::size_t loopIndex, Real absDelta,
                                                       std::vector<DissectedSlice> &result) {
  OffsetLoop<Real> const &offsetLoop = getOffsetLoop(loopIndex);
  std::size_t const parentIndex = offsetLoop.parentLoopIndex;
  Polyline<Real> const &pline = offsetLoop.polyline;
  _loopDissectionPoints.clear();
  for (auto const &setIndex : _slicePointsLookup[loopIndex]) {
    auto const &set = _slicePointSets[setIndex];
    bool isFirstIndex = loopIndex == set.loopIndex1;
    if (isFirstIndex) {
      for (auto const &p : set.slicePoints) {
        _loopDissectionPoints[p.intr.sIndex1].push_back({set.loopIndex2, p.intr.pos});
      }
    } else {
      for (auto const &p : set.slicePoints) {
        _loopDissectionPoints[p.intr.sIndex2].push_back({set.loopIndex1, p.intr.pos});
      }
    }
  }

  // sort points by distance from start vertex
  for (auto &kvp : _loopDissectionPoints) {
    Vector2<Real> startPos = pline[kvp.first].pos();
    auto cmp = [&](DissectionPoint const &p1, DissectionPoint const &p2) {
      return squared_distance(p1.pos, startPos) < squared_distance(p2.pos, startPos);
    };
    std::sort(kvp.second.begin(), kvp.second.end(), cmp);
  }


  for (auto const &kvp : _loopDissectionPoints) {
    // start index for the slice we're about to build
    std::size_t sIndex = kvp.first;
    // self intersect list for this start index
    std::vector<DissectionPoint> const &intrsList = kvp.second;

    const auto &firstSegStartVertex = pline[sIndex];
    std::size_t nextIndex = utils::nextWrappingIndex(sIndex, pline);
    const auto &firstSegEndVertex = pline[nextIndex];

    if (intrsList.size() != 1) {
      // build all the segments between the N intersects in siList (N > 1), skipping the first
      // segment (to be processed at the end)
      SplitResult<Real> firstSplit =
          splitAtPoint(firstSegStartVertex, firstSegEndVertex, intrsList[0].pos);
      auto prevVertex = firstSplit.splitVertex;
      for (std::size_t i = 1; i < intrsList.size(); ++i) {
        std::size_t const sliceStartIndex = intrsList[i - 1].otherLoopIndex;
        std::size_t const sliceEndIndex = intrsList[i].otherLoopIndex;
        SplitResult<Real> split = splitAtPoint(prevVertex, firstSegEndVertex, intrsList[i].pos);
        // update prevVertex for next loop iteration
        prevVertex = split.splitVertex;

        if (fuzzy::equal(split.updatedStart.pos(), split.splitVertex.pos(),
                       utils::realPrecision<Real>())) {
          continue;
        }

        auto sMidpoint = segMidpoint(split.updatedStart, split.splitVertex);
        if (!pointOnOffsetValid(parentIndex, sMidpoint, absDelta)) {
          // skip slice
          continue;
        }

        result.emplace_back();
        result.back().pline.addVertex(split.updatedStart);
        result.back().pline.addVertex(split.splitVertex);
        result.back().sliceParentIndex = loopIndex;
        result.back().startLoopIndex = sliceStartIndex;
        result.back().endLoopIndex = sliceEndIndex;
      }
    }

    // build the segment between the last intersect in instrsList and the next intersect found
    SplitResult<Real> split =
        splitAtPoint(firstSegStartVertex, firstSegEndVertex, intrsList.back().pos);

    DissectedSlice currSlice;
    currSlice.pline.addVertex(split.splitVertex);
    currSlice.sliceParentIndex = loopIndex;
    currSlice.startLoopIndex = intrsList.back().otherLoopIndex;

    std::size_t index = nextIndex;
    std::size_t loopCount = 0;
    const std::size_t maxLoopCount = pline.size();
    while (true) {
      if (loopCount++ > maxLoopCount) {
        PLLIB_ASSERT(false, "Bug detected, should never loop this many times!");
        // break to avoid infinite loop
        break;
      }
      // add vertex
      internal::addOrReplaceIfSamePos(currSlice.pline, pline[index]);

      // check if segment that starts at vertex we just added has an intersect
      auto nextIntr = _loopDissectionPoints.find(index);
      if (nextIntr != _loopDissectionPoints.end()) {
        // there is an intersect, slice is done
        Vector2<Real> const &intersectPos = nextIntr->second[0].pos;

        // trim last added vertex and add final intersect position
        PLVertex<Real> endVertex = PLVertex<Real>(intersectPos, Real(0));
        std::size_t nextIndex = utils::nextWrappingIndex(index, pline);
        SplitResult<Real> split =
            splitAtPoint(currSlice.pline.lastVertex(), pline[nextIndex], intersectPos);
        currSlice.pline.lastVertex() = split.updatedStart;
        internal::addOrReplaceIfSamePos(currSlice.pline, endVertex);
        currSlice.endLoopIndex = nextIntr->second[0].otherLoopIndex;
        break;
      }
      // else there is not an intersect, increment index and continue
      index = utils::nextWrappingIndex(index, pline);
    }

    if (currSlice.pline.size() > 1) {
      auto sMidpoint =
          segMidpoint(currSlice.pline[currSlice.pline.size() - 2], currSlice.pline.lastVertex());
      if (!pointOnOffsetValid(parentIndex, sMidpoint, absDelta)) {
        // skip slice
        continue;
      }
      result.push_back(std::move(currSlice));
    }
  }
}

template <typename Real>
OffsetLoopSet<Real> ParallelOffsetIslands<Real>::compute(const OffsetLoopSet<Real> &input,
                                                         Real offsetDelta) {
  _inputSet = &input;
  OffsetLoopSet<Real> result;
  Real absDelta = std::abs(offsetDelta);
  createOffsetLoops(input, absDelta);
  if (totalOffsetLoopsCount() == 0) {
    return result;
  }

  createOffsetLoopsIndex();
  createSlicePoints();

  std::vector<DissectedSlice> slices;
  std::size_t totalOffsetsCount = totalOffsetLoopsCount();

  std::vector<Polyline<Real>> resultSlices;

  for (std::size_t i = 0; i < totalOffsetsCount; ++i) {
    if (_slicePointsLookup[i].size() == 0) {
      // no intersects but still must test distance of one vertex position since it may be inside
      // another offset (completely eclipsed by island offset)
      auto &loop = getOffsetLoop(i);
      if (!pointOnOffsetValid(loop.parentLoopIndex, loop.polyline[0].pos(), absDelta)) {
        continue;
      }
      if (i < _ccwOffsetLoops.size()) {
        result.ccwLoops.push_back(std::move(loop));
      } else {
        result.cwLoops.push_back(std::move(loop));
      }
      continue;
    }
    createSlicesFromLoop(i, absDelta, slices);
  }

  for (auto &slice : slices) {
    resultSlices.push_back(std::move(slice.pline));
  }

  std::vector<Polyline<Real>> stitched =
      internal::stitchOrderedSlicesIntoClosedPolylines(resultSlices);

  for (auto &r : stitched) {
    Real area = getArea(r);
    if (std::abs(area) < 1e-4) {
      continue;
    }
    auto spatialIndex = createApproxSpatialIndex(r);
    if (area < Real(0)) {
      result.cwLoops.push_back({0, std::move(r), std::move(spatialIndex)});
    } else {
      result.ccwLoops.push_back({0, std::move(r), std::move(spatialIndex)});
    }
  }

  return result;
}

} // namespace cavc

#endif // CAVC_POLYLINEOFFSETISLANDS_HPP
