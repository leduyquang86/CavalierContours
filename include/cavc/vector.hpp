#ifndef CAVC_VECTOR_HPP
#define CAVC_VECTOR_HPP
#include "mathutils.hpp"
#include <array>
#include <cassert>
namespace cavc {
template <typename Real, std::size_t N> class Vector {
public:
  Vector() = default;

  Vector(std::initializer_list<Real> values) {
    if (N == values.size()) {
      std::copy(values.begin(), values.end(), _data.begin());
    } else if (N > values.size()) {
      std::copy(values.begin(), values.end(), _data.begin());
      std::fill(_data.begin() + static_cast<std::ptrdiff_t>(values.size()), _data.end(), Real(0));
    } else {
      std::copy(values.begin(), values.begin() + N, _data.begin());
    }
  }

  Vector(Real x, Real y) {
    static_assert(N == 2, "constructor for Vector2 only");
    _data[0] = x;
    _data[1] = y;
  }

  Vector(Real x, Real y, Real z) {
    static_assert(N == 3, "constructor for Vector3 only");
    _data[0] = x;
    _data[1] = y;
    _data[2] = z;
  }

  inline Real const &operator[](std::size_t i) const { return _data[i]; }

  inline Real &operator[](std::size_t i) { return _data[i]; }

  inline bool operator==(Vector const &vec) const { return _data == vec._data; }

  inline bool operator!=(Vector const &vec) const { return _data != vec._data; }

  inline bool operator<(Vector const &vec) const { return _data < vec._data; }

  inline bool operator<=(Vector const &vec) const { return _data <= vec._data; }

  inline bool operator>(Vector const &vec) const { return _data > vec._data; }

  inline bool operator>=(Vector const &vec) const { return _data >= vec._data; }

  void make_zero() { std::fill(_data.begin(), _data.end(), Real(0)); }

  void make_ones() { std::fill(_data.begin(), _data.end(), Real(1)); }

  void makeUnit(std::size_t d) {
    std::fill(_data.begin(), _data.end(), Real(0));
    _data[d] = Real(1);
  }

  static Vector zero() {
    Vector v;
    v.make_zero();
    return v;
  }

  static Vector ones() {
    Vector v;
    v.make_ones();
    return v;
  }

  static Vector unit(std::size_t d) {
    Vector v;
    v.makeUnit(d);
    return v;
  }

  Real &x() {
    static_assert(N >= 1, "N >= 1 to access x");
    return _data[0];
  }

  Real x() const {
    static_assert(N >= 1, "N >= 1 to access x");
    return _data[0];
  }

  Real &y() {
    static_assert(N >= 2, "N >= 2 to access y");
    return _data[1];
  }

  Real y() const {
    static_assert(N >= 2, "N >= 2 to access y");
    return _data[1];
  }

  Real &z() {
    static_assert(N >= 3, "N >= 3 to access z");
    return _data[2];
  }

  Real z() const {
    static_assert(N >= 3, "N >= 3 to access z");
    return _data[2];
  }

protected:
  std::array<Real, N> _data;
};

template <typename Real, std::size_t N>
bool fuzzy::zero(Vector<Real, N> const &v, Real epsilon = utils::realThreshold<Real>()) {
  bool allCompAreZero = std::abs(v[0]) < epsilon;
  for (std::size_t i = 1; i < N; ++i) {
    allCompAreZero = allCompAreZero && std::abs(v[i]) < epsilon;
  }

  return allCompAreZero;
}

template <typename Real, std::size_t N>
bool fuzzy::equal(Vector<Real, N> const &v1, Vector<Real, N> const &v2,
                Real epsilon = utils::realThreshold<Real>()) {
  for (std::size_t i = 0; i < N; ++i) {
    if (!utils::fuzzy::equal(v1[i], v2[i], epsilon)) {
      return false;
    }
  }

  return true;
}

template <typename Real, std::size_t N> Vector<Real, N> operator+(Vector<Real, N> const &v) {
  return v;
}

template <typename Real, std::size_t N> Vector<Real, N> operator-(Vector<Real, N> const &v) {
  Vector<Real, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = -v[i];
  }
  return result;
}

template <std::size_t N, typename Real>
Vector<Real, N> operator+(Vector<Real, N> const &v0, Vector<Real, N> const &v1) {
  Vector<Real, N> result = v0;
  return result += v1;
}

template <std::size_t N, typename Real>
Vector<Real, N> operator-(Vector<Real, N> const &v0, Vector<Real, N> const &v1) {
  Vector<Real, N> result;
  for (std::size_t i = 0; i < N; i++) {
    result[i] = v0[i] - v1[i];
  }
  return result;
}

template <std::size_t N, typename Real>
Vector<Real, N> operator*(Vector<Real, N> const &v, Real scalar) {
  Vector<Real, N> result = v;
  return result *= scalar;
}

template <std::size_t N, typename Real>
Vector<Real, N> operator*(Real scalar, Vector<Real, N> const &v) {
  Vector<Real, N> result = v;
  return result *= scalar;
}

template <std::size_t N, typename Real>
Vector<Real, N> operator/(Vector<Real, N> const &v, Real scalar) {
  Vector<Real, N> result = v;
  return result /= scalar;
}

template <std::size_t N, typename Real>
Vector<Real, N> &operator+=(Vector<Real, N> &v0, Vector<Real, N> const &v1) {
  for (std::size_t i = 0; i < N; ++i) {
    v0[i] += v1[i];
  }
  return v0;
}

template <std::size_t N, typename Real>
Vector<Real, N> &operator-=(Vector<Real, N> &v0, Vector<Real, N> const &v1) {
  for (std::size_t i = 0; i < N; ++i) {
    v0[i] -= v1[i];
  }
  return v0;
}

template <std::size_t N, typename Real>
Vector<Real, N> &operator*=(Vector<Real, N> &v, Real scalar) {
  for (std::size_t i = 0; i < N; ++i) {
    v[i] *= scalar;
  }
  return v;
}

template <std::size_t N, typename Real>
Vector<Real, N> &operator/=(Vector<Real, N> &v, Real scalar) {
  if (scalar != Real(0)) {
    Real invScalar = Real(1) / scalar;
    for (std::size_t i = 0; i < N; ++i) {
      v[i] *= invScalar;
    }
  } else {
    for (std::size_t i = 0; i < N; ++i) {
      v[i] = Real(0);
    }
  }
  return v;
}

template <std::size_t N, typename Real>
Vector<Real, N> operator*(Vector<Real, N> const &v0, Vector<Real, N> const &v1) {
  Vector<Real, N> result = v0;
  return result *= v1;
}

template <std::size_t N, typename Real>
Vector<Real, N> operator/(Vector<Real, N> const &v0, Vector<Real, N> const &v1) {
  Vector<Real, N> result = v0;
  return result /= v1;
}

template <std::size_t N, typename Real>
Vector<Real, N> &operator*=(Vector<Real, N> &v0, Vector<Real, N> const &v1) {
  for (std::size_t i = 0; i < N; ++i) {
    v0[i] *= v1[i];
  }
  return v0;
}

template <std::size_t N, typename Real>
Vector<Real, N> &operator/=(Vector<Real, N> &v0, Vector<Real, N> const &v1) {
  for (std::size_t i = 0; i < N; ++i) {
    v0[i] /= v1[i];
  }
  return v0;
}

template <std::size_t N, typename Real>
Real dot(Vector<Real, N> const &v0, Vector<Real, N> const &v1) {
  Real dot = v0[0] * v1[0];
  for (std::size_t i = 1; i < N; ++i) {
    dot += v0[i] * v1[i];
  }
  return dot;
}

template <std::size_t N, typename Real> Real length(Vector<Real, N> const &v) {
  return std::sqrt(dot(v, v));
}

template <std::size_t N, typename Real> Real normalize(Vector<Real, N> &v) {
  PLLIB_ASSERT(!fuzzy::zero(v), "normalize not defined for zero vector");
  Real length = std::sqrt(dot(v, v));
  v /= length;
  return length;
}
} // namespace cavc

#endif // CAVC_VECTOR_HPP
