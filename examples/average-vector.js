/*
  Copyright (C) 2013

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  https://github.com/johnmccutchan/ecmascript_simd/blob/master/src/ecmascript_simd.js
*/

// SIMD module.
var SIMD = {};

// private stuff.
var _PRIVATE = {}

_PRIVATE._f32 = new Float32Array(1);
_PRIVATE._i32 = new Int32Array(1);

_PRIVATE._f32x4 = new Float32Array(4);
_PRIVATE._f64x2 = new Float64Array(_PRIVATE._f32x4.buffer);
_PRIVATE._i32x4 = new Int32Array(_PRIVATE._f32x4.buffer);
_PRIVATE._i8x16 = new Int8Array(_PRIVATE._f32x4.buffer);

_PRIVATE._f32x8 = new Float32Array(8);
_PRIVATE._f64x4 = new Float64Array(4);
_PRIVATE._i32x8 = new Int32Array(8);

_PRIVATE.truncatef32 = function(x) {
  _PRIVATE._f32[0] = x;
  return _PRIVATE._f32[0];
}

_PRIVATE.truncatei32 = function(x) {
  _PRIVATE._i32[0] = x;
  return _PRIVATE._i32[0];
}

function checkFloat32x4(t) {
  if (!(t instanceof SIMD.float32x4)) {
    throw new TypeError("argument is not a float32x4.");
  }
}

function checkFloat64x2(t) {
  if (!(t instanceof SIMD.float64x2)) {
    throw new TypeError("argument is not a float64x2.");
  }
}

function checkInt32x4(t) {
  if (!(t instanceof SIMD.int32x4)) {
    throw new TypeError("argument is not a int32x4.");
  }
}

/**
  * Construct a new instance of float32x4 number.
  * @param {double} value used for x lane.
  * @param {double} value used for y lane.
  * @param {double} value used for z lane.
  * @param {double} value used for w lane.
  * @constructor
  */
SIMD.float32x4 = function(x, y, z, w) {
  if (arguments.length == 1) {
    checkFloat32x4(x);
    return x;
  }
  if (!(this instanceof SIMD.float32x4)) {
    return new SIMD.float32x4(x, y, z, w);
  }
  this.x_ = _PRIVATE.truncatef32(x);
  this.y_ = _PRIVATE.truncatef32(y);
  this.z_ = _PRIVATE.truncatef32(z);
  this.w_ = _PRIVATE.truncatef32(w);
}

/**
  * Construct a new instance of float32x4 number with 0.0 in all lanes.
  * @constructor
  */
SIMD.float32x4.zero = function() {
  return SIMD.float32x4(0.0, 0.0, 0.0, 0.0);
}

/**
  * Construct a new instance of float32x4 number with the same value
  * in all lanes.
  * @param {double} value used for all lanes.
  * @constructor
  */
SIMD.float32x4.splat = function(s) {
  return SIMD.float32x4(s, s, s, s);
}

/**
  * @param {float64x2} t An instance of float64x2.
  * @return {float32x4} A float32x4 with .x and .y from t
  */
SIMD.float32x4.fromFloat64x2 = function(t) {
  checkFloat64x2(t);
  var a = SIMD.float32x4.zero();
  a.x_ = t.x_;
  a.y_ = t.y_;
  return a;
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @return {float32x4} A float to integer conversion copy of t.
  */
SIMD.float32x4.fromInt32x4 = function(t) {
  checkInt32x4(t);
  var a = SIMD.float32x4.zero();
  a.x_ = t.x_;
  a.y_ = t.y_;
  a.z_ = t.z_;
  a.w_ = t.w_;
  return a;
}

/**
 * @param {float64x2} t An instance of float64x2.
 * @return {float32x4} a bit-wise copy of t as a float32x4.
 */
SIMD.float32x4.fromFloat64x2Bits = function(t) {
  checkFloat64x2(t);
  _PRIVATE._f64x2[0] = t.x_;
  _PRIVATE._f64x2[1] = t.y_;
  return SIMD.float32x4(_PRIVATE._f32x4[0],
                        _PRIVATE._f32x4[1],
                        _PRIVATE._f32x4[2],
                        _PRIVATE._f32x4[3]);
}

/**
 * @param {int32x4} t An instance of int32x4.
 * @return {float32x4} a bit-wise copy of t as a float32x4.
 */
SIMD.float32x4.fromInt32x4Bits = function(t) {
  checkInt32x4(t);
  _PRIVATE._i32x4[0] = t.x_;
  _PRIVATE._i32x4[1] = t.y_;
  _PRIVATE._i32x4[2] = t.z_;
  _PRIVATE._i32x4[3] = t.w_;
  return SIMD.float32x4(_PRIVATE._f32x4[0],
                        _PRIVATE._f32x4[1],
                        _PRIVATE._f32x4[2],
                        _PRIVATE._f32x4[3]);
}

/**
  * Construct a new instance of float64x2 number.
  * @param {double} value used for x lane.
  * @param {double} value used for y lane.
  * @constructor
  */
SIMD.float64x2 = function(x, y) {
  if (arguments.length == 1) {
    checkFloat64x2(x);
    return x;
  }
  if (!(this instanceof SIMD.float64x2)) {
    return new SIMD.float64x2(x, y);
  }
  this.x_ = x;
  this.y_ = y;
}

/**
  * Construct a new instance of float64x2 number with 0.0 in all lanes.
  * @constructor
  */
SIMD.float64x2.zero = function() {
  return SIMD.float64x2(0.0, 0.0);
}

/**
  * Construct a new instance of float64x2 number with the same value
  * in all lanes.
  * @param {double} value used for all lanes.
  * @constructor
  */
SIMD.float64x2.splat = function(s) {
  return SIMD.float32x4(s, s);
}

/**
  * @param {float32x4} t An instance of float32x4.
  * @return {float64x2} A float64x2 with .x and .y from t
  */
SIMD.float64x2.fromFloat32x4 = function(t) {
  checkFloat32x4(t);
  var a = SIMD.float64x2(t.x_, t.y_);
  return a;
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @return {float64x2} A float64x2 with .x and .y from t
  */
SIMD.float64x2.fromInt32x4 = function(t) {
  checkInt32x4(t);
  var a = SIMD.float64x2.zero();
  a.x_ = t.x_;
  a.y_ = t.y_;
  return a;
}

/**
 * @param {float32x4} t An instance of float32x4.
 * @return {float64x2} a bit-wise copy of t as a float64x2.
 */
SIMD.float64x2.fromFloat32x4Bits = function(t) {
  checkFloat32x4(t);
  _PRIVATE._f32x4[0] = t.x_;
  _PRIVATE._f32x4[1] = t.y_;
  _PRIVATE._f32x4[2] = t.z_;
  _PRIVATE._f32x4[3] = t.w_;
  return SIMD.float64x2(_PRIVATE._f64x2[0], _PRIVATE._f64x2[1]);
}

/**
 * @param {int32x4} t An instance of int32x4.
 * @return {float64x2} a bit-wise copy of t as a float64x2.
 */
SIMD.float64x2.fromInt32x4Bits = function(t) {
  checkInt32x4(t);
  _PRIVATE._i32x4[0] = t.x_;
  _PRIVATE._i32x4[1] = t.y_;
  _PRIVATE._i32x4[2] = t.z_;
  _PRIVATE._i32x4[3] = t.w_;
  return SIMD.float64x2(_PRIVATE._f64x2[0], _PRIVATE._f64x2[1]);
}

/**
  * Construct a new instance of int32x4 number.
  * @param {integer} 32-bit unsigned value used for x lane.
  * @param {integer} 32-bit unsigned value used for y lane.
  * @param {integer} 32-bit unsigned value used for z lane.
  * @param {integer} 32-bit unsigned value used for w lane.
  * @constructor
  */
SIMD.int32x4 = function(x, y, z, w) {
  if (arguments.length == 1) {
    checkInt32x4(x);
    return x;
  }
  if (!(this instanceof SIMD.int32x4)) {
    return new SIMD.int32x4(x, y, z, w);
  }
  this.x_ = _PRIVATE.truncatei32(x);
  this.y_ = _PRIVATE.truncatei32(y);
  this.z_ = _PRIVATE.truncatei32(z);
  this.w_ = _PRIVATE.truncatei32(w);
}

/**
  * Construct a new instance of int32x4 number with 0 in all lanes.
  * @constructor
  */
SIMD.int32x4.zero = function() {
  return SIMD.int32x4(0, 0, 0, 0);
}

/**
  * Construct a new instance of int32x4 number with 0xFFFFFFFF or 0x0 in each
  * lane, depending on the truth value in x, y, z, and w.
  * @param {boolean} flag used for x lane.
  * @param {boolean} flag used for y lane.
  * @param {boolean} flag used for z lane.
  * @param {boolean} flag used for w lane.
  * @constructor
  */
SIMD.int32x4.bool = function(x, y, z, w) {
  return SIMD.int32x4(x ? -1 : 0x0,
                      y ? -1 : 0x0,
                      z ? -1 : 0x0,
                      w ? -1 : 0x0);
}

/**
  * Construct a new instance of int32x4 number with the same value
  * in all lanes.
  * @param {integer} value used for all lanes.
  * @constructor
  */
SIMD.int32x4.splat = function(s) {
  return SIMD.int32x4(s, s, s, s);
}

/**
  * @param {float32x4} t An instance of float32x4.
  * @return {int32x4} with a integer to float conversion of t.
  */
SIMD.int32x4.fromFloat32x4 = function(t) {
  checkFloat32x4(t);
  var a = SIMD.int32x4(Math.floor(t.x_), Math.floor(t.y_),
                       Math.floor(t.z_), Math.floor(t.w_));
  return a;
}

/**
  * @param {float64x2} t An instance of float64x2.
  * @return {int32x4}  An int32x4 with .x and .y from t
  */
SIMD.int32x4.fromFloat64x2 = function(t) {
  checkFloat64x2(t);
  var a = SIMD.int32x4.zero();
  a.x_ = Math.floor(t.x_);
  a.y_ = Math.floor(t.y_);
  return a;
}

/**
  * @param {float32x4} t An instance of float32x4.
  * @return {int32x4} a bit-wise copy of t as a int32x4.
  */
SIMD.int32x4.fromFloat32x4Bits = function(t) {
  checkFloat32x4(t);
  _PRIVATE._f32x4[0] = t.x_;
  _PRIVATE._f32x4[1] = t.y_;
  _PRIVATE._f32x4[2] = t.z_;
  _PRIVATE._f32x4[3] = t.w_;
  var alias = _PRIVATE._i32x4;
  return SIMD.int32x4(alias[0], alias[1], alias[2], alias[3]);
}

/**
 * @param {float64x2} t An instance of float64x2.
 * @return {int32x4} a bit-wise copy of t as an int32x4.
 */
SIMD.int32x4.fromFloat64x2Bits = function(t) {
  checkFloat64x2(t);
  _PRIVATE._f64x2[0] = t.x_;
  _PRIVATE._f64x2[1] = t.y_;
  var alias = _PRIVATE._i32x4;
  var ix4 = SIMD.int32x4(alias[0], alias[1], alias[2], alias[3]);
  return ix4;
}

/**
* @return {float32x4} New instance of float32x4 with absolute values of
* t.
*/
SIMD.float32x4.abs = function(t) {
  checkFloat32x4(t);
  return SIMD.float32x4(Math.abs(t.x), Math.abs(t.y), Math.abs(t.z),
                        Math.abs(t.w));
}

/**
  * @return {float32x4} New instance of float32x4 with negated values of
  * t.
  */
SIMD.float32x4.neg = function(t) {
  checkFloat32x4(t);
  return SIMD.float32x4(-t.x, -t.y, -t.z, -t.w);
}

/**
  * @return {float32x4} New instance of float32x4 with a + b.
  */
SIMD.float32x4.add = function(a, b) {
  checkFloat32x4(a);
  checkFloat32x4(b);
  return SIMD.float32x4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

/**
  * @return {float32x4} New instance of float32x4 with a - b.
  */
SIMD.float32x4.sub = function(a, b) {
  checkFloat32x4(a);
  checkFloat32x4(b);
  return SIMD.float32x4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

/**
  * @return {float32x4} New instance of float32x4 with a * b.
  */
SIMD.float32x4.mul = function(a, b) {
  checkFloat32x4(a);
  checkFloat32x4(b);
  return SIMD.float32x4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

/**
  * @return {float32x4} New instance of float32x4 with a / b.
  */
SIMD.float32x4.div = function(a, b) {
  checkFloat32x4(a);
  checkFloat32x4(b);
  return SIMD.float32x4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}

/**
  * @return {float32x4} New instance of float32x4 with t's values clamped
  * between lowerLimit and upperLimit.
  */
SIMD.float32x4.clamp = function(t, lowerLimit, upperLimit) {
  checkFloat32x4(t);
  var cx = t.x < lowerLimit.x ? lowerLimit.x : t.x;
  var cy = t.y < lowerLimit.y ? lowerLimit.y : t.y;
  var cz = t.z < lowerLimit.z ? lowerLimit.z : t.z;
  var cw = t.w < lowerLimit.w ? lowerLimit.w : t.w;
  cx = cx > upperLimit.x ? upperLimit.x : cx;
  cy = cy > upperLimit.y ? upperLimit.y : cy;
  cz = cz > upperLimit.z ? upperLimit.z : cz;
  cw = cw > upperLimit.w ? upperLimit.w : cw;
  return SIMD.float32x4(cx, cy, cz, cw);
}

/**
  * @return {float32x4} New instance of float32x4 with the minimum value of
  * t and other.
  */
SIMD.float32x4.min = function(t, other) {
  checkFloat32x4(t);
  var cx = t.x > other.x ? other.x : t.x;
  var cy = t.y > other.y ? other.y : t.y;
  var cz = t.z > other.z ? other.z : t.z;
  var cw = t.w > other.w ? other.w : t.w;
  return SIMD.float32x4(cx, cy, cz, cw);
}

/**
  * @return {float32x4} New instance of float32x4 with the maximum value of
  * t and other.
  */
SIMD.float32x4.max = function(t, other) {
  checkFloat32x4(t);
  var cx = t.x < other.x ? other.x : t.x;
  var cy = t.y < other.y ? other.y : t.y;
  var cz = t.z < other.z ? other.z : t.z;
  var cw = t.w < other.w ? other.w : t.w;
  return SIMD.float32x4(cx, cy, cz, cw);
}

/**
  * @return {float32x4} New instance of float32x4 with reciprocal value of
  * t.
  */
SIMD.float32x4.reciprocal = function(t) {
  checkFloat32x4(t);
  return SIMD.float32x4(1.0 / t.x, 1.0 / t.y, 1.0 / t.z, 1.0 / t.w);
}

/**
  * @return {float32x4} New instance of float32x4 with square root of the
  * reciprocal value of t.
  */
SIMD.float32x4.reciprocalSqrt = function(t) {
  checkFloat32x4(t);
  return SIMD.float32x4(Math.sqrt(1.0 / t.x), Math.sqrt(1.0 / t.y),
                        Math.sqrt(1.0 / t.z), Math.sqrt(1.0 / t.w));
}

/**
  * @return {float32x4} New instance of float32x4 with values of t
  * scaled by s.
  */
SIMD.float32x4.scale = function(t, s) {
  checkFloat32x4(t);
  return SIMD.float32x4(s * t.x, s * t.y, s * t.z, s * t.w);
}

/**
  * @return {float32x4} New instance of float32x4 with square root of
  * values of t.
  */
SIMD.float32x4.sqrt = function(t) {
  checkFloat32x4(t);
  return SIMD.float32x4(Math.sqrt(t.x), Math.sqrt(t.y),
                        Math.sqrt(t.z), Math.sqrt(t.w));
}

/**
  * @param {float32x4} t An instance of float32x4 to be shuffled.
  * @param {integer} mask One of the 256 shuffle masks, for example, SIMD.XXXX.
  * @return {float32x4} New instance of float32x4 with lanes shuffled.
  */
SIMD.float32x4.shuffle = function(t, mask) {
  checkFloat32x4(t);
  var _x = (mask) & 0x3;
  var _y = (mask >> 2) & 0x3;
  var _z = (mask >> 4) & 0x3;
  var _w = (mask >> 6) & 0x3;
  _PRIVATE._f32x4[0] = t.x_;
  _PRIVATE._f32x4[1] = t.y_;
  _PRIVATE._f32x4[2] = t.z_;
  _PRIVATE._f32x4[3] = t.w_;
  var storage = _PRIVATE._f32x4;
  return SIMD.float32x4(storage[_x], storage[_y], storage[_z], storage[_w]);
}

/**
  * @param {float32x4} t1 An instance of float32x4 to be shuffled. XY lanes in result
  * @param {float32x4} t2 An instance of float32x4 to be shuffled. ZW lanes in result
  * @param {integer} mask One of the 256 shuffle masks, for example, SIMD.XXXX.
  * @return {float32x4} New instance of float32x4 with lanes shuffled.
  */
SIMD.float32x4.shuffleMix = function(t1, t2, mask) {
  checkFloat32x4(t1);
  checkFloat32x4(t2);
  var _x = (mask) & 0x3;
  var _y = (mask >> 2) & 0x3;
  var _z = (mask >> 4) & 0x3;
  var _w = (mask >> 6) & 0x3;
  var storage = _PRIVATE._f32x8;
  storage[0] = t1.x_;
  storage[1] = t1.y_;
  storage[2] = t1.z_;
  storage[3] = t1.w_;
  storage[4] = t2.x_;
  storage[5] = t2.y_;
  storage[6] = t2.z_;
  storage[7] = t2.w_;
  return SIMD.float32x4(storage[0 + _x], storage[0 + _y],
                        storage[4 + _z], storage[4 + _w]);
}

/**
  * @param {double} value used for x lane.
  * @return {float32x4} New instance of float32x4 with the values in t and
  * x replaced with {x}.
  */
SIMD.float32x4.withX = function(t, x) {
  checkFloat32x4(t);
  return SIMD.float32x4(x, t.y, t.z, t.w);
}

/**
  * @param {double} value used for y lane.
  * @return {float32x4} New instance of float32x4 with the values in t and
  * y replaced with {y}.
  */
SIMD.float32x4.withY = function(t, y) {
  checkFloat32x4(t);
  return SIMD.float32x4(t.x, y, t.z, t.w);
}

/**
  * @param {double} value used for z lane.
  * @return {float32x4} New instance of float32x4 with the values in t and
  * z replaced with {z}.
  */
SIMD.float32x4.withZ = function(t, z) {
  checkFloat32x4(t);
  return SIMD.float32x4(t.x, t.y, z, t.w);
}

/**
  * @param {double} value used for w lane.
  * @return {float32x4} New instance of float32x4 with the values in t and
  * w replaced with {w}.
  */
SIMD.float32x4.withW = function(t, w) {
  checkFloat32x4(t);
  return SIMD.float32x4(t.x, t.y, t.z, w);
}

/**
  * @param {float32x4} t An instance of float32x4.
  * @param {float32x4} other An instance of float32x4.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t < other.
  */
SIMD.float32x4.lessThan = function(t, other) {
  checkFloat32x4(t);
  checkFloat32x4(other);
  var cx = t.x < other.x;
  var cy = t.y < other.y;
  var cz = t.z < other.z;
  var cw = t.w < other.w;
  return SIMD.int32x4.bool(cx, cy, cz, cw);
}

/**
  * @param {float32x4} t An instance of float32x4.
  * @param {float32x4} other An instance of float32x4.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t <= other.
  */
SIMD.float32x4.lessThanOrEqual = function(t, other) {
  checkFloat32x4(t);
  checkFloat32x4(other);
  var cx = t.x <= other.x;
  var cy = t.y <= other.y;
  var cz = t.z <= other.z;
  var cw = t.w <= other.w;
  return SIMD.int32x4.bool(cx, cy, cz, cw);
}

/**
  * @param {float32x4} t An instance of float32x4.
  * @param {float32x4} other An instance of float32x4.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t == other.
  */
SIMD.float32x4.equal = function(t, other) {
  checkFloat32x4(t);
  checkFloat32x4(other);
  var cx = t.x == other.x;
  var cy = t.y == other.y;
  var cz = t.z == other.z;
  var cw = t.w == other.w;
  return SIMD.int32x4.bool(cx, cy, cz, cw);
}

/**
  * @param {float32x4} t An instance of float32x4.
  * @param {float32x4} other An instance of float32x4.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t != other.
  */
SIMD.float32x4.notEqual = function(t, other) {
  checkFloat32x4(t);
  checkFloat32x4(other);
  var cx = t.x != other.x;
  var cy = t.y != other.y;
  var cz = t.z != other.z;
  var cw = t.w != other.w;
  return SIMD.int32x4.bool(cx, cy, cz, cw);
}

/**
  * @param {float32x4} t An instance of float32x4.
  * @param {float32x4} other An instance of float32x4.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t >= other.
  */
SIMD.float32x4.greaterThanOrEqual = function(t, other) {
  checkFloat32x4(t);
  checkFloat32x4(other);
  var cx = t.x >= other.x;
  var cy = t.y >= other.y;
  var cz = t.z >= other.z;
  var cw = t.w >= other.w;
  return SIMD.int32x4.bool(cx, cy, cz, cw);
}

/**
  * @param {float32x4} t An instance of float32x4.
  * @param {float32x4} other An instance of float32x4.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t > other.
  */
SIMD.float32x4.greaterThan = function(t, other) {
  checkFloat32x4(t);
  checkFloat32x4(other);
  var cx = t.x > other.x;
  var cy = t.y > other.y;
  var cz = t.z > other.z;
  var cw = t.w > other.w;
  return SIMD.int32x4.bool(cx, cy, cz, cw);
}

/**
  * @param {int32x4} t Selector mask. An instance of int32x4
  * @param {float32x4} trueValue Pick lane from here if corresponding
  * selector lane is 0xFFFFFFFF
  * @param {float32x4} falseValue Pick lane from here if corresponding
  * selector lane is 0x0
  * @return {float32x4} Mix of lanes from trueValue or falseValue as
  * indicated
  */
SIMD.float32x4.select = function(t, trueValue, falseValue) {
  checkInt32x4(t);
  checkFloat32x4(trueValue);
  checkFloat32x4(falseValue);
  var tv = SIMD.int32x4.fromFloat32x4Bits(trueValue);
  var fv = SIMD.int32x4.fromFloat32x4Bits(falseValue);
  var tr = SIMD.int32x4.and(t, tv);
  var fr = SIMD.int32x4.and(SIMD.int32x4.not(t), fv);
  return SIMD.float32x4.fromInt32x4Bits(SIMD.int32x4.or(tr, fr));
}

/**
  * @param {float32x4} a An instance of float32x4.
  * @param {float32x4} b An instance of float32x4.
  * @return {float32x4} New instance of float32x4 with values of a & b.
  */
SIMD.float32x4.and = function(a, b) {
  checkFloat32x4(a);
  checkFloat32x4(b);
  var aInt = SIMD.int32x4.fromFloat32x4Bits(a);
  var bInt = SIMD.int32x4.fromFloat32x4Bits(b);
  return SIMD.float32x4.fromInt32x4Bits(SIMD.int32x4.and(aInt, bInt));
}

/**
  * @param {float32x4} a An instance of float32x4.
  * @param {float32x4} b An instance of float32x4.
  * @return {float32x4} New instance of float32x4 with values of a | b.
  */
SIMD.float32x4.or = function(a, b) {
  checkFloat32x4(a);
  checkFloat32x4(b);
  var aInt = SIMD.int32x4.fromFloat32x4Bits(a);
  var bInt = SIMD.int32x4.fromFloat32x4Bits(b);
  return SIMD.float32x4.fromInt32x4Bits(SIMD.int32x4.or(aInt, bInt));
}

/**
  * @param {float32x4} a An instance of float32x4.
  * @param {float32x4} b An instance of float32x4.
  * @return {float32x4} New instance of float32x4 with values of a ^ b.
  */
SIMD.float32x4.xor = function(a, b) {
  checkFloat32x4(a);
  checkFloat32x4(b);
  var aInt = SIMD.int32x4.fromFloat32x4Bits(a);
  var bInt = SIMD.int32x4.fromFloat32x4Bits(b);
  return SIMD.float32x4.fromInt32x4Bits(SIMD.int32x4.xor(aInt, bInt));
}

/**
  * @param {float32x4} a An instance of float32x4.
  * @return {float32x4} New instance of float32x4 with values of ~a.
  */
SIMD.float32x4.not = function(a) {
  checkFloat32x4(a);
  var aInt = SIMD.int32x4.fromFloat32x4Bits(a);
  return SIMD.float32x4.fromInt32x4Bits(SIMD.int32x4.not(aInt));
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @return {float32x4} New instance of float32x4.
  */
SIMD.float32x4.load = function(buffer, byteOffset) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 16) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  var i8a = new Int8Array(buffer, byteOffset, 16);
  for (var i = 0; i < 16; i++)
    _PRIVATE._i8x16[i] = i8a[i];
  return SIMD.float32x4(_PRIVATE._f32x4[0], _PRIVATE._f32x4[1],
                        _PRIVATE._f32x4[2], _PRIVATE._f32x4[3]);
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @return {float32x4} New instance of float32x4.
  */
SIMD.float32x4.loadX = function(buffer, byteOffset) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 4) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  var i8a = new Int8Array(buffer, byteOffset, 4);
  for (var i = 0; i < 4; i++)
    _PRIVATE._i8x16[i] = i8a[i];
  return SIMD.float32x4(_PRIVATE._f32x4[0], 0.0, 0.0, 0.0);
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @return {float32x4} New instance of float32x4.
  */
SIMD.float32x4.loadXY = function(buffer, byteOffset) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 8) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  var i8a = new Int8Array(buffer, byteOffset, 8);
  for (var i = 0; i < 8; i++)
    _PRIVATE._i8x16[i] = i8a[i];
  return SIMD.float32x4(_PRIVATE._f32x4[0], _PRIVATE._f32x4[1], 0.0, 0.0);
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @return {float32x4} New instance of float32x4.
  */
SIMD.float32x4.loadXYZ = function(buffer, byteOffset) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 12) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  var i8a = new Int8Array(buffer, byteOffset, 12);
  for (var i = 0; i < 12; i++)
    _PRIVATE._i8x16[i] = i8a[i];
  return SIMD.float32x4(_PRIVATE._f32x4[0], _PRIVATE._f32x4[1],
                        _PRIVATE._f32x4[2], 0.0);
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @param {float32x4} value An instance of float32x4.
  * @return {void}
  */
SIMD.float32x4.store = function(buffer, byteOffset, value) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 16) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkFloat32x4(value);
  _PRIVATE._f32x4[0] = value.x;
  _PRIVATE._f32x4[1] = value.y;
  _PRIVATE._f32x4[2] = value.z;
  _PRIVATE._f32x4[3] = value.w;
  var i8a = new Int8Array(buffer, byteOffset, 16);
  for (var i = 0; i < 16; i++)
    i8a[i] = _PRIVATE._i8x16[i];
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @param {float32x4} value An instance of float32x4.
  * @return {void}
  */
SIMD.float32x4.storeX = function(buffer, byteOffset, value) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 4) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkFloat32x4(value);
  _PRIVATE._f32x4[0] = value.x;
  var i8a = new Int8Array(buffer, byteOffset, 4);
  for (var i = 0; i < 4; i++)
    i8a[i] = _PRIVATE._i8x16[i];
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @param {float32x4} value An instance of float32x4.
  * @return {void}
  */
SIMD.float32x4.storeXY = function(buffer, byteOffset, value) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 8) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkFloat32x4(value);
  _PRIVATE._f32x4[0] = value.x;
  _PRIVATE._f32x4[1] = value.y;
  var i8a = new Int8Array(buffer, byteOffset, 8);
  for (var i = 0; i < 8; i++)
    i8a[i] = _PRIVATE._i8x16[i];
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @param {float32x4} value An instance of float32x4.
  * @return {void}
  */
SIMD.float32x4.storeXYZ = function(buffer, byteOffset, value) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 12) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkFloat32x4(value);
  _PRIVATE._f32x4[0] = value.x;
  _PRIVATE._f32x4[1] = value.y;
  _PRIVATE._f32x4[2] = value.z;
  var i8a = new Int8Array(buffer, byteOffset, 12);
  for (var i = 0; i < 12; i++)
    i8a[i] = _PRIVATE._i8x16[i];
}

/**
* @return {float64x2} New instance of float64x2 with absolute values of
* t.
*/
SIMD.float64x2.abs = function(t) {
  checkFloat64x2(t);
  return SIMD.float64x2(Math.abs(t.x), Math.abs(t.y));
}

/**
  * @return {float64x2} New instance of float64x2 with negated values of
  * t.
  */
SIMD.float64x2.neg = function(t) {
  checkFloat64x2(t);
  return SIMD.float64x2(-t.x, -t.y);
}

/**
  * @return {float64x2} New instance of float64x2 with a + b.
  */
SIMD.float64x2.add = function(a, b) {
  checkFloat64x2(a);
  checkFloat64x2(b);
  return SIMD.float64x2(a.x + b.x, a.y + b.y);
}

/**
  * @return {float64x2} New instance of float64x2 with a - b.
  */
SIMD.float64x2.sub = function(a, b) {
  checkFloat64x2(a);
  checkFloat64x2(b);
  return SIMD.float64x2(a.x - b.x, a.y - b.y);
}

/**
  * @return {float64x2} New instance of float64x2 with a * b.
  */
SIMD.float64x2.mul = function(a, b) {
  checkFloat64x2(a);
  checkFloat64x2(b);
  return SIMD.float64x2(a.x * b.x, a.y * b.y);
}

/**
  * @return {float64x2} New instance of float64x2 with a / b.
  */
SIMD.float64x2.div = function(a, b) {
  checkFloat64x2(a);
  checkFloat64x2(b);
  return SIMD.float64x2(a.x / b.x, a.y / b.y);
}

/**
  * @return {float64x2} New instance of float64x2 with t's values clamped
  * between lowerLimit and upperLimit.
  */
SIMD.float64x2.clamp = function(t, lowerLimit, upperLimit) {
  checkFloat64x2(t);
  var cx = t.x < lowerLimit.x ? lowerLimit.x : t.x;
  var cy = t.y < lowerLimit.y ? lowerLimit.y : t.y;
  cx = cx > upperLimit.x ? upperLimit.x : cx;
  cy = cy > upperLimit.y ? upperLimit.y : cy;
  return SIMD.float64x2(cx, cy);
}

/**
  * @return {float64x2} New instance of float64x2 with the minimum value of
  * t and other.
  */
SIMD.float64x2.min = function(t, other) {
  checkFloat64x2(t);
  checkFloat64x2(other);
  var cx = t.x > other.x ? other.x : t.x;
  var cy = t.y > other.y ? other.y : t.y;
  return SIMD.float64x2(cx, cy);
}

/**
  * @return {float64x2} New instance of float64x2 with the maximum value of
  * t and other.
  */
SIMD.float64x2.max = function(t, other) {
  checkFloat64x2(t);
  checkFloat64x2(other);
  var cx = t.x < other.x ? other.x : t.x;
  var cy = t.y < other.y ? other.y : t.y;
  return SIMD.float64x2(cx, cy);
}

/**
  * @return {float64x2} New instance of float64x2 with reciprocal value of
  * t.
  */
SIMD.float64x2.reciprocal = function(t) {
  checkFloat64x2(t);
  return SIMD.float64x2(1.0 / t.x, 1.0 / t.y);
}

/**
  * @return {float64x2} New instance of float64x2 with square root of the
  * reciprocal value of t.
  */
SIMD.float64x2.reciprocalSqrt = function(t) {
  checkFloat64x2(t);
  return SIMD.float64x2(Math.sqrt(1.0 / t.x), Math.sqrt(1.0 / t.y));
}

/**
  * @return {float64x2} New instance of float32x4 with values of t
  * scaled by s.
  */
SIMD.float64x2.scale = function(t, s) {
  checkFloat64x2(t);
  return SIMD.float64x2(s * t.x, s * t.y);
}

/**
  * @return {float64x2} New instance of float32x4 with square root of
  * values of t.
  */
SIMD.float64x2.sqrt = function(t) {
  checkFloat64x2(t);
  return SIMD.float64x2(Math.sqrt(t.x), Math.sqrt(t.y));
}

/**
  * @param {float64x2} t An instance of float64x2 to be shuffled.
  * @param {integer} mask One of the 4 shuffle masks, for example, SIMD.XY.
  * @return {float64x2} New instance of float64x2 with lanes shuffled.
  */
SIMD.float64x2.shuffle = function(t, mask) {
  checkFloat64x2(t);
  var _x = (mask) & 0x1;
  var _y = (mask >> 1) & 0x1;
  var storage = _PRIVATE._f64x2;
  storage[0] = t.x_;
  storage[1] = t.y_;
  return SIMD.float64x2(storage[_x], storage[_y]);
}

/**
  * @param {float64x2} t1 An instance of float64x2 to be shuffled. X lane in result
  * @param {float64x2} t2 An instance of float64x2 to be shuffled. Y lane in result
  * @param {integer} mask One of the 4 shuffle masks, for example, SIMD.XY.
  * @return {float64x2} New instance of float64x2 with lanes shuffled.
  */
SIMD.float64x2.shuffleMix = function(t1, t2, mask) {
  checkFloat64x2(t1);
  checkFloat64x2(t2);
  var _x = (mask) & 0x1;
  var _y = (mask >> 1) & 0x1;
  var storage = _PRIVATE._f64x4;
  storage[0] = t1.x_;
  storage[1] = t1.y_;
  storage[2] = t2.x_;
  storage[3] = t2.y_;
  return SIMD.float64x2(storage[0 + _x], storage[2 + _y]);
}

/**
  * @param {double} value used for x lane.
  * @return {float64x2} New instance of float64x2 with the values in t and
  * x replaced with {x}.
  */
SIMD.float64x2.withX = function(t, x) {
  checkFloat64x2(t);
  return SIMD.float64x2(x, t.y);
}

/**
  * @param {double} value used for y lane.
  * @return {float64x2} New instance of float64x2 with the values in t and
  * y replaced with {y}.
  */
SIMD.float64x2.withY = function(t, y) {
  checkFloat64x2(t);
  return SIMD.float64x2(t.x, y);
}

/**
  * @param {float64x2} t An instance of float64x2.
  * @param {float64x2} other An instance of float64x2.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t < other.
  */
SIMD.float64x2.lessThan = function(t, other) {
  checkFloat64x2(t);
  checkFloat64x2(other);
  var cx = t.x < other.x;
  var cy = t.y < other.y;
  return SIMD.int32x4.bool(cx, cx, cy, cy);
}

/**
  * @param {float64x2} t An instance of float64x2.
  * @param {float64x2} other An instance of float64x2.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t <= other.
  */
SIMD.float64x2.lessThanOrEqual = function(t, other) {
  checkFloat64x2(t);
  checkFloat64x2(other);
  var cx = t.x <= other.x;
  var cy = t.y <= other.y;
  return SIMD.int32x4.bool(cx, cx, cy, cy);
}

/**
  * @param {float64x2} t An instance of float64x2.
  * @param {float64x2} other An instance of float64x2.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t == other.
  */
SIMD.float64x2.equal = function(t, other) {
  checkFloat64x2(t);
  checkFloat64x2(other);
  var cx = t.x == other.x;
  var cy = t.y == other.y;
  return SIMD.int32x4.bool(cx, cx, cy, cy);
}

/**
  * @param {float64x2} t An instance of float64x2.
  * @param {float64x2} other An instance of float64x2.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t != other.
  */
SIMD.float64x2.notEqual = function(t, other) {
  checkFloat64x2(t);
  checkFloat64x2(other);
  var cx = t.x != other.x;
  var cy = t.y != other.y;
  return SIMD.int32x4.bool(cx, cx, cy, cy);
}

/**
  * @param {float64x2} t An instance of float64x2.
  * @param {float64x2} other An instance of float64x2.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t >= other.
  */
SIMD.float64x2.greaterThanOrEqual = function(t, other) {
  checkFloat64x2(t);
  checkFloat64x2(other);
  var cx = t.x >= other.x;
  var cy = t.y >= other.y;
  return SIMD.int32x4.bool(cx, cx, cy, cy);
}

/**
  * @param {float64x2} t An instance of float64x2.
  * @param {float64x2} other An instance of float64x2.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t > other.
  */
SIMD.float64x2.greaterThan = function(t, other) {
  checkFloat64x2(t);
  checkFloat64x2(other);
  var cx = t.x > other.x;
  var cy = t.y > other.y;
  return SIMD.int32x4.bool(cx, cx, cy, cy);
}

/**
  * @param {int32x4} t Selector mask. An instance of int32x4
  * @param {float64x2} trueValue Pick lane from here if corresponding
  * selector lanes are 0xFFFFFFFF
  * @param {float64x2} falseValue Pick lane from here if corresponding
  * selector lanes are 0x0
  * @return {float64x2} Mix of lanes from trueValue or falseValue as
  * indicated
  */
SIMD.float64x2.select = function(t, trueValue, falseValue) {
  checkInt32x4(t);
  checkFloat64x2(trueValue);
  checkFloat64x2(falseValue);
  var tv = SIMD.int32x4.fromFloat64x2Bits(trueValue);
  var fv = SIMD.int32x4.fromFloat64x2Bits(falseValue);
  var tr = SIMD.int32x4.and(t, tv);
  var fr = SIMD.int32x4.and(SIMD.int32x4.not(t), fv);
  return SIMD.float64x2.fromInt32x4Bits(SIMD.int32x4.or(tr, fr));
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @return {float64x2} New instance of float64x2.
  */
SIMD.float64x2.load = function(buffer, byteOffset) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 16) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  var i8a = new Int8Array(buffer, byteOffset, 16);
  for (var i = 0; i < 16; i++)
    _PRIVATE._i8x16[i] = i8a[i];
  return SIMD.float64x2(_PRIVATE._f64x2[0], _PRIVATE._f64x2[1]);
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @return {float64x2} New instance of float64x2.
  */
SIMD.float64x2.loadX = function(buffer, byteOffset) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 8) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  var i8a = new Int8Array(buffer, byteOffset, 8);
  for (var i = 0; i < 8; i++)
    _PRIVATE._i8x16[i] = i8a[i];
  return SIMD.float64x2(_PRIVATE._f64x2[0], 0.0);
}

/**
  * @param {ArrayBuffer} f64a An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @param {float64x2} value An instance of float64x2.
  * @return {void}
  */
SIMD.float64x2.store = function(buffer, byteOffset, value) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 16) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkFloat64x2(value);
  _PRIVATE._f64x2[0] = value.x;
  _PRIVATE._f64x2[1] = value.y;
  var i8a = new Int8Array(buffer, byteOffset, 16);
  for (var i = 0; i < 16; i++)
    i8a[i] = _PRIVATE._i8x16[i];
}

/**
  * @param {ArrayBuffer} f64a An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @param {float64x2} value An instance of float64x2.
  * @return {void}
  */
SIMD.float64x2.storeX = function(buffer, byteOffset, value) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 8) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkFloat64x2(value);
  _PRIVATE._f64x2[0] = value.x;
  var i8a = new Int8Array(buffer, byteOffset, 8);
  for (var i = 0; i < 8; i++)
    i8a[i] = _PRIVATE._i8x16[i];
}

/**
  * @param {int32x4} a An instance of int32x4.
  * @param {int32x4} b An instance of int32x4.
  * @return {int32x4} New instance of int32x4 with values of a & b.
  */
SIMD.int32x4.and = function(a, b) {
  checkInt32x4(a);
  checkInt32x4(b);
  return SIMD.int32x4(a.x & b.x, a.y & b.y, a.z & b.z, a.w & b.w);
}

/**
  * @param {int32x4} a An instance of int32x4.
  * @param {int32x4} b An instance of int32x4.
  * @return {int32x4} New instance of int32x4 with values of a | b.
  */
SIMD.int32x4.or = function(a, b) {
  checkInt32x4(a);
  checkInt32x4(b);
  return SIMD.int32x4(a.x | b.x, a.y | b.y, a.z | b.z, a.w | b.w);
}

/**
  * @param {int32x4} a An instance of int32x4.
  * @param {int32x4} b An instance of int32x4.
  * @return {int32x4} New instance of int32x4 with values of a ^ b.
  */
SIMD.int32x4.xor = function(a, b) {
  checkInt32x4(a);
  checkInt32x4(b);
  return SIMD.int32x4(a.x ^ b.x, a.y ^ b.y, a.z ^ b.z, a.w ^ b.w);
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @return {int32x4} New instance of int32x4 with values of ~t
  */
SIMD.int32x4.not = function(t) {
  checkInt32x4(t);
  return SIMD.int32x4(~t.x, ~t.y, ~t.z, ~t.w);
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @return {int32x4} New instance of int32x4 with values of -t
  */
SIMD.int32x4.neg = function(t) {
  checkInt32x4(t);
  return SIMD.int32x4(-t.x, -t.y, -t.z, -t.w);
}

/**
  * @param {int32x4} a An instance of int32x4.
  * @param {int32x4} b An instance of int32x4.
  * @return {int32x4} New instance of int32x4 with values of a + b.
  */
SIMD.int32x4.add = function(a, b) {
  checkInt32x4(a);
  checkInt32x4(b);
  return SIMD.int32x4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

/**
  * @param {int32x4} a An instance of int32x4.
  * @param {int32x4} b An instance of int32x4.
  * @return {int32x4} New instance of int32x4 with values of a - b.
  */
SIMD.int32x4.sub = function(a, b) {
  checkInt32x4(a);
  checkInt32x4(b);
  return SIMD.int32x4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

/**
  * @param {int32x4} a An instance of int32x4.
  * @param {int32x4} b An instance of int32x4.
  * @return {int32x4} New instance of int32x4 with values of a * b.
  */
SIMD.int32x4.mul = function(a, b) {
  checkInt32x4(a);
  checkInt32x4(b);
  return SIMD.int32x4(Math.imul(a.x, b.x), Math.imul(a.y, b.y),
                      Math.imul(a.z, b.z), Math.imul(a.w, b.w));
}

/**
  * @param {int32x4} t An instance of float32x4 to be shuffled.
  * @param {integer} mask One of the 256 shuffle masks, for example, SIMD.XXXX.
  * @return {int32x4} New instance of float32x4 with lanes shuffled.
  */
SIMD.int32x4.shuffle = function(t, mask) {
  checkInt32x4(t);
  var _x = (mask) & 0x3;
  var _y = (mask >> 2) & 0x3;
  var _z = (mask >> 4) & 0x3;
  var _w = (mask >> 6) & 0x3;
  var storage = _PRIVATE._i32x4;
  storage[0] = t.x_;
  storage[1] = t.y_;
  storage[2] = t.z_;
  storage[3] = t.w_;
  return SIMD.int32x4(storage[_x], storage[_y], storage[_z], storage[_w]);
}

/**
  * @param {int32x4} t1 An instance of float32x4 to be shuffled. XY lanes in result
  * @param {int32x4} t2 An instance of float32x4 to be shuffled. ZW lanes in result
  * @param {integer} mask One of the 256 shuffle masks, for example, SIMD.XXXX.
  * @return {int32x4} New instance of float32x4 with lanes shuffled.
  */
SIMD.int32x4.shuffleMix = function(t1, t2, mask) {
  checkInt32x4(t1);
  checkInt32x4(t2);
  var _x = (mask) & 0x3;
  var _y = (mask >> 2) & 0x3;
  var _z = (mask >> 4) & 0x3;
  var _w = (mask >> 6) & 0x3;
  var storage = _PRIVATE._i32x8;
  storage[0] = t1.x_;
  storage[1] = t1.y_;
  storage[2] = t1.z_;
  storage[3] = t1.w_;
  storage[4] = t2.x_;
  storage[5] = t2.y_;
  storage[6] = t2.z_;
  storage[7] = t2.w_;
  return SIMD.float32x4(storage[0 + _x], storage[0 + _y],
                        storage[4 + _z], storage[4 + _w]);
}

/**
  * @param {int32x4} t Selector mask. An instance of int32x4
  * @param {int32x4} trueValue Pick lane from here if corresponding
  * selector lane is 0xFFFFFFFF
  * @param {int32x4} falseValue Pick lane from here if corresponding
  * selector lane is 0x0
  * @return {int32x4} Mix of lanes from trueValue or falseValue as
  * indicated
  */
SIMD.int32x4.select = function(t, trueValue, falseValue) {
  checkInt32x4(t);
  checkInt32x4(trueValue);
  checkInt32x4(falseValue);
  var tr = SIMD.int32x4.and(t, trueValue);
  var fr = SIMD.int32x4.and(SIMD.int32x4.not(t), falseValue);
  return SIMD.int32x4.or(tr, fr);
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @param {integer} 32-bit value used for x lane.
  * @return {int32x4} New instance of int32x4 with the values in t and
  * x lane replaced with {x}.
  */
SIMD.int32x4.withX = function(t, x) {
  checkInt32x4(t);
  return SIMD.int32x4(x, t.y, t.z, t.w);
}

/**
  * param {int32x4} t An instance of int32x4.
  * @param {integer} 32-bit value used for y lane.
  * @return {int32x4} New instance of int32x4 with the values in t and
  * y lane replaced with {y}.
  */
SIMD.int32x4.withY = function(t, y) {
  checkInt32x4(t);
  return SIMD.int32x4(t.x, y, t.z, t.w);
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @param {integer} 32-bit value used for z lane.
  * @return {int32x4} New instance of int32x4 with the values in t and
  * z lane replaced with {z}.
  */
SIMD.int32x4.withZ = function(t, z) {
  checkInt32x4(t);
  return SIMD.int32x4(t.x, t.y, z, t.w);
}

/**
  * @param {integer} 32-bit value used for w lane.
  * @return {int32x4} New instance of int32x4 with the values in t and
  * w lane replaced with {w}.
  */
SIMD.int32x4.withW = function(t, w) {
  checkInt32x4(t);
  return SIMD.int32x4(t.x, t.y, t.z, w);
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @param {boolean} x flag used for x lane.
  * @return {int32x4} New instance of int32x4 with the values in t and
  * x lane replaced with {x}.
  */
SIMD.int32x4.withFlagX = function(t, flagX) {
  checkInt32x4(t);
  var x = flagX ? 0xFFFFFFFF : 0x0;
  return SIMD.int32x4(x, t.y, t.z, t.w);
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @param {boolean} y flag used for y lane.
  * @return {int32x4} New instance of int32x4 with the values in t and
  * y lane replaced with {y}.
  */
SIMD.int32x4.withFlagY = function(t, flagY) {
  checkInt32x4(t);
  var y = flagY ? 0xFFFFFFFF : 0x0;
  return SIMD.int32x4(t.x, y, t.z, t.w);
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @param {boolean} z flag used for z lane.
  * @return {int32x4} New instance of int32x4 with the values in t and
  * z lane replaced with {z}.
  */
SIMD.int32x4.withFlagZ = function(t, flagZ) {
  checkInt32x4(t);
  var z = flagZ ? 0xFFFFFFFF : 0x0;
  return SIMD.int32x4(t.x, t.y, z, t.w);
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @param {boolean} w flag used for w lane.
  * @return {int32x4} New instance of int32x4 with the values in t and
  * w lane replaced with {w}.
  */
SIMD.int32x4.withFlagW = function(t, flagW) {
  checkInt32x4(t);
  var w = flagW ? 0xFFFFFFFF : 0x0;
  return SIMD.int32x4(t.x, t.y, t.z, w);
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @param {int32x4} other An instance of int32x4.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t == other.
  */
SIMD.int32x4.equal = function(t, other) {
  checkInt32x4(t);
  checkInt32x4(other);
  var cx = t.x == other.x;
  var cy = t.y == other.y;
  var cz = t.z == other.z;
  var cw = t.w == other.w;
  return SIMD.int32x4.bool(cx, cy, cz, cw);
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @param {int32x4} other An instance of int32x4.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t > other.
  */
SIMD.int32x4.greaterThan = function(t, other) {
  checkInt32x4(t);
  checkInt32x4(other);
  var cx = t.x > other.x;
  var cy = t.y > other.y;
  var cz = t.z > other.z;
  var cw = t.w > other.w;
  return SIMD.int32x4.bool(cx, cy, cz, cw);
}

/**
  * @param {int32x4} t An instance of int32x4.
  * @param {int32x4} other An instance of int32x4.
  * @return {int32x4} 0xFFFFFFFF or 0x0 in each lane depending on
  * the result of t < other.
  */
SIMD.int32x4.lessThan = function(t, other) {
  checkInt32x4(t);
  checkInt32x4(other);
  var cx = t.x < other.x;
  var cy = t.y < other.y;
  var cz = t.z < other.z;
  var cw = t.w < other.w;
  return SIMD.int32x4.bool(cx, cy, cz, cw);
}

/**
  * @param {int32x4} a An instance of int32x4.
  * @param {int} bits Bit count to shift by.
  * @return {int32x4} lanes in a shifted by bits.
  */
SIMD.int32x4.shiftLeft = function(a, bits) {
  checkInt32x4(a);
  var x = a.x << bits;
  var y = a.y << bits;
  var z = a.z << bits;
  var w = a.w << bits;
  return SIMD.int32x4(x, y, z, w);
}

/**
  * @param {int32x4} a An instance of int32x4.
  * @param {int} bits Bit count to shift by.
  * @return {int32x4} lanes in a shifted by bits.
  */
SIMD.int32x4.shiftRightLogical = function(a, bits) {
  checkInt32x4(a);
  var x = a.x >>> bits;
  var y = a.y >>> bits;
  var z = a.z >>> bits;
  var w = a.w >>> bits;
  return SIMD.int32x4(x, y, z, w);
}

/**
  * @param {int32x4} a An instance of int32x4.
  * @param {int} bits Bit count to shift by.
  * @return {int32x4} lanes in a shifted by bits.
  */
SIMD.int32x4.shiftRightArithmetic = function(a, bits) {
  checkInt32x4(a);
  var x = a.x >> bits;
  var y = a.y >> bits;
  var z = a.z >> bits;
  var w = a.w >> bits;
  return SIMD.int32x4(x, y, z, w);
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @return {int32x4} New instance of int32x4.
  */
SIMD.int32x4.load = function(buffer, byteOffset) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 16) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  var i8a = new Int8Array(buffer, byteOffset, 16);
  for (var i = 0; i < 16; i++)
    _PRIVATE._i8x16[i] = i8a[i];
  return SIMD.int32x4(_PRIVATE._i32x4[0], _PRIVATE._i32x4[1],
                      _PRIVATE._i32x4[2], _PRIVATE._i32x4[3]);
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @return {int32x4} New instance of int32x4.
  */
SIMD.int32x4.loadX = function(buffer, byteOffset) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 4) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  var i8a = new Int8Array(buffer, byteOffset, 4);
  for (var i = 0; i < 4; i++)
    _PRIVATE._i8x16[i] = i8a[i];
  return SIMD.int32x4(_PRIVATE._i32x4[0], 0, 0, 0);
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @return {int32x4} New instance of int32x4.
  */
SIMD.int32x4.loadXY = function(buffer, byteOffset) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 8) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  var i8a = new Int8Array(buffer, byteOffset, 8);
  for (var i = 0; i < 8; i++)
    _PRIVATE._i8x16[i] = i8a[i];
  return SIMD.int32x4(_PRIVATE._i32x4[0], _PRIVATE._i32x4[1], 0, 0);
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @return {int32x4} New instance of int32x4.
  */
SIMD.int32x4.loadXYZ = function(buffer, byteOffset) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 12) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  var i8a = new Int8Array(buffer, byteOffset, 12);
  for (var i = 0; i < 12; i++)
    _PRIVATE._i8x16[i] = i8a[i];
  return SIMD.int32x4(_PRIVATE._i32x4[0], _PRIVATE._i32x4[1],
                      _PRIVATE._i32x4[2], 0);
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @param {int32x4} value An instance of int32x4.
  * @return {void}
  */
SIMD.int32x4.store = function(buffer, byteOffset, value) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 16) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkInt32x4(value);
  _PRIVATE._i32x4[0] = value.x;
  _PRIVATE._i32x4[1] = value.y;
  _PRIVATE._i32x4[2] = value.z;
  _PRIVATE._i32x4[3] = value.w;
  var i8a = new Int8Array(buffer, byteOffset, 16);
  for (var i = 0; i < 16; i++)
    i8a[i] = _PRIVATE._i8x16[i];
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @param {int32x4} value An instance of int32x4.
  * @return {void}
  */
SIMD.int32x4.storeX = function(buffer, byteOffset, value) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 4) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkInt32x4(value);
  _PRIVATE._i32x4[0] = value.x;
  var i8a = new Int8Array(buffer, byteOffset, 4);
  for (var i = 0; i < 4; i++)
    i8a[i] = _PRIVATE._i8x16[i];
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @param {int32x4} value An instance of int32x4.
  * @return {void}
  */
SIMD.int32x4.storeXY = function(buffer, byteOffset, value) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 8) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkInt32x4(value);
  _PRIVATE._i32x4[0] = value.x;
  _PRIVATE._i32x4[1] = value.y;
  var i8a = new Int8Array(buffer, byteOffset, 8);
  for (var i = 0; i < 8; i++)
    i8a[i] = _PRIVATE._i8x16[i];
}

/**
  * @param {ArrayBuffer} buffer An instance of ArrayBuffer.
  * @param {Number} byteOffset An instance of Number.
  * @param {int32x4} value An instance of int32x4.
  * @return {void}
  */
SIMD.int32x4.storeXYZ = function(buffer, byteOffset, value) {
  if (!isArrayBuffer(buffer))
    throw new TypeError("The 1st argument must be an ArrayBuffer.");
  if (!isNumber(byteOffset))
    throw new TypeError("The 2nd argument must be a Number.");
  if (byteOffset < 0 || (byteOffset + 12) > buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkInt32x4(value);
  _PRIVATE._i32x4[0] = value.x;
  _PRIVATE._i32x4[1] = value.y;
  _PRIVATE._i32x4[2] = value.z;
  var i8a = new Int8Array(buffer, byteOffset, 12);
  for (var i = 0; i < 12; i++)
    i8a[i] = _PRIVATE._i8x16[i];
}

Object.defineProperty(SIMD, 'XXXX', { get: function() { return 0x0; } });
Object.defineProperty(SIMD, 'XXXY', { get: function() { return 0x40; } });
Object.defineProperty(SIMD, 'XXXZ', { get: function() { return 0x80; } });
Object.defineProperty(SIMD, 'XXXW', { get: function() { return 0xC0; } });
Object.defineProperty(SIMD, 'XXYX', { get: function() { return 0x10; } });
Object.defineProperty(SIMD, 'XXYY', { get: function() { return 0x50; } });
Object.defineProperty(SIMD, 'XXYZ', { get: function() { return 0x90; } });
Object.defineProperty(SIMD, 'XXYW', { get: function() { return 0xD0; } });
Object.defineProperty(SIMD, 'XXZX', { get: function() { return 0x20; } });
Object.defineProperty(SIMD, 'XXZY', { get: function() { return 0x60; } });
Object.defineProperty(SIMD, 'XXZZ', { get: function() { return 0xA0; } });
Object.defineProperty(SIMD, 'XXZW', { get: function() { return 0xE0; } });
Object.defineProperty(SIMD, 'XXWX', { get: function() { return 0x30; } });
Object.defineProperty(SIMD, 'XXWY', { get: function() { return 0x70; } });
Object.defineProperty(SIMD, 'XXWZ', { get: function() { return 0xB0; } });
Object.defineProperty(SIMD, 'XXWW', { get: function() { return 0xF0; } });
Object.defineProperty(SIMD, 'XYXX', { get: function() { return 0x4; } });
Object.defineProperty(SIMD, 'XYXY', { get: function() { return 0x44; } });
Object.defineProperty(SIMD, 'XYXZ', { get: function() { return 0x84; } });
Object.defineProperty(SIMD, 'XYXW', { get: function() { return 0xC4; } });
Object.defineProperty(SIMD, 'XYYX', { get: function() { return 0x14; } });
Object.defineProperty(SIMD, 'XYYY', { get: function() { return 0x54; } });
Object.defineProperty(SIMD, 'XYYZ', { get: function() { return 0x94; } });
Object.defineProperty(SIMD, 'XYYW', { get: function() { return 0xD4; } });
Object.defineProperty(SIMD, 'XYZX', { get: function() { return 0x24; } });
Object.defineProperty(SIMD, 'XYZY', { get: function() { return 0x64; } });
Object.defineProperty(SIMD, 'XYZZ', { get: function() { return 0xA4; } });
Object.defineProperty(SIMD, 'XYZW', { get: function() { return 0xE4; } });
Object.defineProperty(SIMD, 'XYWX', { get: function() { return 0x34; } });
Object.defineProperty(SIMD, 'XYWY', { get: function() { return 0x74; } });
Object.defineProperty(SIMD, 'XYWZ', { get: function() { return 0xB4; } });
Object.defineProperty(SIMD, 'XYWW', { get: function() { return 0xF4; } });
Object.defineProperty(SIMD, 'XZXX', { get: function() { return 0x8; } });
Object.defineProperty(SIMD, 'XZXY', { get: function() { return 0x48; } });
Object.defineProperty(SIMD, 'XZXZ', { get: function() { return 0x88; } });
Object.defineProperty(SIMD, 'XZXW', { get: function() { return 0xC8; } });
Object.defineProperty(SIMD, 'XZYX', { get: function() { return 0x18; } });
Object.defineProperty(SIMD, 'XZYY', { get: function() { return 0x58; } });
Object.defineProperty(SIMD, 'XZYZ', { get: function() { return 0x98; } });
Object.defineProperty(SIMD, 'XZYW', { get: function() { return 0xD8; } });
Object.defineProperty(SIMD, 'XZZX', { get: function() { return 0x28; } });
Object.defineProperty(SIMD, 'XZZY', { get: function() { return 0x68; } });
Object.defineProperty(SIMD, 'XZZZ', { get: function() { return 0xA8; } });
Object.defineProperty(SIMD, 'XZZW', { get: function() { return 0xE8; } });
Object.defineProperty(SIMD, 'XZWX', { get: function() { return 0x38; } });
Object.defineProperty(SIMD, 'XZWY', { get: function() { return 0x78; } });
Object.defineProperty(SIMD, 'XZWZ', { get: function() { return 0xB8; } });
Object.defineProperty(SIMD, 'XZWW', { get: function() { return 0xF8; } });
Object.defineProperty(SIMD, 'XWXX', { get: function() { return 0xC; } });
Object.defineProperty(SIMD, 'XWXY', { get: function() { return 0x4C; } });
Object.defineProperty(SIMD, 'XWXZ', { get: function() { return 0x8C; } });
Object.defineProperty(SIMD, 'XWXW', { get: function() { return 0xCC; } });
Object.defineProperty(SIMD, 'XWYX', { get: function() { return 0x1C; } });
Object.defineProperty(SIMD, 'XWYY', { get: function() { return 0x5C; } });
Object.defineProperty(SIMD, 'XWYZ', { get: function() { return 0x9C; } });
Object.defineProperty(SIMD, 'XWYW', { get: function() { return 0xDC; } });
Object.defineProperty(SIMD, 'XWZX', { get: function() { return 0x2C; } });
Object.defineProperty(SIMD, 'XWZY', { get: function() { return 0x6C; } });
Object.defineProperty(SIMD, 'XWZZ', { get: function() { return 0xAC; } });
Object.defineProperty(SIMD, 'XWZW', { get: function() { return 0xEC; } });
Object.defineProperty(SIMD, 'XWWX', { get: function() { return 0x3C; } });
Object.defineProperty(SIMD, 'XWWY', { get: function() { return 0x7C; } });
Object.defineProperty(SIMD, 'XWWZ', { get: function() { return 0xBC; } });
Object.defineProperty(SIMD, 'XWWW', { get: function() { return 0xFC; } });
Object.defineProperty(SIMD, 'YXXX', { get: function() { return 0x1; } });
Object.defineProperty(SIMD, 'YXXY', { get: function() { return 0x41; } });
Object.defineProperty(SIMD, 'YXXZ', { get: function() { return 0x81; } });
Object.defineProperty(SIMD, 'YXXW', { get: function() { return 0xC1; } });
Object.defineProperty(SIMD, 'YXYX', { get: function() { return 0x11; } });
Object.defineProperty(SIMD, 'YXYY', { get: function() { return 0x51; } });
Object.defineProperty(SIMD, 'YXYZ', { get: function() { return 0x91; } });
Object.defineProperty(SIMD, 'YXYW', { get: function() { return 0xD1; } });
Object.defineProperty(SIMD, 'YXZX', { get: function() { return 0x21; } });
Object.defineProperty(SIMD, 'YXZY', { get: function() { return 0x61; } });
Object.defineProperty(SIMD, 'YXZZ', { get: function() { return 0xA1; } });
Object.defineProperty(SIMD, 'YXZW', { get: function() { return 0xE1; } });
Object.defineProperty(SIMD, 'YXWX', { get: function() { return 0x31; } });
Object.defineProperty(SIMD, 'YXWY', { get: function() { return 0x71; } });
Object.defineProperty(SIMD, 'YXWZ', { get: function() { return 0xB1; } });
Object.defineProperty(SIMD, 'YXWW', { get: function() { return 0xF1; } });
Object.defineProperty(SIMD, 'YYXX', { get: function() { return 0x5; } });
Object.defineProperty(SIMD, 'YYXY', { get: function() { return 0x45; } });
Object.defineProperty(SIMD, 'YYXZ', { get: function() { return 0x85; } });
Object.defineProperty(SIMD, 'YYXW', { get: function() { return 0xC5; } });
Object.defineProperty(SIMD, 'YYYX', { get: function() { return 0x15; } });
Object.defineProperty(SIMD, 'YYYY', { get: function() { return 0x55; } });
Object.defineProperty(SIMD, 'YYYZ', { get: function() { return 0x95; } });
Object.defineProperty(SIMD, 'YYYW', { get: function() { return 0xD5; } });
Object.defineProperty(SIMD, 'YYZX', { get: function() { return 0x25; } });
Object.defineProperty(SIMD, 'YYZY', { get: function() { return 0x65; } });
Object.defineProperty(SIMD, 'YYZZ', { get: function() { return 0xA5; } });
Object.defineProperty(SIMD, 'YYZW', { get: function() { return 0xE5; } });
Object.defineProperty(SIMD, 'YYWX', { get: function() { return 0x35; } });
Object.defineProperty(SIMD, 'YYWY', { get: function() { return 0x75; } });
Object.defineProperty(SIMD, 'YYWZ', { get: function() { return 0xB5; } });
Object.defineProperty(SIMD, 'YYWW', { get: function() { return 0xF5; } });
Object.defineProperty(SIMD, 'YZXX', { get: function() { return 0x9; } });
Object.defineProperty(SIMD, 'YZXY', { get: function() { return 0x49; } });
Object.defineProperty(SIMD, 'YZXZ', { get: function() { return 0x89; } });
Object.defineProperty(SIMD, 'YZXW', { get: function() { return 0xC9; } });
Object.defineProperty(SIMD, 'YZYX', { get: function() { return 0x19; } });
Object.defineProperty(SIMD, 'YZYY', { get: function() { return 0x59; } });
Object.defineProperty(SIMD, 'YZYZ', { get: function() { return 0x99; } });
Object.defineProperty(SIMD, 'YZYW', { get: function() { return 0xD9; } });
Object.defineProperty(SIMD, 'YZZX', { get: function() { return 0x29; } });
Object.defineProperty(SIMD, 'YZZY', { get: function() { return 0x69; } });
Object.defineProperty(SIMD, 'YZZZ', { get: function() { return 0xA9; } });
Object.defineProperty(SIMD, 'YZZW', { get: function() { return 0xE9; } });
Object.defineProperty(SIMD, 'YZWX', { get: function() { return 0x39; } });
Object.defineProperty(SIMD, 'YZWY', { get: function() { return 0x79; } });
Object.defineProperty(SIMD, 'YZWZ', { get: function() { return 0xB9; } });
Object.defineProperty(SIMD, 'YZWW', { get: function() { return 0xF9; } });
Object.defineProperty(SIMD, 'YWXX', { get: function() { return 0xD; } });
Object.defineProperty(SIMD, 'YWXY', { get: function() { return 0x4D; } });
Object.defineProperty(SIMD, 'YWXZ', { get: function() { return 0x8D; } });
Object.defineProperty(SIMD, 'YWXW', { get: function() { return 0xCD; } });
Object.defineProperty(SIMD, 'YWYX', { get: function() { return 0x1D; } });
Object.defineProperty(SIMD, 'YWYY', { get: function() { return 0x5D; } });
Object.defineProperty(SIMD, 'YWYZ', { get: function() { return 0x9D; } });
Object.defineProperty(SIMD, 'YWYW', { get: function() { return 0xDD; } });
Object.defineProperty(SIMD, 'YWZX', { get: function() { return 0x2D; } });
Object.defineProperty(SIMD, 'YWZY', { get: function() { return 0x6D; } });
Object.defineProperty(SIMD, 'YWZZ', { get: function() { return 0xAD; } });
Object.defineProperty(SIMD, 'YWZW', { get: function() { return 0xED; } });
Object.defineProperty(SIMD, 'YWWX', { get: function() { return 0x3D; } });
Object.defineProperty(SIMD, 'YWWY', { get: function() { return 0x7D; } });
Object.defineProperty(SIMD, 'YWWZ', { get: function() { return 0xBD; } });
Object.defineProperty(SIMD, 'YWWW', { get: function() { return 0xFD; } });
Object.defineProperty(SIMD, 'ZXXX', { get: function() { return 0x2; } });
Object.defineProperty(SIMD, 'ZXXY', { get: function() { return 0x42; } });
Object.defineProperty(SIMD, 'ZXXZ', { get: function() { return 0x82; } });
Object.defineProperty(SIMD, 'ZXXW', { get: function() { return 0xC2; } });
Object.defineProperty(SIMD, 'ZXYX', { get: function() { return 0x12; } });
Object.defineProperty(SIMD, 'ZXYY', { get: function() { return 0x52; } });
Object.defineProperty(SIMD, 'ZXYZ', { get: function() { return 0x92; } });
Object.defineProperty(SIMD, 'ZXYW', { get: function() { return 0xD2; } });
Object.defineProperty(SIMD, 'ZXZX', { get: function() { return 0x22; } });
Object.defineProperty(SIMD, 'ZXZY', { get: function() { return 0x62; } });
Object.defineProperty(SIMD, 'ZXZZ', { get: function() { return 0xA2; } });
Object.defineProperty(SIMD, 'ZXZW', { get: function() { return 0xE2; } });
Object.defineProperty(SIMD, 'ZXWX', { get: function() { return 0x32; } });
Object.defineProperty(SIMD, 'ZXWY', { get: function() { return 0x72; } });
Object.defineProperty(SIMD, 'ZXWZ', { get: function() { return 0xB2; } });
Object.defineProperty(SIMD, 'ZXWW', { get: function() { return 0xF2; } });
Object.defineProperty(SIMD, 'ZYXX', { get: function() { return 0x6; } });
Object.defineProperty(SIMD, 'ZYXY', { get: function() { return 0x46; } });
Object.defineProperty(SIMD, 'ZYXZ', { get: function() { return 0x86; } });
Object.defineProperty(SIMD, 'ZYXW', { get: function() { return 0xC6; } });
Object.defineProperty(SIMD, 'ZYYX', { get: function() { return 0x16; } });
Object.defineProperty(SIMD, 'ZYYY', { get: function() { return 0x56; } });
Object.defineProperty(SIMD, 'ZYYZ', { get: function() { return 0x96; } });
Object.defineProperty(SIMD, 'ZYYW', { get: function() { return 0xD6; } });
Object.defineProperty(SIMD, 'ZYZX', { get: function() { return 0x26; } });
Object.defineProperty(SIMD, 'ZYZY', { get: function() { return 0x66; } });
Object.defineProperty(SIMD, 'ZYZZ', { get: function() { return 0xA6; } });
Object.defineProperty(SIMD, 'ZYZW', { get: function() { return 0xE6; } });
Object.defineProperty(SIMD, 'ZYWX', { get: function() { return 0x36; } });
Object.defineProperty(SIMD, 'ZYWY', { get: function() { return 0x76; } });
Object.defineProperty(SIMD, 'ZYWZ', { get: function() { return 0xB6; } });
Object.defineProperty(SIMD, 'ZYWW', { get: function() { return 0xF6; } });
Object.defineProperty(SIMD, 'ZZXX', { get: function() { return 0xA; } });
Object.defineProperty(SIMD, 'ZZXY', { get: function() { return 0x4A; } });
Object.defineProperty(SIMD, 'ZZXZ', { get: function() { return 0x8A; } });
Object.defineProperty(SIMD, 'ZZXW', { get: function() { return 0xCA; } });
Object.defineProperty(SIMD, 'ZZYX', { get: function() { return 0x1A; } });
Object.defineProperty(SIMD, 'ZZYY', { get: function() { return 0x5A; } });
Object.defineProperty(SIMD, 'ZZYZ', { get: function() { return 0x9A; } });
Object.defineProperty(SIMD, 'ZZYW', { get: function() { return 0xDA; } });
Object.defineProperty(SIMD, 'ZZZX', { get: function() { return 0x2A; } });
Object.defineProperty(SIMD, 'ZZZY', { get: function() { return 0x6A; } });
Object.defineProperty(SIMD, 'ZZZZ', { get: function() { return 0xAA; } });
Object.defineProperty(SIMD, 'ZZZW', { get: function() { return 0xEA; } });
Object.defineProperty(SIMD, 'ZZWX', { get: function() { return 0x3A; } });
Object.defineProperty(SIMD, 'ZZWY', { get: function() { return 0x7A; } });
Object.defineProperty(SIMD, 'ZZWZ', { get: function() { return 0xBA; } });
Object.defineProperty(SIMD, 'ZZWW', { get: function() { return 0xFA; } });
Object.defineProperty(SIMD, 'ZWXX', { get: function() { return 0xE; } });
Object.defineProperty(SIMD, 'ZWXY', { get: function() { return 0x4E; } });
Object.defineProperty(SIMD, 'ZWXZ', { get: function() { return 0x8E; } });
Object.defineProperty(SIMD, 'ZWXW', { get: function() { return 0xCE; } });
Object.defineProperty(SIMD, 'ZWYX', { get: function() { return 0x1E; } });
Object.defineProperty(SIMD, 'ZWYY', { get: function() { return 0x5E; } });
Object.defineProperty(SIMD, 'ZWYZ', { get: function() { return 0x9E; } });
Object.defineProperty(SIMD, 'ZWYW', { get: function() { return 0xDE; } });
Object.defineProperty(SIMD, 'ZWZX', { get: function() { return 0x2E; } });
Object.defineProperty(SIMD, 'ZWZY', { get: function() { return 0x6E; } });
Object.defineProperty(SIMD, 'ZWZZ', { get: function() { return 0xAE; } });
Object.defineProperty(SIMD, 'ZWZW', { get: function() { return 0xEE; } });
Object.defineProperty(SIMD, 'ZWWX', { get: function() { return 0x3E; } });
Object.defineProperty(SIMD, 'ZWWY', { get: function() { return 0x7E; } });
Object.defineProperty(SIMD, 'ZWWZ', { get: function() { return 0xBE; } });
Object.defineProperty(SIMD, 'ZWWW', { get: function() { return 0xFE; } });
Object.defineProperty(SIMD, 'WXXX', { get: function() { return 0x3; } });
Object.defineProperty(SIMD, 'WXXY', { get: function() { return 0x43; } });
Object.defineProperty(SIMD, 'WXXZ', { get: function() { return 0x83; } });
Object.defineProperty(SIMD, 'WXXW', { get: function() { return 0xC3; } });
Object.defineProperty(SIMD, 'WXYX', { get: function() { return 0x13; } });
Object.defineProperty(SIMD, 'WXYY', { get: function() { return 0x53; } });
Object.defineProperty(SIMD, 'WXYZ', { get: function() { return 0x93; } });
Object.defineProperty(SIMD, 'WXYW', { get: function() { return 0xD3; } });
Object.defineProperty(SIMD, 'WXZX', { get: function() { return 0x23; } });
Object.defineProperty(SIMD, 'WXZY', { get: function() { return 0x63; } });
Object.defineProperty(SIMD, 'WXZZ', { get: function() { return 0xA3; } });
Object.defineProperty(SIMD, 'WXZW', { get: function() { return 0xE3; } });
Object.defineProperty(SIMD, 'WXWX', { get: function() { return 0x33; } });
Object.defineProperty(SIMD, 'WXWY', { get: function() { return 0x73; } });
Object.defineProperty(SIMD, 'WXWZ', { get: function() { return 0xB3; } });
Object.defineProperty(SIMD, 'WXWW', { get: function() { return 0xF3; } });
Object.defineProperty(SIMD, 'WYXX', { get: function() { return 0x7; } });
Object.defineProperty(SIMD, 'WYXY', { get: function() { return 0x47; } });
Object.defineProperty(SIMD, 'WYXZ', { get: function() { return 0x87; } });
Object.defineProperty(SIMD, 'WYXW', { get: function() { return 0xC7; } });
Object.defineProperty(SIMD, 'WYYX', { get: function() { return 0x17; } });
Object.defineProperty(SIMD, 'WYYY', { get: function() { return 0x57; } });
Object.defineProperty(SIMD, 'WYYZ', { get: function() { return 0x97; } });
Object.defineProperty(SIMD, 'WYYW', { get: function() { return 0xD7; } });
Object.defineProperty(SIMD, 'WYZX', { get: function() { return 0x27; } });
Object.defineProperty(SIMD, 'WYZY', { get: function() { return 0x67; } });
Object.defineProperty(SIMD, 'WYZZ', { get: function() { return 0xA7; } });
Object.defineProperty(SIMD, 'WYZW', { get: function() { return 0xE7; } });
Object.defineProperty(SIMD, 'WYWX', { get: function() { return 0x37; } });
Object.defineProperty(SIMD, 'WYWY', { get: function() { return 0x77; } });
Object.defineProperty(SIMD, 'WYWZ', { get: function() { return 0xB7; } });
Object.defineProperty(SIMD, 'WYWW', { get: function() { return 0xF7; } });
Object.defineProperty(SIMD, 'WZXX', { get: function() { return 0xB; } });
Object.defineProperty(SIMD, 'WZXY', { get: function() { return 0x4B; } });
Object.defineProperty(SIMD, 'WZXZ', { get: function() { return 0x8B; } });
Object.defineProperty(SIMD, 'WZXW', { get: function() { return 0xCB; } });
Object.defineProperty(SIMD, 'WZYX', { get: function() { return 0x1B; } });
Object.defineProperty(SIMD, 'WZYY', { get: function() { return 0x5B; } });
Object.defineProperty(SIMD, 'WZYZ', { get: function() { return 0x9B; } });
Object.defineProperty(SIMD, 'WZYW', { get: function() { return 0xDB; } });
Object.defineProperty(SIMD, 'WZZX', { get: function() { return 0x2B; } });
Object.defineProperty(SIMD, 'WZZY', { get: function() { return 0x6B; } });
Object.defineProperty(SIMD, 'WZZZ', { get: function() { return 0xAB; } });
Object.defineProperty(SIMD, 'WZZW', { get: function() { return 0xEB; } });
Object.defineProperty(SIMD, 'WZWX', { get: function() { return 0x3B; } });
Object.defineProperty(SIMD, 'WZWY', { get: function() { return 0x7B; } });
Object.defineProperty(SIMD, 'WZWZ', { get: function() { return 0xBB; } });
Object.defineProperty(SIMD, 'WZWW', { get: function() { return 0xFB; } });
Object.defineProperty(SIMD, 'WWXX', { get: function() { return 0xF; } });
Object.defineProperty(SIMD, 'WWXY', { get: function() { return 0x4F; } });
Object.defineProperty(SIMD, 'WWXZ', { get: function() { return 0x8F; } });
Object.defineProperty(SIMD, 'WWXW', { get: function() { return 0xCF; } });
Object.defineProperty(SIMD, 'WWYX', { get: function() { return 0x1F; } });
Object.defineProperty(SIMD, 'WWYY', { get: function() { return 0x5F; } });
Object.defineProperty(SIMD, 'WWYZ', { get: function() { return 0x9F; } });
Object.defineProperty(SIMD, 'WWYW', { get: function() { return 0xDF; } });
Object.defineProperty(SIMD, 'WWZX', { get: function() { return 0x2F; } });
Object.defineProperty(SIMD, 'WWZY', { get: function() { return 0x6F; } });
Object.defineProperty(SIMD, 'WWZZ', { get: function() { return 0xAF; } });
Object.defineProperty(SIMD, 'WWZW', { get: function() { return 0xEF; } });
Object.defineProperty(SIMD, 'WWWX', { get: function() { return 0x3F; } });
Object.defineProperty(SIMD, 'WWWY', { get: function() { return 0x7F; } });
Object.defineProperty(SIMD, 'WWWZ', { get: function() { return 0xBF; } });
Object.defineProperty(SIMD, 'WWWW', { get: function() { return 0xFF; } });

Object.defineProperty(SIMD, 'XX',  { get: function() { return 0x0; } });
Object.defineProperty(SIMD, 'XY',  { get: function() { return 0x2; } });
Object.defineProperty(SIMD, 'YX',  { get: function() { return 0x1; } });
Object.defineProperty(SIMD, 'YY',  { get: function() { return 0x3; } });

Object.defineProperty(SIMD.float32x4.prototype, 'x', {
  get: function() { return this.x_; }
});

Object.defineProperty(SIMD.float32x4.prototype, 'y', {
  get: function() { return this.y_; }
});

Object.defineProperty(SIMD.float32x4.prototype, 'z', {
  get: function() { return this.z_; }
});

Object.defineProperty(SIMD.float32x4.prototype, 'w',
  { get: function() { return this.w_; }
});

/**
  * Extract the sign bit from each lane return them in the first 4 bits.
  */
Object.defineProperty(SIMD.float32x4.prototype, 'signMask', {
  get: function() {
    var mx = (this.x < 0.0 || 1/this.x === -Infinity) ? 1 : 0;
    var my = (this.y < 0.0 || 1/this.y === -Infinity) ? 1 : 0;
    var mz = (this.z < 0.0 || 1/this.z === -Infinity) ? 1 : 0;
    var mw = (this.w < 0.0 || 1/this.w === -Infinity) ? 1 : 0;
    return mx | my << 1 | mz << 2 | mw << 3;
  }
});

Object.defineProperty(SIMD.float64x2.prototype, 'x', {
  get: function() { return this.x_; }
});

Object.defineProperty(SIMD.float64x2.prototype, 'y', {
  get: function() { return this.y_; }
});

/**
  * Extract the sign bit from each lane return them in the first 2 bits.
  */
Object.defineProperty(SIMD.float64x2.prototype, 'signMask', {
  get: function() {
    var mx = (this.x < 0.0 || 1/this.x === -Infinity) ? 1 : 0;
    var my = (this.y < 0.0 || 1/this.y === -Infinity) ? 1 : 0;
    return mx | my << 1;
  }
});

Object.defineProperty(SIMD.int32x4.prototype, 'x', {
  get: function() { return this.x_; }
});

Object.defineProperty(SIMD.int32x4.prototype, 'y', {
  get: function() { return this.y_; }
});

Object.defineProperty(SIMD.int32x4.prototype, 'z', {
  get: function() { return this.z_; }
});

Object.defineProperty(SIMD.int32x4.prototype, 'w',
  { get: function() { return this.w_; }
});

Object.defineProperty(SIMD.int32x4.prototype, 'flagX', {
  get: function() { return this.x_ != 0x0; }
});

Object.defineProperty(SIMD.int32x4.prototype, 'flagY', {
  get: function() { return this.y_ != 0x0; }
});

Object.defineProperty(SIMD.int32x4.prototype, 'flagZ', {
  get: function() { return this.z_ != 0x0; }
});

Object.defineProperty(SIMD.int32x4.prototype, 'flagW',
  { get: function() { return this.w_ != 0x0; }
});

/**
  * Extract the sign bit from each lane return them in the first 4 bits.
  */
Object.defineProperty(SIMD.int32x4.prototype, 'signMask', {
  get: function() {
    var mx = (this.x_ & 0x80000000) >>> 31;
    var my = (this.y_ & 0x80000000) >>> 31;
    var mz = (this.z_ & 0x80000000) >>> 31;
    var mw = (this.w_ & 0x80000000) >>> 31;
    return mx | my << 1 | mz << 2 | mw << 3;
  }
});

function isNumber(o) {
    return typeof o == "number" || (typeof o == "object" && o.constructor === Number);
}

function isTypedArray(o) {
  return (o instanceof Int8Array) ||
         (o instanceof Uint8Array) ||
         (o instanceof Uint8ClampedArray) ||
         (o instanceof Int16Array) ||
         (o instanceof Uint16Array) ||
         (o instanceof Int32Array) ||
         (o instanceof Uint32Array) ||
         (o instanceof Float32Array) ||
         (o instanceof Float64Array) ||
         (o instanceof Float32x4Array);
}

function isArrayBuffer(o) {
  return (o instanceof ArrayBuffer);
}

function Float32x4Array(a, b, c) {
  if (isNumber(a)) {
    this.storage_ = new Float32Array(a*4);
    this.length_ = a;
    this.byteOffset_ = 0;
    return;
  } else if (isTypedArray(a)) {
    if (!(a instanceof Float32x4Array)) {
      throw "Copying typed array of non-Float32x4Array is unimplemented.";
    }
    this.storage_ = new Float32Array(a.length * 4);
    this.length_ = a.length;
    this.byteOffset_ = 0;
    // Copy floats.
    for (var i = 0; i < a.length*4; i++) {
      this.storage_[i] = a.storage_[i];
    }
  } else if (isArrayBuffer(a)) {
    if ((b != undefined) && (b % Float32x4Array.BYTES_PER_ELEMENT) != 0) {
      throw "byteOffset must be a multiple of 16.";
    }
    if (c != undefined) {
      c *= 4;
      this.storage_ = new Float32Array(a, b, c);
    }
    else {
      // Note = new Float32Array(a, b) is NOT equivalent to new Float32Array(a, b, undefined)
      this.storage_ = new Float32Array(a, b);
    }
    this.length_ = this.storage_.length / 4;
    this.byteOffset_ = b != undefined ? b : 0;
  } else {
    throw "Unknown type of first argument.";
  }
}

Object.defineProperty(Float32x4Array.prototype, 'length',
  { get: function() { return this.length_; }
});

Object.defineProperty(Float32x4Array.prototype, 'byteLength',
  { get: function() { return this.length_ * Float32x4Array.BYTES_PER_ELEMENT; }
});

Object.defineProperty(Float32x4Array, 'BYTES_PER_ELEMENT',
  { get: function() { return 16; }
});

Object.defineProperty(Float32x4Array.prototype, 'BYTES_PER_ELEMENT',
  { get: function() { return 16; }
});

Object.defineProperty(Float32x4Array.prototype, 'byteOffset',
  { get: function() { return this.byteOffset_; }
});

Object.defineProperty(Float32x4Array.prototype, 'buffer',
  { get: function() { return this.storage_.buffer; }
});

Float32x4Array.prototype.getAt = function(i) {
  if (i < 0) {
    throw "Index must be >= 0.";
  }
  if (i >= this.length) {
    throw "Index out of bounds.";
  }
  var x = this.storage_[i*4+0];
  var y = this.storage_[i*4+1];
  var z = this.storage_[i*4+2];
  var w = this.storage_[i*4+3];
  return SIMD.float32x4(x, y, z, w);
}

Float32x4Array.prototype.setAt = function(i, v) {
  if (i < 0) {
    throw "Index must be >= 0.";
  }
  if (i >= this.length) {
    throw "Index out of bounds.";
  }
  if (!(v instanceof SIMD.float32x4)) {
    throw "Value is not a float32x4.";
  }
  this.storage_[i*4+0] = v.x;
  this.storage_[i*4+1] = v.y;
  this.storage_[i*4+2] = v.z;
  this.storage_[i*4+3] = v.w;
}


function Int32x4Array(a, b, c) {

  function isNumber(o) {
      return typeof o == "number" || (typeof o == "object" && o.constructor === Number);
  }

  function isTypedArray(o) {
    return (o instanceof Int8Array) ||
           (o instanceof Uint8Array) ||
           (o instanceof Uint8ClampedArray) ||
           (o instanceof Int16Array) ||
           (o instanceof Uint16Array) ||
           (o instanceof Int32Array) ||
           (o instanceof Uint32Array) ||
           (o instanceof Float32Array) ||
           (o instanceof Float64Array) ||
           (o instanceof Int32x4Array) ||
           (o instanceof Float32x4Array);
  }

  function isArrayBuffer(o) {
    return (o instanceof ArrayBuffer);
  }

  if (isNumber(a)) {
    this.storage_ = new Int32Array(a*4);
    this.length_ = a;
    this.byteOffset_ = 0;
    return;
  } else if (isTypedArray(a)) {
    if (!(a instanceof Int32x4Array)) {
      throw "Copying typed array of non-Int32x4Array is unimplemented.";
    }
    this.storage_ = new Int32Array(a.length * 4);
    this.length_ = a.length;
    this.byteOffset_ = 0;
    // Copy ints.
    for (var i = 0; i < a.length*4; i++) {
      this.storage_[i] = a.storage_[i];
    }
  } else if (isArrayBuffer(a)) {
    if ((b != undefined) && (b % Int32x4Array.BYTES_PER_ELEMENT) != 0) {
      throw "byteOffset must be a multiple of 16.";
    }
    if (c != undefined) {
      c *= 4;
      this.storage_ = new Int32Array(a, b, c);
    }
    else {
      // Note = new Int32Array(a, b) is NOT equivalent to new Float32Array(a, b, undefined)
      this.storage_ = new Int32Array(a, b);
    }
    this.length_ = this.storage_.length / 4;
    this.byteOffset_ = b != undefined ? b : 0;
  } else {
    throw "Unknown type of first argument.";
  }
}

Object.defineProperty(Int32x4Array.prototype, 'length',
  { get: function() { return this.length_; }
});

Object.defineProperty(Int32x4Array.prototype, 'byteLength',
  { get: function() { return this.length_ * Int32x4Array.BYTES_PER_ELEMENT; }
});

Object.defineProperty(Int32x4Array, 'BYTES_PER_ELEMENT',
  { get: function() { return 16; }
});

Object.defineProperty(Int32x4Array.prototype, 'BYTES_PER_ELEMENT',
  { get: function() { return 16; }
});

Object.defineProperty(Int32x4Array.prototype, 'byteOffset',
  { get: function() { return this.byteOffset_; }
});

Object.defineProperty(Int32x4Array.prototype, 'buffer',
  { get: function() { return this.storage_.buffer; }
});

Int32x4Array.prototype.getAt = function(i) {
  if (i < 0) {
    throw "Index must be >= 0.";
  }
  if (i >= this.length) {
    throw "Index out of bounds.";
  }
  var x = this.storage_[i*4+0];
  var y = this.storage_[i*4+1];
  var z = this.storage_[i*4+2];
  var w = this.storage_[i*4+3];
  return SIMD.int32x4(x, y, z, w);
}

Int32x4Array.prototype.setAt = function(i, v) {
  if (i < 0) {
    throw "Index must be >= 0.";
  }
  if (i >= this.length) {
    throw "Index out of bounds.";
  }
  if (!(v instanceof SIMD.int32x4)) {
    throw "Value is not a int32x4.";
  }
  this.storage_[i*4+0] = v.x;
  this.storage_[i*4+1] = v.y;
  this.storage_[i*4+2] = v.z;
  this.storage_[i*4+3] = v.w;
}

function isDataView(v) {
  return v instanceof DataView;
}

DataView.prototype.getFloat32x4 = function(byteOffset, littleEndian) {
  if (!isDataView(this))
    throw new TypeError("This is not a DataView.");
  if (byteOffset < 0 || (byteOffset + 16) > this.buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  if (typeof littleEndian == 'undefined')
    littleEndian = false;
  return SIMD.float32x4(this.getFloat32(byteOffset, littleEndian),
                        this.getFloat32(byteOffset + 4, littleEndian),
                        this.getFloat32(byteOffset + 8, littleEndian),
                        this.getFloat32(byteOffset + 12, littleEndian));
}

DataView.prototype.getFloat64x2 = function(byteOffset, littleEndian) {
  if (!isDataView(this))
    throw new TypeError("This is not a DataView.");
  if (byteOffset < 0 || (byteOffset + 16) > this.buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  if (typeof littleEndian == 'undefined')
    littleEndian = false;
  return SIMD.float64x2(this.getFloat64(byteOffset, littleEndian),
                        this.getFloat64(byteOffset + 8, littleEndian));
}

DataView.prototype.getInt32x4 = function(byteOffset, littleEndian) {
  if (!isDataView(this))
    throw new TypeError("This is not a DataView.");
  if (byteOffset < 0 || (byteOffset + 16) > this.buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  if (typeof littleEndian == 'undefined')
    littleEndian = false;
  return SIMD.int32x4(this.getInt32(byteOffset, littleEndian),
                      this.getInt32(byteOffset + 4, littleEndian),
                      this.getInt32(byteOffset + 8, littleEndian),
                      this.getInt32(byteOffset + 12, littleEndian));
}

DataView.prototype.setFloat32x4 = function(byteOffset, value, littleEndian) {
  if (!isDataView(this))
    throw new TypeError("This is not a DataView.");
  if (byteOffset < 0 || (byteOffset + 16) > this.buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkFloat32x4(value);
  if (typeof littleEndian == 'undefined')
    littleEndian = false;
  this.setFloat32(byteOffset, value.x, littleEndian);
  this.setFloat32(byteOffset + 4, value.y, littleEndian);
  this.setFloat32(byteOffset + 8, value.z, littleEndian);
  this.setFloat32(byteOffset + 12, value.w, littleEndian);
}

DataView.prototype.setFloat64x2 = function(byteOffset, value, littleEndian) {
  if (!isDataView(this))
    throw new TypeError("This is not a DataView.");
  if (byteOffset < 0 || (byteOffset + 16) > this.buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkFloat64x2(value);
  if (typeof littleEndian == 'undefined')
    littleEndian = false;
  this.setFloat64(byteOffset, value.x, littleEndian);
  this.setFloat64(byteOffset + 8, value.y, littleEndian);
}

DataView.prototype.setInt32x4 = function(byteOffset, value, littleEndian) {
  if (!isDataView(this))
    throw new TypeError("This is not a DataView.");
  if (byteOffset < 0 || (byteOffset + 16) > this.buffer.byteLength)
    throw new RangeError("The value of byteOffset is invalid.");
  checkInt32x4(value);
  if (typeof littleEndian == 'undefined')
    littleEndian = false;
  this.setInt32(byteOffset, value.x, littleEndian);
  this.setInt32(byteOffset + 4, value.y, littleEndian);
  this.setInt32(byteOffset + 8, value.z, littleEndian);
  this.setInt32(byteOffset + 12, value.w, littleEndian);
}


// The Module object: Our interface to the outside world. We import
// and export values on it, and do the work to get that through
// closure compiler if necessary. There are various ways Module can be used:
// 1. Not defined. We create it here
// 2. A function parameter, function(Module) { ..generated code.. }
// 3. pre-run appended it, var Module = {}; ..generated code..
// 4. External script tag defines var Module.
// We need to do an eval in order to handle the closure compiler
// case, where this code here is minified but Module was defined
// elsewhere (e.g. case 4 above). We also need to check if Module
// already exists (e.g. case 3 above).
// Note that if you want to run closure, and also to use Module
// after the generated code, you will need to define   var Module = {};
// before the code. Then that object will be used in the code, and you
// can continue to use Module afterwards as well.
var Module;
if (!Module) Module = (typeof Module !== 'undefined' ? Module : null) || {};

// Sometimes an existing Module object exists with properties
// meant to overwrite the default module functionality. Here
// we collect those properties and reapply _after_ we configure
// the current environment's defaults to avoid having to be so
// defensive during initialization.
var moduleOverrides = {};
for (var key in Module) {
  if (Module.hasOwnProperty(key)) {
    moduleOverrides[key] = Module[key];
  }
}

// The environment setup code below is customized to use Module.
// *** Environment setup code ***
var ENVIRONMENT_IS_NODE = typeof process === 'object' && typeof require === 'function';
var ENVIRONMENT_IS_WEB = typeof window === 'object';
var ENVIRONMENT_IS_WORKER = typeof importScripts === 'function';
var ENVIRONMENT_IS_SHELL = !ENVIRONMENT_IS_WEB && !ENVIRONMENT_IS_NODE && !ENVIRONMENT_IS_WORKER;

if (ENVIRONMENT_IS_NODE) {
  // Expose functionality in the same simple way that the shells work
  // Note that we pollute the global namespace here, otherwise we break in node
  if (!Module['print']) Module['print'] = function print(x) {
    process['stdout'].write(x + '\n');
  };
  if (!Module['printErr']) Module['printErr'] = function printErr(x) {
    process['stderr'].write(x + '\n');
  };

  var nodeFS = require('fs');
  var nodePath = require('path');

  Module['read'] = function read(filename, binary) {
    filename = nodePath['normalize'](filename);
    var ret = nodeFS['readFileSync'](filename);
    // The path is absolute if the normalized version is the same as the resolved.
    if (!ret && filename != nodePath['resolve'](filename)) {
      filename = path.join(__dirname, '..', 'src', filename);
      ret = nodeFS['readFileSync'](filename);
    }
    if (ret && !binary) ret = ret.toString();
    return ret;
  };

  Module['readBinary'] = function readBinary(filename) { return Module['read'](filename, true) };

  Module['load'] = function load(f) {
    globalEval(read(f));
  };

  Module['thisProgram'] = process['argv'][1].replace(/\\/g, '/');
  Module['arguments'] = process['argv'].slice(2);

  if (typeof module !== 'undefined') {
    module['exports'] = Module;
  }

  process['on']('uncaughtException', function(ex) {
    // suppress ExitStatus exceptions from showing an error
    if (!(ex instanceof ExitStatus)) {
      throw ex;
    }
  });
}
else if (ENVIRONMENT_IS_SHELL) {
  if (!Module['print']) Module['print'] = print;
  if (typeof printErr != 'undefined') Module['printErr'] = printErr; // not present in v8 or older sm

  if (typeof read != 'undefined') {
    Module['read'] = read;
  } else {
    Module['read'] = function read() { throw 'no read() available (jsc?)' };
  }

  Module['readBinary'] = function readBinary(f) {
    if (typeof readbuffer === 'function') {
      return new Uint8Array(readbuffer(f));
    }
    var data = read(f, 'binary');
    assert(typeof data === 'object');
    return data;
  };

  if (typeof scriptArgs != 'undefined') {
    Module['arguments'] = scriptArgs;
  } else if (typeof arguments != 'undefined') {
    Module['arguments'] = arguments;
  }

  this['Module'] = Module;

}
else if (ENVIRONMENT_IS_WEB || ENVIRONMENT_IS_WORKER) {
  Module['read'] = function read(url) {
    var xhr = new XMLHttpRequest();
    xhr.open('GET', url, false);
    xhr.send(null);
    return xhr.responseText;
  };

  if (typeof arguments != 'undefined') {
    Module['arguments'] = arguments;
  }

  if (typeof console !== 'undefined') {
    if (!Module['print']) Module['print'] = function print(x) {
      console.log(x);
    };
    if (!Module['printErr']) Module['printErr'] = function printErr(x) {
      console.log(x);
    };
  } else {
    // Probably a worker, and without console.log. We can do very little here...
    var TRY_USE_DUMP = false;
    if (!Module['print']) Module['print'] = (TRY_USE_DUMP && (typeof(dump) !== "undefined") ? (function(x) {
      dump(x);
    }) : (function(x) {
      // self.postMessage(x); // enable this if you want stdout to be sent as messages
    }));
  }

  if (ENVIRONMENT_IS_WEB) {
    window['Module'] = Module;
  } else {
    Module['load'] = importScripts;
  }
}
else {
  // Unreachable because SHELL is dependant on the others
  throw 'Unknown runtime environment. Where are we?';
}

function globalEval(x) {
  eval.call(null, x);
}
if (!Module['load'] && Module['read']) {
  Module['load'] = function load(f) {
    globalEval(Module['read'](f));
  };
}
if (!Module['print']) {
  Module['print'] = function(){};
}
if (!Module['printErr']) {
  Module['printErr'] = Module['print'];
}
if (!Module['arguments']) {
  Module['arguments'] = [];
}
if (!Module['thisProgram']) {
  Module['thisProgram'] = './this.program';
}

// *** Environment setup code ***

// Closure helpers
Module.print = Module['print'];
Module.printErr = Module['printErr'];

// Callbacks
Module['preRun'] = [];
Module['postRun'] = [];

// Merge back in the overrides
for (var key in moduleOverrides) {
  if (moduleOverrides.hasOwnProperty(key)) {
    Module[key] = moduleOverrides[key];
  }
}



// === Preamble library stuff ===

// Documentation for the public APIs defined in this file must be updated in: 
//    site/source/docs/api_reference/preamble.js.rst
// A prebuilt local version of the documentation is available at: 
//    site/build/text/docs/api_reference/preamble.js.txt
// You can also build docs locally as HTML or other formats in site/
// An online HTML version (which may be of a different version of Emscripten)
//    is up at http://kripken.github.io/emscripten-site/docs/api_reference/preamble.js.html

//========================================
// Runtime code shared with compiler
//========================================

var Runtime = {
  setTempRet0: function (value) {
    tempRet0 = value;
  },
  getTempRet0: function () {
    return tempRet0;
  },
  stackSave: function () {
    return STACKTOP;
  },
  stackRestore: function (stackTop) {
    STACKTOP = stackTop;
  },
  getNativeTypeSize: function (type) {
    switch (type) {
      case 'i1': case 'i8': return 1;
      case 'i16': return 2;
      case 'i32': return 4;
      case 'i64': return 8;
      case 'float': return 4;
      case 'double': return 8;
      default: {
        if (type[type.length-1] === '*') {
          return Runtime.QUANTUM_SIZE; // A pointer
        } else if (type[0] === 'i') {
          var bits = parseInt(type.substr(1));
          assert(bits % 8 === 0);
          return bits/8;
        } else {
          return 0;
        }
      }
    }
  },
  getNativeFieldSize: function (type) {
    return Math.max(Runtime.getNativeTypeSize(type), Runtime.QUANTUM_SIZE);
  },
  STACK_ALIGN: 16,
  getAlignSize: function (type, size, vararg) {
    // we align i64s and doubles on 64-bit boundaries, unlike x86
    if (!vararg && (type == 'i64' || type == 'double')) return 8;
    if (!type) return Math.min(size, 8); // align structures internally to 64 bits
    return Math.min(size || (type ? Runtime.getNativeFieldSize(type) : 0), Runtime.QUANTUM_SIZE);
  },
  dynCall: function (sig, ptr, args) {
    if (args && args.length) {
      if (!args.splice) args = Array.prototype.slice.call(args);
      args.splice(0, 0, ptr);
      return Module['dynCall_' + sig].apply(null, args);
    } else {
      return Module['dynCall_' + sig].call(null, ptr);
    }
  },
  functionPointers: [],
  addFunction: function (func) {
    for (var i = 0; i < Runtime.functionPointers.length; i++) {
      if (!Runtime.functionPointers[i]) {
        Runtime.functionPointers[i] = func;
        return 2*(1 + i);
      }
    }
    throw 'Finished up all reserved function pointers. Use a higher value for RESERVED_FUNCTION_POINTERS.';
  },
  removeFunction: function (index) {
    Runtime.functionPointers[(index-2)/2] = null;
  },
  getAsmConst: function (code, numArgs) {
    // code is a constant string on the heap, so we can cache these
    if (!Runtime.asmConstCache) Runtime.asmConstCache = {};
    var func = Runtime.asmConstCache[code];
    if (func) return func;
    var args = [];
    for (var i = 0; i < numArgs; i++) {
      args.push(String.fromCharCode(36) + i); // $0, $1 etc
    }
    var source = Pointer_stringify(code);
    if (source[0] === '"') {
      // tolerate EM_ASM("..code..") even though EM_ASM(..code..) is correct
      if (source.indexOf('"', 1) === source.length-1) {
        source = source.substr(1, source.length-2);
      } else {
        // something invalid happened, e.g. EM_ASM("..code($0)..", input)
        abort('invalid EM_ASM input |' + source + '|. Please use EM_ASM(..code..) (no quotes) or EM_ASM({ ..code($0).. }, input) (to input values)');
      }
    }
    try {
      // Module is the only 'upvar', which we provide directly. We also provide FS for legacy support.
      var evalled = eval('(function(Module, FS) { return function(' + args.join(',') + '){ ' + source + ' } })')(Module, typeof FS !== 'undefined' ? FS : null);
    } catch(e) {
      Module.printErr('error in executing inline EM_ASM code: ' + e + ' on: \n\n' + source + '\n\nwith args |' + args + '| (make sure to use the right one out of EM_ASM, EM_ASM_ARGS, etc.)');
      throw e;
    }
    return Runtime.asmConstCache[code] = evalled;
  },
  warnOnce: function (text) {
    if (!Runtime.warnOnce.shown) Runtime.warnOnce.shown = {};
    if (!Runtime.warnOnce.shown[text]) {
      Runtime.warnOnce.shown[text] = 1;
      Module.printErr(text);
    }
  },
  funcWrappers: {},
  getFuncWrapper: function (func, sig) {
    assert(sig);
    if (!Runtime.funcWrappers[sig]) {
      Runtime.funcWrappers[sig] = {};
    }
    var sigCache = Runtime.funcWrappers[sig];
    if (!sigCache[func]) {
      sigCache[func] = function dynCall_wrapper() {
        return Runtime.dynCall(sig, func, arguments);
      };
    }
    return sigCache[func];
  },
  UTF8Processor: function () {
    var buffer = [];
    var needed = 0;
    this.processCChar = function (code) {
      code = code & 0xFF;

      if (buffer.length == 0) {
        if ((code & 0x80) == 0x00) {        // 0xxxxxxx
          return String.fromCharCode(code);
        }
        buffer.push(code);
        if ((code & 0xE0) == 0xC0) {        // 110xxxxx
          needed = 1;
        } else if ((code & 0xF0) == 0xE0) { // 1110xxxx
          needed = 2;
        } else {                            // 11110xxx
          needed = 3;
        }
        return '';
      }

      if (needed) {
        buffer.push(code);
        needed--;
        if (needed > 0) return '';
      }

      var c1 = buffer[0];
      var c2 = buffer[1];
      var c3 = buffer[2];
      var c4 = buffer[3];
      var ret;
      if (buffer.length == 2) {
        ret = String.fromCharCode(((c1 & 0x1F) << 6)  | (c2 & 0x3F));
      } else if (buffer.length == 3) {
        ret = String.fromCharCode(((c1 & 0x0F) << 12) | ((c2 & 0x3F) << 6)  | (c3 & 0x3F));
      } else {
        // http://mathiasbynens.be/notes/javascript-encoding#surrogate-formulae
        var codePoint = ((c1 & 0x07) << 18) | ((c2 & 0x3F) << 12) |
                        ((c3 & 0x3F) << 6)  | (c4 & 0x3F);
        ret = String.fromCharCode(
          (((codePoint - 0x10000) / 0x400)|0) + 0xD800,
          (codePoint - 0x10000) % 0x400 + 0xDC00);
      }
      buffer.length = 0;
      return ret;
    }
    this.processJSString = function processJSString(string) {
      /* TODO: use TextEncoder when present,
        var encoder = new TextEncoder();
        encoder['encoding'] = "utf-8";
        var utf8Array = encoder['encode'](aMsg.data);
      */
      string = unescape(encodeURIComponent(string));
      var ret = [];
      for (var i = 0; i < string.length; i++) {
        ret.push(string.charCodeAt(i));
      }
      return ret;
    }
  },
  getCompilerSetting: function (name) {
    throw 'You must build with -s RETAIN_COMPILER_SETTINGS=1 for Runtime.getCompilerSetting or emscripten_get_compiler_setting to work';
  },
  stackAlloc: function (size) { var ret = STACKTOP;STACKTOP = (STACKTOP + size)|0;STACKTOP = (((STACKTOP)+15)&-16); return ret; },
  staticAlloc: function (size) { var ret = STATICTOP;STATICTOP = (STATICTOP + size)|0;STATICTOP = (((STATICTOP)+15)&-16); return ret; },
  dynamicAlloc: function (size) { var ret = DYNAMICTOP;DYNAMICTOP = (DYNAMICTOP + size)|0;DYNAMICTOP = (((DYNAMICTOP)+15)&-16); if (DYNAMICTOP >= TOTAL_MEMORY) enlargeMemory();; return ret; },
  alignMemory: function (size,quantum) { var ret = size = Math.ceil((size)/(quantum ? quantum : 16))*(quantum ? quantum : 16); return ret; },
  makeBigInt: function (low,high,unsigned) { var ret = (unsigned ? ((+((low>>>0)))+((+((high>>>0)))*(+4294967296))) : ((+((low>>>0)))+((+((high|0)))*(+4294967296)))); return ret; },
  GLOBAL_BASE: 8,
  QUANTUM_SIZE: 4,
  __dummy__: 0
}


Module['Runtime'] = Runtime;









//========================================
// Runtime essentials
//========================================

var __THREW__ = 0; // Used in checking for thrown exceptions.

var ABORT = false; // whether we are quitting the application. no code should run after this. set in exit() and abort()
var EXITSTATUS = 0;

var undef = 0;
// tempInt is used for 32-bit signed values or smaller. tempBigInt is used
// for 32-bit unsigned values or more than 32 bits. TODO: audit all uses of tempInt
var tempValue, tempInt, tempBigInt, tempInt2, tempBigInt2, tempPair, tempBigIntI, tempBigIntR, tempBigIntS, tempBigIntP, tempBigIntD, tempDouble, tempFloat;
var tempI64, tempI64b;
var tempRet0, tempRet1, tempRet2, tempRet3, tempRet4, tempRet5, tempRet6, tempRet7, tempRet8, tempRet9;

function assert(condition, text) {
  if (!condition) {
    abort('Assertion failed: ' + text);
  }
}

var globalScope = this;

// Returns the C function with a specified identifier (for C++, you need to do manual name mangling)
function getCFunc(ident) {
  var func = Module['_' + ident]; // closure exported function
  if (!func) {
    try {
      func = eval('_' + ident); // explicit lookup
    } catch(e) {}
  }
  assert(func, 'Cannot call unknown function ' + ident + ' (perhaps LLVM optimizations or closure removed it?)');
  return func;
}

var cwrap, ccall;
(function(){
  var stack = 0;
  var JSfuncs = {
    'stackSave' : function() {
      stack = Runtime.stackSave();
    },
    'stackRestore' : function() {
      Runtime.stackRestore(stack);
    },
    // type conversion from js to c
    'arrayToC' : function(arr) {
      var ret = Runtime.stackAlloc(arr.length);
      writeArrayToMemory(arr, ret);
      return ret;
    },
    'stringToC' : function(str) {
      var ret = 0;
      if (str !== null && str !== undefined && str !== 0) { // null string
        // at most 4 bytes per UTF-8 code point, +1 for the trailing '\0'
        ret = Runtime.stackAlloc((str.length << 2) + 1);
        writeStringToMemory(str, ret);
      }
      return ret;
    }
  };
  // For fast lookup of conversion functions
  var toC = {'string' : JSfuncs['stringToC'], 'array' : JSfuncs['arrayToC']};

  // C calling interface. 
  ccall = function ccallFunc(ident, returnType, argTypes, args) {
    var func = getCFunc(ident);
    var cArgs = [];
    if (args) {
      for (var i = 0; i < args.length; i++) {
        var converter = toC[argTypes[i]];
        if (converter) {
          if (stack === 0) stack = Runtime.stackSave();
          cArgs[i] = converter(args[i]);
        } else {
          cArgs[i] = args[i];
        }
      }
    }
    var ret = func.apply(null, cArgs);
    if (returnType === 'string') ret = Pointer_stringify(ret);
    if (stack !== 0) JSfuncs['stackRestore']();
    return ret;
  }

  var sourceRegex = /^function\s*\(([^)]*)\)\s*{\s*([^*]*?)[\s;]*(?:return\s*(.*?)[;\s]*)?}$/;
  function parseJSFunc(jsfunc) {
    // Match the body and the return value of a javascript function source
    var parsed = jsfunc.toString().match(sourceRegex).slice(1);
    return {arguments : parsed[0], body : parsed[1], returnValue: parsed[2]}
  }
  var JSsource = {};
  for (var fun in JSfuncs) {
    if (JSfuncs.hasOwnProperty(fun)) {
      // Elements of toCsource are arrays of three items:
      // the code, and the return value
      JSsource[fun] = parseJSFunc(JSfuncs[fun]);
    }
  }

  
  cwrap = function cwrap(ident, returnType, argTypes) {
    argTypes = argTypes || [];
    var cfunc = getCFunc(ident);
    // When the function takes numbers and returns a number, we can just return
    // the original function
    var numericArgs = argTypes.every(function(type){ return type === 'number'});
    var numericRet = (returnType !== 'string');
    if ( numericRet && numericArgs) {
      return cfunc;
    }
    // Creation of the arguments list (["$1","$2",...,"$nargs"])
    var argNames = argTypes.map(function(x,i){return '$'+i});
    var funcstr = "(function(" + argNames.join(',') + ") {";
    var nargs = argTypes.length;
    if (!numericArgs) {
      // Generate the code needed to convert the arguments from javascript
      // values to pointers
      funcstr += JSsource['stackSave'].body + ';';
      for (var i = 0; i < nargs; i++) {
        var arg = argNames[i], type = argTypes[i];
        if (type === 'number') continue;
        var convertCode = JSsource[type + 'ToC']; // [code, return]
        funcstr += 'var ' + convertCode.arguments + ' = ' + arg + ';';
        funcstr += convertCode.body + ';';
        funcstr += arg + '=' + convertCode.returnValue + ';';
      }
    }

    // When the code is compressed, the name of cfunc is not literally 'cfunc' anymore
    var cfuncname = parseJSFunc(function(){return cfunc}).returnValue;
    // Call the function
    funcstr += 'var ret = ' + cfuncname + '(' + argNames.join(',') + ');';
    if (!numericRet) { // Return type can only by 'string' or 'number'
      // Convert the result to a string
      var strgfy = parseJSFunc(function(){return Pointer_stringify}).returnValue;
      funcstr += 'ret = ' + strgfy + '(ret);';
    }
    if (!numericArgs) {
      // If we had a stack, restore it
      funcstr += JSsource['stackRestore'].body + ';';
    }
    funcstr += 'return ret})';
    return eval(funcstr);
  };
})();
Module["cwrap"] = cwrap;
Module["ccall"] = ccall;


function setValue(ptr, value, type, noSafe) {
  type = type || 'i8';
  if (type.charAt(type.length-1) === '*') type = 'i32'; // pointers are 32-bit
    switch(type) {
      case 'i1': HEAP8[((ptr)>>0)]=value; break;
      case 'i8': HEAP8[((ptr)>>0)]=value; break;
      case 'i16': HEAP16[((ptr)>>1)]=value; break;
      case 'i32': HEAP32[((ptr)>>2)]=value; break;
      case 'i64': (tempI64 = [value>>>0,(tempDouble=value,(+(Math_abs(tempDouble))) >= (+1) ? (tempDouble > (+0) ? ((Math_min((+(Math_floor((tempDouble)/(+4294967296)))), (+4294967295)))|0)>>>0 : (~~((+(Math_ceil((tempDouble - +(((~~(tempDouble)))>>>0))/(+4294967296))))))>>>0) : 0)],HEAP32[((ptr)>>2)]=tempI64[0],HEAP32[(((ptr)+(4))>>2)]=tempI64[1]); break;
      case 'float': HEAPF32[((ptr)>>2)]=value; break;
      case 'double': HEAPF64[((ptr)>>3)]=value; break;
      default: abort('invalid type for setValue: ' + type);
    }
}
Module['setValue'] = setValue;


function getValue(ptr, type, noSafe) {
  type = type || 'i8';
  if (type.charAt(type.length-1) === '*') type = 'i32'; // pointers are 32-bit
    switch(type) {
      case 'i1': return HEAP8[((ptr)>>0)];
      case 'i8': return HEAP8[((ptr)>>0)];
      case 'i16': return HEAP16[((ptr)>>1)];
      case 'i32': return HEAP32[((ptr)>>2)];
      case 'i64': return HEAP32[((ptr)>>2)];
      case 'float': return HEAPF32[((ptr)>>2)];
      case 'double': return HEAPF64[((ptr)>>3)];
      default: abort('invalid type for setValue: ' + type);
    }
  return null;
}
Module['getValue'] = getValue;

var ALLOC_NORMAL = 0; // Tries to use _malloc()
var ALLOC_STACK = 1; // Lives for the duration of the current function call
var ALLOC_STATIC = 2; // Cannot be freed
var ALLOC_DYNAMIC = 3; // Cannot be freed except through sbrk
var ALLOC_NONE = 4; // Do not allocate
Module['ALLOC_NORMAL'] = ALLOC_NORMAL;
Module['ALLOC_STACK'] = ALLOC_STACK;
Module['ALLOC_STATIC'] = ALLOC_STATIC;
Module['ALLOC_DYNAMIC'] = ALLOC_DYNAMIC;
Module['ALLOC_NONE'] = ALLOC_NONE;

// allocate(): This is for internal use. You can use it yourself as well, but the interface
//             is a little tricky (see docs right below). The reason is that it is optimized
//             for multiple syntaxes to save space in generated code. So you should
//             normally not use allocate(), and instead allocate memory using _malloc(),
//             initialize it with setValue(), and so forth.
// @slab: An array of data, or a number. If a number, then the size of the block to allocate,
//        in *bytes* (note that this is sometimes confusing: the next parameter does not
//        affect this!)
// @types: Either an array of types, one for each byte (or 0 if no type at that position),
//         or a single type which is used for the entire block. This only matters if there
//         is initial data - if @slab is a number, then this does not matter at all and is
//         ignored.
// @allocator: How to allocate memory, see ALLOC_*
function allocate(slab, types, allocator, ptr) {
  var zeroinit, size;
  if (typeof slab === 'number') {
    zeroinit = true;
    size = slab;
  } else {
    zeroinit = false;
    size = slab.length;
  }

  var singleType = typeof types === 'string' ? types : null;

  var ret;
  if (allocator == ALLOC_NONE) {
    ret = ptr;
  } else {
    ret = [_malloc, Runtime.stackAlloc, Runtime.staticAlloc, Runtime.dynamicAlloc][allocator === undefined ? ALLOC_STATIC : allocator](Math.max(size, singleType ? 1 : types.length));
  }

  if (zeroinit) {
    var ptr = ret, stop;
    assert((ret & 3) == 0);
    stop = ret + (size & ~3);
    for (; ptr < stop; ptr += 4) {
      HEAP32[((ptr)>>2)]=0;
    }
    stop = ret + size;
    while (ptr < stop) {
      HEAP8[((ptr++)>>0)]=0;
    }
    return ret;
  }

  if (singleType === 'i8') {
    if (slab.subarray || slab.slice) {
      HEAPU8.set(slab, ret);
    } else {
      HEAPU8.set(new Uint8Array(slab), ret);
    }
    return ret;
  }

  var i = 0, type, typeSize, previousType;
  while (i < size) {
    var curr = slab[i];

    if (typeof curr === 'function') {
      curr = Runtime.getFunctionIndex(curr);
    }

    type = singleType || types[i];
    if (type === 0) {
      i++;
      continue;
    }

    if (type == 'i64') type = 'i32'; // special case: we have one i32 here, and one i32 later

    setValue(ret+i, curr, type);

    // no need to look up size unless type changes, so cache it
    if (previousType !== type) {
      typeSize = Runtime.getNativeTypeSize(type);
      previousType = type;
    }
    i += typeSize;
  }

  return ret;
}
Module['allocate'] = allocate;

function Pointer_stringify(ptr, /* optional */ length) {
  // TODO: use TextDecoder
  // Find the length, and check for UTF while doing so
  var hasUtf = false;
  var t;
  var i = 0;
  while (1) {
    t = HEAPU8[(((ptr)+(i))>>0)];
    if (t >= 128) hasUtf = true;
    else if (t == 0 && !length) break;
    i++;
    if (length && i == length) break;
  }
  if (!length) length = i;

  var ret = '';

  if (!hasUtf) {
    var MAX_CHUNK = 1024; // split up into chunks, because .apply on a huge string can overflow the stack
    var curr;
    while (length > 0) {
      curr = String.fromCharCode.apply(String, HEAPU8.subarray(ptr, ptr + Math.min(length, MAX_CHUNK)));
      ret = ret ? ret + curr : curr;
      ptr += MAX_CHUNK;
      length -= MAX_CHUNK;
    }
    return ret;
  }

  var utf8 = new Runtime.UTF8Processor();
  for (i = 0; i < length; i++) {
    t = HEAPU8[(((ptr)+(i))>>0)];
    ret += utf8.processCChar(t);
  }
  return ret;
}
Module['Pointer_stringify'] = Pointer_stringify;

function UTF16ToString(ptr) {
  var i = 0;

  var str = '';
  while (1) {
    var codeUnit = HEAP16[(((ptr)+(i*2))>>1)];
    if (codeUnit == 0)
      return str;
    ++i;
    // fromCharCode constructs a character from a UTF-16 code unit, so we can pass the UTF16 string right through.
    str += String.fromCharCode(codeUnit);
  }
}
Module['UTF16ToString'] = UTF16ToString;


function stringToUTF16(str, outPtr) {
  for(var i = 0; i < str.length; ++i) {
    // charCodeAt returns a UTF-16 encoded code unit, so it can be directly written to the HEAP.
    var codeUnit = str.charCodeAt(i); // possibly a lead surrogate
    HEAP16[(((outPtr)+(i*2))>>1)]=codeUnit;
  }
  // Null-terminate the pointer to the HEAP.
  HEAP16[(((outPtr)+(str.length*2))>>1)]=0;
}
Module['stringToUTF16'] = stringToUTF16;


function UTF32ToString(ptr) {
  var i = 0;

  var str = '';
  while (1) {
    var utf32 = HEAP32[(((ptr)+(i*4))>>2)];
    if (utf32 == 0)
      return str;
    ++i;
    // Gotcha: fromCharCode constructs a character from a UTF-16 encoded code (pair), not from a Unicode code point! So encode the code point to UTF-16 for constructing.
    if (utf32 >= 0x10000) {
      var ch = utf32 - 0x10000;
      str += String.fromCharCode(0xD800 | (ch >> 10), 0xDC00 | (ch & 0x3FF));
    } else {
      str += String.fromCharCode(utf32);
    }
  }
}
Module['UTF32ToString'] = UTF32ToString;


function stringToUTF32(str, outPtr) {
  var iChar = 0;
  for(var iCodeUnit = 0; iCodeUnit < str.length; ++iCodeUnit) {
    // Gotcha: charCodeAt returns a 16-bit word that is a UTF-16 encoded code unit, not a Unicode code point of the character! We must decode the string to UTF-32 to the heap.
    var codeUnit = str.charCodeAt(iCodeUnit); // possibly a lead surrogate
    if (codeUnit >= 0xD800 && codeUnit <= 0xDFFF) {
      var trailSurrogate = str.charCodeAt(++iCodeUnit);
      codeUnit = 0x10000 + ((codeUnit & 0x3FF) << 10) | (trailSurrogate & 0x3FF);
    }
    HEAP32[(((outPtr)+(iChar*4))>>2)]=codeUnit;
    ++iChar;
  }
  // Null-terminate the pointer to the HEAP.
  HEAP32[(((outPtr)+(iChar*4))>>2)]=0;
}
Module['stringToUTF32'] = stringToUTF32;

function demangle(func) {
  var hasLibcxxabi = !!Module['___cxa_demangle'];
  if (hasLibcxxabi) {
    try {
      var buf = _malloc(func.length);
      writeStringToMemory(func.substr(1), buf);
      var status = _malloc(4);
      var ret = Module['___cxa_demangle'](buf, 0, 0, status);
      if (getValue(status, 'i32') === 0 && ret) {
        return Pointer_stringify(ret);
      }
      // otherwise, libcxxabi failed, we can try ours which may return a partial result
    } catch(e) {
      // failure when using libcxxabi, we can try ours which may return a partial result
    } finally {
      if (buf) _free(buf);
      if (status) _free(status);
      if (ret) _free(ret);
    }
  }
  var i = 3;
  // params, etc.
  var basicTypes = {
    'v': 'void',
    'b': 'bool',
    'c': 'char',
    's': 'short',
    'i': 'int',
    'l': 'long',
    'f': 'float',
    'd': 'double',
    'w': 'wchar_t',
    'a': 'signed char',
    'h': 'unsigned char',
    't': 'unsigned short',
    'j': 'unsigned int',
    'm': 'unsigned long',
    'x': 'long long',
    'y': 'unsigned long long',
    'z': '...'
  };
  var subs = [];
  var first = true;
  function dump(x) {
    //return;
    if (x) Module.print(x);
    Module.print(func);
    var pre = '';
    for (var a = 0; a < i; a++) pre += ' ';
    Module.print (pre + '^');
  }
  function parseNested() {
    i++;
    if (func[i] === 'K') i++; // ignore const
    var parts = [];
    while (func[i] !== 'E') {
      if (func[i] === 'S') { // substitution
        i++;
        var next = func.indexOf('_', i);
        var num = func.substring(i, next) || 0;
        parts.push(subs[num] || '?');
        i = next+1;
        continue;
      }
      if (func[i] === 'C') { // constructor
        parts.push(parts[parts.length-1]);
        i += 2;
        continue;
      }
      var size = parseInt(func.substr(i));
      var pre = size.toString().length;
      if (!size || !pre) { i--; break; } // counter i++ below us
      var curr = func.substr(i + pre, size);
      parts.push(curr);
      subs.push(curr);
      i += pre + size;
    }
    i++; // skip E
    return parts;
  }
  function parse(rawList, limit, allowVoid) { // main parser
    limit = limit || Infinity;
    var ret = '', list = [];
    function flushList() {
      return '(' + list.join(', ') + ')';
    }
    var name;
    if (func[i] === 'N') {
      // namespaced N-E
      name = parseNested().join('::');
      limit--;
      if (limit === 0) return rawList ? [name] : name;
    } else {
      // not namespaced
      if (func[i] === 'K' || (first && func[i] === 'L')) i++; // ignore const and first 'L'
      var size = parseInt(func.substr(i));
      if (size) {
        var pre = size.toString().length;
        name = func.substr(i + pre, size);
        i += pre + size;
      }
    }
    first = false;
    if (func[i] === 'I') {
      i++;
      var iList = parse(true);
      var iRet = parse(true, 1, true);
      ret += iRet[0] + ' ' + name + '<' + iList.join(', ') + '>';
    } else {
      ret = name;
    }
    paramLoop: while (i < func.length && limit-- > 0) {
      //dump('paramLoop');
      var c = func[i++];
      if (c in basicTypes) {
        list.push(basicTypes[c]);
      } else {
        switch (c) {
          case 'P': list.push(parse(true, 1, true)[0] + '*'); break; // pointer
          case 'R': list.push(parse(true, 1, true)[0] + '&'); break; // reference
          case 'L': { // literal
            i++; // skip basic type
            var end = func.indexOf('E', i);
            var size = end - i;
            list.push(func.substr(i, size));
            i += size + 2; // size + 'EE'
            break;
          }
          case 'A': { // array
            var size = parseInt(func.substr(i));
            i += size.toString().length;
            if (func[i] !== '_') throw '?';
            i++; // skip _
            list.push(parse(true, 1, true)[0] + ' [' + size + ']');
            break;
          }
          case 'E': break paramLoop;
          default: ret += '?' + c; break paramLoop;
        }
      }
    }
    if (!allowVoid && list.length === 1 && list[0] === 'void') list = []; // avoid (void)
    if (rawList) {
      if (ret) {
        list.push(ret + '?');
      }
      return list;
    } else {
      return ret + flushList();
    }
  }
  var final = func;
  try {
    // Special-case the entry point, since its name differs from other name mangling.
    if (func == 'Object._main' || func == '_main') {
      return 'main()';
    }
    if (typeof func === 'number') func = Pointer_stringify(func);
    if (func[0] !== '_') return func;
    if (func[1] !== '_') return func; // C function
    if (func[2] !== 'Z') return func;
    switch (func[3]) {
      case 'n': return 'operator new()';
      case 'd': return 'operator delete()';
    }
    final = parse();
  } catch(e) {
    final += '?';
  }
  if (final.indexOf('?') >= 0 && !hasLibcxxabi) {
    Runtime.warnOnce('warning: a problem occurred in builtin C++ name demangling; build with  -s DEMANGLE_SUPPORT=1  to link in libcxxabi demangling');
  }
  return final;
}

function demangleAll(text) {
  return text.replace(/__Z[\w\d_]+/g, function(x) { var y = demangle(x); return x === y ? x : (x + ' [' + y + ']') });
}

function jsStackTrace() {
  var err = new Error();
  if (!err.stack) {
    // IE10+ special cases: It does have callstack info, but it is only populated if an Error object is thrown,
    // so try that as a special-case.
    try {
      throw new Error(0);
    } catch(e) {
      err = e;
    }
    if (!err.stack) {
      return '(no stack trace available)';
    }
  }
  return err.stack.toString();
}

function stackTrace() {
  return demangleAll(jsStackTrace());
}
Module['stackTrace'] = stackTrace;

// Memory management

var PAGE_SIZE = 4096;
function alignMemoryPage(x) {
  return (x+4095)&-4096;
}

var HEAP;
var HEAP8, HEAPU8, HEAP16, HEAPU16, HEAP32, HEAPU32, HEAPF32, HEAPF64;

var STATIC_BASE = 0, STATICTOP = 0, staticSealed = false; // static area
var STACK_BASE = 0, STACKTOP = 0, STACK_MAX = 0; // stack area
var DYNAMIC_BASE = 0, DYNAMICTOP = 0; // dynamic area handled by sbrk

function enlargeMemory() {
  abort('Cannot enlarge memory arrays. Either (1) compile with -s TOTAL_MEMORY=X with X higher than the current value ' + TOTAL_MEMORY + ', (2) compile with ALLOW_MEMORY_GROWTH which adjusts the size at runtime but prevents some optimizations, or (3) set Module.TOTAL_MEMORY before the program runs.');
}


var TOTAL_STACK = Module['TOTAL_STACK'] || 5242880;
var TOTAL_MEMORY = Module['TOTAL_MEMORY'] || 16777216;
var FAST_MEMORY = Module['FAST_MEMORY'] || 2097152;

var totalMemory = 64*1024;
while (totalMemory < TOTAL_MEMORY || totalMemory < 2*TOTAL_STACK) {
  if (totalMemory < 16*1024*1024) {
    totalMemory *= 2;
  } else {
    totalMemory += 16*1024*1024
  }
}
if (totalMemory !== TOTAL_MEMORY) {
  Module.printErr('increasing TOTAL_MEMORY to ' + totalMemory + ' to be more reasonable');
  TOTAL_MEMORY = totalMemory;
}

// Initialize the runtime's memory
// check for full engine support (use string 'subarray' to avoid closure compiler confusion)
assert(typeof Int32Array !== 'undefined' && typeof Float64Array !== 'undefined' && !!(new Int32Array(1)['subarray']) && !!(new Int32Array(1)['set']),
       'JS engine does not provide full typed array support');

var buffer = new ArrayBuffer(TOTAL_MEMORY);
HEAP8 = new Int8Array(buffer);
HEAP16 = new Int16Array(buffer);
HEAP32 = new Int32Array(buffer);
HEAPU8 = new Uint8Array(buffer);
HEAPU16 = new Uint16Array(buffer);
HEAPU32 = new Uint32Array(buffer);
HEAPF32 = new Float32Array(buffer);
HEAPF64 = new Float64Array(buffer);

// Endianness check (note: assumes compiler arch was little-endian)
HEAP32[0] = 255;
assert(HEAPU8[0] === 255 && HEAPU8[3] === 0, 'Typed arrays 2 must be run on a little-endian system');

Module['HEAP'] = HEAP;
Module['HEAP8'] = HEAP8;
Module['HEAP16'] = HEAP16;
Module['HEAP32'] = HEAP32;
Module['HEAPU8'] = HEAPU8;
Module['HEAPU16'] = HEAPU16;
Module['HEAPU32'] = HEAPU32;
Module['HEAPF32'] = HEAPF32;
Module['HEAPF64'] = HEAPF64;

function callRuntimeCallbacks(callbacks) {
  while(callbacks.length > 0) {
    var callback = callbacks.shift();
    if (typeof callback == 'function') {
      callback();
      continue;
    }
    var func = callback.func;
    if (typeof func === 'number') {
      if (callback.arg === undefined) {
        Runtime.dynCall('v', func);
      } else {
        Runtime.dynCall('vi', func, [callback.arg]);
      }
    } else {
      func(callback.arg === undefined ? null : callback.arg);
    }
  }
}

var __ATPRERUN__  = []; // functions called before the runtime is initialized
var __ATINIT__    = []; // functions called during startup
var __ATMAIN__    = []; // functions called when main() is to be run
var __ATEXIT__    = []; // functions called during shutdown
var __ATPOSTRUN__ = []; // functions called after the runtime has exited

var runtimeInitialized = false;
var runtimeExited = false;

function preRun() {
  // compatibility - merge in anything from Module['preRun'] at this time
  if (Module['preRun']) {
    if (typeof Module['preRun'] == 'function') Module['preRun'] = [Module['preRun']];
    while (Module['preRun'].length) {
      addOnPreRun(Module['preRun'].shift());
    }
  }
  callRuntimeCallbacks(__ATPRERUN__);
}

function ensureInitRuntime() {
  if (runtimeInitialized) return;
  runtimeInitialized = true;
  callRuntimeCallbacks(__ATINIT__);
}

function preMain() {
  callRuntimeCallbacks(__ATMAIN__);
}

function exitRuntime() {
  callRuntimeCallbacks(__ATEXIT__);
  runtimeExited = true;
}

function postRun() {
  // compatibility - merge in anything from Module['postRun'] at this time
  if (Module['postRun']) {
    if (typeof Module['postRun'] == 'function') Module['postRun'] = [Module['postRun']];
    while (Module['postRun'].length) {
      addOnPostRun(Module['postRun'].shift());
    }
  }
  callRuntimeCallbacks(__ATPOSTRUN__);
}

function addOnPreRun(cb) {
  __ATPRERUN__.unshift(cb);
}
Module['addOnPreRun'] = Module.addOnPreRun = addOnPreRun;

function addOnInit(cb) {
  __ATINIT__.unshift(cb);
}
Module['addOnInit'] = Module.addOnInit = addOnInit;

function addOnPreMain(cb) {
  __ATMAIN__.unshift(cb);
}
Module['addOnPreMain'] = Module.addOnPreMain = addOnPreMain;

function addOnExit(cb) {
  __ATEXIT__.unshift(cb);
}
Module['addOnExit'] = Module.addOnExit = addOnExit;

function addOnPostRun(cb) {
  __ATPOSTRUN__.unshift(cb);
}
Module['addOnPostRun'] = Module.addOnPostRun = addOnPostRun;

// Tools


function intArrayFromString(stringy, dontAddNull, length /* optional */) {
  var ret = (new Runtime.UTF8Processor()).processJSString(stringy);
  if (length) {
    ret.length = length;
  }
  if (!dontAddNull) {
    ret.push(0);
  }
  return ret;
}
Module['intArrayFromString'] = intArrayFromString;

function intArrayToString(array) {
  var ret = [];
  for (var i = 0; i < array.length; i++) {
    var chr = array[i];
    if (chr > 0xFF) {
      chr &= 0xFF;
    }
    ret.push(String.fromCharCode(chr));
  }
  return ret.join('');
}
Module['intArrayToString'] = intArrayToString;

function writeStringToMemory(string, buffer, dontAddNull) {
  var array = intArrayFromString(string, dontAddNull);
  var i = 0;
  while (i < array.length) {
    var chr = array[i];
    HEAP8[(((buffer)+(i))>>0)]=chr;
    i = i + 1;
  }
}
Module['writeStringToMemory'] = writeStringToMemory;

function writeArrayToMemory(array, buffer) {
  for (var i = 0; i < array.length; i++) {
    HEAP8[(((buffer)+(i))>>0)]=array[i];
  }
}
Module['writeArrayToMemory'] = writeArrayToMemory;

function writeAsciiToMemory(str, buffer, dontAddNull) {
  for (var i = 0; i < str.length; i++) {
    HEAP8[(((buffer)+(i))>>0)]=str.charCodeAt(i);
  }
  if (!dontAddNull) HEAP8[(((buffer)+(str.length))>>0)]=0;
}
Module['writeAsciiToMemory'] = writeAsciiToMemory;

function unSign(value, bits, ignore) {
  if (value >= 0) {
    return value;
  }
  return bits <= 32 ? 2*Math.abs(1 << (bits-1)) + value // Need some trickery, since if bits == 32, we are right at the limit of the bits JS uses in bitshifts
                    : Math.pow(2, bits)         + value;
}
function reSign(value, bits, ignore) {
  if (value <= 0) {
    return value;
  }
  var half = bits <= 32 ? Math.abs(1 << (bits-1)) // abs is needed if bits == 32
                        : Math.pow(2, bits-1);
  if (value >= half && (bits <= 32 || value > half)) { // for huge values, we can hit the precision limit and always get true here. so don't do that
                                                       // but, in general there is no perfect solution here. With 64-bit ints, we get rounding and errors
                                                       // TODO: In i64 mode 1, resign the two parts separately and safely
    value = -2*half + value; // Cannot bitshift half, as it may be at the limit of the bits JS uses in bitshifts
  }
  return value;
}

// check for imul support, and also for correctness ( https://bugs.webkit.org/show_bug.cgi?id=126345 )
if (!Math['imul'] || Math['imul'](0xffffffff, 5) !== -5) Math['imul'] = function imul(a, b) {
  var ah  = a >>> 16;
  var al = a & 0xffff;
  var bh  = b >>> 16;
  var bl = b & 0xffff;
  return (al*bl + ((ah*bl + al*bh) << 16))|0;
};
Math.imul = Math['imul'];


var Math_abs = Math.abs;
var Math_cos = Math.cos;
var Math_sin = Math.sin;
var Math_tan = Math.tan;
var Math_acos = Math.acos;
var Math_asin = Math.asin;
var Math_atan = Math.atan;
var Math_atan2 = Math.atan2;
var Math_exp = Math.exp;
var Math_log = Math.log;
var Math_sqrt = Math.sqrt;
var Math_ceil = Math.ceil;
var Math_floor = Math.floor;
var Math_pow = Math.pow;
var Math_imul = Math.imul;
var Math_fround = Math.fround;
var Math_min = Math.min;

// A counter of dependencies for calling run(). If we need to
// do asynchronous work before running, increment this and
// decrement it. Incrementing must happen in a place like
// PRE_RUN_ADDITIONS (used by emcc to add file preloading).
// Note that you can add dependencies in preRun, even though
// it happens right before run - run will be postponed until
// the dependencies are met.
var runDependencies = 0;
var runDependencyWatcher = null;
var dependenciesFulfilled = null; // overridden to take different actions when all run dependencies are fulfilled

function addRunDependency(id) {
  runDependencies++;
  if (Module['monitorRunDependencies']) {
    Module['monitorRunDependencies'](runDependencies);
  }
}
Module['addRunDependency'] = addRunDependency;
function removeRunDependency(id) {
  runDependencies--;
  if (Module['monitorRunDependencies']) {
    Module['monitorRunDependencies'](runDependencies);
  }
  if (runDependencies == 0) {
    if (runDependencyWatcher !== null) {
      clearInterval(runDependencyWatcher);
      runDependencyWatcher = null;
    }
    if (dependenciesFulfilled) {
      var callback = dependenciesFulfilled;
      dependenciesFulfilled = null;
      callback(); // can add another dependenciesFulfilled
    }
  }
}
Module['removeRunDependency'] = removeRunDependency;

Module["preloadedImages"] = {}; // maps url to image data
Module["preloadedAudios"] = {}; // maps url to audio data


var memoryInitializer = null;

// === Body ===





STATIC_BASE = 8;

STATICTOP = STATIC_BASE + Runtime.alignMemory(4499);
  /* global initializers */ __ATINIT__.push();
  

var memoryInitializer = "a.out.js.mem";




var tempDoublePtr = Runtime.alignMemory(allocate(12, "i8", ALLOC_STATIC), 8);

assert(tempDoublePtr % 8 == 0);

function copyTempFloat(ptr) { // functions, because inlining this code increases code size too much

  HEAP8[tempDoublePtr] = HEAP8[ptr];

  HEAP8[tempDoublePtr+1] = HEAP8[ptr+1];

  HEAP8[tempDoublePtr+2] = HEAP8[ptr+2];

  HEAP8[tempDoublePtr+3] = HEAP8[ptr+3];

}

function copyTempDouble(ptr) {

  HEAP8[tempDoublePtr] = HEAP8[ptr];

  HEAP8[tempDoublePtr+1] = HEAP8[ptr+1];

  HEAP8[tempDoublePtr+2] = HEAP8[ptr+2];

  HEAP8[tempDoublePtr+3] = HEAP8[ptr+3];

  HEAP8[tempDoublePtr+4] = HEAP8[ptr+4];

  HEAP8[tempDoublePtr+5] = HEAP8[ptr+5];

  HEAP8[tempDoublePtr+6] = HEAP8[ptr+6];

  HEAP8[tempDoublePtr+7] = HEAP8[ptr+7];

}


  function _sbrk(bytes) {
      // Implement a Linux-like 'memory area' for our 'process'.
      // Changes the size of the memory area by |bytes|; returns the
      // address of the previous top ('break') of the memory area
      // We control the "dynamic" memory - DYNAMIC_BASE to DYNAMICTOP
      var self = _sbrk;
      if (!self.called) {
        DYNAMICTOP = alignMemoryPage(DYNAMICTOP); // make sure we start out aligned
        self.called = true;
        assert(Runtime.dynamicAlloc);
        self.alloc = Runtime.dynamicAlloc;
        Runtime.dynamicAlloc = function() { abort('cannot dynamically allocate, sbrk now has control') };
      }
      var ret = DYNAMICTOP;
      if (bytes != 0) self.alloc(bytes);
      return ret;  // Previous break location.
    }

  
  
  var ___errno_state=0;function ___setErrNo(value) {
      // For convenient setting and returning of errno.
      HEAP32[((___errno_state)>>2)]=value;
      return value;
    }
  
  var ERRNO_CODES={EPERM:1,ENOENT:2,ESRCH:3,EINTR:4,EIO:5,ENXIO:6,E2BIG:7,ENOEXEC:8,EBADF:9,ECHILD:10,EAGAIN:11,EWOULDBLOCK:11,ENOMEM:12,EACCES:13,EFAULT:14,ENOTBLK:15,EBUSY:16,EEXIST:17,EXDEV:18,ENODEV:19,ENOTDIR:20,EISDIR:21,EINVAL:22,ENFILE:23,EMFILE:24,ENOTTY:25,ETXTBSY:26,EFBIG:27,ENOSPC:28,ESPIPE:29,EROFS:30,EMLINK:31,EPIPE:32,EDOM:33,ERANGE:34,ENOMSG:42,EIDRM:43,ECHRNG:44,EL2NSYNC:45,EL3HLT:46,EL3RST:47,ELNRNG:48,EUNATCH:49,ENOCSI:50,EL2HLT:51,EDEADLK:35,ENOLCK:37,EBADE:52,EBADR:53,EXFULL:54,ENOANO:55,EBADRQC:56,EBADSLT:57,EDEADLOCK:35,EBFONT:59,ENOSTR:60,ENODATA:61,ETIME:62,ENOSR:63,ENONET:64,ENOPKG:65,EREMOTE:66,ENOLINK:67,EADV:68,ESRMNT:69,ECOMM:70,EPROTO:71,EMULTIHOP:72,EDOTDOT:73,EBADMSG:74,ENOTUNIQ:76,EBADFD:77,EREMCHG:78,ELIBACC:79,ELIBBAD:80,ELIBSCN:81,ELIBMAX:82,ELIBEXEC:83,ENOSYS:38,ENOTEMPTY:39,ENAMETOOLONG:36,ELOOP:40,EOPNOTSUPP:95,EPFNOSUPPORT:96,ECONNRESET:104,ENOBUFS:105,EAFNOSUPPORT:97,EPROTOTYPE:91,ENOTSOCK:88,ENOPROTOOPT:92,ESHUTDOWN:108,ECONNREFUSED:111,EADDRINUSE:98,ECONNABORTED:103,ENETUNREACH:101,ENETDOWN:100,ETIMEDOUT:110,EHOSTDOWN:112,EHOSTUNREACH:113,EINPROGRESS:115,EALREADY:114,EDESTADDRREQ:89,EMSGSIZE:90,EPROTONOSUPPORT:93,ESOCKTNOSUPPORT:94,EADDRNOTAVAIL:99,ENETRESET:102,EISCONN:106,ENOTCONN:107,ETOOMANYREFS:109,EUSERS:87,EDQUOT:122,ESTALE:116,ENOTSUP:95,ENOMEDIUM:123,EILSEQ:84,EOVERFLOW:75,ECANCELED:125,ENOTRECOVERABLE:131,EOWNERDEAD:130,ESTRPIPE:86};function _sysconf(name) {
      // long sysconf(int name);
      // http://pubs.opengroup.org/onlinepubs/009695399/functions/sysconf.html
      switch(name) {
        case 30: return PAGE_SIZE;
        case 132:
        case 133:
        case 12:
        case 137:
        case 138:
        case 15:
        case 235:
        case 16:
        case 17:
        case 18:
        case 19:
        case 20:
        case 149:
        case 13:
        case 10:
        case 236:
        case 153:
        case 9:
        case 21:
        case 22:
        case 159:
        case 154:
        case 14:
        case 77:
        case 78:
        case 139:
        case 80:
        case 81:
        case 79:
        case 82:
        case 68:
        case 67:
        case 164:
        case 11:
        case 29:
        case 47:
        case 48:
        case 95:
        case 52:
        case 51:
        case 46:
          return 200809;
        case 27:
        case 246:
        case 127:
        case 128:
        case 23:
        case 24:
        case 160:
        case 161:
        case 181:
        case 182:
        case 242:
        case 183:
        case 184:
        case 243:
        case 244:
        case 245:
        case 165:
        case 178:
        case 179:
        case 49:
        case 50:
        case 168:
        case 169:
        case 175:
        case 170:
        case 171:
        case 172:
        case 97:
        case 76:
        case 32:
        case 173:
        case 35:
          return -1;
        case 176:
        case 177:
        case 7:
        case 155:
        case 8:
        case 157:
        case 125:
        case 126:
        case 92:
        case 93:
        case 129:
        case 130:
        case 131:
        case 94:
        case 91:
          return 1;
        case 74:
        case 60:
        case 69:
        case 70:
        case 4:
          return 1024;
        case 31:
        case 42:
        case 72:
          return 32;
        case 87:
        case 26:
        case 33:
          return 2147483647;
        case 34:
        case 1:
          return 47839;
        case 38:
        case 36:
          return 99;
        case 43:
        case 37:
          return 2048;
        case 0: return 2097152;
        case 3: return 65536;
        case 28: return 32768;
        case 44: return 32767;
        case 75: return 16384;
        case 39: return 1000;
        case 89: return 700;
        case 71: return 256;
        case 40: return 255;
        case 2: return 100;
        case 180: return 64;
        case 25: return 20;
        case 5: return 16;
        case 6: return 6;
        case 73: return 4;
        case 84: {
          if (typeof navigator === 'object') return navigator['hardwareConcurrency'] || 1;
          return 1;
        }
      }
      ___setErrNo(ERRNO_CODES.EINVAL);
      return -1;
    }

   
  Module["_memset"] = _memset;

  function ___errno_location() {
      return ___errno_state;
    }

  function _abort() {
      Module['abort']();
    }

  
  
  
  var ERRNO_MESSAGES={0:"Success",1:"Not super-user",2:"No such file or directory",3:"No such process",4:"Interrupted system call",5:"I/O error",6:"No such device or address",7:"Arg list too long",8:"Exec format error",9:"Bad file number",10:"No children",11:"No more processes",12:"Not enough core",13:"Permission denied",14:"Bad address",15:"Block device required",16:"Mount device busy",17:"File exists",18:"Cross-device link",19:"No such device",20:"Not a directory",21:"Is a directory",22:"Invalid argument",23:"Too many open files in system",24:"Too many open files",25:"Not a typewriter",26:"Text file busy",27:"File too large",28:"No space left on device",29:"Illegal seek",30:"Read only file system",31:"Too many links",32:"Broken pipe",33:"Math arg out of domain of func",34:"Math result not representable",35:"File locking deadlock error",36:"File or path name too long",37:"No record locks available",38:"Function not implemented",39:"Directory not empty",40:"Too many symbolic links",42:"No message of desired type",43:"Identifier removed",44:"Channel number out of range",45:"Level 2 not synchronized",46:"Level 3 halted",47:"Level 3 reset",48:"Link number out of range",49:"Protocol driver not attached",50:"No CSI structure available",51:"Level 2 halted",52:"Invalid exchange",53:"Invalid request descriptor",54:"Exchange full",55:"No anode",56:"Invalid request code",57:"Invalid slot",59:"Bad font file fmt",60:"Device not a stream",61:"No data (for no delay io)",62:"Timer expired",63:"Out of streams resources",64:"Machine is not on the network",65:"Package not installed",66:"The object is remote",67:"The link has been severed",68:"Advertise error",69:"Srmount error",70:"Communication error on send",71:"Protocol error",72:"Multihop attempted",73:"Cross mount point (not really error)",74:"Trying to read unreadable message",75:"Value too large for defined data type",76:"Given log. name not unique",77:"f.d. invalid for this operation",78:"Remote address changed",79:"Can   access a needed shared lib",80:"Accessing a corrupted shared lib",81:".lib section in a.out corrupted",82:"Attempting to link in too many libs",83:"Attempting to exec a shared library",84:"Illegal byte sequence",86:"Streams pipe error",87:"Too many users",88:"Socket operation on non-socket",89:"Destination address required",90:"Message too long",91:"Protocol wrong type for socket",92:"Protocol not available",93:"Unknown protocol",94:"Socket type not supported",95:"Not supported",96:"Protocol family not supported",97:"Address family not supported by protocol family",98:"Address already in use",99:"Address not available",100:"Network interface is not configured",101:"Network is unreachable",102:"Connection reset by network",103:"Connection aborted",104:"Connection reset by peer",105:"No buffer space available",106:"Socket is already connected",107:"Socket is not connected",108:"Can't send after socket shutdown",109:"Too many references",110:"Connection timed out",111:"Connection refused",112:"Host is down",113:"Host is unreachable",114:"Socket already connected",115:"Connection already in progress",116:"Stale file handle",122:"Quota exceeded",123:"No medium (in tape drive)",125:"Operation canceled",130:"Previous owner died",131:"State not recoverable"};
  
  var TTY={ttys:[],init:function () {
        // https://github.com/kripken/emscripten/pull/1555
        // if (ENVIRONMENT_IS_NODE) {
        //   // currently, FS.init does not distinguish if process.stdin is a file or TTY
        //   // device, it always assumes it's a TTY device. because of this, we're forcing
        //   // process.stdin to UTF8 encoding to at least make stdin reading compatible
        //   // with text files until FS.init can be refactored.
        //   process['stdin']['setEncoding']('utf8');
        // }
      },shutdown:function () {
        // https://github.com/kripken/emscripten/pull/1555
        // if (ENVIRONMENT_IS_NODE) {
        //   // inolen: any idea as to why node -e 'process.stdin.read()' wouldn't exit immediately (with process.stdin being a tty)?
        //   // isaacs: because now it's reading from the stream, you've expressed interest in it, so that read() kicks off a _read() which creates a ReadReq operation
        //   // inolen: I thought read() in that case was a synchronous operation that just grabbed some amount of buffered data if it exists?
        //   // isaacs: it is. but it also triggers a _read() call, which calls readStart() on the handle
        //   // isaacs: do process.stdin.pause() and i'd think it'd probably close the pending call
        //   process['stdin']['pause']();
        // }
      },register:function (dev, ops) {
        TTY.ttys[dev] = { input: [], output: [], ops: ops };
        FS.registerDevice(dev, TTY.stream_ops);
      },stream_ops:{open:function (stream) {
          var tty = TTY.ttys[stream.node.rdev];
          if (!tty) {
            throw new FS.ErrnoError(ERRNO_CODES.ENODEV);
          }
          stream.tty = tty;
          stream.seekable = false;
        },close:function (stream) {
          // flush any pending line data
          if (stream.tty.output.length) {
            stream.tty.ops.put_char(stream.tty, 10);
          }
        },read:function (stream, buffer, offset, length, pos /* ignored */) {
          if (!stream.tty || !stream.tty.ops.get_char) {
            throw new FS.ErrnoError(ERRNO_CODES.ENXIO);
          }
          var bytesRead = 0;
          for (var i = 0; i < length; i++) {
            var result;
            try {
              result = stream.tty.ops.get_char(stream.tty);
            } catch (e) {
              throw new FS.ErrnoError(ERRNO_CODES.EIO);
            }
            if (result === undefined && bytesRead === 0) {
              throw new FS.ErrnoError(ERRNO_CODES.EAGAIN);
            }
            if (result === null || result === undefined) break;
            bytesRead++;
            buffer[offset+i] = result;
          }
          if (bytesRead) {
            stream.node.timestamp = Date.now();
          }
          return bytesRead;
        },write:function (stream, buffer, offset, length, pos) {
          if (!stream.tty || !stream.tty.ops.put_char) {
            throw new FS.ErrnoError(ERRNO_CODES.ENXIO);
          }
          for (var i = 0; i < length; i++) {
            try {
              stream.tty.ops.put_char(stream.tty, buffer[offset+i]);
            } catch (e) {
              throw new FS.ErrnoError(ERRNO_CODES.EIO);
            }
          }
          if (length) {
            stream.node.timestamp = Date.now();
          }
          return i;
        }},default_tty_ops:{get_char:function (tty) {
          if (!tty.input.length) {
            var result = null;
            if (ENVIRONMENT_IS_NODE) {
              result = process['stdin']['read']();
              if (!result) {
                if (process['stdin']['_readableState'] && process['stdin']['_readableState']['ended']) {
                  return null;  // EOF
                }
                return undefined;  // no data available
              }
            } else if (typeof window != 'undefined' &&
              typeof window.prompt == 'function') {
              // Browser.
              result = window.prompt('Input: ');  // returns null on cancel
              if (result !== null) {
                result += '\n';
              }
            } else if (typeof readline == 'function') {
              // Command line.
              result = readline();
              if (result !== null) {
                result += '\n';
              }
            }
            if (!result) {
              return null;
            }
            tty.input = intArrayFromString(result, true);
          }
          return tty.input.shift();
        },put_char:function (tty, val) {
          if (val === null || val === 10) {
            Module['print'](tty.output.join(''));
            tty.output = [];
          } else {
            tty.output.push(TTY.utf8.processCChar(val));
          }
        }},default_tty1_ops:{put_char:function (tty, val) {
          if (val === null || val === 10) {
            Module['printErr'](tty.output.join(''));
            tty.output = [];
          } else {
            tty.output.push(TTY.utf8.processCChar(val));
          }
        }}};
  
  var MEMFS={ops_table:null,mount:function (mount) {
        return MEMFS.createNode(null, '/', 16384 | 511 /* 0777 */, 0);
      },createNode:function (parent, name, mode, dev) {
        if (FS.isBlkdev(mode) || FS.isFIFO(mode)) {
          // no supported
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        if (!MEMFS.ops_table) {
          MEMFS.ops_table = {
            dir: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr,
                lookup: MEMFS.node_ops.lookup,
                mknod: MEMFS.node_ops.mknod,
                rename: MEMFS.node_ops.rename,
                unlink: MEMFS.node_ops.unlink,
                rmdir: MEMFS.node_ops.rmdir,
                readdir: MEMFS.node_ops.readdir,
                symlink: MEMFS.node_ops.symlink
              },
              stream: {
                llseek: MEMFS.stream_ops.llseek
              }
            },
            file: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr
              },
              stream: {
                llseek: MEMFS.stream_ops.llseek,
                read: MEMFS.stream_ops.read,
                write: MEMFS.stream_ops.write,
                allocate: MEMFS.stream_ops.allocate,
                mmap: MEMFS.stream_ops.mmap
              }
            },
            link: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr,
                readlink: MEMFS.node_ops.readlink
              },
              stream: {}
            },
            chrdev: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr
              },
              stream: FS.chrdev_stream_ops
            },
          };
        }
        var node = FS.createNode(parent, name, mode, dev);
        if (FS.isDir(node.mode)) {
          node.node_ops = MEMFS.ops_table.dir.node;
          node.stream_ops = MEMFS.ops_table.dir.stream;
          node.contents = {};
        } else if (FS.isFile(node.mode)) {
          node.node_ops = MEMFS.ops_table.file.node;
          node.stream_ops = MEMFS.ops_table.file.stream;
          node.usedBytes = 0; // The actual number of bytes used in the typed array, as opposed to contents.buffer.byteLength which gives the whole capacity.
          // When the byte data of the file is populated, this will point to either a typed array, or a normal JS array. Typed arrays are preferred
          // for performance, and used by default. However, typed arrays are not resizable like normal JS arrays are, so there is a small disk size
          // penalty involved for appending file writes that continuously grow a file similar to std::vector capacity vs used -scheme.
          node.contents = null; 
        } else if (FS.isLink(node.mode)) {
          node.node_ops = MEMFS.ops_table.link.node;
          node.stream_ops = MEMFS.ops_table.link.stream;
        } else if (FS.isChrdev(node.mode)) {
          node.node_ops = MEMFS.ops_table.chrdev.node;
          node.stream_ops = MEMFS.ops_table.chrdev.stream;
        }
        node.timestamp = Date.now();
        // add the new node to the parent
        if (parent) {
          parent.contents[name] = node;
        }
        return node;
      },getFileDataAsRegularArray:function (node) {
        if (node.contents && node.contents.subarray) {
          var arr = [];
          for (var i = 0; i < node.usedBytes; ++i) arr.push(node.contents[i]);
          return arr; // Returns a copy of the original data.
        }
        return node.contents; // No-op, the file contents are already in a JS array. Return as-is.
      },getFileDataAsTypedArray:function (node) {
        if (node.contents && node.contents.subarray) return node.contents.subarray(0, node.usedBytes); // Make sure to not return excess unused bytes.
        return new Uint8Array(node.contents);
      },expandFileStorage:function (node, newCapacity) {
  
        // If we are asked to expand the size of a file that already exists, revert to using a standard JS array to store the file
        // instead of a typed array. This makes resizing the array more flexible because we can just .push() elements at the back to
        // increase the size.
        if (node.contents && node.contents.subarray && newCapacity > node.contents.length) {
          node.contents = MEMFS.getFileDataAsRegularArray(node);
          node.usedBytes = node.contents.length; // We might be writing to a lazy-loaded file which had overridden this property, so force-reset it.
        }
  
        if (!node.contents || node.contents.subarray) { // Keep using a typed array if creating a new storage, or if old one was a typed array as well.
          var prevCapacity = node.contents ? node.contents.buffer.byteLength : 0;
          if (prevCapacity >= newCapacity) return; // No need to expand, the storage was already large enough.
          // Don't expand strictly to the given requested limit if it's only a very small increase, but instead geometrically grow capacity.
          // For small filesizes (<1MB), perform size*2 geometric increase, but for large sizes, do a much more conservative size*1.125 increase to
          // avoid overshooting the allocation cap by a very large margin.
          var CAPACITY_DOUBLING_MAX = 1024 * 1024;
          newCapacity = Math.max(newCapacity, (prevCapacity * (prevCapacity < CAPACITY_DOUBLING_MAX ? 2.0 : 1.125)) | 0);
          if (prevCapacity != 0) newCapacity = Math.max(newCapacity, 256); // At minimum allocate 256b for each file when expanding.
          var oldContents = node.contents;
          node.contents = new Uint8Array(newCapacity); // Allocate new storage.
          if (node.usedBytes > 0) node.contents.set(oldContents.subarray(0, node.usedBytes), 0); // Copy old data over to the new storage.
          return;
        }
        // Not using a typed array to back the file storage. Use a standard JS array instead.
        if (!node.contents && newCapacity > 0) node.contents = [];
        while (node.contents.length < newCapacity) node.contents.push(0);
      },resizeFileStorage:function (node, newSize) {
        if (node.usedBytes == newSize) return;
        if (newSize == 0) {
          node.contents = null; // Fully decommit when requesting a resize to zero.
          node.usedBytes = 0;
          return;
        }
  
        if (!node.contents || node.contents.subarray) { // Resize a typed array if that is being used as the backing store.
          var oldContents = node.contents;
          node.contents = new Uint8Array(new ArrayBuffer(newSize)); // Allocate new storage.
          if (oldContents) {
            node.contents.set(oldContents.subarray(0, Math.min(newSize, node.usedBytes))); // Copy old data over to the new storage.
          }
          node.usedBytes = newSize;
          return;
        }
        // Backing with a JS array.
        if (!node.contents) node.contents = [];
        if (node.contents.length > newSize) node.contents.length = newSize;
        else while (node.contents.length < newSize) node.contents.push(0);
        node.usedBytes = newSize;
      },node_ops:{getattr:function (node) {
          var attr = {};
          // device numbers reuse inode numbers.
          attr.dev = FS.isChrdev(node.mode) ? node.id : 1;
          attr.ino = node.id;
          attr.mode = node.mode;
          attr.nlink = 1;
          attr.uid = 0;
          attr.gid = 0;
          attr.rdev = node.rdev;
          if (FS.isDir(node.mode)) {
            attr.size = 4096;
          } else if (FS.isFile(node.mode)) {
            attr.size = node.usedBytes;
          } else if (FS.isLink(node.mode)) {
            attr.size = node.link.length;
          } else {
            attr.size = 0;
          }
          attr.atime = new Date(node.timestamp);
          attr.mtime = new Date(node.timestamp);
          attr.ctime = new Date(node.timestamp);
          // NOTE: In our implementation, st_blocks = Math.ceil(st_size/st_blksize),
          //       but this is not required by the standard.
          attr.blksize = 4096;
          attr.blocks = Math.ceil(attr.size / attr.blksize);
          return attr;
        },setattr:function (node, attr) {
          if (attr.mode !== undefined) {
            node.mode = attr.mode;
          }
          if (attr.timestamp !== undefined) {
            node.timestamp = attr.timestamp;
          }
          if (attr.size !== undefined) {
            MEMFS.resizeFileStorage(node, attr.size);
          }
        },lookup:function (parent, name) {
          throw FS.genericErrors[ERRNO_CODES.ENOENT];
        },mknod:function (parent, name, mode, dev) {
          return MEMFS.createNode(parent, name, mode, dev);
        },rename:function (old_node, new_dir, new_name) {
          // if we're overwriting a directory at new_name, make sure it's empty.
          if (FS.isDir(old_node.mode)) {
            var new_node;
            try {
              new_node = FS.lookupNode(new_dir, new_name);
            } catch (e) {
            }
            if (new_node) {
              for (var i in new_node.contents) {
                throw new FS.ErrnoError(ERRNO_CODES.ENOTEMPTY);
              }
            }
          }
          // do the internal rewiring
          delete old_node.parent.contents[old_node.name];
          old_node.name = new_name;
          new_dir.contents[new_name] = old_node;
          old_node.parent = new_dir;
        },unlink:function (parent, name) {
          delete parent.contents[name];
        },rmdir:function (parent, name) {
          var node = FS.lookupNode(parent, name);
          for (var i in node.contents) {
            throw new FS.ErrnoError(ERRNO_CODES.ENOTEMPTY);
          }
          delete parent.contents[name];
        },readdir:function (node) {
          var entries = ['.', '..']
          for (var key in node.contents) {
            if (!node.contents.hasOwnProperty(key)) {
              continue;
            }
            entries.push(key);
          }
          return entries;
        },symlink:function (parent, newname, oldpath) {
          var node = MEMFS.createNode(parent, newname, 511 /* 0777 */ | 40960, 0);
          node.link = oldpath;
          return node;
        },readlink:function (node) {
          if (!FS.isLink(node.mode)) {
            throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
          }
          return node.link;
        }},stream_ops:{read:function (stream, buffer, offset, length, position) {
          var contents = stream.node.contents;
          if (position >= stream.node.usedBytes) return 0;
          var size = Math.min(stream.node.usedBytes - position, length);
          assert(size >= 0);
          if (size > 8 && contents.subarray) { // non-trivial, and typed array
            buffer.set(contents.subarray(position, position + size), offset);
          } else
          {
            for (var i = 0; i < size; i++) buffer[offset + i] = contents[position + i];
          }
          return size;
        },write:function (stream, buffer, offset, length, position, canOwn) {
          if (!length) return 0;
          var node = stream.node;
          node.timestamp = Date.now();
  
          if (buffer.subarray && (!node.contents || node.contents.subarray)) { // This write is from a typed array to a typed array?
            if (canOwn) { // Can we just reuse the buffer we are given?
              node.contents = buffer.subarray(offset, offset + length);
              node.usedBytes = length;
              return length;
            } else if (node.usedBytes === 0 && position === 0) { // If this is a simple first write to an empty file, do a fast set since we don't need to care about old data.
              node.contents = new Uint8Array(buffer.subarray(offset, offset + length));
              node.usedBytes = length;
              return length;
            } else if (position + length <= node.usedBytes) { // Writing to an already allocated and used subrange of the file?
              node.contents.set(buffer.subarray(offset, offset + length), position);
              return length;
            }
          }
          // Appending to an existing file and we need to reallocate, or source data did not come as a typed array.
          MEMFS.expandFileStorage(node, position+length);
          if (node.contents.subarray && buffer.subarray) node.contents.set(buffer.subarray(offset, offset + length), position); // Use typed array write if available.
          else
            for (var i = 0; i < length; i++) {
             node.contents[position + i] = buffer[offset + i]; // Or fall back to manual write if not.
            }
          node.usedBytes = Math.max(node.usedBytes, position+length);
          return length;
        },llseek:function (stream, offset, whence) {
          var position = offset;
          if (whence === 1) {  // SEEK_CUR.
            position += stream.position;
          } else if (whence === 2) {  // SEEK_END.
            if (FS.isFile(stream.node.mode)) {
              position += stream.node.usedBytes;
            }
          }
          if (position < 0) {
            throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
          }
          stream.ungotten = [];
          stream.position = position;
          return position;
        },allocate:function (stream, offset, length) {
          MEMFS.expandFileStorage(stream.node, offset + length);
          stream.node.usedBytes = Math.max(stream.node.usedBytes, offset + length);
        },mmap:function (stream, buffer, offset, length, position, prot, flags) {
          if (!FS.isFile(stream.node.mode)) {
            throw new FS.ErrnoError(ERRNO_CODES.ENODEV);
          }
          var ptr;
          var allocated;
          var contents = stream.node.contents;
          // Only make a new copy when MAP_PRIVATE is specified.
          if ( !(flags & 2) &&
                (contents.buffer === buffer || contents.buffer === buffer.buffer) ) {
            // We can't emulate MAP_SHARED when the file is not backed by the buffer
            // we're mapping to (e.g. the HEAP buffer).
            allocated = false;
            ptr = contents.byteOffset;
          } else {
            // Try to avoid unnecessary slices.
            if (position > 0 || position + length < stream.node.usedBytes) {
              if (contents.subarray) {
                contents = contents.subarray(position, position + length);
              } else {
                contents = Array.prototype.slice.call(contents, position, position + length);
              }
            }
            allocated = true;
            ptr = _malloc(length);
            if (!ptr) {
              throw new FS.ErrnoError(ERRNO_CODES.ENOMEM);
            }
            buffer.set(contents, ptr);
          }
          return { ptr: ptr, allocated: allocated };
        }}};
  
  var IDBFS={dbs:{},indexedDB:function () {
        if (typeof indexedDB !== 'undefined') return indexedDB;
        var ret = null;
        if (typeof window === 'object') ret = window.indexedDB || window.mozIndexedDB || window.webkitIndexedDB || window.msIndexedDB;
        assert(ret, 'IDBFS used, but indexedDB not supported');
        return ret;
      },DB_VERSION:21,DB_STORE_NAME:"FILE_DATA",mount:function (mount) {
        // reuse all of the core MEMFS functionality
        return MEMFS.mount.apply(null, arguments);
      },syncfs:function (mount, populate, callback) {
        IDBFS.getLocalSet(mount, function(err, local) {
          if (err) return callback(err);
  
          IDBFS.getRemoteSet(mount, function(err, remote) {
            if (err) return callback(err);
  
            var src = populate ? remote : local;
            var dst = populate ? local : remote;
  
            IDBFS.reconcile(src, dst, callback);
          });
        });
      },getDB:function (name, callback) {
        // check the cache first
        var db = IDBFS.dbs[name];
        if (db) {
          return callback(null, db);
        }
  
        var req;
        try {
          req = IDBFS.indexedDB().open(name, IDBFS.DB_VERSION);
        } catch (e) {
          return callback(e);
        }
        req.onupgradeneeded = function(e) {
          var db = e.target.result;
          var transaction = e.target.transaction;
  
          var fileStore;
  
          if (db.objectStoreNames.contains(IDBFS.DB_STORE_NAME)) {
            fileStore = transaction.objectStore(IDBFS.DB_STORE_NAME);
          } else {
            fileStore = db.createObjectStore(IDBFS.DB_STORE_NAME);
          }
  
          fileStore.createIndex('timestamp', 'timestamp', { unique: false });
        };
        req.onsuccess = function() {
          db = req.result;
  
          // add to the cache
          IDBFS.dbs[name] = db;
          callback(null, db);
        };
        req.onerror = function() {
          callback(this.error);
        };
      },getLocalSet:function (mount, callback) {
        var entries = {};
  
        function isRealDir(p) {
          return p !== '.' && p !== '..';
        };
        function toAbsolute(root) {
          return function(p) {
            return PATH.join2(root, p);
          }
        };
  
        var check = FS.readdir(mount.mountpoint).filter(isRealDir).map(toAbsolute(mount.mountpoint));
  
        while (check.length) {
          var path = check.pop();
          var stat;
  
          try {
            stat = FS.stat(path);
          } catch (e) {
            return callback(e);
          }
  
          if (FS.isDir(stat.mode)) {
            check.push.apply(check, FS.readdir(path).filter(isRealDir).map(toAbsolute(path)));
          }
  
          entries[path] = { timestamp: stat.mtime };
        }
  
        return callback(null, { type: 'local', entries: entries });
      },getRemoteSet:function (mount, callback) {
        var entries = {};
  
        IDBFS.getDB(mount.mountpoint, function(err, db) {
          if (err) return callback(err);
  
          var transaction = db.transaction([IDBFS.DB_STORE_NAME], 'readonly');
          transaction.onerror = function() { callback(this.error); };
  
          var store = transaction.objectStore(IDBFS.DB_STORE_NAME);
          var index = store.index('timestamp');
  
          index.openKeyCursor().onsuccess = function(event) {
            var cursor = event.target.result;
  
            if (!cursor) {
              return callback(null, { type: 'remote', db: db, entries: entries });
            }
  
            entries[cursor.primaryKey] = { timestamp: cursor.key };
  
            cursor.continue();
          };
        });
      },loadLocalEntry:function (path, callback) {
        var stat, node;
  
        try {
          var lookup = FS.lookupPath(path);
          node = lookup.node;
          stat = FS.stat(path);
        } catch (e) {
          return callback(e);
        }
  
        if (FS.isDir(stat.mode)) {
          return callback(null, { timestamp: stat.mtime, mode: stat.mode });
        } else if (FS.isFile(stat.mode)) {
          // Performance consideration: storing a normal JavaScript array to a IndexedDB is much slower than storing a typed array.
          // Therefore always convert the file contents to a typed array first before writing the data to IndexedDB.
          node.contents = MEMFS.getFileDataAsTypedArray(node);
          return callback(null, { timestamp: stat.mtime, mode: stat.mode, contents: node.contents });
        } else {
          return callback(new Error('node type not supported'));
        }
      },storeLocalEntry:function (path, entry, callback) {
        try {
          if (FS.isDir(entry.mode)) {
            FS.mkdir(path, entry.mode);
          } else if (FS.isFile(entry.mode)) {
            FS.writeFile(path, entry.contents, { encoding: 'binary', canOwn: true });
          } else {
            return callback(new Error('node type not supported'));
          }
  
          FS.chmod(path, entry.mode);
          FS.utime(path, entry.timestamp, entry.timestamp);
        } catch (e) {
          return callback(e);
        }
  
        callback(null);
      },removeLocalEntry:function (path, callback) {
        try {
          var lookup = FS.lookupPath(path);
          var stat = FS.stat(path);
  
          if (FS.isDir(stat.mode)) {
            FS.rmdir(path);
          } else if (FS.isFile(stat.mode)) {
            FS.unlink(path);
          }
        } catch (e) {
          return callback(e);
        }
  
        callback(null);
      },loadRemoteEntry:function (store, path, callback) {
        var req = store.get(path);
        req.onsuccess = function(event) { callback(null, event.target.result); };
        req.onerror = function() { callback(this.error); };
      },storeRemoteEntry:function (store, path, entry, callback) {
        var req = store.put(entry, path);
        req.onsuccess = function() { callback(null); };
        req.onerror = function() { callback(this.error); };
      },removeRemoteEntry:function (store, path, callback) {
        var req = store.delete(path);
        req.onsuccess = function() { callback(null); };
        req.onerror = function() { callback(this.error); };
      },reconcile:function (src, dst, callback) {
        var total = 0;
  
        var create = [];
        Object.keys(src.entries).forEach(function (key) {
          var e = src.entries[key];
          var e2 = dst.entries[key];
          if (!e2 || e.timestamp > e2.timestamp) {
            create.push(key);
            total++;
          }
        });
  
        var remove = [];
        Object.keys(dst.entries).forEach(function (key) {
          var e = dst.entries[key];
          var e2 = src.entries[key];
          if (!e2) {
            remove.push(key);
            total++;
          }
        });
  
        if (!total) {
          return callback(null);
        }
  
        var errored = false;
        var completed = 0;
        var db = src.type === 'remote' ? src.db : dst.db;
        var transaction = db.transaction([IDBFS.DB_STORE_NAME], 'readwrite');
        var store = transaction.objectStore(IDBFS.DB_STORE_NAME);
  
        function done(err) {
          if (err) {
            if (!done.errored) {
              done.errored = true;
              return callback(err);
            }
            return;
          }
          if (++completed >= total) {
            return callback(null);
          }
        };
  
        transaction.onerror = function() { done(this.error); };
  
        // sort paths in ascending order so directory entries are created
        // before the files inside them
        create.sort().forEach(function (path) {
          if (dst.type === 'local') {
            IDBFS.loadRemoteEntry(store, path, function (err, entry) {
              if (err) return done(err);
              IDBFS.storeLocalEntry(path, entry, done);
            });
          } else {
            IDBFS.loadLocalEntry(path, function (err, entry) {
              if (err) return done(err);
              IDBFS.storeRemoteEntry(store, path, entry, done);
            });
          }
        });
  
        // sort paths in descending order so files are deleted before their
        // parent directories
        remove.sort().reverse().forEach(function(path) {
          if (dst.type === 'local') {
            IDBFS.removeLocalEntry(path, done);
          } else {
            IDBFS.removeRemoteEntry(store, path, done);
          }
        });
      }};
  
  var NODEFS={isWindows:false,staticInit:function () {
        NODEFS.isWindows = !!process.platform.match(/^win/);
      },mount:function (mount) {
        assert(ENVIRONMENT_IS_NODE);
        return NODEFS.createNode(null, '/', NODEFS.getMode(mount.opts.root), 0);
      },createNode:function (parent, name, mode, dev) {
        if (!FS.isDir(mode) && !FS.isFile(mode) && !FS.isLink(mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        var node = FS.createNode(parent, name, mode);
        node.node_ops = NODEFS.node_ops;
        node.stream_ops = NODEFS.stream_ops;
        return node;
      },getMode:function (path) {
        var stat;
        try {
          stat = fs.lstatSync(path);
          if (NODEFS.isWindows) {
            // On Windows, directories return permission bits 'rw-rw-rw-', even though they have 'rwxrwxrwx', so 
            // propagate write bits to execute bits.
            stat.mode = stat.mode | ((stat.mode & 146) >> 1);
          }
        } catch (e) {
          if (!e.code) throw e;
          throw new FS.ErrnoError(ERRNO_CODES[e.code]);
        }
        return stat.mode;
      },realPath:function (node) {
        var parts = [];
        while (node.parent !== node) {
          parts.push(node.name);
          node = node.parent;
        }
        parts.push(node.mount.opts.root);
        parts.reverse();
        return PATH.join.apply(null, parts);
      },flagsToPermissionStringMap:{0:"r",1:"r+",2:"r+",64:"r",65:"r+",66:"r+",129:"rx+",193:"rx+",514:"w+",577:"w",578:"w+",705:"wx",706:"wx+",1024:"a",1025:"a",1026:"a+",1089:"a",1090:"a+",1153:"ax",1154:"ax+",1217:"ax",1218:"ax+",4096:"rs",4098:"rs+"},flagsToPermissionString:function (flags) {
        if (flags in NODEFS.flagsToPermissionStringMap) {
          return NODEFS.flagsToPermissionStringMap[flags];
        } else {
          return flags;
        }
      },node_ops:{getattr:function (node) {
          var path = NODEFS.realPath(node);
          var stat;
          try {
            stat = fs.lstatSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
          // node.js v0.10.20 doesn't report blksize and blocks on Windows. Fake them with default blksize of 4096.
          // See http://support.microsoft.com/kb/140365
          if (NODEFS.isWindows && !stat.blksize) {
            stat.blksize = 4096;
          }
          if (NODEFS.isWindows && !stat.blocks) {
            stat.blocks = (stat.size+stat.blksize-1)/stat.blksize|0;
          }
          return {
            dev: stat.dev,
            ino: stat.ino,
            mode: stat.mode,
            nlink: stat.nlink,
            uid: stat.uid,
            gid: stat.gid,
            rdev: stat.rdev,
            size: stat.size,
            atime: stat.atime,
            mtime: stat.mtime,
            ctime: stat.ctime,
            blksize: stat.blksize,
            blocks: stat.blocks
          };
        },setattr:function (node, attr) {
          var path = NODEFS.realPath(node);
          try {
            if (attr.mode !== undefined) {
              fs.chmodSync(path, attr.mode);
              // update the common node structure mode as well
              node.mode = attr.mode;
            }
            if (attr.timestamp !== undefined) {
              var date = new Date(attr.timestamp);
              fs.utimesSync(path, date, date);
            }
            if (attr.size !== undefined) {
              fs.truncateSync(path, attr.size);
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },lookup:function (parent, name) {
          var path = PATH.join2(NODEFS.realPath(parent), name);
          var mode = NODEFS.getMode(path);
          return NODEFS.createNode(parent, name, mode);
        },mknod:function (parent, name, mode, dev) {
          var node = NODEFS.createNode(parent, name, mode, dev);
          // create the backing node for this in the fs root as well
          var path = NODEFS.realPath(node);
          try {
            if (FS.isDir(node.mode)) {
              fs.mkdirSync(path, node.mode);
            } else {
              fs.writeFileSync(path, '', { mode: node.mode });
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
          return node;
        },rename:function (oldNode, newDir, newName) {
          var oldPath = NODEFS.realPath(oldNode);
          var newPath = PATH.join2(NODEFS.realPath(newDir), newName);
          try {
            fs.renameSync(oldPath, newPath);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },unlink:function (parent, name) {
          var path = PATH.join2(NODEFS.realPath(parent), name);
          try {
            fs.unlinkSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },rmdir:function (parent, name) {
          var path = PATH.join2(NODEFS.realPath(parent), name);
          try {
            fs.rmdirSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },readdir:function (node) {
          var path = NODEFS.realPath(node);
          try {
            return fs.readdirSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },symlink:function (parent, newName, oldPath) {
          var newPath = PATH.join2(NODEFS.realPath(parent), newName);
          try {
            fs.symlinkSync(oldPath, newPath);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },readlink:function (node) {
          var path = NODEFS.realPath(node);
          try {
            return fs.readlinkSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        }},stream_ops:{open:function (stream) {
          var path = NODEFS.realPath(stream.node);
          try {
            if (FS.isFile(stream.node.mode)) {
              stream.nfd = fs.openSync(path, NODEFS.flagsToPermissionString(stream.flags));
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },close:function (stream) {
          try {
            if (FS.isFile(stream.node.mode) && stream.nfd) {
              fs.closeSync(stream.nfd);
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },read:function (stream, buffer, offset, length, position) {
          // FIXME this is terrible.
          var nbuffer = new Buffer(length);
          var res;
          try {
            res = fs.readSync(stream.nfd, nbuffer, 0, length, position);
          } catch (e) {
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
          if (res > 0) {
            for (var i = 0; i < res; i++) {
              buffer[offset + i] = nbuffer[i];
            }
          }
          return res;
        },write:function (stream, buffer, offset, length, position) {
          // FIXME this is terrible.
          var nbuffer = new Buffer(buffer.subarray(offset, offset + length));
          var res;
          try {
            res = fs.writeSync(stream.nfd, nbuffer, 0, length, position);
          } catch (e) {
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
          return res;
        },llseek:function (stream, offset, whence) {
          var position = offset;
          if (whence === 1) {  // SEEK_CUR.
            position += stream.position;
          } else if (whence === 2) {  // SEEK_END.
            if (FS.isFile(stream.node.mode)) {
              try {
                var stat = fs.fstatSync(stream.nfd);
                position += stat.size;
              } catch (e) {
                throw new FS.ErrnoError(ERRNO_CODES[e.code]);
              }
            }
          }
  
          if (position < 0) {
            throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
          }
  
          stream.position = position;
          return position;
        }}};
  
  var _stdin=allocate(1, "i32*", ALLOC_STATIC);
  
  var _stdout=allocate(1, "i32*", ALLOC_STATIC);
  
  var _stderr=allocate(1, "i32*", ALLOC_STATIC);
  
  function _fflush(stream) {
      // int fflush(FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fflush.html
      // we don't currently perform any user-space buffering of data
    }var FS={root:null,mounts:[],devices:[null],streams:[],nextInode:1,nameTable:null,currentPath:"/",initialized:false,ignorePermissions:true,trackingDelegate:{},tracking:{openFlags:{READ:1,WRITE:2}},ErrnoError:null,genericErrors:{},handleFSError:function (e) {
        if (!(e instanceof FS.ErrnoError)) throw e + ' : ' + stackTrace();
        return ___setErrNo(e.errno);
      },lookupPath:function (path, opts) {
        path = PATH.resolve(FS.cwd(), path);
        opts = opts || {};
  
        if (!path) return { path: '', node: null };
  
        var defaults = {
          follow_mount: true,
          recurse_count: 0
        };
        for (var key in defaults) {
          if (opts[key] === undefined) {
            opts[key] = defaults[key];
          }
        }
  
        if (opts.recurse_count > 8) {  // max recursive lookup of 8
          throw new FS.ErrnoError(ERRNO_CODES.ELOOP);
        }
  
        // split the path
        var parts = PATH.normalizeArray(path.split('/').filter(function(p) {
          return !!p;
        }), false);
  
        // start at the root
        var current = FS.root;
        var current_path = '/';
  
        for (var i = 0; i < parts.length; i++) {
          var islast = (i === parts.length-1);
          if (islast && opts.parent) {
            // stop resolving
            break;
          }
  
          current = FS.lookupNode(current, parts[i]);
          current_path = PATH.join2(current_path, parts[i]);
  
          // jump to the mount's root node if this is a mountpoint
          if (FS.isMountpoint(current)) {
            if (!islast || (islast && opts.follow_mount)) {
              current = current.mounted.root;
            }
          }
  
          // by default, lookupPath will not follow a symlink if it is the final path component.
          // setting opts.follow = true will override this behavior.
          if (!islast || opts.follow) {
            var count = 0;
            while (FS.isLink(current.mode)) {
              var link = FS.readlink(current_path);
              current_path = PATH.resolve(PATH.dirname(current_path), link);
              
              var lookup = FS.lookupPath(current_path, { recurse_count: opts.recurse_count });
              current = lookup.node;
  
              if (count++ > 40) {  // limit max consecutive symlinks to 40 (SYMLOOP_MAX).
                throw new FS.ErrnoError(ERRNO_CODES.ELOOP);
              }
            }
          }
        }
  
        return { path: current_path, node: current };
      },getPath:function (node) {
        var path;
        while (true) {
          if (FS.isRoot(node)) {
            var mount = node.mount.mountpoint;
            if (!path) return mount;
            return mount[mount.length-1] !== '/' ? mount + '/' + path : mount + path;
          }
          path = path ? node.name + '/' + path : node.name;
          node = node.parent;
        }
      },hashName:function (parentid, name) {
        var hash = 0;
  
  
        for (var i = 0; i < name.length; i++) {
          hash = ((hash << 5) - hash + name.charCodeAt(i)) | 0;
        }
        return ((parentid + hash) >>> 0) % FS.nameTable.length;
      },hashAddNode:function (node) {
        var hash = FS.hashName(node.parent.id, node.name);
        node.name_next = FS.nameTable[hash];
        FS.nameTable[hash] = node;
      },hashRemoveNode:function (node) {
        var hash = FS.hashName(node.parent.id, node.name);
        if (FS.nameTable[hash] === node) {
          FS.nameTable[hash] = node.name_next;
        } else {
          var current = FS.nameTable[hash];
          while (current) {
            if (current.name_next === node) {
              current.name_next = node.name_next;
              break;
            }
            current = current.name_next;
          }
        }
      },lookupNode:function (parent, name) {
        var err = FS.mayLookup(parent);
        if (err) {
          throw new FS.ErrnoError(err, parent);
        }
        var hash = FS.hashName(parent.id, name);
        for (var node = FS.nameTable[hash]; node; node = node.name_next) {
          var nodeName = node.name;
          if (node.parent.id === parent.id && nodeName === name) {
            return node;
          }
        }
        // if we failed to find it in the cache, call into the VFS
        return FS.lookup(parent, name);
      },createNode:function (parent, name, mode, rdev) {
        if (!FS.FSNode) {
          FS.FSNode = function(parent, name, mode, rdev) {
            if (!parent) {
              parent = this;  // root node sets parent to itself
            }
            this.parent = parent;
            this.mount = parent.mount;
            this.mounted = null;
            this.id = FS.nextInode++;
            this.name = name;
            this.mode = mode;
            this.node_ops = {};
            this.stream_ops = {};
            this.rdev = rdev;
          };
  
          FS.FSNode.prototype = {};
  
          // compatibility
          var readMode = 292 | 73;
          var writeMode = 146;
  
          // NOTE we must use Object.defineProperties instead of individual calls to
          // Object.defineProperty in order to make closure compiler happy
          Object.defineProperties(FS.FSNode.prototype, {
            read: {
              get: function() { return (this.mode & readMode) === readMode; },
              set: function(val) { val ? this.mode |= readMode : this.mode &= ~readMode; }
            },
            write: {
              get: function() { return (this.mode & writeMode) === writeMode; },
              set: function(val) { val ? this.mode |= writeMode : this.mode &= ~writeMode; }
            },
            isFolder: {
              get: function() { return FS.isDir(this.mode); },
            },
            isDevice: {
              get: function() { return FS.isChrdev(this.mode); },
            },
          });
        }
  
        var node = new FS.FSNode(parent, name, mode, rdev);
  
        FS.hashAddNode(node);
  
        return node;
      },destroyNode:function (node) {
        FS.hashRemoveNode(node);
      },isRoot:function (node) {
        return node === node.parent;
      },isMountpoint:function (node) {
        return !!node.mounted;
      },isFile:function (mode) {
        return (mode & 61440) === 32768;
      },isDir:function (mode) {
        return (mode & 61440) === 16384;
      },isLink:function (mode) {
        return (mode & 61440) === 40960;
      },isChrdev:function (mode) {
        return (mode & 61440) === 8192;
      },isBlkdev:function (mode) {
        return (mode & 61440) === 24576;
      },isFIFO:function (mode) {
        return (mode & 61440) === 4096;
      },isSocket:function (mode) {
        return (mode & 49152) === 49152;
      },flagModes:{"r":0,"rs":1052672,"r+":2,"w":577,"wx":705,"xw":705,"w+":578,"wx+":706,"xw+":706,"a":1089,"ax":1217,"xa":1217,"a+":1090,"ax+":1218,"xa+":1218},modeStringToFlags:function (str) {
        var flags = FS.flagModes[str];
        if (typeof flags === 'undefined') {
          throw new Error('Unknown file open mode: ' + str);
        }
        return flags;
      },flagsToPermissionString:function (flag) {
        var accmode = flag & 2097155;
        var perms = ['r', 'w', 'rw'][accmode];
        if ((flag & 512)) {
          perms += 'w';
        }
        return perms;
      },nodePermissions:function (node, perms) {
        if (FS.ignorePermissions) {
          return 0;
        }
        // return 0 if any user, group or owner bits are set.
        if (perms.indexOf('r') !== -1 && !(node.mode & 292)) {
          return ERRNO_CODES.EACCES;
        } else if (perms.indexOf('w') !== -1 && !(node.mode & 146)) {
          return ERRNO_CODES.EACCES;
        } else if (perms.indexOf('x') !== -1 && !(node.mode & 73)) {
          return ERRNO_CODES.EACCES;
        }
        return 0;
      },mayLookup:function (dir) {
        var err = FS.nodePermissions(dir, 'x');
        if (err) return err;
        if (!dir.node_ops.lookup) return ERRNO_CODES.EACCES;
        return 0;
      },mayCreate:function (dir, name) {
        try {
          var node = FS.lookupNode(dir, name);
          return ERRNO_CODES.EEXIST;
        } catch (e) {
        }
        return FS.nodePermissions(dir, 'wx');
      },mayDelete:function (dir, name, isdir) {
        var node;
        try {
          node = FS.lookupNode(dir, name);
        } catch (e) {
          return e.errno;
        }
        var err = FS.nodePermissions(dir, 'wx');
        if (err) {
          return err;
        }
        if (isdir) {
          if (!FS.isDir(node.mode)) {
            return ERRNO_CODES.ENOTDIR;
          }
          if (FS.isRoot(node) || FS.getPath(node) === FS.cwd()) {
            return ERRNO_CODES.EBUSY;
          }
        } else {
          if (FS.isDir(node.mode)) {
            return ERRNO_CODES.EISDIR;
          }
        }
        return 0;
      },mayOpen:function (node, flags) {
        if (!node) {
          return ERRNO_CODES.ENOENT;
        }
        if (FS.isLink(node.mode)) {
          return ERRNO_CODES.ELOOP;
        } else if (FS.isDir(node.mode)) {
          if ((flags & 2097155) !== 0 ||  // opening for write
              (flags & 512)) {
            return ERRNO_CODES.EISDIR;
          }
        }
        return FS.nodePermissions(node, FS.flagsToPermissionString(flags));
      },MAX_OPEN_FDS:4096,nextfd:function (fd_start, fd_end) {
        fd_start = fd_start || 0;
        fd_end = fd_end || FS.MAX_OPEN_FDS;
        for (var fd = fd_start; fd <= fd_end; fd++) {
          if (!FS.streams[fd]) {
            return fd;
          }
        }
        throw new FS.ErrnoError(ERRNO_CODES.EMFILE);
      },getStream:function (fd) {
        return FS.streams[fd];
      },createStream:function (stream, fd_start, fd_end) {
        if (!FS.FSStream) {
          FS.FSStream = function(){};
          FS.FSStream.prototype = {};
          // compatibility
          Object.defineProperties(FS.FSStream.prototype, {
            object: {
              get: function() { return this.node; },
              set: function(val) { this.node = val; }
            },
            isRead: {
              get: function() { return (this.flags & 2097155) !== 1; }
            },
            isWrite: {
              get: function() { return (this.flags & 2097155) !== 0; }
            },
            isAppend: {
              get: function() { return (this.flags & 1024); }
            }
          });
        }
        // clone it, so we can return an instance of FSStream
        var newStream = new FS.FSStream();
        for (var p in stream) {
          newStream[p] = stream[p];
        }
        stream = newStream;
        var fd = FS.nextfd(fd_start, fd_end);
        stream.fd = fd;
        FS.streams[fd] = stream;
        return stream;
      },closeStream:function (fd) {
        FS.streams[fd] = null;
      },getStreamFromPtr:function (ptr) {
        return FS.streams[ptr - 1];
      },getPtrForStream:function (stream) {
        return stream ? stream.fd + 1 : 0;
      },chrdev_stream_ops:{open:function (stream) {
          var device = FS.getDevice(stream.node.rdev);
          // override node's stream ops with the device's
          stream.stream_ops = device.stream_ops;
          // forward the open call
          if (stream.stream_ops.open) {
            stream.stream_ops.open(stream);
          }
        },llseek:function () {
          throw new FS.ErrnoError(ERRNO_CODES.ESPIPE);
        }},major:function (dev) {
        return ((dev) >> 8);
      },minor:function (dev) {
        return ((dev) & 0xff);
      },makedev:function (ma, mi) {
        return ((ma) << 8 | (mi));
      },registerDevice:function (dev, ops) {
        FS.devices[dev] = { stream_ops: ops };
      },getDevice:function (dev) {
        return FS.devices[dev];
      },getMounts:function (mount) {
        var mounts = [];
        var check = [mount];
  
        while (check.length) {
          var m = check.pop();
  
          mounts.push(m);
  
          check.push.apply(check, m.mounts);
        }
  
        return mounts;
      },syncfs:function (populate, callback) {
        if (typeof(populate) === 'function') {
          callback = populate;
          populate = false;
        }
  
        var mounts = FS.getMounts(FS.root.mount);
        var completed = 0;
  
        function done(err) {
          if (err) {
            if (!done.errored) {
              done.errored = true;
              return callback(err);
            }
            return;
          }
          if (++completed >= mounts.length) {
            callback(null);
          }
        };
  
        // sync all mounts
        mounts.forEach(function (mount) {
          if (!mount.type.syncfs) {
            return done(null);
          }
          mount.type.syncfs(mount, populate, done);
        });
      },mount:function (type, opts, mountpoint) {
        var root = mountpoint === '/';
        var pseudo = !mountpoint;
        var node;
  
        if (root && FS.root) {
          throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
        } else if (!root && !pseudo) {
          var lookup = FS.lookupPath(mountpoint, { follow_mount: false });
  
          mountpoint = lookup.path;  // use the absolute path
          node = lookup.node;
  
          if (FS.isMountpoint(node)) {
            throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
          }
  
          if (!FS.isDir(node.mode)) {
            throw new FS.ErrnoError(ERRNO_CODES.ENOTDIR);
          }
        }
  
        var mount = {
          type: type,
          opts: opts,
          mountpoint: mountpoint,
          mounts: []
        };
  
        // create a root node for the fs
        var mountRoot = type.mount(mount);
        mountRoot.mount = mount;
        mount.root = mountRoot;
  
        if (root) {
          FS.root = mountRoot;
        } else if (node) {
          // set as a mountpoint
          node.mounted = mount;
  
          // add the new mount to the current mount's children
          if (node.mount) {
            node.mount.mounts.push(mount);
          }
        }
  
        return mountRoot;
      },unmount:function (mountpoint) {
        var lookup = FS.lookupPath(mountpoint, { follow_mount: false });
  
        if (!FS.isMountpoint(lookup.node)) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
  
        // destroy the nodes for this mount, and all its child mounts
        var node = lookup.node;
        var mount = node.mounted;
        var mounts = FS.getMounts(mount);
  
        Object.keys(FS.nameTable).forEach(function (hash) {
          var current = FS.nameTable[hash];
  
          while (current) {
            var next = current.name_next;
  
            if (mounts.indexOf(current.mount) !== -1) {
              FS.destroyNode(current);
            }
  
            current = next;
          }
        });
  
        // no longer a mountpoint
        node.mounted = null;
  
        // remove this mount from the child mounts
        var idx = node.mount.mounts.indexOf(mount);
        assert(idx !== -1);
        node.mount.mounts.splice(idx, 1);
      },lookup:function (parent, name) {
        return parent.node_ops.lookup(parent, name);
      },mknod:function (path, mode, dev) {
        var lookup = FS.lookupPath(path, { parent: true });
        var parent = lookup.node;
        var name = PATH.basename(path);
        if (!name || name === '.' || name === '..') {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        var err = FS.mayCreate(parent, name);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        if (!parent.node_ops.mknod) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        return parent.node_ops.mknod(parent, name, mode, dev);
      },create:function (path, mode) {
        mode = mode !== undefined ? mode : 438 /* 0666 */;
        mode &= 4095;
        mode |= 32768;
        return FS.mknod(path, mode, 0);
      },mkdir:function (path, mode) {
        mode = mode !== undefined ? mode : 511 /* 0777 */;
        mode &= 511 | 512;
        mode |= 16384;
        return FS.mknod(path, mode, 0);
      },mkdev:function (path, mode, dev) {
        if (typeof(dev) === 'undefined') {
          dev = mode;
          mode = 438 /* 0666 */;
        }
        mode |= 8192;
        return FS.mknod(path, mode, dev);
      },symlink:function (oldpath, newpath) {
        if (!PATH.resolve(oldpath)) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOENT);
        }
        var lookup = FS.lookupPath(newpath, { parent: true });
        var parent = lookup.node;
        if (!parent) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOENT);
        }
        var newname = PATH.basename(newpath);
        var err = FS.mayCreate(parent, newname);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        if (!parent.node_ops.symlink) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        return parent.node_ops.symlink(parent, newname, oldpath);
      },rename:function (old_path, new_path) {
        var old_dirname = PATH.dirname(old_path);
        var new_dirname = PATH.dirname(new_path);
        var old_name = PATH.basename(old_path);
        var new_name = PATH.basename(new_path);
        // parents must exist
        var lookup, old_dir, new_dir;
        try {
          lookup = FS.lookupPath(old_path, { parent: true });
          old_dir = lookup.node;
          lookup = FS.lookupPath(new_path, { parent: true });
          new_dir = lookup.node;
        } catch (e) {
          throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
        }
        if (!old_dir || !new_dir) throw new FS.ErrnoError(ERRNO_CODES.ENOENT);
        // need to be part of the same mount
        if (old_dir.mount !== new_dir.mount) {
          throw new FS.ErrnoError(ERRNO_CODES.EXDEV);
        }
        // source must exist
        var old_node = FS.lookupNode(old_dir, old_name);
        // old path should not be an ancestor of the new path
        var relative = PATH.relative(old_path, new_dirname);
        if (relative.charAt(0) !== '.') {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        // new path should not be an ancestor of the old path
        relative = PATH.relative(new_path, old_dirname);
        if (relative.charAt(0) !== '.') {
          throw new FS.ErrnoError(ERRNO_CODES.ENOTEMPTY);
        }
        // see if the new path already exists
        var new_node;
        try {
          new_node = FS.lookupNode(new_dir, new_name);
        } catch (e) {
          // not fatal
        }
        // early out if nothing needs to change
        if (old_node === new_node) {
          return;
        }
        // we'll need to delete the old entry
        var isdir = FS.isDir(old_node.mode);
        var err = FS.mayDelete(old_dir, old_name, isdir);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        // need delete permissions if we'll be overwriting.
        // need create permissions if new doesn't already exist.
        err = new_node ?
          FS.mayDelete(new_dir, new_name, isdir) :
          FS.mayCreate(new_dir, new_name);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        if (!old_dir.node_ops.rename) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        if (FS.isMountpoint(old_node) || (new_node && FS.isMountpoint(new_node))) {
          throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
        }
        // if we are going to change the parent, check write permissions
        if (new_dir !== old_dir) {
          err = FS.nodePermissions(old_dir, 'w');
          if (err) {
            throw new FS.ErrnoError(err);
          }
        }
        try {
          if (FS.trackingDelegate['willMovePath']) {
            FS.trackingDelegate['willMovePath'](old_path, new_path);
          }
        } catch(e) {
          console.log("FS.trackingDelegate['willMovePath']('"+old_path+"', '"+new_path+"') threw an exception: " + e.message);
        }
        // remove the node from the lookup hash
        FS.hashRemoveNode(old_node);
        // do the underlying fs rename
        try {
          old_dir.node_ops.rename(old_node, new_dir, new_name);
        } catch (e) {
          throw e;
        } finally {
          // add the node back to the hash (in case node_ops.rename
          // changed its name)
          FS.hashAddNode(old_node);
        }
        try {
          if (FS.trackingDelegate['onMovePath']) FS.trackingDelegate['onMovePath'](old_path, new_path);
        } catch(e) {
          console.log("FS.trackingDelegate['onMovePath']('"+old_path+"', '"+new_path+"') threw an exception: " + e.message);
        }
      },rmdir:function (path) {
        var lookup = FS.lookupPath(path, { parent: true });
        var parent = lookup.node;
        var name = PATH.basename(path);
        var node = FS.lookupNode(parent, name);
        var err = FS.mayDelete(parent, name, true);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        if (!parent.node_ops.rmdir) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        if (FS.isMountpoint(node)) {
          throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
        }
        try {
          if (FS.trackingDelegate['willDeletePath']) {
            FS.trackingDelegate['willDeletePath'](path);
          }
        } catch(e) {
          console.log("FS.trackingDelegate['willDeletePath']('"+path+"') threw an exception: " + e.message);
        }
        parent.node_ops.rmdir(parent, name);
        FS.destroyNode(node);
        try {
          if (FS.trackingDelegate['onDeletePath']) FS.trackingDelegate['onDeletePath'](path);
        } catch(e) {
          console.log("FS.trackingDelegate['onDeletePath']('"+path+"') threw an exception: " + e.message);
        }
      },readdir:function (path) {
        var lookup = FS.lookupPath(path, { follow: true });
        var node = lookup.node;
        if (!node.node_ops.readdir) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOTDIR);
        }
        return node.node_ops.readdir(node);
      },unlink:function (path) {
        var lookup = FS.lookupPath(path, { parent: true });
        var parent = lookup.node;
        var name = PATH.basename(path);
        var node = FS.lookupNode(parent, name);
        var err = FS.mayDelete(parent, name, false);
        if (err) {
          // POSIX says unlink should set EPERM, not EISDIR
          if (err === ERRNO_CODES.EISDIR) err = ERRNO_CODES.EPERM;
          throw new FS.ErrnoError(err);
        }
        if (!parent.node_ops.unlink) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        if (FS.isMountpoint(node)) {
          throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
        }
        try {
          if (FS.trackingDelegate['willDeletePath']) {
            FS.trackingDelegate['willDeletePath'](path);
          }
        } catch(e) {
          console.log("FS.trackingDelegate['willDeletePath']('"+path+"') threw an exception: " + e.message);
        }
        parent.node_ops.unlink(parent, name);
        FS.destroyNode(node);
        try {
          if (FS.trackingDelegate['onDeletePath']) FS.trackingDelegate['onDeletePath'](path);
        } catch(e) {
          console.log("FS.trackingDelegate['onDeletePath']('"+path+"') threw an exception: " + e.message);
        }
      },readlink:function (path) {
        var lookup = FS.lookupPath(path);
        var link = lookup.node;
        if (!link) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOENT);
        }
        if (!link.node_ops.readlink) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        return link.node_ops.readlink(link);
      },stat:function (path, dontFollow) {
        var lookup = FS.lookupPath(path, { follow: !dontFollow });
        var node = lookup.node;
        if (!node) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOENT);
        }
        if (!node.node_ops.getattr) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        return node.node_ops.getattr(node);
      },lstat:function (path) {
        return FS.stat(path, true);
      },chmod:function (path, mode, dontFollow) {
        var node;
        if (typeof path === 'string') {
          var lookup = FS.lookupPath(path, { follow: !dontFollow });
          node = lookup.node;
        } else {
          node = path;
        }
        if (!node.node_ops.setattr) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        node.node_ops.setattr(node, {
          mode: (mode & 4095) | (node.mode & ~4095),
          timestamp: Date.now()
        });
      },lchmod:function (path, mode) {
        FS.chmod(path, mode, true);
      },fchmod:function (fd, mode) {
        var stream = FS.getStream(fd);
        if (!stream) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        FS.chmod(stream.node, mode);
      },chown:function (path, uid, gid, dontFollow) {
        var node;
        if (typeof path === 'string') {
          var lookup = FS.lookupPath(path, { follow: !dontFollow });
          node = lookup.node;
        } else {
          node = path;
        }
        if (!node.node_ops.setattr) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        node.node_ops.setattr(node, {
          timestamp: Date.now()
          // we ignore the uid / gid for now
        });
      },lchown:function (path, uid, gid) {
        FS.chown(path, uid, gid, true);
      },fchown:function (fd, uid, gid) {
        var stream = FS.getStream(fd);
        if (!stream) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        FS.chown(stream.node, uid, gid);
      },truncate:function (path, len) {
        if (len < 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        var node;
        if (typeof path === 'string') {
          var lookup = FS.lookupPath(path, { follow: true });
          node = lookup.node;
        } else {
          node = path;
        }
        if (!node.node_ops.setattr) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        if (FS.isDir(node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.EISDIR);
        }
        if (!FS.isFile(node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        var err = FS.nodePermissions(node, 'w');
        if (err) {
          throw new FS.ErrnoError(err);
        }
        node.node_ops.setattr(node, {
          size: len,
          timestamp: Date.now()
        });
      },ftruncate:function (fd, len) {
        var stream = FS.getStream(fd);
        if (!stream) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        if ((stream.flags & 2097155) === 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        FS.truncate(stream.node, len);
      },utime:function (path, atime, mtime) {
        var lookup = FS.lookupPath(path, { follow: true });
        var node = lookup.node;
        node.node_ops.setattr(node, {
          timestamp: Math.max(atime, mtime)
        });
      },open:function (path, flags, mode, fd_start, fd_end) {
        if (path === "") {
          throw new FS.ErrnoError(ERRNO_CODES.ENOENT);
        }
        flags = typeof flags === 'string' ? FS.modeStringToFlags(flags) : flags;
        mode = typeof mode === 'undefined' ? 438 /* 0666 */ : mode;
        if ((flags & 64)) {
          mode = (mode & 4095) | 32768;
        } else {
          mode = 0;
        }
        var node;
        if (typeof path === 'object') {
          node = path;
        } else {
          path = PATH.normalize(path);
          try {
            var lookup = FS.lookupPath(path, {
              follow: !(flags & 131072)
            });
            node = lookup.node;
          } catch (e) {
            // ignore
          }
        }
        // perhaps we need to create the node
        var created = false;
        if ((flags & 64)) {
          if (node) {
            // if O_CREAT and O_EXCL are set, error out if the node already exists
            if ((flags & 128)) {
              throw new FS.ErrnoError(ERRNO_CODES.EEXIST);
            }
          } else {
            // node doesn't exist, try to create it
            node = FS.mknod(path, mode, 0);
            created = true;
          }
        }
        if (!node) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOENT);
        }
        // can't truncate a device
        if (FS.isChrdev(node.mode)) {
          flags &= ~512;
        }
        // check permissions, if this is not a file we just created now (it is ok to
        // create and write to a file with read-only permissions; it is read-only
        // for later use)
        if (!created) {
          var err = FS.mayOpen(node, flags);
          if (err) {
            throw new FS.ErrnoError(err);
          }
        }
        // do truncation if necessary
        if ((flags & 512)) {
          FS.truncate(node, 0);
        }
        // we've already handled these, don't pass down to the underlying vfs
        flags &= ~(128 | 512);
  
        // register the stream with the filesystem
        var stream = FS.createStream({
          node: node,
          path: FS.getPath(node),  // we want the absolute path to the node
          flags: flags,
          seekable: true,
          position: 0,
          stream_ops: node.stream_ops,
          // used by the file family libc calls (fopen, fwrite, ferror, etc.)
          ungotten: [],
          error: false
        }, fd_start, fd_end);
        // call the new stream's open function
        if (stream.stream_ops.open) {
          stream.stream_ops.open(stream);
        }
        if (Module['logReadFiles'] && !(flags & 1)) {
          if (!FS.readFiles) FS.readFiles = {};
          if (!(path in FS.readFiles)) {
            FS.readFiles[path] = 1;
            Module['printErr']('read file: ' + path);
          }
        }
        try {
          if (FS.trackingDelegate['onOpenFile']) {
            var trackingFlags = 0;
            if ((flags & 2097155) !== 1) {
              trackingFlags |= FS.tracking.openFlags.READ;
            }
            if ((flags & 2097155) !== 0) {
              trackingFlags |= FS.tracking.openFlags.WRITE;
            }
            FS.trackingDelegate['onOpenFile'](path, trackingFlags);
          }
        } catch(e) {
          console.log("FS.trackingDelegate['onOpenFile']('"+path+"', flags) threw an exception: " + e.message);
        }
        return stream;
      },close:function (stream) {
        try {
          if (stream.stream_ops.close) {
            stream.stream_ops.close(stream);
          }
        } catch (e) {
          throw e;
        } finally {
          FS.closeStream(stream.fd);
        }
      },llseek:function (stream, offset, whence) {
        if (!stream.seekable || !stream.stream_ops.llseek) {
          throw new FS.ErrnoError(ERRNO_CODES.ESPIPE);
        }
        return stream.stream_ops.llseek(stream, offset, whence);
      },read:function (stream, buffer, offset, length, position) {
        if (length < 0 || position < 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        if ((stream.flags & 2097155) === 1) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        if (FS.isDir(stream.node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.EISDIR);
        }
        if (!stream.stream_ops.read) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        var seeking = true;
        if (typeof position === 'undefined') {
          position = stream.position;
          seeking = false;
        } else if (!stream.seekable) {
          throw new FS.ErrnoError(ERRNO_CODES.ESPIPE);
        }
        var bytesRead = stream.stream_ops.read(stream, buffer, offset, length, position);
        if (!seeking) stream.position += bytesRead;
        return bytesRead;
      },write:function (stream, buffer, offset, length, position, canOwn) {
        if (length < 0 || position < 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        if ((stream.flags & 2097155) === 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        if (FS.isDir(stream.node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.EISDIR);
        }
        if (!stream.stream_ops.write) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        if (stream.flags & 1024) {
          // seek to the end before writing in append mode
          FS.llseek(stream, 0, 2);
        }
        var seeking = true;
        if (typeof position === 'undefined') {
          position = stream.position;
          seeking = false;
        } else if (!stream.seekable) {
          throw new FS.ErrnoError(ERRNO_CODES.ESPIPE);
        }
        var bytesWritten = stream.stream_ops.write(stream, buffer, offset, length, position, canOwn);
        if (!seeking) stream.position += bytesWritten;
        try {
          if (stream.path && FS.trackingDelegate['onWriteToFile']) FS.trackingDelegate['onWriteToFile'](stream.path);
        } catch(e) {
          console.log("FS.trackingDelegate['onWriteToFile']('"+path+"') threw an exception: " + e.message);
        }
        return bytesWritten;
      },allocate:function (stream, offset, length) {
        if (offset < 0 || length <= 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        if ((stream.flags & 2097155) === 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        if (!FS.isFile(stream.node.mode) && !FS.isDir(node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.ENODEV);
        }
        if (!stream.stream_ops.allocate) {
          throw new FS.ErrnoError(ERRNO_CODES.EOPNOTSUPP);
        }
        stream.stream_ops.allocate(stream, offset, length);
      },mmap:function (stream, buffer, offset, length, position, prot, flags) {
        // TODO if PROT is PROT_WRITE, make sure we have write access
        if ((stream.flags & 2097155) === 1) {
          throw new FS.ErrnoError(ERRNO_CODES.EACCES);
        }
        if (!stream.stream_ops.mmap) {
          throw new FS.ErrnoError(ERRNO_CODES.ENODEV);
        }
        return stream.stream_ops.mmap(stream, buffer, offset, length, position, prot, flags);
      },ioctl:function (stream, cmd, arg) {
        if (!stream.stream_ops.ioctl) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOTTY);
        }
        return stream.stream_ops.ioctl(stream, cmd, arg);
      },readFile:function (path, opts) {
        opts = opts || {};
        opts.flags = opts.flags || 'r';
        opts.encoding = opts.encoding || 'binary';
        if (opts.encoding !== 'utf8' && opts.encoding !== 'binary') {
          throw new Error('Invalid encoding type "' + opts.encoding + '"');
        }
        var ret;
        var stream = FS.open(path, opts.flags);
        var stat = FS.stat(path);
        var length = stat.size;
        var buf = new Uint8Array(length);
        FS.read(stream, buf, 0, length, 0);
        if (opts.encoding === 'utf8') {
          ret = '';
          var utf8 = new Runtime.UTF8Processor();
          for (var i = 0; i < length; i++) {
            ret += utf8.processCChar(buf[i]);
          }
        } else if (opts.encoding === 'binary') {
          ret = buf;
        }
        FS.close(stream);
        return ret;
      },writeFile:function (path, data, opts) {
        opts = opts || {};
        opts.flags = opts.flags || 'w';
        opts.encoding = opts.encoding || 'utf8';
        if (opts.encoding !== 'utf8' && opts.encoding !== 'binary') {
          throw new Error('Invalid encoding type "' + opts.encoding + '"');
        }
        var stream = FS.open(path, opts.flags, opts.mode);
        if (opts.encoding === 'utf8') {
          var utf8 = new Runtime.UTF8Processor();
          var buf = new Uint8Array(utf8.processJSString(data));
          FS.write(stream, buf, 0, buf.length, 0, opts.canOwn);
        } else if (opts.encoding === 'binary') {
          FS.write(stream, data, 0, data.length, 0, opts.canOwn);
        }
        FS.close(stream);
      },cwd:function () {
        return FS.currentPath;
      },chdir:function (path) {
        var lookup = FS.lookupPath(path, { follow: true });
        if (!FS.isDir(lookup.node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOTDIR);
        }
        var err = FS.nodePermissions(lookup.node, 'x');
        if (err) {
          throw new FS.ErrnoError(err);
        }
        FS.currentPath = lookup.path;
      },createDefaultDirectories:function () {
        FS.mkdir('/tmp');
        FS.mkdir('/home');
        FS.mkdir('/home/web_user');
      },createDefaultDevices:function () {
        // create /dev
        FS.mkdir('/dev');
        // setup /dev/null
        FS.registerDevice(FS.makedev(1, 3), {
          read: function() { return 0; },
          write: function() { return 0; }
        });
        FS.mkdev('/dev/null', FS.makedev(1, 3));
        // setup /dev/tty and /dev/tty1
        // stderr needs to print output using Module['printErr']
        // so we register a second tty just for it.
        TTY.register(FS.makedev(5, 0), TTY.default_tty_ops);
        TTY.register(FS.makedev(6, 0), TTY.default_tty1_ops);
        FS.mkdev('/dev/tty', FS.makedev(5, 0));
        FS.mkdev('/dev/tty1', FS.makedev(6, 0));
        // setup /dev/[u]random
        var random_device;
        if (typeof crypto !== 'undefined') {
          // for modern web browsers
          var randomBuffer = new Uint8Array(1);
          random_device = function() { crypto.getRandomValues(randomBuffer); return randomBuffer[0]; };
        } else if (ENVIRONMENT_IS_NODE) {
          // for nodejs
          random_device = function() { return require('crypto').randomBytes(1)[0]; };
        } else {
          // default for ES5 platforms
          random_device = function() { return (Math.random()*256)|0; };
        }
        FS.createDevice('/dev', 'random', random_device);
        FS.createDevice('/dev', 'urandom', random_device);
        // we're not going to emulate the actual shm device,
        // just create the tmp dirs that reside in it commonly
        FS.mkdir('/dev/shm');
        FS.mkdir('/dev/shm/tmp');
      },createStandardStreams:function () {
        // TODO deprecate the old functionality of a single
        // input / output callback and that utilizes FS.createDevice
        // and instead require a unique set of stream ops
  
        // by default, we symlink the standard streams to the
        // default tty devices. however, if the standard streams
        // have been overwritten we create a unique device for
        // them instead.
        if (Module['stdin']) {
          FS.createDevice('/dev', 'stdin', Module['stdin']);
        } else {
          FS.symlink('/dev/tty', '/dev/stdin');
        }
        if (Module['stdout']) {
          FS.createDevice('/dev', 'stdout', null, Module['stdout']);
        } else {
          FS.symlink('/dev/tty', '/dev/stdout');
        }
        if (Module['stderr']) {
          FS.createDevice('/dev', 'stderr', null, Module['stderr']);
        } else {
          FS.symlink('/dev/tty1', '/dev/stderr');
        }
  
        // open default streams for the stdin, stdout and stderr devices
        var stdin = FS.open('/dev/stdin', 'r');
        HEAP32[((_stdin)>>2)]=FS.getPtrForStream(stdin);
        assert(stdin.fd === 0, 'invalid handle for stdin (' + stdin.fd + ')');
  
        var stdout = FS.open('/dev/stdout', 'w');
        HEAP32[((_stdout)>>2)]=FS.getPtrForStream(stdout);
        assert(stdout.fd === 1, 'invalid handle for stdout (' + stdout.fd + ')');
  
        var stderr = FS.open('/dev/stderr', 'w');
        HEAP32[((_stderr)>>2)]=FS.getPtrForStream(stderr);
        assert(stderr.fd === 2, 'invalid handle for stderr (' + stderr.fd + ')');
      },ensureErrnoError:function () {
        if (FS.ErrnoError) return;
        FS.ErrnoError = function ErrnoError(errno, node) {
          this.node = node;
          this.setErrno = function(errno) {
            this.errno = errno;
            for (var key in ERRNO_CODES) {
              if (ERRNO_CODES[key] === errno) {
                this.code = key;
                break;
              }
            }
          };
          this.setErrno(errno);
          this.message = ERRNO_MESSAGES[errno];
        };
        FS.ErrnoError.prototype = new Error();
        FS.ErrnoError.prototype.constructor = FS.ErrnoError;
        // Some errors may happen quite a bit, to avoid overhead we reuse them (and suffer a lack of stack info)
        [ERRNO_CODES.ENOENT].forEach(function(code) {
          FS.genericErrors[code] = new FS.ErrnoError(code);
          FS.genericErrors[code].stack = '<generic error, no stack>';
        });
      },staticInit:function () {
        FS.ensureErrnoError();
  
        FS.nameTable = new Array(4096);
  
        FS.mount(MEMFS, {}, '/');
  
        FS.createDefaultDirectories();
        FS.createDefaultDevices();
      },init:function (input, output, error) {
        assert(!FS.init.initialized, 'FS.init was previously called. If you want to initialize later with custom parameters, remove any earlier calls (note that one is automatically added to the generated code)');
        FS.init.initialized = true;
  
        FS.ensureErrnoError();
  
        // Allow Module.stdin etc. to provide defaults, if none explicitly passed to us here
        Module['stdin'] = input || Module['stdin'];
        Module['stdout'] = output || Module['stdout'];
        Module['stderr'] = error || Module['stderr'];
  
        FS.createStandardStreams();
      },quit:function () {
        FS.init.initialized = false;
        for (var i = 0; i < FS.streams.length; i++) {
          var stream = FS.streams[i];
          if (!stream) {
            continue;
          }
          FS.close(stream);
        }
      },getMode:function (canRead, canWrite) {
        var mode = 0;
        if (canRead) mode |= 292 | 73;
        if (canWrite) mode |= 146;
        return mode;
      },joinPath:function (parts, forceRelative) {
        var path = PATH.join.apply(null, parts);
        if (forceRelative && path[0] == '/') path = path.substr(1);
        return path;
      },absolutePath:function (relative, base) {
        return PATH.resolve(base, relative);
      },standardizePath:function (path) {
        return PATH.normalize(path);
      },findObject:function (path, dontResolveLastLink) {
        var ret = FS.analyzePath(path, dontResolveLastLink);
        if (ret.exists) {
          return ret.object;
        } else {
          ___setErrNo(ret.error);
          return null;
        }
      },analyzePath:function (path, dontResolveLastLink) {
        // operate from within the context of the symlink's target
        try {
          var lookup = FS.lookupPath(path, { follow: !dontResolveLastLink });
          path = lookup.path;
        } catch (e) {
        }
        var ret = {
          isRoot: false, exists: false, error: 0, name: null, path: null, object: null,
          parentExists: false, parentPath: null, parentObject: null
        };
        try {
          var lookup = FS.lookupPath(path, { parent: true });
          ret.parentExists = true;
          ret.parentPath = lookup.path;
          ret.parentObject = lookup.node;
          ret.name = PATH.basename(path);
          lookup = FS.lookupPath(path, { follow: !dontResolveLastLink });
          ret.exists = true;
          ret.path = lookup.path;
          ret.object = lookup.node;
          ret.name = lookup.node.name;
          ret.isRoot = lookup.path === '/';
        } catch (e) {
          ret.error = e.errno;
        };
        return ret;
      },createFolder:function (parent, name, canRead, canWrite) {
        var path = PATH.join2(typeof parent === 'string' ? parent : FS.getPath(parent), name);
        var mode = FS.getMode(canRead, canWrite);
        return FS.mkdir(path, mode);
      },createPath:function (parent, path, canRead, canWrite) {
        parent = typeof parent === 'string' ? parent : FS.getPath(parent);
        var parts = path.split('/').reverse();
        while (parts.length) {
          var part = parts.pop();
          if (!part) continue;
          var current = PATH.join2(parent, part);
          try {
            FS.mkdir(current);
          } catch (e) {
            // ignore EEXIST
          }
          parent = current;
        }
        return current;
      },createFile:function (parent, name, properties, canRead, canWrite) {
        var path = PATH.join2(typeof parent === 'string' ? parent : FS.getPath(parent), name);
        var mode = FS.getMode(canRead, canWrite);
        return FS.create(path, mode);
      },createDataFile:function (parent, name, data, canRead, canWrite, canOwn) {
        var path = name ? PATH.join2(typeof parent === 'string' ? parent : FS.getPath(parent), name) : parent;
        var mode = FS.getMode(canRead, canWrite);
        var node = FS.create(path, mode);
        if (data) {
          if (typeof data === 'string') {
            var arr = new Array(data.length);
            for (var i = 0, len = data.length; i < len; ++i) arr[i] = data.charCodeAt(i);
            data = arr;
          }
          // make sure we can write to the file
          FS.chmod(node, mode | 146);
          var stream = FS.open(node, 'w');
          FS.write(stream, data, 0, data.length, 0, canOwn);
          FS.close(stream);
          FS.chmod(node, mode);
        }
        return node;
      },createDevice:function (parent, name, input, output) {
        var path = PATH.join2(typeof parent === 'string' ? parent : FS.getPath(parent), name);
        var mode = FS.getMode(!!input, !!output);
        if (!FS.createDevice.major) FS.createDevice.major = 64;
        var dev = FS.makedev(FS.createDevice.major++, 0);
        // Create a fake device that a set of stream ops to emulate
        // the old behavior.
        FS.registerDevice(dev, {
          open: function(stream) {
            stream.seekable = false;
          },
          close: function(stream) {
            // flush any pending line data
            if (output && output.buffer && output.buffer.length) {
              output(10);
            }
          },
          read: function(stream, buffer, offset, length, pos /* ignored */) {
            var bytesRead = 0;
            for (var i = 0; i < length; i++) {
              var result;
              try {
                result = input();
              } catch (e) {
                throw new FS.ErrnoError(ERRNO_CODES.EIO);
              }
              if (result === undefined && bytesRead === 0) {
                throw new FS.ErrnoError(ERRNO_CODES.EAGAIN);
              }
              if (result === null || result === undefined) break;
              bytesRead++;
              buffer[offset+i] = result;
            }
            if (bytesRead) {
              stream.node.timestamp = Date.now();
            }
            return bytesRead;
          },
          write: function(stream, buffer, offset, length, pos) {
            for (var i = 0; i < length; i++) {
              try {
                output(buffer[offset+i]);
              } catch (e) {
                throw new FS.ErrnoError(ERRNO_CODES.EIO);
              }
            }
            if (length) {
              stream.node.timestamp = Date.now();
            }
            return i;
          }
        });
        return FS.mkdev(path, mode, dev);
      },createLink:function (parent, name, target, canRead, canWrite) {
        var path = PATH.join2(typeof parent === 'string' ? parent : FS.getPath(parent), name);
        return FS.symlink(target, path);
      },forceLoadFile:function (obj) {
        if (obj.isDevice || obj.isFolder || obj.link || obj.contents) return true;
        var success = true;
        if (typeof XMLHttpRequest !== 'undefined') {
          throw new Error("Lazy loading should have been performed (contents set) in createLazyFile, but it was not. Lazy loading only works in web workers. Use --embed-file or --preload-file in emcc on the main thread.");
        } else if (Module['read']) {
          // Command-line.
          try {
            // WARNING: Can't read binary files in V8's d8 or tracemonkey's js, as
            //          read() will try to parse UTF8.
            obj.contents = intArrayFromString(Module['read'](obj.url), true);
            obj.usedBytes = obj.contents.length;
          } catch (e) {
            success = false;
          }
        } else {
          throw new Error('Cannot load without read() or XMLHttpRequest.');
        }
        if (!success) ___setErrNo(ERRNO_CODES.EIO);
        return success;
      },createLazyFile:function (parent, name, url, canRead, canWrite) {
        // Lazy chunked Uint8Array (implements get and length from Uint8Array). Actual getting is abstracted away for eventual reuse.
        function LazyUint8Array() {
          this.lengthKnown = false;
          this.chunks = []; // Loaded chunks. Index is the chunk number
        }
        LazyUint8Array.prototype.get = function LazyUint8Array_get(idx) {
          if (idx > this.length-1 || idx < 0) {
            return undefined;
          }
          var chunkOffset = idx % this.chunkSize;
          var chunkNum = (idx / this.chunkSize)|0;
          return this.getter(chunkNum)[chunkOffset];
        }
        LazyUint8Array.prototype.setDataGetter = function LazyUint8Array_setDataGetter(getter) {
          this.getter = getter;
        }
        LazyUint8Array.prototype.cacheLength = function LazyUint8Array_cacheLength() {
          // Find length
          var xhr = new XMLHttpRequest();
          xhr.open('HEAD', url, false);
          xhr.send(null);
          if (!(xhr.status >= 200 && xhr.status < 300 || xhr.status === 304)) throw new Error("Couldn't load " + url + ". Status: " + xhr.status);
          var datalength = Number(xhr.getResponseHeader("Content-length"));
          var header;
          var hasByteServing = (header = xhr.getResponseHeader("Accept-Ranges")) && header === "bytes";
          var chunkSize = 1024*1024; // Chunk size in bytes
  
          if (!hasByteServing) chunkSize = datalength;
  
          // Function to get a range from the remote URL.
          var doXHR = (function(from, to) {
            if (from > to) throw new Error("invalid range (" + from + ", " + to + ") or no bytes requested!");
            if (to > datalength-1) throw new Error("only " + datalength + " bytes available! programmer error!");
  
            // TODO: Use mozResponseArrayBuffer, responseStream, etc. if available.
            var xhr = new XMLHttpRequest();
            xhr.open('GET', url, false);
            if (datalength !== chunkSize) xhr.setRequestHeader("Range", "bytes=" + from + "-" + to);
  
            // Some hints to the browser that we want binary data.
            if (typeof Uint8Array != 'undefined') xhr.responseType = 'arraybuffer';
            if (xhr.overrideMimeType) {
              xhr.overrideMimeType('text/plain; charset=x-user-defined');
            }
  
            xhr.send(null);
            if (!(xhr.status >= 200 && xhr.status < 300 || xhr.status === 304)) throw new Error("Couldn't load " + url + ". Status: " + xhr.status);
            if (xhr.response !== undefined) {
              return new Uint8Array(xhr.response || []);
            } else {
              return intArrayFromString(xhr.responseText || '', true);
            }
          });
          var lazyArray = this;
          lazyArray.setDataGetter(function(chunkNum) {
            var start = chunkNum * chunkSize;
            var end = (chunkNum+1) * chunkSize - 1; // including this byte
            end = Math.min(end, datalength-1); // if datalength-1 is selected, this is the last block
            if (typeof(lazyArray.chunks[chunkNum]) === "undefined") {
              lazyArray.chunks[chunkNum] = doXHR(start, end);
            }
            if (typeof(lazyArray.chunks[chunkNum]) === "undefined") throw new Error("doXHR failed!");
            return lazyArray.chunks[chunkNum];
          });
  
          this._length = datalength;
          this._chunkSize = chunkSize;
          this.lengthKnown = true;
        }
        if (typeof XMLHttpRequest !== 'undefined') {
          if (!ENVIRONMENT_IS_WORKER) throw 'Cannot do synchronous binary XHRs outside webworkers in modern browsers. Use --embed-file or --preload-file in emcc';
          var lazyArray = new LazyUint8Array();
          Object.defineProperty(lazyArray, "length", {
              get: function() {
                  if(!this.lengthKnown) {
                      this.cacheLength();
                  }
                  return this._length;
              }
          });
          Object.defineProperty(lazyArray, "chunkSize", {
              get: function() {
                  if(!this.lengthKnown) {
                      this.cacheLength();
                  }
                  return this._chunkSize;
              }
          });
  
          var properties = { isDevice: false, contents: lazyArray };
        } else {
          var properties = { isDevice: false, url: url };
        }
  
        var node = FS.createFile(parent, name, properties, canRead, canWrite);
        // This is a total hack, but I want to get this lazy file code out of the
        // core of MEMFS. If we want to keep this lazy file concept I feel it should
        // be its own thin LAZYFS proxying calls to MEMFS.
        if (properties.contents) {
          node.contents = properties.contents;
        } else if (properties.url) {
          node.contents = null;
          node.url = properties.url;
        }
        // Add a function that defers querying the file size until it is asked the first time.
        Object.defineProperty(node, "usedBytes", {
            get: function() { return this.contents.length; }
        });
        // override each stream op with one that tries to force load the lazy file first
        var stream_ops = {};
        var keys = Object.keys(node.stream_ops);
        keys.forEach(function(key) {
          var fn = node.stream_ops[key];
          stream_ops[key] = function forceLoadLazyFile() {
            if (!FS.forceLoadFile(node)) {
              throw new FS.ErrnoError(ERRNO_CODES.EIO);
            }
            return fn.apply(null, arguments);
          };
        });
        // use a custom read function
        stream_ops.read = function stream_ops_read(stream, buffer, offset, length, position) {
          if (!FS.forceLoadFile(node)) {
            throw new FS.ErrnoError(ERRNO_CODES.EIO);
          }
          var contents = stream.node.contents;
          if (position >= contents.length)
            return 0;
          var size = Math.min(contents.length - position, length);
          assert(size >= 0);
          if (contents.slice) { // normal array
            for (var i = 0; i < size; i++) {
              buffer[offset + i] = contents[position + i];
            }
          } else {
            for (var i = 0; i < size; i++) { // LazyUint8Array from sync binary XHR
              buffer[offset + i] = contents.get(position + i);
            }
          }
          return size;
        };
        node.stream_ops = stream_ops;
        return node;
      },createPreloadedFile:function (parent, name, url, canRead, canWrite, onload, onerror, dontCreateFile, canOwn) {
        Browser.init();
        // TODO we should allow people to just pass in a complete filename instead
        // of parent and name being that we just join them anyways
        var fullname = name ? PATH.resolve(PATH.join2(parent, name)) : parent;
        function processData(byteArray) {
          function finish(byteArray) {
            if (!dontCreateFile) {
              FS.createDataFile(parent, name, byteArray, canRead, canWrite, canOwn);
            }
            if (onload) onload();
            removeRunDependency('cp ' + fullname);
          }
          var handled = false;
          Module['preloadPlugins'].forEach(function(plugin) {
            if (handled) return;
            if (plugin['canHandle'](fullname)) {
              plugin['handle'](byteArray, fullname, finish, function() {
                if (onerror) onerror();
                removeRunDependency('cp ' + fullname);
              });
              handled = true;
            }
          });
          if (!handled) finish(byteArray);
        }
        addRunDependency('cp ' + fullname);
        if (typeof url == 'string') {
          Browser.asyncLoad(url, function(byteArray) {
            processData(byteArray);
          }, onerror);
        } else {
          processData(url);
        }
      },indexedDB:function () {
        return window.indexedDB || window.mozIndexedDB || window.webkitIndexedDB || window.msIndexedDB;
      },DB_NAME:function () {
        return 'EM_FS_' + window.location.pathname;
      },DB_VERSION:20,DB_STORE_NAME:"FILE_DATA",saveFilesToDB:function (paths, onload, onerror) {
        onload = onload || function(){};
        onerror = onerror || function(){};
        var indexedDB = FS.indexedDB();
        try {
          var openRequest = indexedDB.open(FS.DB_NAME(), FS.DB_VERSION);
        } catch (e) {
          return onerror(e);
        }
        openRequest.onupgradeneeded = function openRequest_onupgradeneeded() {
          console.log('creating db');
          var db = openRequest.result;
          db.createObjectStore(FS.DB_STORE_NAME);
        };
        openRequest.onsuccess = function openRequest_onsuccess() {
          var db = openRequest.result;
          var transaction = db.transaction([FS.DB_STORE_NAME], 'readwrite');
          var files = transaction.objectStore(FS.DB_STORE_NAME);
          var ok = 0, fail = 0, total = paths.length;
          function finish() {
            if (fail == 0) onload(); else onerror();
          }
          paths.forEach(function(path) {
            var putRequest = files.put(FS.analyzePath(path).object.contents, path);
            putRequest.onsuccess = function putRequest_onsuccess() { ok++; if (ok + fail == total) finish() };
            putRequest.onerror = function putRequest_onerror() { fail++; if (ok + fail == total) finish() };
          });
          transaction.onerror = onerror;
        };
        openRequest.onerror = onerror;
      },loadFilesFromDB:function (paths, onload, onerror) {
        onload = onload || function(){};
        onerror = onerror || function(){};
        var indexedDB = FS.indexedDB();
        try {
          var openRequest = indexedDB.open(FS.DB_NAME(), FS.DB_VERSION);
        } catch (e) {
          return onerror(e);
        }
        openRequest.onupgradeneeded = onerror; // no database to load from
        openRequest.onsuccess = function openRequest_onsuccess() {
          var db = openRequest.result;
          try {
            var transaction = db.transaction([FS.DB_STORE_NAME], 'readonly');
          } catch(e) {
            onerror(e);
            return;
          }
          var files = transaction.objectStore(FS.DB_STORE_NAME);
          var ok = 0, fail = 0, total = paths.length;
          function finish() {
            if (fail == 0) onload(); else onerror();
          }
          paths.forEach(function(path) {
            var getRequest = files.get(path);
            getRequest.onsuccess = function getRequest_onsuccess() {
              if (FS.analyzePath(path).exists) {
                FS.unlink(path);
              }
              FS.createDataFile(PATH.dirname(path), PATH.basename(path), getRequest.result, true, true, true);
              ok++;
              if (ok + fail == total) finish();
            };
            getRequest.onerror = function getRequest_onerror() { fail++; if (ok + fail == total) finish() };
          });
          transaction.onerror = onerror;
        };
        openRequest.onerror = onerror;
      }};var PATH={splitPath:function (filename) {
        var splitPathRe = /^(\/?|)([\s\S]*?)((?:\.{1,2}|[^\/]+?|)(\.[^.\/]*|))(?:[\/]*)$/;
        return splitPathRe.exec(filename).slice(1);
      },normalizeArray:function (parts, allowAboveRoot) {
        // if the path tries to go above the root, `up` ends up > 0
        var up = 0;
        for (var i = parts.length - 1; i >= 0; i--) {
          var last = parts[i];
          if (last === '.') {
            parts.splice(i, 1);
          } else if (last === '..') {
            parts.splice(i, 1);
            up++;
          } else if (up) {
            parts.splice(i, 1);
            up--;
          }
        }
        // if the path is allowed to go above the root, restore leading ..s
        if (allowAboveRoot) {
          for (; up--; up) {
            parts.unshift('..');
          }
        }
        return parts;
      },normalize:function (path) {
        var isAbsolute = path.charAt(0) === '/',
            trailingSlash = path.substr(-1) === '/';
        // Normalize the path
        path = PATH.normalizeArray(path.split('/').filter(function(p) {
          return !!p;
        }), !isAbsolute).join('/');
        if (!path && !isAbsolute) {
          path = '.';
        }
        if (path && trailingSlash) {
          path += '/';
        }
        return (isAbsolute ? '/' : '') + path;
      },dirname:function (path) {
        var result = PATH.splitPath(path),
            root = result[0],
            dir = result[1];
        if (!root && !dir) {
          // No dirname whatsoever
          return '.';
        }
        if (dir) {
          // It has a dirname, strip trailing slash
          dir = dir.substr(0, dir.length - 1);
        }
        return root + dir;
      },basename:function (path) {
        // EMSCRIPTEN return '/'' for '/', not an empty string
        if (path === '/') return '/';
        var lastSlash = path.lastIndexOf('/');
        if (lastSlash === -1) return path;
        return path.substr(lastSlash+1);
      },extname:function (path) {
        return PATH.splitPath(path)[3];
      },join:function () {
        var paths = Array.prototype.slice.call(arguments, 0);
        return PATH.normalize(paths.join('/'));
      },join2:function (l, r) {
        return PATH.normalize(l + '/' + r);
      },resolve:function () {
        var resolvedPath = '',
          resolvedAbsolute = false;
        for (var i = arguments.length - 1; i >= -1 && !resolvedAbsolute; i--) {
          var path = (i >= 0) ? arguments[i] : FS.cwd();
          // Skip empty and invalid entries
          if (typeof path !== 'string') {
            throw new TypeError('Arguments to path.resolve must be strings');
          } else if (!path) {
            return ''; // an invalid portion invalidates the whole thing
          }
          resolvedPath = path + '/' + resolvedPath;
          resolvedAbsolute = path.charAt(0) === '/';
        }
        // At this point the path should be resolved to a full absolute path, but
        // handle relative paths to be safe (might happen when process.cwd() fails)
        resolvedPath = PATH.normalizeArray(resolvedPath.split('/').filter(function(p) {
          return !!p;
        }), !resolvedAbsolute).join('/');
        return ((resolvedAbsolute ? '/' : '') + resolvedPath) || '.';
      },relative:function (from, to) {
        from = PATH.resolve(from).substr(1);
        to = PATH.resolve(to).substr(1);
        function trim(arr) {
          var start = 0;
          for (; start < arr.length; start++) {
            if (arr[start] !== '') break;
          }
          var end = arr.length - 1;
          for (; end >= 0; end--) {
            if (arr[end] !== '') break;
          }
          if (start > end) return [];
          return arr.slice(start, end - start + 1);
        }
        var fromParts = trim(from.split('/'));
        var toParts = trim(to.split('/'));
        var length = Math.min(fromParts.length, toParts.length);
        var samePartsLength = length;
        for (var i = 0; i < length; i++) {
          if (fromParts[i] !== toParts[i]) {
            samePartsLength = i;
            break;
          }
        }
        var outputParts = [];
        for (var i = samePartsLength; i < fromParts.length; i++) {
          outputParts.push('..');
        }
        outputParts = outputParts.concat(toParts.slice(samePartsLength));
        return outputParts.join('/');
      }};var Browser={mainLoop:{scheduler:null,method:"",shouldPause:false,paused:false,queue:[],pause:function () {
          Browser.mainLoop.shouldPause = true;
        },resume:function () {
          if (Browser.mainLoop.paused) {
            Browser.mainLoop.paused = false;
            Browser.mainLoop.scheduler();
          }
          Browser.mainLoop.shouldPause = false;
        },updateStatus:function () {
          if (Module['setStatus']) {
            var message = Module['statusMessage'] || 'Please wait...';
            var remaining = Browser.mainLoop.remainingBlockers;
            var expected = Browser.mainLoop.expectedBlockers;
            if (remaining) {
              if (remaining < expected) {
                Module['setStatus'](message + ' (' + (expected - remaining) + '/' + expected + ')');
              } else {
                Module['setStatus'](message);
              }
            } else {
              Module['setStatus']('');
            }
          }
        },runIter:function (func) {
          if (ABORT) return;
          if (Module['preMainLoop']) {
            var preRet = Module['preMainLoop']();
            if (preRet === false) {
              return; // |return false| skips a frame
            }
          }
          try {
            func();
          } catch (e) {
            if (e instanceof ExitStatus) {
              return;
            } else {
              if (e && typeof e === 'object' && e.stack) Module.printErr('exception thrown: ' + [e, e.stack]);
              throw e;
            }
          }
          if (Module['postMainLoop']) Module['postMainLoop']();
        }},isFullScreen:false,pointerLock:false,moduleContextCreatedCallbacks:[],workers:[],init:function () {
        if (!Module["preloadPlugins"]) Module["preloadPlugins"] = []; // needs to exist even in workers
  
        if (Browser.initted) return;
        Browser.initted = true;
  
        try {
          new Blob();
          Browser.hasBlobConstructor = true;
        } catch(e) {
          Browser.hasBlobConstructor = false;
          console.log("warning: no blob constructor, cannot create blobs with mimetypes");
        }
        Browser.BlobBuilder = typeof MozBlobBuilder != "undefined" ? MozBlobBuilder : (typeof WebKitBlobBuilder != "undefined" ? WebKitBlobBuilder : (!Browser.hasBlobConstructor ? console.log("warning: no BlobBuilder") : null));
        Browser.URLObject = typeof window != "undefined" ? (window.URL ? window.URL : window.webkitURL) : undefined;
        if (!Module.noImageDecoding && typeof Browser.URLObject === 'undefined') {
          console.log("warning: Browser does not support creating object URLs. Built-in browser image decoding will not be available.");
          Module.noImageDecoding = true;
        }
  
        // Support for plugins that can process preloaded files. You can add more of these to
        // your app by creating and appending to Module.preloadPlugins.
        //
        // Each plugin is asked if it can handle a file based on the file's name. If it can,
        // it is given the file's raw data. When it is done, it calls a callback with the file's
        // (possibly modified) data. For example, a plugin might decompress a file, or it
        // might create some side data structure for use later (like an Image element, etc.).
  
        var imagePlugin = {};
        imagePlugin['canHandle'] = function imagePlugin_canHandle(name) {
          return !Module.noImageDecoding && /\.(jpg|jpeg|png|bmp)$/i.test(name);
        };
        imagePlugin['handle'] = function imagePlugin_handle(byteArray, name, onload, onerror) {
          var b = null;
          if (Browser.hasBlobConstructor) {
            try {
              b = new Blob([byteArray], { type: Browser.getMimetype(name) });
              if (b.size !== byteArray.length) { // Safari bug #118630
                // Safari's Blob can only take an ArrayBuffer
                b = new Blob([(new Uint8Array(byteArray)).buffer], { type: Browser.getMimetype(name) });
              }
            } catch(e) {
              Runtime.warnOnce('Blob constructor present but fails: ' + e + '; falling back to blob builder');
            }
          }
          if (!b) {
            var bb = new Browser.BlobBuilder();
            bb.append((new Uint8Array(byteArray)).buffer); // we need to pass a buffer, and must copy the array to get the right data range
            b = bb.getBlob();
          }
          var url = Browser.URLObject.createObjectURL(b);
          var img = new Image();
          img.onload = function img_onload() {
            assert(img.complete, 'Image ' + name + ' could not be decoded');
            var canvas = document.createElement('canvas');
            canvas.width = img.width;
            canvas.height = img.height;
            var ctx = canvas.getContext('2d');
            ctx.drawImage(img, 0, 0);
            Module["preloadedImages"][name] = canvas;
            Browser.URLObject.revokeObjectURL(url);
            if (onload) onload(byteArray);
          };
          img.onerror = function img_onerror(event) {
            console.log('Image ' + url + ' could not be decoded');
            if (onerror) onerror();
          };
          img.src = url;
        };
        Module['preloadPlugins'].push(imagePlugin);
  
        var audioPlugin = {};
        audioPlugin['canHandle'] = function audioPlugin_canHandle(name) {
          return !Module.noAudioDecoding && name.substr(-4) in { '.ogg': 1, '.wav': 1, '.mp3': 1 };
        };
        audioPlugin['handle'] = function audioPlugin_handle(byteArray, name, onload, onerror) {
          var done = false;
          function finish(audio) {
            if (done) return;
            done = true;
            Module["preloadedAudios"][name] = audio;
            if (onload) onload(byteArray);
          }
          function fail() {
            if (done) return;
            done = true;
            Module["preloadedAudios"][name] = new Audio(); // empty shim
            if (onerror) onerror();
          }
          if (Browser.hasBlobConstructor) {
            try {
              var b = new Blob([byteArray], { type: Browser.getMimetype(name) });
            } catch(e) {
              return fail();
            }
            var url = Browser.URLObject.createObjectURL(b); // XXX we never revoke this!
            var audio = new Audio();
            audio.addEventListener('canplaythrough', function() { finish(audio) }, false); // use addEventListener due to chromium bug 124926
            audio.onerror = function audio_onerror(event) {
              if (done) return;
              console.log('warning: browser could not fully decode audio ' + name + ', trying slower base64 approach');
              function encode64(data) {
                var BASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
                var PAD = '=';
                var ret = '';
                var leftchar = 0;
                var leftbits = 0;
                for (var i = 0; i < data.length; i++) {
                  leftchar = (leftchar << 8) | data[i];
                  leftbits += 8;
                  while (leftbits >= 6) {
                    var curr = (leftchar >> (leftbits-6)) & 0x3f;
                    leftbits -= 6;
                    ret += BASE[curr];
                  }
                }
                if (leftbits == 2) {
                  ret += BASE[(leftchar&3) << 4];
                  ret += PAD + PAD;
                } else if (leftbits == 4) {
                  ret += BASE[(leftchar&0xf) << 2];
                  ret += PAD;
                }
                return ret;
              }
              audio.src = 'data:audio/x-' + name.substr(-3) + ';base64,' + encode64(byteArray);
              finish(audio); // we don't wait for confirmation this worked - but it's worth trying
            };
            audio.src = url;
            // workaround for chrome bug 124926 - we do not always get oncanplaythrough or onerror
            Browser.safeSetTimeout(function() {
              finish(audio); // try to use it even though it is not necessarily ready to play
            }, 10000);
          } else {
            return fail();
          }
        };
        Module['preloadPlugins'].push(audioPlugin);
  
        // Canvas event setup
  
        var canvas = Module['canvas'];
        function pointerLockChange() {
          Browser.pointerLock = document['pointerLockElement'] === canvas ||
                                document['mozPointerLockElement'] === canvas ||
                                document['webkitPointerLockElement'] === canvas ||
                                document['msPointerLockElement'] === canvas;
        }
        if (canvas) {
          // forced aspect ratio can be enabled by defining 'forcedAspectRatio' on Module
          // Module['forcedAspectRatio'] = 4 / 3;
          
          canvas.requestPointerLock = canvas['requestPointerLock'] ||
                                      canvas['mozRequestPointerLock'] ||
                                      canvas['webkitRequestPointerLock'] ||
                                      canvas['msRequestPointerLock'] ||
                                      function(){};
          canvas.exitPointerLock = document['exitPointerLock'] ||
                                   document['mozExitPointerLock'] ||
                                   document['webkitExitPointerLock'] ||
                                   document['msExitPointerLock'] ||
                                   function(){}; // no-op if function does not exist
          canvas.exitPointerLock = canvas.exitPointerLock.bind(document);
  
  
          document.addEventListener('pointerlockchange', pointerLockChange, false);
          document.addEventListener('mozpointerlockchange', pointerLockChange, false);
          document.addEventListener('webkitpointerlockchange', pointerLockChange, false);
          document.addEventListener('mspointerlockchange', pointerLockChange, false);
  
          if (Module['elementPointerLock']) {
            canvas.addEventListener("click", function(ev) {
              if (!Browser.pointerLock && canvas.requestPointerLock) {
                canvas.requestPointerLock();
                ev.preventDefault();
              }
            }, false);
          }
        }
      },createContext:function (canvas, useWebGL, setInModule, webGLContextAttributes) {
        if (useWebGL && Module.ctx && canvas == Module.canvas) return Module.ctx; // no need to recreate GL context if it's already been created for this canvas.
  
        var ctx;
        var contextHandle;
        if (useWebGL) {
          // For GLES2/desktop GL compatibility, adjust a few defaults to be different to WebGL defaults, so that they align better with the desktop defaults.
          var contextAttributes = {
            antialias: false,
            alpha: false
          };
  
          if (webGLContextAttributes) {
            for (var attribute in webGLContextAttributes) {
              contextAttributes[attribute] = webGLContextAttributes[attribute];
            }
          }
  
          contextHandle = GL.createContext(canvas, contextAttributes);
          if (contextHandle) {
            ctx = GL.getContext(contextHandle).GLctx;
          }
          // Set the background of the WebGL canvas to black
          canvas.style.backgroundColor = "black";
        } else {
          ctx = canvas.getContext('2d');
        }
  
        if (!ctx) return null;
  
        if (setInModule) {
          if (!useWebGL) assert(typeof GLctx === 'undefined', 'cannot set in module if GLctx is used, but we are a non-GL context that would replace it');
  
          Module.ctx = ctx;
          if (useWebGL) GL.makeContextCurrent(contextHandle);
          Module.useWebGL = useWebGL;
          Browser.moduleContextCreatedCallbacks.forEach(function(callback) { callback() });
          Browser.init();
        }
        return ctx;
      },destroyContext:function (canvas, useWebGL, setInModule) {},fullScreenHandlersInstalled:false,lockPointer:undefined,resizeCanvas:undefined,requestFullScreen:function (lockPointer, resizeCanvas) {
        Browser.lockPointer = lockPointer;
        Browser.resizeCanvas = resizeCanvas;
        if (typeof Browser.lockPointer === 'undefined') Browser.lockPointer = true;
        if (typeof Browser.resizeCanvas === 'undefined') Browser.resizeCanvas = false;
  
        var canvas = Module['canvas'];
        function fullScreenChange() {
          Browser.isFullScreen = false;
          var canvasContainer = canvas.parentNode;
          if ((document['webkitFullScreenElement'] || document['webkitFullscreenElement'] ||
               document['mozFullScreenElement'] || document['mozFullscreenElement'] ||
               document['fullScreenElement'] || document['fullscreenElement'] ||
               document['msFullScreenElement'] || document['msFullscreenElement'] ||
               document['webkitCurrentFullScreenElement']) === canvasContainer) {
            canvas.cancelFullScreen = document['cancelFullScreen'] ||
                                      document['mozCancelFullScreen'] ||
                                      document['webkitCancelFullScreen'] ||
                                      document['msExitFullscreen'] ||
                                      document['exitFullscreen'] ||
                                      function() {};
            canvas.cancelFullScreen = canvas.cancelFullScreen.bind(document);
            if (Browser.lockPointer) canvas.requestPointerLock();
            Browser.isFullScreen = true;
            if (Browser.resizeCanvas) Browser.setFullScreenCanvasSize();
          } else {
            
            // remove the full screen specific parent of the canvas again to restore the HTML structure from before going full screen
            canvasContainer.parentNode.insertBefore(canvas, canvasContainer);
            canvasContainer.parentNode.removeChild(canvasContainer);
            
            if (Browser.resizeCanvas) Browser.setWindowedCanvasSize();
          }
          if (Module['onFullScreen']) Module['onFullScreen'](Browser.isFullScreen);
          Browser.updateCanvasDimensions(canvas);
        }
  
        if (!Browser.fullScreenHandlersInstalled) {
          Browser.fullScreenHandlersInstalled = true;
          document.addEventListener('fullscreenchange', fullScreenChange, false);
          document.addEventListener('mozfullscreenchange', fullScreenChange, false);
          document.addEventListener('webkitfullscreenchange', fullScreenChange, false);
          document.addEventListener('MSFullscreenChange', fullScreenChange, false);
        }
  
        // create a new parent to ensure the canvas has no siblings. this allows browsers to optimize full screen performance when its parent is the full screen root
        var canvasContainer = document.createElement("div");
        canvas.parentNode.insertBefore(canvasContainer, canvas);
        canvasContainer.appendChild(canvas);
        
        // use parent of canvas as full screen root to allow aspect ratio correction (Firefox stretches the root to screen size)
        canvasContainer.requestFullScreen = canvasContainer['requestFullScreen'] ||
                                            canvasContainer['mozRequestFullScreen'] ||
                                            canvasContainer['msRequestFullscreen'] ||
                                           (canvasContainer['webkitRequestFullScreen'] ? function() { canvasContainer['webkitRequestFullScreen'](Element['ALLOW_KEYBOARD_INPUT']) } : null);
        canvasContainer.requestFullScreen();
      },nextRAF:0,fakeRequestAnimationFrame:function (func) {
        // try to keep 60fps between calls to here
        var now = Date.now();
        if (Browser.nextRAF === 0) {
          Browser.nextRAF = now + 1000/60;
        } else {
          while (now + 2 >= Browser.nextRAF) { // fudge a little, to avoid timer jitter causing us to do lots of delay:0
            Browser.nextRAF += 1000/60;
          }
        }
        var delay = Math.max(Browser.nextRAF - now, 0);
        setTimeout(func, delay);
      },requestAnimationFrame:function requestAnimationFrame(func) {
        if (typeof window === 'undefined') { // Provide fallback to setTimeout if window is undefined (e.g. in Node.js)
          Browser.fakeRequestAnimationFrame(func);
        } else {
          if (!window.requestAnimationFrame) {
            window.requestAnimationFrame = window['requestAnimationFrame'] ||
                                           window['mozRequestAnimationFrame'] ||
                                           window['webkitRequestAnimationFrame'] ||
                                           window['msRequestAnimationFrame'] ||
                                           window['oRequestAnimationFrame'] ||
                                           Browser.fakeRequestAnimationFrame;
          }
          window.requestAnimationFrame(func);
        }
      },safeCallback:function (func) {
        return function() {
          if (!ABORT) return func.apply(null, arguments);
        };
      },safeRequestAnimationFrame:function (func) {
        return Browser.requestAnimationFrame(function() {
          if (!ABORT) func();
        });
      },safeSetTimeout:function (func, timeout) {
        Module['noExitRuntime'] = true;
        return setTimeout(function() {
          if (!ABORT) func();
        }, timeout);
      },safeSetInterval:function (func, timeout) {
        Module['noExitRuntime'] = true;
        return setInterval(function() {
          if (!ABORT) func();
        }, timeout);
      },getMimetype:function (name) {
        return {
          'jpg': 'image/jpeg',
          'jpeg': 'image/jpeg',
          'png': 'image/png',
          'bmp': 'image/bmp',
          'ogg': 'audio/ogg',
          'wav': 'audio/wav',
          'mp3': 'audio/mpeg'
        }[name.substr(name.lastIndexOf('.')+1)];
      },getUserMedia:function (func) {
        if(!window.getUserMedia) {
          window.getUserMedia = navigator['getUserMedia'] ||
                                navigator['mozGetUserMedia'];
        }
        window.getUserMedia(func);
      },getMovementX:function (event) {
        return event['movementX'] ||
               event['mozMovementX'] ||
               event['webkitMovementX'] ||
               0;
      },getMovementY:function (event) {
        return event['movementY'] ||
               event['mozMovementY'] ||
               event['webkitMovementY'] ||
               0;
      },getMouseWheelDelta:function (event) {
        var delta = 0;
        switch (event.type) {
          case 'DOMMouseScroll': 
            delta = event.detail;
            break;
          case 'mousewheel': 
            delta = event.wheelDelta;
            break;
          case 'wheel': 
            delta = event['deltaY'];
            break;
          default:
            throw 'unrecognized mouse wheel event: ' + event.type;
        }
        return delta;
      },mouseX:0,mouseY:0,mouseMovementX:0,mouseMovementY:0,touches:{},lastTouches:{},calculateMouseEvent:function (event) { // event should be mousemove, mousedown or mouseup
        if (Browser.pointerLock) {
          // When the pointer is locked, calculate the coordinates
          // based on the movement of the mouse.
          // Workaround for Firefox bug 764498
          if (event.type != 'mousemove' &&
              ('mozMovementX' in event)) {
            Browser.mouseMovementX = Browser.mouseMovementY = 0;
          } else {
            Browser.mouseMovementX = Browser.getMovementX(event);
            Browser.mouseMovementY = Browser.getMovementY(event);
          }
          
          // check if SDL is available
          if (typeof SDL != "undefined") {
          	Browser.mouseX = SDL.mouseX + Browser.mouseMovementX;
          	Browser.mouseY = SDL.mouseY + Browser.mouseMovementY;
          } else {
          	// just add the mouse delta to the current absolut mouse position
          	// FIXME: ideally this should be clamped against the canvas size and zero
          	Browser.mouseX += Browser.mouseMovementX;
          	Browser.mouseY += Browser.mouseMovementY;
          }        
        } else {
          // Otherwise, calculate the movement based on the changes
          // in the coordinates.
          var rect = Module["canvas"].getBoundingClientRect();
          var cw = Module["canvas"].width;
          var ch = Module["canvas"].height;
  
          // Neither .scrollX or .pageXOffset are defined in a spec, but
          // we prefer .scrollX because it is currently in a spec draft.
          // (see: http://www.w3.org/TR/2013/WD-cssom-view-20131217/)
          var scrollX = ((typeof window.scrollX !== 'undefined') ? window.scrollX : window.pageXOffset);
          var scrollY = ((typeof window.scrollY !== 'undefined') ? window.scrollY : window.pageYOffset);
  
          if (event.type === 'touchstart' || event.type === 'touchend' || event.type === 'touchmove') {
            var touch = event.touch;
            if (touch === undefined) {
              return; // the "touch" property is only defined in SDL
  
            }
            var adjustedX = touch.pageX - (scrollX + rect.left);
            var adjustedY = touch.pageY - (scrollY + rect.top);
  
            adjustedX = adjustedX * (cw / rect.width);
            adjustedY = adjustedY * (ch / rect.height);
  
            var coords = { x: adjustedX, y: adjustedY };
            
            if (event.type === 'touchstart') {
              Browser.lastTouches[touch.identifier] = coords;
              Browser.touches[touch.identifier] = coords;
            } else if (event.type === 'touchend' || event.type === 'touchmove') {
              Browser.lastTouches[touch.identifier] = Browser.touches[touch.identifier];
              Browser.touches[touch.identifier] = { x: adjustedX, y: adjustedY };
            } 
            return;
          }
  
          var x = event.pageX - (scrollX + rect.left);
          var y = event.pageY - (scrollY + rect.top);
  
          // the canvas might be CSS-scaled compared to its backbuffer;
          // SDL-using content will want mouse coordinates in terms
          // of backbuffer units.
          x = x * (cw / rect.width);
          y = y * (ch / rect.height);
  
          Browser.mouseMovementX = x - Browser.mouseX;
          Browser.mouseMovementY = y - Browser.mouseY;
          Browser.mouseX = x;
          Browser.mouseY = y;
        }
      },xhrLoad:function (url, onload, onerror) {
        var xhr = new XMLHttpRequest();
        xhr.open('GET', url, true);
        xhr.responseType = 'arraybuffer';
        xhr.onload = function xhr_onload() {
          if (xhr.status == 200 || (xhr.status == 0 && xhr.response)) { // file URLs can return 0
            onload(xhr.response);
          } else {
            onerror();
          }
        };
        xhr.onerror = onerror;
        xhr.send(null);
      },asyncLoad:function (url, onload, onerror, noRunDep) {
        Browser.xhrLoad(url, function(arrayBuffer) {
          assert(arrayBuffer, 'Loading data file "' + url + '" failed (no arrayBuffer).');
          onload(new Uint8Array(arrayBuffer));
          if (!noRunDep) removeRunDependency('al ' + url);
        }, function(event) {
          if (onerror) {
            onerror();
          } else {
            throw 'Loading data file "' + url + '" failed.';
          }
        });
        if (!noRunDep) addRunDependency('al ' + url);
      },resizeListeners:[],updateResizeListeners:function () {
        var canvas = Module['canvas'];
        Browser.resizeListeners.forEach(function(listener) {
          listener(canvas.width, canvas.height);
        });
      },setCanvasSize:function (width, height, noUpdates) {
        var canvas = Module['canvas'];
        Browser.updateCanvasDimensions(canvas, width, height);
        if (!noUpdates) Browser.updateResizeListeners();
      },windowedWidth:0,windowedHeight:0,setFullScreenCanvasSize:function () {
        // check if SDL is available   
        if (typeof SDL != "undefined") {
        	var flags = HEAPU32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)];
        	flags = flags | 0x00800000; // set SDL_FULLSCREEN flag
        	HEAP32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)]=flags
        }
        Browser.updateResizeListeners();
      },setWindowedCanvasSize:function () {
        // check if SDL is available       
        if (typeof SDL != "undefined") {
        	var flags = HEAPU32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)];
        	flags = flags & ~0x00800000; // clear SDL_FULLSCREEN flag
        	HEAP32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)]=flags
        }
        Browser.updateResizeListeners();
      },updateCanvasDimensions:function (canvas, wNative, hNative) {
        if (wNative && hNative) {
          canvas.widthNative = wNative;
          canvas.heightNative = hNative;
        } else {
          wNative = canvas.widthNative;
          hNative = canvas.heightNative;
        }
        var w = wNative;
        var h = hNative;
        if (Module['forcedAspectRatio'] && Module['forcedAspectRatio'] > 0) {
          if (w/h < Module['forcedAspectRatio']) {
            w = Math.round(h * Module['forcedAspectRatio']);
          } else {
            h = Math.round(w / Module['forcedAspectRatio']);
          }
        }
        if (((document['webkitFullScreenElement'] || document['webkitFullscreenElement'] ||
             document['mozFullScreenElement'] || document['mozFullscreenElement'] ||
             document['fullScreenElement'] || document['fullscreenElement'] ||
             document['msFullScreenElement'] || document['msFullscreenElement'] ||
             document['webkitCurrentFullScreenElement']) === canvas.parentNode) && (typeof screen != 'undefined')) {
           var factor = Math.min(screen.width / w, screen.height / h);
           w = Math.round(w * factor);
           h = Math.round(h * factor);
        }
        if (Browser.resizeCanvas) {
          if (canvas.width  != w) canvas.width  = w;
          if (canvas.height != h) canvas.height = h;
          if (typeof canvas.style != 'undefined') {
            canvas.style.removeProperty( "width");
            canvas.style.removeProperty("height");
          }
        } else {
          if (canvas.width  != wNative) canvas.width  = wNative;
          if (canvas.height != hNative) canvas.height = hNative;
          if (typeof canvas.style != 'undefined') {
            if (w != wNative || h != hNative) {
              canvas.style.setProperty( "width", w + "px", "important");
              canvas.style.setProperty("height", h + "px", "important");
            } else {
              canvas.style.removeProperty( "width");
              canvas.style.removeProperty("height");
            }
          }
        }
      },wgetRequests:{},nextWgetRequestHandle:0,getNextWgetRequestHandle:function () {
        var handle = Browser.nextWgetRequestHandle;
        Browser.nextWgetRequestHandle++;
        return handle;
      }};

  function _time(ptr) {
      var ret = (Date.now()/1000)|0;
      if (ptr) {
        HEAP32[((ptr)>>2)]=ret;
      }
      return ret;
    }

   
  Module["_strlen"] = _strlen;

  
  function _emscripten_memcpy_big(dest, src, num) {
      HEAPU8.set(HEAPU8.subarray(src, src+num), dest);
      return dest;
    } 
  Module["_memcpy"] = _memcpy;
___errno_state = Runtime.staticAlloc(4); HEAP32[((___errno_state)>>2)]=0;
Module["requestFullScreen"] = function Module_requestFullScreen(lockPointer, resizeCanvas) { Browser.requestFullScreen(lockPointer, resizeCanvas) };
  Module["requestAnimationFrame"] = function Module_requestAnimationFrame(func) { Browser.requestAnimationFrame(func) };
  Module["setCanvasSize"] = function Module_setCanvasSize(width, height, noUpdates) { Browser.setCanvasSize(width, height, noUpdates) };
  Module["pauseMainLoop"] = function Module_pauseMainLoop() { Browser.mainLoop.pause() };
  Module["resumeMainLoop"] = function Module_resumeMainLoop() { Browser.mainLoop.resume() };
  Module["getUserMedia"] = function Module_getUserMedia() { Browser.getUserMedia() }
FS.staticInit();__ATINIT__.unshift({ func: function() { if (!Module["noFSInit"] && !FS.init.initialized) FS.init() } });__ATMAIN__.push({ func: function() { FS.ignorePermissions = false } });__ATEXIT__.push({ func: function() { FS.quit() } });Module["FS_createFolder"] = FS.createFolder;Module["FS_createPath"] = FS.createPath;Module["FS_createDataFile"] = FS.createDataFile;Module["FS_createPreloadedFile"] = FS.createPreloadedFile;Module["FS_createLazyFile"] = FS.createLazyFile;Module["FS_createLink"] = FS.createLink;Module["FS_createDevice"] = FS.createDevice;
__ATINIT__.unshift({ func: function() { TTY.init() } });__ATEXIT__.push({ func: function() { TTY.shutdown() } });TTY.utf8 = new Runtime.UTF8Processor();
if (ENVIRONMENT_IS_NODE) { var fs = require("fs"); NODEFS.staticInit(); }
STACK_BASE = STACKTOP = Runtime.alignMemory(STATICTOP);

staticSealed = true; // seal the static portion of memory

STACK_MAX = STACK_BASE + TOTAL_STACK;

DYNAMIC_BASE = DYNAMICTOP = Runtime.alignMemory(STACK_MAX);

assert(DYNAMIC_BASE < TOTAL_MEMORY, "TOTAL_MEMORY not big enough for stack");


  var Math_min = Math.min;
  // EMSCRIPTEN_START_ASM
  var asm = (function(global, env, buffer) {
    'almost asm';
    
    var HEAP8 = new global.Int8Array(buffer);
    var HEAP16 = new global.Int16Array(buffer);
    var HEAP32 = new global.Int32Array(buffer);
    var HEAPU8 = new global.Uint8Array(buffer);
    var HEAPU16 = new global.Uint16Array(buffer);
    var HEAPU32 = new global.Uint32Array(buffer);
    var HEAPF32 = new global.Float32Array(buffer);
    var HEAPF64 = new global.Float64Array(buffer);

  
  var STACKTOP=env.STACKTOP|0;
  var STACK_MAX=env.STACK_MAX|0;
  var tempDoublePtr=env.tempDoublePtr|0;
  var ABORT=env.ABORT|0;

    var __THREW__ = 0;
    var threwValue = 0;
    var setjmpId = 0;
    var undef = 0;
    var nan = +env.NaN, inf = +env.Infinity;
    var tempInt = 0, tempBigInt = 0, tempBigIntP = 0, tempBigIntS = 0, tempBigIntR = 0.0, tempBigIntI = 0, tempBigIntD = 0, tempValue = 0, tempDouble = 0.0;
  
    var tempRet0 = 0;
    var tempRet1 = 0;
    var tempRet2 = 0;
    var tempRet3 = 0;
    var tempRet4 = 0;
    var tempRet5 = 0;
    var tempRet6 = 0;
    var tempRet7 = 0;
    var tempRet8 = 0;
    var tempRet9 = 0;
  var Math_floor=global.Math.floor;
  var Math_abs=global.Math.abs;
  var Math_sqrt=global.Math.sqrt;
  var Math_pow=global.Math.pow;
  var Math_cos=global.Math.cos;
  var Math_sin=global.Math.sin;
  var Math_tan=global.Math.tan;
  var Math_acos=global.Math.acos;
  var Math_asin=global.Math.asin;
  var Math_atan=global.Math.atan;
  var Math_atan2=global.Math.atan2;
  var Math_exp=global.Math.exp;
  var Math_log=global.Math.log;
  var Math_ceil=global.Math.ceil;
  var Math_imul=global.Math.imul;
  var abort=env.abort;
  var assert=env.assert;
  var Math_min=env.min;
  var _fflush=env._fflush;
  var _abort=env._abort;
  var ___setErrNo=env.___setErrNo;
  var _sbrk=env._sbrk;
  var _time=env._time;
  var _emscripten_memcpy_big=env._emscripten_memcpy_big;
  var _sysconf=env._sysconf;
  var ___errno_location=env.___errno_location;
  var SIMD_float32x4=global.SIMD.float32x4;
  var SIMD_int32x4=global.SIMD.int32x4;
  var SIMD_int32x4_add=SIMD_int32x4.add;
  var SIMD_int32x4_sub=SIMD_int32x4.sub;
  var SIMD_int32x4_equal=SIMD_int32x4.equal;
  var SIMD_int32x4_notEqual=SIMD_int32x4.notEqual;
  var SIMD_int32x4_lessThan=SIMD_int32x4.lessThan;
  var SIMD_int32x4_lessThanOrEqual=SIMD_int32x4.lessThanOrEqual;
  var SIMD_int32x4_greaterThan=SIMD_int32x4.greaterThan;
  var SIMD_int32x4_greaterThanOrEqual=SIMD_int32x4.greaterThanOrEqual;
  var SIMD_int32x4_select=SIMD_int32x4.select;
  var SIMD_int32x4_and=SIMD_int32x4.and;
  var SIMD_int32x4_or=SIMD_int32x4.or;
  var SIMD_int32x4_xor=SIMD_int32x4.xor;
  var SIMD_int32x4_not=SIMD_int32x4.not;
  var SIMD_int32x4_splat=SIMD_int32x4.splat;
  var SIMD_int32x4_shuffle=SIMD_int32x4.shuffle;
  var SIMD_int32x4_shuffleMix=SIMD_int32x4.shuffleMix;
  var SIMD_int32x4_withX=SIMD_int32x4.withX;
  var SIMD_int32x4_withY=SIMD_int32x4.withY;
  var SIMD_int32x4_withZ=SIMD_int32x4.withZ;
  var SIMD_int32x4_withW=SIMD_int32x4.withW;
  var SIMD_int32x4_load=SIMD_int32x4.load;
  var SIMD_int32x4_loadX=SIMD_int32x4.loadX;
  var SIMD_int32x4_loadXY=SIMD_int32x4.loadXY;
  var SIMD_int32x4_loadXYZ=SIMD_int32x4.loadXYZ;
  var SIMD_int32x4_store=SIMD_int32x4.store;
  var SIMD_int32x4_storeX=SIMD_int32x4.storeX;
  var SIMD_int32x4_storeXY=SIMD_int32x4.storeXY;
  var SIMD_int32x4_storeXYZ=SIMD_int32x4.storeXYZ;
  var SIMD_int32x4_fromFloat32x4=SIMD_int32x4.fromFloat32x4;
  var SIMD_int32x4_fromFloat32x4Bits=SIMD_int32x4.fromFloat32x4Bits;
  var SIMD_float32x4_add=SIMD_float32x4.add;
  var SIMD_float32x4_sub=SIMD_float32x4.sub;
  var SIMD_float32x4_equal=SIMD_float32x4.equal;
  var SIMD_float32x4_notEqual=SIMD_float32x4.notEqual;
  var SIMD_float32x4_lessThan=SIMD_float32x4.lessThan;
  var SIMD_float32x4_lessThanOrEqual=SIMD_float32x4.lessThanOrEqual;
  var SIMD_float32x4_greaterThan=SIMD_float32x4.greaterThan;
  var SIMD_float32x4_greaterThanOrEqual=SIMD_float32x4.greaterThanOrEqual;
  var SIMD_float32x4_select=SIMD_float32x4.select;
  var SIMD_float32x4_and=SIMD_float32x4.and;
  var SIMD_float32x4_or=SIMD_float32x4.or;
  var SIMD_float32x4_xor=SIMD_float32x4.xor;
  var SIMD_float32x4_not=SIMD_float32x4.not;
  var SIMD_float32x4_splat=SIMD_float32x4.splat;
  var SIMD_float32x4_shuffle=SIMD_float32x4.shuffle;
  var SIMD_float32x4_shuffleMix=SIMD_float32x4.shuffleMix;
  var SIMD_float32x4_withX=SIMD_float32x4.withX;
  var SIMD_float32x4_withY=SIMD_float32x4.withY;
  var SIMD_float32x4_withZ=SIMD_float32x4.withZ;
  var SIMD_float32x4_withW=SIMD_float32x4.withW;
  var SIMD_float32x4_load=SIMD_float32x4.load;
  var SIMD_float32x4_loadX=SIMD_float32x4.loadX;
  var SIMD_float32x4_loadXY=SIMD_float32x4.loadXY;
  var SIMD_float32x4_loadXYZ=SIMD_float32x4.loadXYZ;
  var SIMD_float32x4_store=SIMD_float32x4.store;
  var SIMD_float32x4_storeX=SIMD_float32x4.storeX;
  var SIMD_float32x4_storeXY=SIMD_float32x4.storeXY;
  var SIMD_float32x4_storeXYZ=SIMD_float32x4.storeXYZ;
  var SIMD_float32x4_mul=SIMD_float32x4.mul;
  var SIMD_float32x4_div=SIMD_float32x4.div;
  var SIMD_float32x4_min=SIMD_float32x4.min;
  var SIMD_float32x4_max=SIMD_float32x4.max;
  var SIMD_float32x4_sqrt=SIMD_float32x4.sqrt;
  var SIMD_float32x4_fromInt32x4=SIMD_float32x4.fromInt32x4;
  var SIMD_float32x4_fromInt32x44Bits=SIMD_float32x4.fromInt32x44Bits;
  var tempFloat = 0.0;

  // EMSCRIPTEN_START_FUNCS
function _malloc($bytes) {
 $bytes = $bytes | 0;
 var $$pre$phi$i$iZ2D = 0, $$pre$phi$i145Z2D = 0, $$pre$phi$i67$iZ2D = 0, $$pre$phi$iZ2D = 0, $$pre$phiZ2D = 0, $0 = 0, $1 = 0, $10 = 0, $101 = 0, $102 = 0, $103 = 0, $105 = 0, $106 = 0, $108 = 0, $109 = 0, $110 = 0, $111 = 0, $112 = 0, $116 = 0, $120 = 0, $121 = 0, $125 = 0, $128 = 0, $129 = 0, $13 = 0, $130 = 0, $133 = 0, $138 = 0, $14 = 0, $141 = 0, $143 = 0, $149 = 0, $15 = 0, $150 = 0, $151 = 0, $157 = 0, $158 = 0, $159 = 0, $16 = 0, $163 = 0, $164 = 0, $165 = 0, $166 = 0, $168 = 0, $17 = 0, $174 = 0, $176 = 0, $179 = 0, $180 = 0, $181 = 0, $183 = 0, $184 = 0, $186 = 0, $189 = 0, $19 = 0, $190 = 0, $191 = 0, $192 = 0, $194 = 0, $196 = 0, $2 = 0, $20 = 0, $203 = 0, $204 = 0, $205 = 0, $208 = 0, $209 = 0, $211 = 0, $214 = 0, $215 = 0, $216 = 0, $217 = 0, $22 = 0, $23 = 0, $25 = 0, $26 = 0, $27 = 0, $28 = 0, $3 = 0, $31 = 0, $32 = 0, $33 = 0, $34 = 0, $35 = 0, $41 = 0, $43 = 0, $46 = 0, $47 = 0, $48 = 0, $49 = 0, $50 = 0, $52 = 0, $53 = 0, $55 = 0, $59 = 0, $62 = 0, $63 = 0, $64 = 0, $65 = 0, $68 = 0, $69 = 0, $70 = 0, $71 = 0, $72 = 0, $78 = 0, $8 = 0, $80 = 0, $83 = 0, $84 = 0, $85 = 0, $87 = 0, $88 = 0, $9 = 0, $90 = 0, $93 = 0, $94 = 0, $95 = 0, $96 = 0, $98 = 0, $99 = 0, $F$0$i$i = 0, $F104$0 = 0, $F197$0$i = 0, $F224$0$i$i = 0, $F289$0$i = 0, $I252$0$i$i = 0, $I315$0$i = 0, $I57$0$i$i = 0, $K105$017$i$i = 0, $K305$043$i$i = 0, $K372$024$i = 0, $R$0$i = 0, $R$0$i$i = 0, $R$0$i135 = 0, $R$1$i = 0, $R$1$i$i = 0, $R$1$i137 = 0, $RP$0$i = 0, $RP$0$i$i = 0, $RP$0$i134 = 0, $T$0$lcssa$i = 0, $T$0$lcssa$i$i = 0, $T$0$lcssa$i69$i = 0, $T$016$i$i = 0, $T$023$i = 0, $T$042$i$i = 0, $add$i$i = 0, $add$i147 = 0, $add$ptr$i = 0, $add$ptr$i$i$i = 0, $add$ptr$i126 = 0, $add$ptr14$i$i = 0, $add$ptr16$i$i = 0, $add$ptr16$sum$i$i = 0, $add$ptr16$sum2627$i$i = 0, $add$ptr16$sum56$i$i = 0, $add$ptr17$i$i = 0, $add$ptr224$i = 0, $add$ptr2418$i$i = 0, $add$ptr2420$i$i = 0, $add$ptr4$sum$i49$i = 0, $add$ptr7$i$i = 0, $add$ptr95 = 0, $add143 = 0, $add147$i = 0, $add17$i = 0, $add17$i150 = 0, $add177$i = 0, $add212$i = 0, $add26$i$i = 0, $add267$i = 0, $add278$i$i = 0, $add345$i = 0, $add51$i = 0, $add64 = 0, $add8 = 0, $add83$i$i = 0, $add9$i = 0, $and$i110 = 0, $and101$i = 0, $and11$i = 0, $and12$i = 0, $and13$i = 0, $and144 = 0, $and17$i = 0, $and264$i$i = 0, $and268$i$i = 0, $and273$i$i = 0, $and3$i = 0, $and32$i = 0, $and330$i = 0, $and335$i = 0, $and340$i = 0, $and37$i$i = 0, $and41 = 0, $and46 = 0, $and49 = 0, $and53 = 0, $and57 = 0, $and6$i = 0, $and61 = 0, $and63$i = 0, $and69$i$i = 0, $and72$i = 0, $and73$i$i = 0, $and76$i = 0, $and77$i = 0, $and78$i$i = 0, $and8$i = 0, $and80$i = 0, $and84$i = 0, $and88$i = 0, $and9$i = 0, $arrayidx = 0, $arrayidx$i$i = 0, $arrayidx$i21$i = 0, $arrayidx$i57$i = 0, $arrayidx103 = 0, $arrayidx103$i$i = 0, $arrayidx107$i$i = 0, $arrayidx113$i = 0, $arrayidx123$i$i = 0, $arrayidx126$i$i = 0, $arrayidx143$i$i = 0, $arrayidx150$i = 0, $arrayidx154$i131 = 0, $arrayidx160$i = 0, $arrayidx164$i = 0, $arrayidx183$i = 0, $arrayidx196$i = 0, $arrayidx203$i = 0, $arrayidx223$i$i = 0, $arrayidx287$i$i = 0, $arrayidx288$i = 0, $arrayidx325$i$i = 0, $arrayidx354$i = 0, $arrayidx393$i = 0, $arrayidx61$i = 0, $arrayidx65$i = 0, $arrayidx66 = 0, $arrayidx71$i = 0, $arrayidx75$i = 0, $arrayidx91$i$i = 0, $arrayidx94$i = 0, $arrayidx96$i$i = 0, $bk = 0, $bk135$i = 0, $bk47$i = 0, $bk78 = 0, $bk82$i$i = 0, $br$0$i = 0, $call$i$i = 0, $call128$i = 0, $call129$i = 0, $call34$i = 0, $call65$i = 0, $call80$i = 0, $child$i$i = 0, $cmp101$i = 0, $cmp138$i164 = 0, $cmp32$i = 0, $cmp66$i158 = 0, $cmp82$i = 0, $cond = 0, $cond$i = 0, $cond$i$i = 0, $cond$i$i$i = 0, $cond$i17$i = 0, $cond$i27$i = 0, $cond$i42$i = 0, $cond115$i$i = 0, $cond13$i$i = 0, $cond15$i$i = 0, $cond315$i$i = 0, $cond382$i = 0, $cond6$i = 0, $fd138$i = 0, $fd145$i$i = 0, $fd344$i$i = 0, $fd412$i = 0, $fd50$i = 0, $fd59$i$i = 0, $fd68$pre$phi$i$iZ2D = 0, $fd69 = 0, $fd85$i$i = 0, $fd9 = 0, $head178 = 0, $head182$i = 0, $head208$i$i = 0, $head25 = 0, $head273$i = 0, $head31$i$i = 0, $i$02$i$i = 0, $idx$0$i = 0, $mem$0 = 0, $nb$0 = 0, $neg$i149 = 0, $oldfirst$0$i$i = 0, $qsize$0$i$i = 0, $rsize$0$i = 0, $rsize$0$i120 = 0, $rsize$1$i = 0, $rsize$2$i = 0, $rsize$3$lcssa$i = 0, $rsize$328$i = 0, $rst$0$i = 0, $rst$1$i = 0, $shl = 0, $shl$i$i = 0, $shl$i111 = 0, $shl$i20$i = 0, $shl102 = 0, $shl105 = 0, $shl195$i = 0, $shl198$i = 0, $shl22 = 0, $shl221$i$i = 0, $shl226$i$i = 0, $shl265$i$i = 0, $shl270$i$i = 0, $shl287$i = 0, $shl290$i = 0, $shl294$i$i = 0, $shl332$i = 0, $shl337$i = 0, $shl361$i = 0, $shl37 = 0, $shl39$i$i = 0, $shl59$i = 0, $shl65 = 0, $shl70$i$i = 0, $shl75$i$i = 0, $shl9$i = 0, $shl90 = 0, $shl95$i$i = 0, $shr = 0, $shr$i$i = 0, $shr$i106 = 0, $shr$i54$i = 0, $shr101 = 0, $shr11$i = 0, $shr15$i = 0, $shr194$i = 0, $shr214$i$i = 0, $shr253$i$i = 0, $shr282$i = 0, $shr3 = 0, $shr317$i = 0, $shr4$i = 0, $shr47 = 0, $shr51 = 0, $shr55 = 0, $shr58$i$i = 0, $shr59 = 0, $shr7$i = 0, $shr74$i = 0, $shr78$i = 0, $shr82$i = 0, $shr86$i = 0, $size$i$i = 0, $size185$i = 0, $size242$i = 0, $sizebits$0$i = 0, $sp$0$i$i = 0, $sp$0$i$i$i = 0, $sp$0109$i = 0, $sp$1105$i = 0, $ssize$0$i = 0, $ssize$1$i = 0, $ssize$2$i = 0, $sub$i105 = 0, $sub$i148 = 0, $sub$ptr$sub$i = 0, $sub$ptr$sub$i$i = 0, $sub100$i = 0, $sub100$rsize$3$i = 0, $sub109$i = 0, $sub159 = 0, $sub18$i$i = 0, $sub187 = 0, $sub2$i = 0, $sub253$i = 0, $sub31$i = 0, $sub33$i = 0, $sub38$i = 0, $sub44 = 0, $sub5$i$i = 0, $sub5$i$i$i = 0, $sub5$i29$i = 0, $sub69$i = 0, $sub91 = 0, $t$0$i = 0, $t$0$i119 = 0, $t$1$i = 0, $t$2$ph$i = 0, $t$2$v$3$i = 0, $t$227$i = 0, $tbase$0$i = 0, $tbase$291$i = 0, $tsize$0$i = 0, $tsize$0748284$i = 0, $tsize$1$i = 0, $tsize$290$i = 0, $v$0$i = 0, $v$0$i121 = 0, $v$1$i = 0, $v$2$i = 0, $v$3$lcssa$i = 0, $v$329$i = 0, label = 0, sp = 0, $add$ptr2420$i$i$looptemp = 0;
 sp = STACKTOP;
 do if ($bytes >>> 0 < 245) {
  if ($bytes >>> 0 < 11) $cond = 16; else $cond = $bytes + 11 & -8;
  $shr = $cond >>> 3;
  $0 = HEAP32[1002] | 0;
  $shr3 = $0 >>> $shr;
  if (($shr3 & 3 | 0) != 0) {
   $add8 = ($shr3 & 1 ^ 1) + $shr | 0;
   $shl = $add8 << 1;
   $arrayidx = 4048 + ($shl << 2) | 0;
   $1 = 4048 + ($shl + 2 << 2) | 0;
   $2 = HEAP32[$1 >> 2] | 0;
   $fd9 = $2 + 8 | 0;
   $3 = HEAP32[$fd9 >> 2] | 0;
   do if (($arrayidx | 0) == ($3 | 0)) HEAP32[1002] = $0 & ~(1 << $add8); else {
    if ($3 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
    $bk = $3 + 12 | 0;
    if ((HEAP32[$bk >> 2] | 0) == ($2 | 0)) {
     HEAP32[$bk >> 2] = $arrayidx;
     HEAP32[$1 >> 2] = $3;
     break;
    } else _abort();
   } while (0);
   $shl22 = $add8 << 3;
   HEAP32[$2 + 4 >> 2] = $shl22 | 3;
   $head25 = $2 + ($shl22 | 4) | 0;
   HEAP32[$head25 >> 2] = HEAP32[$head25 >> 2] | 1;
   $mem$0 = $fd9;
   STACKTOP = sp;
   return $mem$0 | 0;
  }
  if ($cond >>> 0 > (HEAP32[1004] | 0) >>> 0) {
   if (($shr3 | 0) != 0) {
    $shl37 = 2 << $shr;
    $and41 = $shr3 << $shr & ($shl37 | 0 - $shl37);
    $sub44 = ($and41 & 0 - $and41) + -1 | 0;
    $and46 = $sub44 >>> 12 & 16;
    $shr47 = $sub44 >>> $and46;
    $and49 = $shr47 >>> 5 & 8;
    $shr51 = $shr47 >>> $and49;
    $and53 = $shr51 >>> 2 & 4;
    $shr55 = $shr51 >>> $and53;
    $and57 = $shr55 >>> 1 & 2;
    $shr59 = $shr55 >>> $and57;
    $and61 = $shr59 >>> 1 & 1;
    $add64 = ($and49 | $and46 | $and53 | $and57 | $and61) + ($shr59 >>> $and61) | 0;
    $shl65 = $add64 << 1;
    $arrayidx66 = 4048 + ($shl65 << 2) | 0;
    $8 = 4048 + ($shl65 + 2 << 2) | 0;
    $9 = HEAP32[$8 >> 2] | 0;
    $fd69 = $9 + 8 | 0;
    $10 = HEAP32[$fd69 >> 2] | 0;
    do if (($arrayidx66 | 0) == ($10 | 0)) HEAP32[1002] = $0 & ~(1 << $add64); else {
     if ($10 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
     $bk78 = $10 + 12 | 0;
     if ((HEAP32[$bk78 >> 2] | 0) == ($9 | 0)) {
      HEAP32[$bk78 >> 2] = $arrayidx66;
      HEAP32[$8 >> 2] = $10;
      break;
     } else _abort();
    } while (0);
    $shl90 = $add64 << 3;
    $sub91 = $shl90 - $cond | 0;
    HEAP32[$9 + 4 >> 2] = $cond | 3;
    $add$ptr95 = $9 + $cond | 0;
    HEAP32[$9 + ($cond | 4) >> 2] = $sub91 | 1;
    HEAP32[$9 + $shl90 >> 2] = $sub91;
    $13 = HEAP32[1004] | 0;
    if (($13 | 0) != 0) {
     $14 = HEAP32[1007] | 0;
     $shr101 = $13 >>> 3;
     $shl102 = $shr101 << 1;
     $arrayidx103 = 4048 + ($shl102 << 2) | 0;
     $15 = HEAP32[1002] | 0;
     $shl105 = 1 << $shr101;
     if (($15 & $shl105 | 0) == 0) {
      HEAP32[1002] = $15 | $shl105;
      $$pre$phiZ2D = 4048 + ($shl102 + 2 << 2) | 0;
      $F104$0 = $arrayidx103;
     } else {
      $16 = 4048 + ($shl102 + 2 << 2) | 0;
      $17 = HEAP32[$16 >> 2] | 0;
      if ($17 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
       $$pre$phiZ2D = $16;
       $F104$0 = $17;
      }
     }
     HEAP32[$$pre$phiZ2D >> 2] = $14;
     HEAP32[$F104$0 + 12 >> 2] = $14;
     HEAP32[$14 + 8 >> 2] = $F104$0;
     HEAP32[$14 + 12 >> 2] = $arrayidx103;
    }
    HEAP32[1004] = $sub91;
    HEAP32[1007] = $add$ptr95;
    $mem$0 = $fd69;
    STACKTOP = sp;
    return $mem$0 | 0;
   }
   $19 = HEAP32[1003] | 0;
   if (($19 | 0) == 0) $nb$0 = $cond; else {
    $sub2$i = ($19 & 0 - $19) + -1 | 0;
    $and3$i = $sub2$i >>> 12 & 16;
    $shr4$i = $sub2$i >>> $and3$i;
    $and6$i = $shr4$i >>> 5 & 8;
    $shr7$i = $shr4$i >>> $and6$i;
    $and9$i = $shr7$i >>> 2 & 4;
    $shr11$i = $shr7$i >>> $and9$i;
    $and13$i = $shr11$i >>> 1 & 2;
    $shr15$i = $shr11$i >>> $and13$i;
    $and17$i = $shr15$i >>> 1 & 1;
    $20 = HEAP32[4312 + (($and6$i | $and3$i | $and9$i | $and13$i | $and17$i) + ($shr15$i >>> $and17$i) << 2) >> 2] | 0;
    $rsize$0$i = (HEAP32[$20 + 4 >> 2] & -8) - $cond | 0;
    $t$0$i = $20;
    $v$0$i = $20;
    while (1) {
     $22 = HEAP32[$t$0$i + 16 >> 2] | 0;
     if (($22 | 0) == 0) {
      $23 = HEAP32[$t$0$i + 20 >> 2] | 0;
      if (($23 | 0) == 0) break; else $cond6$i = $23;
     } else $cond6$i = $22;
     $sub31$i = (HEAP32[$cond6$i + 4 >> 2] & -8) - $cond | 0;
     $cmp32$i = $sub31$i >>> 0 < $rsize$0$i >>> 0;
     $rsize$0$i = $cmp32$i ? $sub31$i : $rsize$0$i;
     $t$0$i = $cond6$i;
     $v$0$i = $cmp32$i ? $cond6$i : $v$0$i;
    }
    $25 = HEAP32[1006] | 0;
    if ($v$0$i >>> 0 < $25 >>> 0) _abort();
    $add$ptr$i = $v$0$i + $cond | 0;
    if (!($v$0$i >>> 0 < $add$ptr$i >>> 0)) _abort();
    $26 = HEAP32[$v$0$i + 24 >> 2] | 0;
    $27 = HEAP32[$v$0$i + 12 >> 2] | 0;
    do if (($27 | 0) == ($v$0$i | 0)) {
     $arrayidx61$i = $v$0$i + 20 | 0;
     $31 = HEAP32[$arrayidx61$i >> 2] | 0;
     if (($31 | 0) == 0) {
      $arrayidx65$i = $v$0$i + 16 | 0;
      $32 = HEAP32[$arrayidx65$i >> 2] | 0;
      if (($32 | 0) == 0) {
       $R$1$i = 0;
       break;
      } else {
       $R$0$i = $32;
       $RP$0$i = $arrayidx65$i;
      }
     } else {
      $R$0$i = $31;
      $RP$0$i = $arrayidx61$i;
     }
     while (1) {
      $arrayidx71$i = $R$0$i + 20 | 0;
      $33 = HEAP32[$arrayidx71$i >> 2] | 0;
      if (($33 | 0) != 0) {
       $R$0$i = $33;
       $RP$0$i = $arrayidx71$i;
       continue;
      }
      $arrayidx75$i = $R$0$i + 16 | 0;
      $34 = HEAP32[$arrayidx75$i >> 2] | 0;
      if (($34 | 0) == 0) break; else {
       $R$0$i = $34;
       $RP$0$i = $arrayidx75$i;
      }
     }
     if ($RP$0$i >>> 0 < $25 >>> 0) _abort(); else {
      HEAP32[$RP$0$i >> 2] = 0;
      $R$1$i = $R$0$i;
      break;
     }
    } else {
     $28 = HEAP32[$v$0$i + 8 >> 2] | 0;
     if ($28 >>> 0 < $25 >>> 0) _abort();
     $bk47$i = $28 + 12 | 0;
     if ((HEAP32[$bk47$i >> 2] | 0) != ($v$0$i | 0)) _abort();
     $fd50$i = $27 + 8 | 0;
     if ((HEAP32[$fd50$i >> 2] | 0) == ($v$0$i | 0)) {
      HEAP32[$bk47$i >> 2] = $27;
      HEAP32[$fd50$i >> 2] = $28;
      $R$1$i = $27;
      break;
     } else _abort();
    } while (0);
    do if (($26 | 0) != 0) {
     $35 = HEAP32[$v$0$i + 28 >> 2] | 0;
     $arrayidx94$i = 4312 + ($35 << 2) | 0;
     if (($v$0$i | 0) == (HEAP32[$arrayidx94$i >> 2] | 0)) {
      HEAP32[$arrayidx94$i >> 2] = $R$1$i;
      if (($R$1$i | 0) == 0) {
       HEAP32[1003] = HEAP32[1003] & ~(1 << $35);
       break;
      }
     } else {
      if ($26 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
      $arrayidx113$i = $26 + 16 | 0;
      if ((HEAP32[$arrayidx113$i >> 2] | 0) == ($v$0$i | 0)) HEAP32[$arrayidx113$i >> 2] = $R$1$i; else HEAP32[$26 + 20 >> 2] = $R$1$i;
      if (($R$1$i | 0) == 0) break;
     }
     if ($R$1$i >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
     HEAP32[$R$1$i + 24 >> 2] = $26;
     $41 = HEAP32[$v$0$i + 16 >> 2] | 0;
     do if (($41 | 0) != 0) if ($41 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
      HEAP32[$R$1$i + 16 >> 2] = $41;
      HEAP32[$41 + 24 >> 2] = $R$1$i;
      break;
     } while (0);
     $43 = HEAP32[$v$0$i + 20 >> 2] | 0;
     if (($43 | 0) != 0) if ($43 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
      HEAP32[$R$1$i + 20 >> 2] = $43;
      HEAP32[$43 + 24 >> 2] = $R$1$i;
      break;
     }
    } while (0);
    if ($rsize$0$i >>> 0 < 16) {
     $add177$i = $rsize$0$i + $cond | 0;
     HEAP32[$v$0$i + 4 >> 2] = $add177$i | 3;
     $head182$i = $v$0$i + ($add177$i + 4) | 0;
     HEAP32[$head182$i >> 2] = HEAP32[$head182$i >> 2] | 1;
    } else {
     HEAP32[$v$0$i + 4 >> 2] = $cond | 3;
     HEAP32[$v$0$i + ($cond | 4) >> 2] = $rsize$0$i | 1;
     HEAP32[$v$0$i + ($rsize$0$i + $cond) >> 2] = $rsize$0$i;
     $46 = HEAP32[1004] | 0;
     if (($46 | 0) != 0) {
      $47 = HEAP32[1007] | 0;
      $shr194$i = $46 >>> 3;
      $shl195$i = $shr194$i << 1;
      $arrayidx196$i = 4048 + ($shl195$i << 2) | 0;
      $48 = HEAP32[1002] | 0;
      $shl198$i = 1 << $shr194$i;
      if (($48 & $shl198$i | 0) == 0) {
       HEAP32[1002] = $48 | $shl198$i;
       $$pre$phi$iZ2D = 4048 + ($shl195$i + 2 << 2) | 0;
       $F197$0$i = $arrayidx196$i;
      } else {
       $49 = 4048 + ($shl195$i + 2 << 2) | 0;
       $50 = HEAP32[$49 >> 2] | 0;
       if ($50 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
        $$pre$phi$iZ2D = $49;
        $F197$0$i = $50;
       }
      }
      HEAP32[$$pre$phi$iZ2D >> 2] = $47;
      HEAP32[$F197$0$i + 12 >> 2] = $47;
      HEAP32[$47 + 8 >> 2] = $F197$0$i;
      HEAP32[$47 + 12 >> 2] = $arrayidx196$i;
     }
     HEAP32[1004] = $rsize$0$i;
     HEAP32[1007] = $add$ptr$i;
    }
    $mem$0 = $v$0$i + 8 | 0;
    STACKTOP = sp;
    return $mem$0 | 0;
   }
  } else $nb$0 = $cond;
 } else if ($bytes >>> 0 > 4294967231) $nb$0 = -1; else {
  $add143 = $bytes + 11 | 0;
  $and144 = $add143 & -8;
  $52 = HEAP32[1003] | 0;
  if (($52 | 0) == 0) $nb$0 = $and144; else {
   $sub$i105 = 0 - $and144 | 0;
   $shr$i106 = $add143 >>> 8;
   if (($shr$i106 | 0) == 0) $idx$0$i = 0; else if ($and144 >>> 0 > 16777215) $idx$0$i = 31; else {
    $and$i110 = ($shr$i106 + 1048320 | 0) >>> 16 & 8;
    $shl$i111 = $shr$i106 << $and$i110;
    $and8$i = ($shl$i111 + 520192 | 0) >>> 16 & 4;
    $shl9$i = $shl$i111 << $and8$i;
    $and12$i = ($shl9$i + 245760 | 0) >>> 16 & 2;
    $add17$i = 14 - ($and8$i | $and$i110 | $and12$i) + ($shl9$i << $and12$i >>> 15) | 0;
    $idx$0$i = $and144 >>> ($add17$i + 7 | 0) & 1 | $add17$i << 1;
   }
   $53 = HEAP32[4312 + ($idx$0$i << 2) >> 2] | 0;
   L126 : do if (($53 | 0) == 0) {
    $rsize$2$i = $sub$i105;
    $t$1$i = 0;
    $v$2$i = 0;
   } else {
    if (($idx$0$i | 0) == 31) $cond$i = 0; else $cond$i = 25 - ($idx$0$i >>> 1) | 0;
    $rsize$0$i120 = $sub$i105;
    $rst$0$i = 0;
    $sizebits$0$i = $and144 << $cond$i;
    $t$0$i119 = $53;
    $v$0$i121 = 0;
    while (1) {
     $and32$i = HEAP32[$t$0$i119 + 4 >> 2] & -8;
     $sub33$i = $and32$i - $and144 | 0;
     if ($sub33$i >>> 0 < $rsize$0$i120 >>> 0) if (($and32$i | 0) == ($and144 | 0)) {
      $rsize$2$i = $sub33$i;
      $t$1$i = $t$0$i119;
      $v$2$i = $t$0$i119;
      break L126;
     } else {
      $rsize$1$i = $sub33$i;
      $v$1$i = $t$0$i119;
     } else {
      $rsize$1$i = $rsize$0$i120;
      $v$1$i = $v$0$i121;
     }
     $55 = HEAP32[$t$0$i119 + 20 >> 2] | 0;
     $t$0$i119 = HEAP32[$t$0$i119 + ($sizebits$0$i >>> 31 << 2) + 16 >> 2] | 0;
     $rst$1$i = ($55 | 0) == 0 | ($55 | 0) == ($t$0$i119 | 0) ? $rst$0$i : $55;
     if (($t$0$i119 | 0) == 0) {
      $rsize$2$i = $rsize$1$i;
      $t$1$i = $rst$1$i;
      $v$2$i = $v$1$i;
      break;
     } else {
      $rsize$0$i120 = $rsize$1$i;
      $rst$0$i = $rst$1$i;
      $sizebits$0$i = $sizebits$0$i << 1;
      $v$0$i121 = $v$1$i;
     }
    }
   } while (0);
   if (($t$1$i | 0) == 0 & ($v$2$i | 0) == 0) {
    $shl59$i = 2 << $idx$0$i;
    $and63$i = $52 & ($shl59$i | 0 - $shl59$i);
    if (($and63$i | 0) == 0) {
     $nb$0 = $and144;
     break;
    }
    $sub69$i = ($and63$i & 0 - $and63$i) + -1 | 0;
    $and72$i = $sub69$i >>> 12 & 16;
    $shr74$i = $sub69$i >>> $and72$i;
    $and76$i = $shr74$i >>> 5 & 8;
    $shr78$i = $shr74$i >>> $and76$i;
    $and80$i = $shr78$i >>> 2 & 4;
    $shr82$i = $shr78$i >>> $and80$i;
    $and84$i = $shr82$i >>> 1 & 2;
    $shr86$i = $shr82$i >>> $and84$i;
    $and88$i = $shr86$i >>> 1 & 1;
    $t$2$ph$i = HEAP32[4312 + (($and76$i | $and72$i | $and80$i | $and84$i | $and88$i) + ($shr86$i >>> $and88$i) << 2) >> 2] | 0;
   } else $t$2$ph$i = $t$1$i;
   if (($t$2$ph$i | 0) == 0) {
    $rsize$3$lcssa$i = $rsize$2$i;
    $v$3$lcssa$i = $v$2$i;
   } else {
    $rsize$328$i = $rsize$2$i;
    $t$227$i = $t$2$ph$i;
    $v$329$i = $v$2$i;
    while (1) {
     $sub100$i = (HEAP32[$t$227$i + 4 >> 2] & -8) - $and144 | 0;
     $cmp101$i = $sub100$i >>> 0 < $rsize$328$i >>> 0;
     $sub100$rsize$3$i = $cmp101$i ? $sub100$i : $rsize$328$i;
     $t$2$v$3$i = $cmp101$i ? $t$227$i : $v$329$i;
     $59 = HEAP32[$t$227$i + 16 >> 2] | 0;
     if (($59 | 0) != 0) {
      $rsize$328$i = $sub100$rsize$3$i;
      $t$227$i = $59;
      $v$329$i = $t$2$v$3$i;
      continue;
     }
     $t$227$i = HEAP32[$t$227$i + 20 >> 2] | 0;
     if (($t$227$i | 0) == 0) {
      $rsize$3$lcssa$i = $sub100$rsize$3$i;
      $v$3$lcssa$i = $t$2$v$3$i;
      break;
     } else {
      $rsize$328$i = $sub100$rsize$3$i;
      $v$329$i = $t$2$v$3$i;
     }
    }
   }
   if (($v$3$lcssa$i | 0) == 0) $nb$0 = $and144; else if ($rsize$3$lcssa$i >>> 0 < ((HEAP32[1004] | 0) - $and144 | 0) >>> 0) {
    $62 = HEAP32[1006] | 0;
    if ($v$3$lcssa$i >>> 0 < $62 >>> 0) _abort();
    $add$ptr$i126 = $v$3$lcssa$i + $and144 | 0;
    if (!($v$3$lcssa$i >>> 0 < $add$ptr$i126 >>> 0)) _abort();
    $63 = HEAP32[$v$3$lcssa$i + 24 >> 2] | 0;
    $64 = HEAP32[$v$3$lcssa$i + 12 >> 2] | 0;
    do if (($64 | 0) == ($v$3$lcssa$i | 0)) {
     $arrayidx150$i = $v$3$lcssa$i + 20 | 0;
     $68 = HEAP32[$arrayidx150$i >> 2] | 0;
     if (($68 | 0) == 0) {
      $arrayidx154$i131 = $v$3$lcssa$i + 16 | 0;
      $69 = HEAP32[$arrayidx154$i131 >> 2] | 0;
      if (($69 | 0) == 0) {
       $R$1$i137 = 0;
       break;
      } else {
       $R$0$i135 = $69;
       $RP$0$i134 = $arrayidx154$i131;
      }
     } else {
      $R$0$i135 = $68;
      $RP$0$i134 = $arrayidx150$i;
     }
     while (1) {
      $arrayidx160$i = $R$0$i135 + 20 | 0;
      $70 = HEAP32[$arrayidx160$i >> 2] | 0;
      if (($70 | 0) != 0) {
       $R$0$i135 = $70;
       $RP$0$i134 = $arrayidx160$i;
       continue;
      }
      $arrayidx164$i = $R$0$i135 + 16 | 0;
      $71 = HEAP32[$arrayidx164$i >> 2] | 0;
      if (($71 | 0) == 0) break; else {
       $R$0$i135 = $71;
       $RP$0$i134 = $arrayidx164$i;
      }
     }
     if ($RP$0$i134 >>> 0 < $62 >>> 0) _abort(); else {
      HEAP32[$RP$0$i134 >> 2] = 0;
      $R$1$i137 = $R$0$i135;
      break;
     }
    } else {
     $65 = HEAP32[$v$3$lcssa$i + 8 >> 2] | 0;
     if ($65 >>> 0 < $62 >>> 0) _abort();
     $bk135$i = $65 + 12 | 0;
     if ((HEAP32[$bk135$i >> 2] | 0) != ($v$3$lcssa$i | 0)) _abort();
     $fd138$i = $64 + 8 | 0;
     if ((HEAP32[$fd138$i >> 2] | 0) == ($v$3$lcssa$i | 0)) {
      HEAP32[$bk135$i >> 2] = $64;
      HEAP32[$fd138$i >> 2] = $65;
      $R$1$i137 = $64;
      break;
     } else _abort();
    } while (0);
    do if (($63 | 0) != 0) {
     $72 = HEAP32[$v$3$lcssa$i + 28 >> 2] | 0;
     $arrayidx183$i = 4312 + ($72 << 2) | 0;
     if (($v$3$lcssa$i | 0) == (HEAP32[$arrayidx183$i >> 2] | 0)) {
      HEAP32[$arrayidx183$i >> 2] = $R$1$i137;
      if (($R$1$i137 | 0) == 0) {
       HEAP32[1003] = HEAP32[1003] & ~(1 << $72);
       break;
      }
     } else {
      if ($63 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
      $arrayidx203$i = $63 + 16 | 0;
      if ((HEAP32[$arrayidx203$i >> 2] | 0) == ($v$3$lcssa$i | 0)) HEAP32[$arrayidx203$i >> 2] = $R$1$i137; else HEAP32[$63 + 20 >> 2] = $R$1$i137;
      if (($R$1$i137 | 0) == 0) break;
     }
     if ($R$1$i137 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
     HEAP32[$R$1$i137 + 24 >> 2] = $63;
     $78 = HEAP32[$v$3$lcssa$i + 16 >> 2] | 0;
     do if (($78 | 0) != 0) if ($78 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
      HEAP32[$R$1$i137 + 16 >> 2] = $78;
      HEAP32[$78 + 24 >> 2] = $R$1$i137;
      break;
     } while (0);
     $80 = HEAP32[$v$3$lcssa$i + 20 >> 2] | 0;
     if (($80 | 0) != 0) if ($80 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
      HEAP32[$R$1$i137 + 20 >> 2] = $80;
      HEAP32[$80 + 24 >> 2] = $R$1$i137;
      break;
     }
    } while (0);
    L204 : do if ($rsize$3$lcssa$i >>> 0 < 16) {
     $add267$i = $rsize$3$lcssa$i + $and144 | 0;
     HEAP32[$v$3$lcssa$i + 4 >> 2] = $add267$i | 3;
     $head273$i = $v$3$lcssa$i + ($add267$i + 4) | 0;
     HEAP32[$head273$i >> 2] = HEAP32[$head273$i >> 2] | 1;
    } else {
     HEAP32[$v$3$lcssa$i + 4 >> 2] = $and144 | 3;
     HEAP32[$v$3$lcssa$i + ($and144 | 4) >> 2] = $rsize$3$lcssa$i | 1;
     HEAP32[$v$3$lcssa$i + ($rsize$3$lcssa$i + $and144) >> 2] = $rsize$3$lcssa$i;
     $shr282$i = $rsize$3$lcssa$i >>> 3;
     if ($rsize$3$lcssa$i >>> 0 < 256) {
      $shl287$i = $shr282$i << 1;
      $arrayidx288$i = 4048 + ($shl287$i << 2) | 0;
      $83 = HEAP32[1002] | 0;
      $shl290$i = 1 << $shr282$i;
      do if (($83 & $shl290$i | 0) == 0) {
       HEAP32[1002] = $83 | $shl290$i;
       $$pre$phi$i145Z2D = 4048 + ($shl287$i + 2 << 2) | 0;
       $F289$0$i = $arrayidx288$i;
      } else {
       $84 = 4048 + ($shl287$i + 2 << 2) | 0;
       $85 = HEAP32[$84 >> 2] | 0;
       if (!($85 >>> 0 < (HEAP32[1006] | 0) >>> 0)) {
        $$pre$phi$i145Z2D = $84;
        $F289$0$i = $85;
        break;
       }
       _abort();
      } while (0);
      HEAP32[$$pre$phi$i145Z2D >> 2] = $add$ptr$i126;
      HEAP32[$F289$0$i + 12 >> 2] = $add$ptr$i126;
      HEAP32[$v$3$lcssa$i + ($and144 + 8) >> 2] = $F289$0$i;
      HEAP32[$v$3$lcssa$i + ($and144 + 12) >> 2] = $arrayidx288$i;
      break;
     }
     $shr317$i = $rsize$3$lcssa$i >>> 8;
     if (($shr317$i | 0) == 0) $I315$0$i = 0; else if ($rsize$3$lcssa$i >>> 0 > 16777215) $I315$0$i = 31; else {
      $and330$i = ($shr317$i + 1048320 | 0) >>> 16 & 8;
      $shl332$i = $shr317$i << $and330$i;
      $and335$i = ($shl332$i + 520192 | 0) >>> 16 & 4;
      $shl337$i = $shl332$i << $and335$i;
      $and340$i = ($shl337$i + 245760 | 0) >>> 16 & 2;
      $add345$i = 14 - ($and335$i | $and330$i | $and340$i) + ($shl337$i << $and340$i >>> 15) | 0;
      $I315$0$i = $rsize$3$lcssa$i >>> ($add345$i + 7 | 0) & 1 | $add345$i << 1;
     }
     $arrayidx354$i = 4312 + ($I315$0$i << 2) | 0;
     HEAP32[$v$3$lcssa$i + ($and144 + 28) >> 2] = $I315$0$i;
     HEAP32[$v$3$lcssa$i + ($and144 + 20) >> 2] = 0;
     HEAP32[$v$3$lcssa$i + ($and144 + 16) >> 2] = 0;
     $87 = HEAP32[1003] | 0;
     $shl361$i = 1 << $I315$0$i;
     if (($87 & $shl361$i | 0) == 0) {
      HEAP32[1003] = $87 | $shl361$i;
      HEAP32[$arrayidx354$i >> 2] = $add$ptr$i126;
      HEAP32[$v$3$lcssa$i + ($and144 + 24) >> 2] = $arrayidx354$i;
      HEAP32[$v$3$lcssa$i + ($and144 + 12) >> 2] = $add$ptr$i126;
      HEAP32[$v$3$lcssa$i + ($and144 + 8) >> 2] = $add$ptr$i126;
      break;
     }
     $88 = HEAP32[$arrayidx354$i >> 2] | 0;
     if (($I315$0$i | 0) == 31) $cond382$i = 0; else $cond382$i = 25 - ($I315$0$i >>> 1) | 0;
     L225 : do if ((HEAP32[$88 + 4 >> 2] & -8 | 0) == ($rsize$3$lcssa$i | 0)) $T$0$lcssa$i = $88; else {
      $K372$024$i = $rsize$3$lcssa$i << $cond382$i;
      $T$023$i = $88;
      while (1) {
       $arrayidx393$i = $T$023$i + ($K372$024$i >>> 31 << 2) + 16 | 0;
       $90 = HEAP32[$arrayidx393$i >> 2] | 0;
       if (($90 | 0) == 0) break;
       if ((HEAP32[$90 + 4 >> 2] & -8 | 0) == ($rsize$3$lcssa$i | 0)) {
        $T$0$lcssa$i = $90;
        break L225;
       } else {
        $K372$024$i = $K372$024$i << 1;
        $T$023$i = $90;
       }
      }
      if ($arrayidx393$i >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
       HEAP32[$arrayidx393$i >> 2] = $add$ptr$i126;
       HEAP32[$v$3$lcssa$i + ($and144 + 24) >> 2] = $T$023$i;
       HEAP32[$v$3$lcssa$i + ($and144 + 12) >> 2] = $add$ptr$i126;
       HEAP32[$v$3$lcssa$i + ($and144 + 8) >> 2] = $add$ptr$i126;
       break L204;
      }
     } while (0);
     $fd412$i = $T$0$lcssa$i + 8 | 0;
     $93 = HEAP32[$fd412$i >> 2] | 0;
     $94 = HEAP32[1006] | 0;
     if ($T$0$lcssa$i >>> 0 < $94 >>> 0) _abort();
     if ($93 >>> 0 < $94 >>> 0) _abort(); else {
      HEAP32[$93 + 12 >> 2] = $add$ptr$i126;
      HEAP32[$fd412$i >> 2] = $add$ptr$i126;
      HEAP32[$v$3$lcssa$i + ($and144 + 8) >> 2] = $93;
      HEAP32[$v$3$lcssa$i + ($and144 + 12) >> 2] = $T$0$lcssa$i;
      HEAP32[$v$3$lcssa$i + ($and144 + 24) >> 2] = 0;
      break;
     }
    } while (0);
    $mem$0 = $v$3$lcssa$i + 8 | 0;
    STACKTOP = sp;
    return $mem$0 | 0;
   } else $nb$0 = $and144;
  }
 } while (0);
 $95 = HEAP32[1004] | 0;
 if (!($nb$0 >>> 0 > $95 >>> 0)) {
  $sub159 = $95 - $nb$0 | 0;
  $96 = HEAP32[1007] | 0;
  if ($sub159 >>> 0 > 15) {
   HEAP32[1007] = $96 + $nb$0;
   HEAP32[1004] = $sub159;
   HEAP32[$96 + ($nb$0 + 4) >> 2] = $sub159 | 1;
   HEAP32[$96 + $95 >> 2] = $sub159;
   HEAP32[$96 + 4 >> 2] = $nb$0 | 3;
  } else {
   HEAP32[1004] = 0;
   HEAP32[1007] = 0;
   HEAP32[$96 + 4 >> 2] = $95 | 3;
   $head178 = $96 + ($95 + 4) | 0;
   HEAP32[$head178 >> 2] = HEAP32[$head178 >> 2] | 1;
  }
  $mem$0 = $96 + 8 | 0;
  STACKTOP = sp;
  return $mem$0 | 0;
 }
 $98 = HEAP32[1005] | 0;
 if ($nb$0 >>> 0 < $98 >>> 0) {
  $sub187 = $98 - $nb$0 | 0;
  HEAP32[1005] = $sub187;
  $99 = HEAP32[1008] | 0;
  HEAP32[1008] = $99 + $nb$0;
  HEAP32[$99 + ($nb$0 + 4) >> 2] = $sub187 | 1;
  HEAP32[$99 + 4 >> 2] = $nb$0 | 3;
  $mem$0 = $99 + 8 | 0;
  STACKTOP = sp;
  return $mem$0 | 0;
 }
 do if ((HEAP32[1120] | 0) == 0) {
  $call$i$i = _sysconf(30) | 0;
  if (($call$i$i + -1 & $call$i$i | 0) == 0) {
   HEAP32[1122] = $call$i$i;
   HEAP32[1121] = $call$i$i;
   HEAP32[1123] = -1;
   HEAP32[1124] = -1;
   HEAP32[1125] = 0;
   HEAP32[1113] = 0;
   HEAP32[1120] = (_time(0) | 0) & -16 ^ 1431655768;
   break;
  } else _abort();
 } while (0);
 $add$i147 = $nb$0 + 48 | 0;
 $101 = HEAP32[1122] | 0;
 $sub$i148 = $nb$0 + 47 | 0;
 $add9$i = $101 + $sub$i148 | 0;
 $neg$i149 = 0 - $101 | 0;
 $and11$i = $add9$i & $neg$i149;
 if (!($and11$i >>> 0 > $nb$0 >>> 0)) {
  $mem$0 = 0;
  STACKTOP = sp;
  return $mem$0 | 0;
 }
 $102 = HEAP32[1112] | 0;
 if (($102 | 0) != 0) {
  $103 = HEAP32[1110] | 0;
  $add17$i150 = $103 + $and11$i | 0;
  if ($add17$i150 >>> 0 <= $103 >>> 0 | $add17$i150 >>> 0 > $102 >>> 0) {
   $mem$0 = 0;
   STACKTOP = sp;
   return $mem$0 | 0;
  }
 }
 L269 : do if ((HEAP32[1113] & 4 | 0) == 0) {
  $105 = HEAP32[1008] | 0;
  L271 : do if (($105 | 0) == 0) label = 182; else {
   $sp$0$i$i = 4456 | 0;
   while (1) {
    $106 = HEAP32[$sp$0$i$i >> 2] | 0;
    if (!($106 >>> 0 > $105 >>> 0)) {
     $size$i$i = $sp$0$i$i + 4 | 0;
     if (($106 + (HEAP32[$size$i$i >> 2] | 0) | 0) >>> 0 > $105 >>> 0) break;
    }
    $108 = HEAP32[$sp$0$i$i + 8 >> 2] | 0;
    if (($108 | 0) == 0) {
     label = 182;
     break L271;
    } else $sp$0$i$i = $108;
   }
   if (($sp$0$i$i | 0) == 0) label = 182; else {
    $and77$i = $add9$i - (HEAP32[1005] | 0) & $neg$i149;
    if ($and77$i >>> 0 < 2147483647) {
     $call80$i = _sbrk($and77$i | 0) | 0;
     $cmp82$i = ($call80$i | 0) == ((HEAP32[$sp$0$i$i >> 2] | 0) + (HEAP32[$size$i$i >> 2] | 0) | 0);
     $br$0$i = $call80$i;
     $ssize$1$i = $and77$i;
     $tbase$0$i = $cmp82$i ? $call80$i : -1;
     $tsize$0$i = $cmp82$i ? $and77$i : 0;
     label = 191;
    } else $tsize$0748284$i = 0;
   }
  } while (0);
  do if ((label | 0) == 182) {
   $call34$i = _sbrk(0) | 0;
   if (($call34$i | 0) == (-1 | 0)) $tsize$0748284$i = 0; else {
    $109 = $call34$i;
    $110 = HEAP32[1121] | 0;
    $sub38$i = $110 + -1 | 0;
    if (($sub38$i & $109 | 0) == 0) $ssize$0$i = $and11$i; else $ssize$0$i = $and11$i - $109 + ($sub38$i + $109 & 0 - $110) | 0;
    $111 = HEAP32[1110] | 0;
    $add51$i = $111 + $ssize$0$i | 0;
    if ($ssize$0$i >>> 0 > $nb$0 >>> 0 & $ssize$0$i >>> 0 < 2147483647) {
     $112 = HEAP32[1112] | 0;
     if (($112 | 0) != 0) if ($add51$i >>> 0 <= $111 >>> 0 | $add51$i >>> 0 > $112 >>> 0) {
      $tsize$0748284$i = 0;
      break;
     }
     $call65$i = _sbrk($ssize$0$i | 0) | 0;
     $cmp66$i158 = ($call65$i | 0) == ($call34$i | 0);
     $br$0$i = $call65$i;
     $ssize$1$i = $ssize$0$i;
     $tbase$0$i = $cmp66$i158 ? $call34$i : -1;
     $tsize$0$i = $cmp66$i158 ? $ssize$0$i : 0;
     label = 191;
    } else $tsize$0748284$i = 0;
   }
  } while (0);
  L291 : do if ((label | 0) == 191) {
   $sub109$i = 0 - $ssize$1$i | 0;
   if (($tbase$0$i | 0) != (-1 | 0)) {
    $tbase$291$i = $tbase$0$i;
    $tsize$290$i = $tsize$0$i;
    label = 202;
    break L269;
   }
   do if (($br$0$i | 0) != (-1 | 0) & $ssize$1$i >>> 0 < 2147483647 & $ssize$1$i >>> 0 < $add$i147 >>> 0) {
    $116 = HEAP32[1122] | 0;
    $and101$i = $sub$i148 - $ssize$1$i + $116 & 0 - $116;
    if ($and101$i >>> 0 < 2147483647) if ((_sbrk($and101$i | 0) | 0) == (-1 | 0)) {
     _sbrk($sub109$i | 0) | 0;
     $tsize$0748284$i = $tsize$0$i;
     break L291;
    } else {
     $ssize$2$i = $and101$i + $ssize$1$i | 0;
     break;
    } else $ssize$2$i = $ssize$1$i;
   } else $ssize$2$i = $ssize$1$i; while (0);
   if (($br$0$i | 0) == (-1 | 0)) $tsize$0748284$i = $tsize$0$i; else {
    $tbase$291$i = $br$0$i;
    $tsize$290$i = $ssize$2$i;
    label = 202;
    break L269;
   }
  } while (0);
  HEAP32[1113] = HEAP32[1113] | 4;
  $tsize$1$i = $tsize$0748284$i;
  label = 199;
 } else {
  $tsize$1$i = 0;
  label = 199;
 } while (0);
 if ((label | 0) == 199) if ($and11$i >>> 0 < 2147483647) {
  $call128$i = _sbrk($and11$i | 0) | 0;
  $call129$i = _sbrk(0) | 0;
  if (($call129$i | 0) != (-1 | 0) & ($call128$i | 0) != (-1 | 0) & $call128$i >>> 0 < $call129$i >>> 0) {
   $sub$ptr$sub$i = $call129$i - $call128$i | 0;
   $cmp138$i164 = $sub$ptr$sub$i >>> 0 > ($nb$0 + 40 | 0) >>> 0;
   if ($cmp138$i164) {
    $tbase$291$i = $call128$i;
    $tsize$290$i = $cmp138$i164 ? $sub$ptr$sub$i : $tsize$1$i;
    label = 202;
   }
  }
 }
 if ((label | 0) == 202) {
  $add147$i = (HEAP32[1110] | 0) + $tsize$290$i | 0;
  HEAP32[1110] = $add147$i;
  if ($add147$i >>> 0 > (HEAP32[1111] | 0) >>> 0) HEAP32[1111] = $add147$i;
  $120 = HEAP32[1008] | 0;
  L311 : do if (($120 | 0) == 0) {
   $121 = HEAP32[1006] | 0;
   if (($121 | 0) == 0 | $tbase$291$i >>> 0 < $121 >>> 0) HEAP32[1006] = $tbase$291$i;
   HEAP32[1114] = $tbase$291$i;
   HEAP32[1115] = $tsize$290$i;
   HEAP32[1117] = 0;
   HEAP32[1011] = HEAP32[1120];
   HEAP32[1010] = -1;
   $i$02$i$i = 0;
   do {
    $shl$i$i = $i$02$i$i << 1;
    $arrayidx$i$i = 4048 + ($shl$i$i << 2) | 0;
    HEAP32[4048 + ($shl$i$i + 3 << 2) >> 2] = $arrayidx$i$i;
    HEAP32[4048 + ($shl$i$i + 2 << 2) >> 2] = $arrayidx$i$i;
    $i$02$i$i = $i$02$i$i + 1 | 0;
   } while (($i$02$i$i | 0) != 32);
   $125 = $tbase$291$i + 8 | 0;
   if (($125 & 7 | 0) == 0) $cond$i$i = 0; else $cond$i$i = 0 - $125 & 7;
   $sub5$i$i = $tsize$290$i + -40 - $cond$i$i | 0;
   HEAP32[1008] = $tbase$291$i + $cond$i$i;
   HEAP32[1005] = $sub5$i$i;
   HEAP32[$tbase$291$i + ($cond$i$i + 4) >> 2] = $sub5$i$i | 1;
   HEAP32[$tbase$291$i + ($tsize$290$i + -36) >> 2] = 40;
   HEAP32[1009] = HEAP32[1124];
  } else {
   $sp$0109$i = 4456 | 0;
   while (1) {
    $128 = HEAP32[$sp$0109$i >> 2] | 0;
    $size185$i = $sp$0109$i + 4 | 0;
    $129 = HEAP32[$size185$i >> 2] | 0;
    if (($tbase$291$i | 0) == ($128 + $129 | 0)) {
     label = 214;
     break;
    }
    $130 = HEAP32[$sp$0109$i + 8 >> 2] | 0;
    if (($130 | 0) == 0) break; else $sp$0109$i = $130;
   }
   if ((label | 0) == 214) if ((HEAP32[$sp$0109$i + 12 >> 2] & 8 | 0) == 0) if ($120 >>> 0 >= $128 >>> 0 & $120 >>> 0 < $tbase$291$i >>> 0) {
    HEAP32[$size185$i >> 2] = $129 + $tsize$290$i;
    $add212$i = (HEAP32[1005] | 0) + $tsize$290$i | 0;
    $133 = $120 + 8 | 0;
    if (($133 & 7 | 0) == 0) $cond$i27$i = 0; else $cond$i27$i = 0 - $133 & 7;
    $sub5$i29$i = $add212$i - $cond$i27$i | 0;
    HEAP32[1008] = $120 + $cond$i27$i;
    HEAP32[1005] = $sub5$i29$i;
    HEAP32[$120 + ($cond$i27$i + 4) >> 2] = $sub5$i29$i | 1;
    HEAP32[$120 + ($add212$i + 4) >> 2] = 40;
    HEAP32[1009] = HEAP32[1124];
    break;
   }
   if ($tbase$291$i >>> 0 < (HEAP32[1006] | 0) >>> 0) HEAP32[1006] = $tbase$291$i;
   $add$ptr224$i = $tbase$291$i + $tsize$290$i | 0;
   $sp$1105$i = 4456 | 0;
   while (1) {
    if ((HEAP32[$sp$1105$i >> 2] | 0) == ($add$ptr224$i | 0)) {
     label = 224;
     break;
    }
    $138 = HEAP32[$sp$1105$i + 8 >> 2] | 0;
    if (($138 | 0) == 0) break; else $sp$1105$i = $138;
   }
   if ((label | 0) == 224) if ((HEAP32[$sp$1105$i + 12 >> 2] & 8 | 0) == 0) {
    HEAP32[$sp$1105$i >> 2] = $tbase$291$i;
    $size242$i = $sp$1105$i + 4 | 0;
    HEAP32[$size242$i >> 2] = (HEAP32[$size242$i >> 2] | 0) + $tsize$290$i;
    $141 = $tbase$291$i + 8 | 0;
    if (($141 & 7 | 0) == 0) $cond$i42$i = 0; else $cond$i42$i = 0 - $141 & 7;
    $143 = $tbase$291$i + ($tsize$290$i + 8) | 0;
    if (($143 & 7 | 0) == 0) $cond15$i$i = 0; else $cond15$i$i = 0 - $143 & 7;
    $add$ptr16$i$i = $tbase$291$i + ($cond15$i$i + $tsize$290$i) | 0;
    $add$ptr4$sum$i49$i = $cond$i42$i + $nb$0 | 0;
    $add$ptr17$i$i = $tbase$291$i + $add$ptr4$sum$i49$i | 0;
    $sub18$i$i = $add$ptr16$i$i - ($tbase$291$i + $cond$i42$i) - $nb$0 | 0;
    HEAP32[$tbase$291$i + ($cond$i42$i + 4) >> 2] = $nb$0 | 3;
    L338 : do if (($add$ptr16$i$i | 0) == (HEAP32[1008] | 0)) {
     $add$i$i = (HEAP32[1005] | 0) + $sub18$i$i | 0;
     HEAP32[1005] = $add$i$i;
     HEAP32[1008] = $add$ptr17$i$i;
     HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 4) >> 2] = $add$i$i | 1;
    } else {
     if (($add$ptr16$i$i | 0) == (HEAP32[1007] | 0)) {
      $add26$i$i = (HEAP32[1004] | 0) + $sub18$i$i | 0;
      HEAP32[1004] = $add26$i$i;
      HEAP32[1007] = $add$ptr17$i$i;
      HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 4) >> 2] = $add26$i$i | 1;
      HEAP32[$tbase$291$i + ($add26$i$i + $add$ptr4$sum$i49$i) >> 2] = $add26$i$i;
      break;
     }
     $add$ptr16$sum$i$i = $tsize$290$i + 4 | 0;
     $149 = HEAP32[$tbase$291$i + ($add$ptr16$sum$i$i + $cond15$i$i) >> 2] | 0;
     if (($149 & 3 | 0) == 1) {
      $and37$i$i = $149 & -8;
      $shr$i54$i = $149 >>> 3;
      L346 : do if ($149 >>> 0 < 256) {
       $150 = HEAP32[$tbase$291$i + (($cond15$i$i | 8) + $tsize$290$i) >> 2] | 0;
       $151 = HEAP32[$tbase$291$i + ($tsize$290$i + 12 + $cond15$i$i) >> 2] | 0;
       $arrayidx$i57$i = 4048 + ($shr$i54$i << 1 << 2) | 0;
       do if (($150 | 0) != ($arrayidx$i57$i | 0)) {
        if ($150 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
        if ((HEAP32[$150 + 12 >> 2] | 0) == ($add$ptr16$i$i | 0)) break;
        _abort();
       } while (0);
       if (($151 | 0) == ($150 | 0)) {
        HEAP32[1002] = HEAP32[1002] & ~(1 << $shr$i54$i);
        break;
       }
       do if (($151 | 0) == ($arrayidx$i57$i | 0)) $fd68$pre$phi$i$iZ2D = $151 + 8 | 0; else {
        if ($151 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
        $fd59$i$i = $151 + 8 | 0;
        if ((HEAP32[$fd59$i$i >> 2] | 0) == ($add$ptr16$i$i | 0)) {
         $fd68$pre$phi$i$iZ2D = $fd59$i$i;
         break;
        }
        _abort();
       } while (0);
       HEAP32[$150 + 12 >> 2] = $151;
       HEAP32[$fd68$pre$phi$i$iZ2D >> 2] = $150;
      } else {
       $157 = HEAP32[$tbase$291$i + (($cond15$i$i | 24) + $tsize$290$i) >> 2] | 0;
       $158 = HEAP32[$tbase$291$i + ($tsize$290$i + 12 + $cond15$i$i) >> 2] | 0;
       do if (($158 | 0) == ($add$ptr16$i$i | 0)) {
        $add$ptr16$sum56$i$i = $cond15$i$i | 16;
        $arrayidx96$i$i = $tbase$291$i + ($add$ptr16$sum$i$i + $add$ptr16$sum56$i$i) | 0;
        $163 = HEAP32[$arrayidx96$i$i >> 2] | 0;
        if (($163 | 0) == 0) {
         $child$i$i = $tbase$291$i + ($add$ptr16$sum56$i$i + $tsize$290$i) | 0;
         $164 = HEAP32[$child$i$i >> 2] | 0;
         if (($164 | 0) == 0) {
          $R$1$i$i = 0;
          break;
         } else {
          $R$0$i$i = $164;
          $RP$0$i$i = $child$i$i;
         }
        } else {
         $R$0$i$i = $163;
         $RP$0$i$i = $arrayidx96$i$i;
        }
        while (1) {
         $arrayidx103$i$i = $R$0$i$i + 20 | 0;
         $165 = HEAP32[$arrayidx103$i$i >> 2] | 0;
         if (($165 | 0) != 0) {
          $R$0$i$i = $165;
          $RP$0$i$i = $arrayidx103$i$i;
          continue;
         }
         $arrayidx107$i$i = $R$0$i$i + 16 | 0;
         $166 = HEAP32[$arrayidx107$i$i >> 2] | 0;
         if (($166 | 0) == 0) break; else {
          $R$0$i$i = $166;
          $RP$0$i$i = $arrayidx107$i$i;
         }
        }
        if ($RP$0$i$i >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
         HEAP32[$RP$0$i$i >> 2] = 0;
         $R$1$i$i = $R$0$i$i;
         break;
        }
       } else {
        $159 = HEAP32[$tbase$291$i + (($cond15$i$i | 8) + $tsize$290$i) >> 2] | 0;
        if ($159 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
        $bk82$i$i = $159 + 12 | 0;
        if ((HEAP32[$bk82$i$i >> 2] | 0) != ($add$ptr16$i$i | 0)) _abort();
        $fd85$i$i = $158 + 8 | 0;
        if ((HEAP32[$fd85$i$i >> 2] | 0) == ($add$ptr16$i$i | 0)) {
         HEAP32[$bk82$i$i >> 2] = $158;
         HEAP32[$fd85$i$i >> 2] = $159;
         $R$1$i$i = $158;
         break;
        } else _abort();
       } while (0);
       if (($157 | 0) == 0) break;
       $168 = HEAP32[$tbase$291$i + ($tsize$290$i + 28 + $cond15$i$i) >> 2] | 0;
       $arrayidx123$i$i = 4312 + ($168 << 2) | 0;
       do if (($add$ptr16$i$i | 0) == (HEAP32[$arrayidx123$i$i >> 2] | 0)) {
        HEAP32[$arrayidx123$i$i >> 2] = $R$1$i$i;
        if (($R$1$i$i | 0) != 0) break;
        HEAP32[1003] = HEAP32[1003] & ~(1 << $168);
        break L346;
       } else {
        if ($157 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
        $arrayidx143$i$i = $157 + 16 | 0;
        if ((HEAP32[$arrayidx143$i$i >> 2] | 0) == ($add$ptr16$i$i | 0)) HEAP32[$arrayidx143$i$i >> 2] = $R$1$i$i; else HEAP32[$157 + 20 >> 2] = $R$1$i$i;
        if (($R$1$i$i | 0) == 0) break L346;
       } while (0);
       if ($R$1$i$i >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
       HEAP32[$R$1$i$i + 24 >> 2] = $157;
       $add$ptr16$sum2627$i$i = $cond15$i$i | 16;
       $174 = HEAP32[$tbase$291$i + ($add$ptr16$sum2627$i$i + $tsize$290$i) >> 2] | 0;
       do if (($174 | 0) != 0) if ($174 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
        HEAP32[$R$1$i$i + 16 >> 2] = $174;
        HEAP32[$174 + 24 >> 2] = $R$1$i$i;
        break;
       } while (0);
       $176 = HEAP32[$tbase$291$i + ($add$ptr16$sum$i$i + $add$ptr16$sum2627$i$i) >> 2] | 0;
       if (($176 | 0) == 0) break;
       if ($176 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
        HEAP32[$R$1$i$i + 20 >> 2] = $176;
        HEAP32[$176 + 24 >> 2] = $R$1$i$i;
        break;
       }
      } while (0);
      $oldfirst$0$i$i = $tbase$291$i + (($and37$i$i | $cond15$i$i) + $tsize$290$i) | 0;
      $qsize$0$i$i = $and37$i$i + $sub18$i$i | 0;
     } else {
      $oldfirst$0$i$i = $add$ptr16$i$i;
      $qsize$0$i$i = $sub18$i$i;
     }
     $head208$i$i = $oldfirst$0$i$i + 4 | 0;
     HEAP32[$head208$i$i >> 2] = HEAP32[$head208$i$i >> 2] & -2;
     HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 4) >> 2] = $qsize$0$i$i | 1;
     HEAP32[$tbase$291$i + ($qsize$0$i$i + $add$ptr4$sum$i49$i) >> 2] = $qsize$0$i$i;
     $shr214$i$i = $qsize$0$i$i >>> 3;
     if ($qsize$0$i$i >>> 0 < 256) {
      $shl221$i$i = $shr214$i$i << 1;
      $arrayidx223$i$i = 4048 + ($shl221$i$i << 2) | 0;
      $179 = HEAP32[1002] | 0;
      $shl226$i$i = 1 << $shr214$i$i;
      do if (($179 & $shl226$i$i | 0) == 0) {
       HEAP32[1002] = $179 | $shl226$i$i;
       $$pre$phi$i67$iZ2D = 4048 + ($shl221$i$i + 2 << 2) | 0;
       $F224$0$i$i = $arrayidx223$i$i;
      } else {
       $180 = 4048 + ($shl221$i$i + 2 << 2) | 0;
       $181 = HEAP32[$180 >> 2] | 0;
       if (!($181 >>> 0 < (HEAP32[1006] | 0) >>> 0)) {
        $$pre$phi$i67$iZ2D = $180;
        $F224$0$i$i = $181;
        break;
       }
       _abort();
      } while (0);
      HEAP32[$$pre$phi$i67$iZ2D >> 2] = $add$ptr17$i$i;
      HEAP32[$F224$0$i$i + 12 >> 2] = $add$ptr17$i$i;
      HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 8) >> 2] = $F224$0$i$i;
      HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 12) >> 2] = $arrayidx223$i$i;
      break;
     }
     $shr253$i$i = $qsize$0$i$i >>> 8;
     do if (($shr253$i$i | 0) == 0) $I252$0$i$i = 0; else {
      if ($qsize$0$i$i >>> 0 > 16777215) {
       $I252$0$i$i = 31;
       break;
      }
      $and264$i$i = ($shr253$i$i + 1048320 | 0) >>> 16 & 8;
      $shl265$i$i = $shr253$i$i << $and264$i$i;
      $and268$i$i = ($shl265$i$i + 520192 | 0) >>> 16 & 4;
      $shl270$i$i = $shl265$i$i << $and268$i$i;
      $and273$i$i = ($shl270$i$i + 245760 | 0) >>> 16 & 2;
      $add278$i$i = 14 - ($and268$i$i | $and264$i$i | $and273$i$i) + ($shl270$i$i << $and273$i$i >>> 15) | 0;
      $I252$0$i$i = $qsize$0$i$i >>> ($add278$i$i + 7 | 0) & 1 | $add278$i$i << 1;
     } while (0);
     $arrayidx287$i$i = 4312 + ($I252$0$i$i << 2) | 0;
     HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 28) >> 2] = $I252$0$i$i;
     HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 20) >> 2] = 0;
     HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 16) >> 2] = 0;
     $183 = HEAP32[1003] | 0;
     $shl294$i$i = 1 << $I252$0$i$i;
     if (($183 & $shl294$i$i | 0) == 0) {
      HEAP32[1003] = $183 | $shl294$i$i;
      HEAP32[$arrayidx287$i$i >> 2] = $add$ptr17$i$i;
      HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 24) >> 2] = $arrayidx287$i$i;
      HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 12) >> 2] = $add$ptr17$i$i;
      HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 8) >> 2] = $add$ptr17$i$i;
      break;
     }
     $184 = HEAP32[$arrayidx287$i$i >> 2] | 0;
     if (($I252$0$i$i | 0) == 31) $cond315$i$i = 0; else $cond315$i$i = 25 - ($I252$0$i$i >>> 1) | 0;
     L435 : do if ((HEAP32[$184 + 4 >> 2] & -8 | 0) == ($qsize$0$i$i | 0)) $T$0$lcssa$i69$i = $184; else {
      $K305$043$i$i = $qsize$0$i$i << $cond315$i$i;
      $T$042$i$i = $184;
      while (1) {
       $arrayidx325$i$i = $T$042$i$i + ($K305$043$i$i >>> 31 << 2) + 16 | 0;
       $186 = HEAP32[$arrayidx325$i$i >> 2] | 0;
       if (($186 | 0) == 0) break;
       if ((HEAP32[$186 + 4 >> 2] & -8 | 0) == ($qsize$0$i$i | 0)) {
        $T$0$lcssa$i69$i = $186;
        break L435;
       } else {
        $K305$043$i$i = $K305$043$i$i << 1;
        $T$042$i$i = $186;
       }
      }
      if ($arrayidx325$i$i >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
       HEAP32[$arrayidx325$i$i >> 2] = $add$ptr17$i$i;
       HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 24) >> 2] = $T$042$i$i;
       HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 12) >> 2] = $add$ptr17$i$i;
       HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 8) >> 2] = $add$ptr17$i$i;
       break L338;
      }
     } while (0);
     $fd344$i$i = $T$0$lcssa$i69$i + 8 | 0;
     $189 = HEAP32[$fd344$i$i >> 2] | 0;
     $190 = HEAP32[1006] | 0;
     if ($T$0$lcssa$i69$i >>> 0 < $190 >>> 0) _abort();
     if ($189 >>> 0 < $190 >>> 0) _abort(); else {
      HEAP32[$189 + 12 >> 2] = $add$ptr17$i$i;
      HEAP32[$fd344$i$i >> 2] = $add$ptr17$i$i;
      HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 8) >> 2] = $189;
      HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 12) >> 2] = $T$0$lcssa$i69$i;
      HEAP32[$tbase$291$i + ($add$ptr4$sum$i49$i + 24) >> 2] = 0;
      break;
     }
    } while (0);
    $mem$0 = $tbase$291$i + ($cond$i42$i | 8) | 0;
    STACKTOP = sp;
    return $mem$0 | 0;
   }
   $sp$0$i$i$i = 4456 | 0;
   while (1) {
    $191 = HEAP32[$sp$0$i$i$i >> 2] | 0;
    if (!($191 >>> 0 > $120 >>> 0)) {
     $192 = HEAP32[$sp$0$i$i$i + 4 >> 2] | 0;
     $add$ptr$i$i$i = $191 + $192 | 0;
     if ($add$ptr$i$i$i >>> 0 > $120 >>> 0) break;
    }
    $sp$0$i$i$i = HEAP32[$sp$0$i$i$i + 8 >> 2] | 0;
   }
   $194 = $191 + ($192 + -39) | 0;
   if (($194 & 7 | 0) == 0) $cond$i17$i = 0; else $cond$i17$i = 0 - $194 & 7;
   $add$ptr7$i$i = $191 + ($192 + -47 + $cond$i17$i) | 0;
   $cond13$i$i = $add$ptr7$i$i >>> 0 < ($120 + 16 | 0) >>> 0 ? $120 : $add$ptr7$i$i;
   $add$ptr14$i$i = $cond13$i$i + 8 | 0;
   $196 = $tbase$291$i + 8 | 0;
   if (($196 & 7 | 0) == 0) $cond$i$i$i = 0; else $cond$i$i$i = 0 - $196 & 7;
   $sub5$i$i$i = $tsize$290$i + -40 - $cond$i$i$i | 0;
   HEAP32[1008] = $tbase$291$i + $cond$i$i$i;
   HEAP32[1005] = $sub5$i$i$i;
   HEAP32[$tbase$291$i + ($cond$i$i$i + 4) >> 2] = $sub5$i$i$i | 1;
   HEAP32[$tbase$291$i + ($tsize$290$i + -36) >> 2] = 40;
   HEAP32[1009] = HEAP32[1124];
   HEAP32[$cond13$i$i + 4 >> 2] = 27;
   HEAP32[$add$ptr14$i$i + 0 >> 2] = HEAP32[1114];
   HEAP32[$add$ptr14$i$i + 4 >> 2] = HEAP32[1115];
   HEAP32[$add$ptr14$i$i + 8 >> 2] = HEAP32[1116];
   HEAP32[$add$ptr14$i$i + 12 >> 2] = HEAP32[1117];
   HEAP32[1114] = $tbase$291$i;
   HEAP32[1115] = $tsize$290$i;
   HEAP32[1117] = 0;
   HEAP32[1116] = $add$ptr14$i$i;
   $add$ptr2418$i$i = $cond13$i$i + 28 | 0;
   HEAP32[$add$ptr2418$i$i >> 2] = 7;
   if (($cond13$i$i + 32 | 0) >>> 0 < $add$ptr$i$i$i >>> 0) {
    $add$ptr2420$i$i = $add$ptr2418$i$i;
    do {
     $add$ptr2420$i$i$looptemp = $add$ptr2420$i$i;
     $add$ptr2420$i$i = $add$ptr2420$i$i + 4 | 0;
     HEAP32[$add$ptr2420$i$i >> 2] = 7;
    } while (($add$ptr2420$i$i$looptemp + 8 | 0) >>> 0 < $add$ptr$i$i$i >>> 0);
   }
   if (($cond13$i$i | 0) != ($120 | 0)) {
    $sub$ptr$sub$i$i = $cond13$i$i - $120 | 0;
    $head31$i$i = $120 + ($sub$ptr$sub$i$i + 4) | 0;
    HEAP32[$head31$i$i >> 2] = HEAP32[$head31$i$i >> 2] & -2;
    HEAP32[$120 + 4 >> 2] = $sub$ptr$sub$i$i | 1;
    HEAP32[$120 + $sub$ptr$sub$i$i >> 2] = $sub$ptr$sub$i$i;
    $shr$i$i = $sub$ptr$sub$i$i >>> 3;
    if ($sub$ptr$sub$i$i >>> 0 < 256) {
     $shl$i20$i = $shr$i$i << 1;
     $arrayidx$i21$i = 4048 + ($shl$i20$i << 2) | 0;
     $203 = HEAP32[1002] | 0;
     $shl39$i$i = 1 << $shr$i$i;
     do if (($203 & $shl39$i$i | 0) == 0) {
      HEAP32[1002] = $203 | $shl39$i$i;
      $$pre$phi$i$iZ2D = 4048 + ($shl$i20$i + 2 << 2) | 0;
      $F$0$i$i = $arrayidx$i21$i;
     } else {
      $204 = 4048 + ($shl$i20$i + 2 << 2) | 0;
      $205 = HEAP32[$204 >> 2] | 0;
      if (!($205 >>> 0 < (HEAP32[1006] | 0) >>> 0)) {
       $$pre$phi$i$iZ2D = $204;
       $F$0$i$i = $205;
       break;
      }
      _abort();
     } while (0);
     HEAP32[$$pre$phi$i$iZ2D >> 2] = $120;
     HEAP32[$F$0$i$i + 12 >> 2] = $120;
     HEAP32[$120 + 8 >> 2] = $F$0$i$i;
     HEAP32[$120 + 12 >> 2] = $arrayidx$i21$i;
     break;
    }
    $shr58$i$i = $sub$ptr$sub$i$i >>> 8;
    if (($shr58$i$i | 0) == 0) $I57$0$i$i = 0; else if ($sub$ptr$sub$i$i >>> 0 > 16777215) $I57$0$i$i = 31; else {
     $and69$i$i = ($shr58$i$i + 1048320 | 0) >>> 16 & 8;
     $shl70$i$i = $shr58$i$i << $and69$i$i;
     $and73$i$i = ($shl70$i$i + 520192 | 0) >>> 16 & 4;
     $shl75$i$i = $shl70$i$i << $and73$i$i;
     $and78$i$i = ($shl75$i$i + 245760 | 0) >>> 16 & 2;
     $add83$i$i = 14 - ($and73$i$i | $and69$i$i | $and78$i$i) + ($shl75$i$i << $and78$i$i >>> 15) | 0;
     $I57$0$i$i = $sub$ptr$sub$i$i >>> ($add83$i$i + 7 | 0) & 1 | $add83$i$i << 1;
    }
    $arrayidx91$i$i = 4312 + ($I57$0$i$i << 2) | 0;
    HEAP32[$120 + 28 >> 2] = $I57$0$i$i;
    HEAP32[$120 + 20 >> 2] = 0;
    HEAP32[$120 + 16 >> 2] = 0;
    $208 = HEAP32[1003] | 0;
    $shl95$i$i = 1 << $I57$0$i$i;
    if (($208 & $shl95$i$i | 0) == 0) {
     HEAP32[1003] = $208 | $shl95$i$i;
     HEAP32[$arrayidx91$i$i >> 2] = $120;
     HEAP32[$120 + 24 >> 2] = $arrayidx91$i$i;
     HEAP32[$120 + 12 >> 2] = $120;
     HEAP32[$120 + 8 >> 2] = $120;
     break;
    }
    $209 = HEAP32[$arrayidx91$i$i >> 2] | 0;
    if (($I57$0$i$i | 0) == 31) $cond115$i$i = 0; else $cond115$i$i = 25 - ($I57$0$i$i >>> 1) | 0;
    L489 : do if ((HEAP32[$209 + 4 >> 2] & -8 | 0) == ($sub$ptr$sub$i$i | 0)) $T$0$lcssa$i$i = $209; else {
     $K105$017$i$i = $sub$ptr$sub$i$i << $cond115$i$i;
     $T$016$i$i = $209;
     while (1) {
      $arrayidx126$i$i = $T$016$i$i + ($K105$017$i$i >>> 31 << 2) + 16 | 0;
      $211 = HEAP32[$arrayidx126$i$i >> 2] | 0;
      if (($211 | 0) == 0) break;
      if ((HEAP32[$211 + 4 >> 2] & -8 | 0) == ($sub$ptr$sub$i$i | 0)) {
       $T$0$lcssa$i$i = $211;
       break L489;
      } else {
       $K105$017$i$i = $K105$017$i$i << 1;
       $T$016$i$i = $211;
      }
     }
     if ($arrayidx126$i$i >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
      HEAP32[$arrayidx126$i$i >> 2] = $120;
      HEAP32[$120 + 24 >> 2] = $T$016$i$i;
      HEAP32[$120 + 12 >> 2] = $120;
      HEAP32[$120 + 8 >> 2] = $120;
      break L311;
     }
    } while (0);
    $fd145$i$i = $T$0$lcssa$i$i + 8 | 0;
    $214 = HEAP32[$fd145$i$i >> 2] | 0;
    $215 = HEAP32[1006] | 0;
    if ($T$0$lcssa$i$i >>> 0 < $215 >>> 0) _abort();
    if ($214 >>> 0 < $215 >>> 0) _abort(); else {
     HEAP32[$214 + 12 >> 2] = $120;
     HEAP32[$fd145$i$i >> 2] = $120;
     HEAP32[$120 + 8 >> 2] = $214;
     HEAP32[$120 + 12 >> 2] = $T$0$lcssa$i$i;
     HEAP32[$120 + 24 >> 2] = 0;
     break;
    }
   }
  } while (0);
  $216 = HEAP32[1005] | 0;
  if ($216 >>> 0 > $nb$0 >>> 0) {
   $sub253$i = $216 - $nb$0 | 0;
   HEAP32[1005] = $sub253$i;
   $217 = HEAP32[1008] | 0;
   HEAP32[1008] = $217 + $nb$0;
   HEAP32[$217 + ($nb$0 + 4) >> 2] = $sub253$i | 1;
   HEAP32[$217 + 4 >> 2] = $nb$0 | 3;
   $mem$0 = $217 + 8 | 0;
   STACKTOP = sp;
   return $mem$0 | 0;
  }
 }
 HEAP32[(___errno_location() | 0) >> 2] = 12;
 $mem$0 = 0;
 STACKTOP = sp;
 return $mem$0 | 0;
}
function _free($mem) {
 $mem = $mem | 0;
 var $$pre$phiZ2D = 0, $0 = 0, $1 = 0, $10 = 0, $11 = 0, $14 = 0, $15 = 0, $16 = 0, $17 = 0, $18 = 0, $2 = 0, $24 = 0, $26 = 0, $30 = 0, $36 = 0, $37 = 0, $4 = 0, $43 = 0, $44 = 0, $45 = 0, $49 = 0, $5 = 0, $50 = 0, $51 = 0, $52 = 0, $54 = 0, $60 = 0, $62 = 0, $65 = 0, $66 = 0, $67 = 0, $70 = 0, $71 = 0, $73 = 0, $76 = 0, $77 = 0, $9 = 0, $F502$0 = 0, $I526$0 = 0, $K575$0270 = 0, $R$0 = 0, $R$1 = 0, $R327$0 = 0, $R327$1 = 0, $RP$0 = 0, $RP355$0 = 0, $T$0$lcssa = 0, $T$0269 = 0, $add$ptr = 0, $add$ptr$sum230 = 0, $add$ptr16 = 0, $add$ptr6 = 0, $add17 = 0, $add243 = 0, $add254 = 0, $add262 = 0, $add551 = 0, $and = 0, $and5 = 0, $and537 = 0, $and541 = 0, $and546 = 0, $arrayidx = 0, $arrayidx108 = 0, $arrayidx113 = 0, $arrayidx130 = 0, $arrayidx149 = 0, $arrayidx274 = 0, $arrayidx357 = 0, $arrayidx369 = 0, $arrayidx374 = 0, $arrayidx395 = 0, $arrayidx414 = 0, $arrayidx501 = 0, $arrayidx559 = 0, $arrayidx591 = 0, $arrayidx99 = 0, $bk338 = 0, $bk82 = 0, $child = 0, $child356 = 0, $cond = 0, $dec = 0, $fd306 = 0, $fd317$pre$phiZ2D = 0, $fd342 = 0, $fd56 = 0, $fd609 = 0, $fd67$pre$phiZ2D = 0, $fd86 = 0, $head209 = 0, $head228 = 0, $p$0 = 0, $psize$0 = 0, $psize$1 = 0, $shl500 = 0, $shl503 = 0, $shl538 = 0, $shl543 = 0, $shl565 = 0, $shr = 0, $shr263 = 0, $shr493 = 0, $shr527 = 0, $sp$0$i = 0, $sp$0$in$i = 0, sp = 0;
 sp = STACKTOP;
 if (($mem | 0) == 0) {
  STACKTOP = sp;
  return;
 }
 $add$ptr = $mem + -8 | 0;
 $0 = HEAP32[1006] | 0;
 if ($add$ptr >>> 0 < $0 >>> 0) _abort();
 $1 = HEAP32[$mem + -4 >> 2] | 0;
 $and = $1 & 3;
 if (($and | 0) == 1) _abort();
 $and5 = $1 & -8;
 $add$ptr6 = $mem + ($and5 + -8) | 0;
 do if (($1 & 1 | 0) == 0) {
  $2 = HEAP32[$add$ptr >> 2] | 0;
  if (($and | 0) == 0) {
   STACKTOP = sp;
   return;
  }
  $add$ptr$sum230 = -8 - $2 | 0;
  $add$ptr16 = $mem + $add$ptr$sum230 | 0;
  $add17 = $2 + $and5 | 0;
  if ($add$ptr16 >>> 0 < $0 >>> 0) _abort();
  if (($add$ptr16 | 0) == (HEAP32[1007] | 0)) {
   $head209 = $mem + ($and5 + -4) | 0;
   if ((HEAP32[$head209 >> 2] & 3 | 0) != 3) {
    $p$0 = $add$ptr16;
    $psize$0 = $add17;
    break;
   }
   HEAP32[1004] = $add17;
   HEAP32[$head209 >> 2] = HEAP32[$head209 >> 2] & -2;
   HEAP32[$mem + ($add$ptr$sum230 + 4) >> 2] = $add17 | 1;
   HEAP32[$add$ptr6 >> 2] = $add17;
   STACKTOP = sp;
   return;
  }
  $shr = $2 >>> 3;
  if ($2 >>> 0 < 256) {
   $4 = HEAP32[$mem + ($add$ptr$sum230 + 8) >> 2] | 0;
   $5 = HEAP32[$mem + ($add$ptr$sum230 + 12) >> 2] | 0;
   $arrayidx = 4048 + ($shr << 1 << 2) | 0;
   if (($4 | 0) != ($arrayidx | 0)) {
    if ($4 >>> 0 < $0 >>> 0) _abort();
    if ((HEAP32[$4 + 12 >> 2] | 0) != ($add$ptr16 | 0)) _abort();
   }
   if (($5 | 0) == ($4 | 0)) {
    HEAP32[1002] = HEAP32[1002] & ~(1 << $shr);
    $p$0 = $add$ptr16;
    $psize$0 = $add17;
    break;
   }
   if (($5 | 0) == ($arrayidx | 0)) $fd67$pre$phiZ2D = $5 + 8 | 0; else {
    if ($5 >>> 0 < $0 >>> 0) _abort();
    $fd56 = $5 + 8 | 0;
    if ((HEAP32[$fd56 >> 2] | 0) == ($add$ptr16 | 0)) $fd67$pre$phiZ2D = $fd56; else _abort();
   }
   HEAP32[$4 + 12 >> 2] = $5;
   HEAP32[$fd67$pre$phiZ2D >> 2] = $4;
   $p$0 = $add$ptr16;
   $psize$0 = $add17;
   break;
  }
  $9 = HEAP32[$mem + ($add$ptr$sum230 + 24) >> 2] | 0;
  $10 = HEAP32[$mem + ($add$ptr$sum230 + 12) >> 2] | 0;
  do if (($10 | 0) == ($add$ptr16 | 0)) {
   $arrayidx99 = $mem + ($add$ptr$sum230 + 20) | 0;
   $14 = HEAP32[$arrayidx99 >> 2] | 0;
   if (($14 | 0) == 0) {
    $child = $mem + ($add$ptr$sum230 + 16) | 0;
    $15 = HEAP32[$child >> 2] | 0;
    if (($15 | 0) == 0) {
     $R$1 = 0;
     break;
    } else {
     $R$0 = $15;
     $RP$0 = $child;
    }
   } else {
    $R$0 = $14;
    $RP$0 = $arrayidx99;
   }
   while (1) {
    $arrayidx108 = $R$0 + 20 | 0;
    $16 = HEAP32[$arrayidx108 >> 2] | 0;
    if (($16 | 0) != 0) {
     $R$0 = $16;
     $RP$0 = $arrayidx108;
     continue;
    }
    $arrayidx113 = $R$0 + 16 | 0;
    $17 = HEAP32[$arrayidx113 >> 2] | 0;
    if (($17 | 0) == 0) break; else {
     $R$0 = $17;
     $RP$0 = $arrayidx113;
    }
   }
   if ($RP$0 >>> 0 < $0 >>> 0) _abort(); else {
    HEAP32[$RP$0 >> 2] = 0;
    $R$1 = $R$0;
    break;
   }
  } else {
   $11 = HEAP32[$mem + ($add$ptr$sum230 + 8) >> 2] | 0;
   if ($11 >>> 0 < $0 >>> 0) _abort();
   $bk82 = $11 + 12 | 0;
   if ((HEAP32[$bk82 >> 2] | 0) != ($add$ptr16 | 0)) _abort();
   $fd86 = $10 + 8 | 0;
   if ((HEAP32[$fd86 >> 2] | 0) == ($add$ptr16 | 0)) {
    HEAP32[$bk82 >> 2] = $10;
    HEAP32[$fd86 >> 2] = $11;
    $R$1 = $10;
    break;
   } else _abort();
  } while (0);
  if (($9 | 0) == 0) {
   $p$0 = $add$ptr16;
   $psize$0 = $add17;
  } else {
   $18 = HEAP32[$mem + ($add$ptr$sum230 + 28) >> 2] | 0;
   $arrayidx130 = 4312 + ($18 << 2) | 0;
   if (($add$ptr16 | 0) == (HEAP32[$arrayidx130 >> 2] | 0)) {
    HEAP32[$arrayidx130 >> 2] = $R$1;
    if (($R$1 | 0) == 0) {
     HEAP32[1003] = HEAP32[1003] & ~(1 << $18);
     $p$0 = $add$ptr16;
     $psize$0 = $add17;
     break;
    }
   } else {
    if ($9 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
    $arrayidx149 = $9 + 16 | 0;
    if ((HEAP32[$arrayidx149 >> 2] | 0) == ($add$ptr16 | 0)) HEAP32[$arrayidx149 >> 2] = $R$1; else HEAP32[$9 + 20 >> 2] = $R$1;
    if (($R$1 | 0) == 0) {
     $p$0 = $add$ptr16;
     $psize$0 = $add17;
     break;
    }
   }
   if ($R$1 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
   HEAP32[$R$1 + 24 >> 2] = $9;
   $24 = HEAP32[$mem + ($add$ptr$sum230 + 16) >> 2] | 0;
   do if (($24 | 0) != 0) if ($24 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
    HEAP32[$R$1 + 16 >> 2] = $24;
    HEAP32[$24 + 24 >> 2] = $R$1;
    break;
   } while (0);
   $26 = HEAP32[$mem + ($add$ptr$sum230 + 20) >> 2] | 0;
   if (($26 | 0) == 0) {
    $p$0 = $add$ptr16;
    $psize$0 = $add17;
   } else if ($26 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
    HEAP32[$R$1 + 20 >> 2] = $26;
    HEAP32[$26 + 24 >> 2] = $R$1;
    $p$0 = $add$ptr16;
    $psize$0 = $add17;
    break;
   }
  }
 } else {
  $p$0 = $add$ptr;
  $psize$0 = $and5;
 } while (0);
 if (!($p$0 >>> 0 < $add$ptr6 >>> 0)) _abort();
 $head228 = $mem + ($and5 + -4) | 0;
 $30 = HEAP32[$head228 >> 2] | 0;
 if (($30 & 1 | 0) == 0) _abort();
 if (($30 & 2 | 0) == 0) {
  if (($add$ptr6 | 0) == (HEAP32[1008] | 0)) {
   $add243 = (HEAP32[1005] | 0) + $psize$0 | 0;
   HEAP32[1005] = $add243;
   HEAP32[1008] = $p$0;
   HEAP32[$p$0 + 4 >> 2] = $add243 | 1;
   if (($p$0 | 0) != (HEAP32[1007] | 0)) {
    STACKTOP = sp;
    return;
   }
   HEAP32[1007] = 0;
   HEAP32[1004] = 0;
   STACKTOP = sp;
   return;
  }
  if (($add$ptr6 | 0) == (HEAP32[1007] | 0)) {
   $add254 = (HEAP32[1004] | 0) + $psize$0 | 0;
   HEAP32[1004] = $add254;
   HEAP32[1007] = $p$0;
   HEAP32[$p$0 + 4 >> 2] = $add254 | 1;
   HEAP32[$p$0 + $add254 >> 2] = $add254;
   STACKTOP = sp;
   return;
  }
  $add262 = ($30 & -8) + $psize$0 | 0;
  $shr263 = $30 >>> 3;
  do if ($30 >>> 0 < 256) {
   $36 = HEAP32[$mem + $and5 >> 2] | 0;
   $37 = HEAP32[$mem + ($and5 | 4) >> 2] | 0;
   $arrayidx274 = 4048 + ($shr263 << 1 << 2) | 0;
   if (($36 | 0) != ($arrayidx274 | 0)) {
    if ($36 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
    if ((HEAP32[$36 + 12 >> 2] | 0) != ($add$ptr6 | 0)) _abort();
   }
   if (($37 | 0) == ($36 | 0)) {
    HEAP32[1002] = HEAP32[1002] & ~(1 << $shr263);
    break;
   }
   if (($37 | 0) == ($arrayidx274 | 0)) $fd317$pre$phiZ2D = $37 + 8 | 0; else {
    if ($37 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
    $fd306 = $37 + 8 | 0;
    if ((HEAP32[$fd306 >> 2] | 0) == ($add$ptr6 | 0)) $fd317$pre$phiZ2D = $fd306; else _abort();
   }
   HEAP32[$36 + 12 >> 2] = $37;
   HEAP32[$fd317$pre$phiZ2D >> 2] = $36;
  } else {
   $43 = HEAP32[$mem + ($and5 + 16) >> 2] | 0;
   $44 = HEAP32[$mem + ($and5 | 4) >> 2] | 0;
   do if (($44 | 0) == ($add$ptr6 | 0)) {
    $arrayidx357 = $mem + ($and5 + 12) | 0;
    $49 = HEAP32[$arrayidx357 >> 2] | 0;
    if (($49 | 0) == 0) {
     $child356 = $mem + ($and5 + 8) | 0;
     $50 = HEAP32[$child356 >> 2] | 0;
     if (($50 | 0) == 0) {
      $R327$1 = 0;
      break;
     } else {
      $R327$0 = $50;
      $RP355$0 = $child356;
     }
    } else {
     $R327$0 = $49;
     $RP355$0 = $arrayidx357;
    }
    while (1) {
     $arrayidx369 = $R327$0 + 20 | 0;
     $51 = HEAP32[$arrayidx369 >> 2] | 0;
     if (($51 | 0) != 0) {
      $R327$0 = $51;
      $RP355$0 = $arrayidx369;
      continue;
     }
     $arrayidx374 = $R327$0 + 16 | 0;
     $52 = HEAP32[$arrayidx374 >> 2] | 0;
     if (($52 | 0) == 0) break; else {
      $R327$0 = $52;
      $RP355$0 = $arrayidx374;
     }
    }
    if ($RP355$0 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
     HEAP32[$RP355$0 >> 2] = 0;
     $R327$1 = $R327$0;
     break;
    }
   } else {
    $45 = HEAP32[$mem + $and5 >> 2] | 0;
    if ($45 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
    $bk338 = $45 + 12 | 0;
    if ((HEAP32[$bk338 >> 2] | 0) != ($add$ptr6 | 0)) _abort();
    $fd342 = $44 + 8 | 0;
    if ((HEAP32[$fd342 >> 2] | 0) == ($add$ptr6 | 0)) {
     HEAP32[$bk338 >> 2] = $44;
     HEAP32[$fd342 >> 2] = $45;
     $R327$1 = $44;
     break;
    } else _abort();
   } while (0);
   if (($43 | 0) != 0) {
    $54 = HEAP32[$mem + ($and5 + 20) >> 2] | 0;
    $arrayidx395 = 4312 + ($54 << 2) | 0;
    if (($add$ptr6 | 0) == (HEAP32[$arrayidx395 >> 2] | 0)) {
     HEAP32[$arrayidx395 >> 2] = $R327$1;
     if (($R327$1 | 0) == 0) {
      HEAP32[1003] = HEAP32[1003] & ~(1 << $54);
      break;
     }
    } else {
     if ($43 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
     $arrayidx414 = $43 + 16 | 0;
     if ((HEAP32[$arrayidx414 >> 2] | 0) == ($add$ptr6 | 0)) HEAP32[$arrayidx414 >> 2] = $R327$1; else HEAP32[$43 + 20 >> 2] = $R327$1;
     if (($R327$1 | 0) == 0) break;
    }
    if ($R327$1 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort();
    HEAP32[$R327$1 + 24 >> 2] = $43;
    $60 = HEAP32[$mem + ($and5 + 8) >> 2] | 0;
    do if (($60 | 0) != 0) if ($60 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
     HEAP32[$R327$1 + 16 >> 2] = $60;
     HEAP32[$60 + 24 >> 2] = $R327$1;
     break;
    } while (0);
    $62 = HEAP32[$mem + ($and5 + 12) >> 2] | 0;
    if (($62 | 0) != 0) if ($62 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
     HEAP32[$R327$1 + 20 >> 2] = $62;
     HEAP32[$62 + 24 >> 2] = $R327$1;
     break;
    }
   }
  } while (0);
  HEAP32[$p$0 + 4 >> 2] = $add262 | 1;
  HEAP32[$p$0 + $add262 >> 2] = $add262;
  if (($p$0 | 0) == (HEAP32[1007] | 0)) {
   HEAP32[1004] = $add262;
   STACKTOP = sp;
   return;
  } else $psize$1 = $add262;
 } else {
  HEAP32[$head228 >> 2] = $30 & -2;
  HEAP32[$p$0 + 4 >> 2] = $psize$0 | 1;
  HEAP32[$p$0 + $psize$0 >> 2] = $psize$0;
  $psize$1 = $psize$0;
 }
 $shr493 = $psize$1 >>> 3;
 if ($psize$1 >>> 0 < 256) {
  $shl500 = $shr493 << 1;
  $arrayidx501 = 4048 + ($shl500 << 2) | 0;
  $65 = HEAP32[1002] | 0;
  $shl503 = 1 << $shr493;
  if (($65 & $shl503 | 0) == 0) {
   HEAP32[1002] = $65 | $shl503;
   $$pre$phiZ2D = 4048 + ($shl500 + 2 << 2) | 0;
   $F502$0 = $arrayidx501;
  } else {
   $66 = 4048 + ($shl500 + 2 << 2) | 0;
   $67 = HEAP32[$66 >> 2] | 0;
   if ($67 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
    $$pre$phiZ2D = $66;
    $F502$0 = $67;
   }
  }
  HEAP32[$$pre$phiZ2D >> 2] = $p$0;
  HEAP32[$F502$0 + 12 >> 2] = $p$0;
  HEAP32[$p$0 + 8 >> 2] = $F502$0;
  HEAP32[$p$0 + 12 >> 2] = $arrayidx501;
  STACKTOP = sp;
  return;
 }
 $shr527 = $psize$1 >>> 8;
 if (($shr527 | 0) == 0) $I526$0 = 0; else if ($psize$1 >>> 0 > 16777215) $I526$0 = 31; else {
  $and537 = ($shr527 + 1048320 | 0) >>> 16 & 8;
  $shl538 = $shr527 << $and537;
  $and541 = ($shl538 + 520192 | 0) >>> 16 & 4;
  $shl543 = $shl538 << $and541;
  $and546 = ($shl543 + 245760 | 0) >>> 16 & 2;
  $add551 = 14 - ($and541 | $and537 | $and546) + ($shl543 << $and546 >>> 15) | 0;
  $I526$0 = $psize$1 >>> ($add551 + 7 | 0) & 1 | $add551 << 1;
 }
 $arrayidx559 = 4312 + ($I526$0 << 2) | 0;
 HEAP32[$p$0 + 28 >> 2] = $I526$0;
 HEAP32[$p$0 + 20 >> 2] = 0;
 HEAP32[$p$0 + 16 >> 2] = 0;
 $70 = HEAP32[1003] | 0;
 $shl565 = 1 << $I526$0;
 L199 : do if (($70 & $shl565 | 0) == 0) {
  HEAP32[1003] = $70 | $shl565;
  HEAP32[$arrayidx559 >> 2] = $p$0;
  HEAP32[$p$0 + 24 >> 2] = $arrayidx559;
  HEAP32[$p$0 + 12 >> 2] = $p$0;
  HEAP32[$p$0 + 8 >> 2] = $p$0;
 } else {
  $71 = HEAP32[$arrayidx559 >> 2] | 0;
  if (($I526$0 | 0) == 31) $cond = 0; else $cond = 25 - ($I526$0 >>> 1) | 0;
  L205 : do if ((HEAP32[$71 + 4 >> 2] & -8 | 0) == ($psize$1 | 0)) $T$0$lcssa = $71; else {
   $K575$0270 = $psize$1 << $cond;
   $T$0269 = $71;
   while (1) {
    $arrayidx591 = $T$0269 + ($K575$0270 >>> 31 << 2) + 16 | 0;
    $73 = HEAP32[$arrayidx591 >> 2] | 0;
    if (($73 | 0) == 0) break;
    if ((HEAP32[$73 + 4 >> 2] & -8 | 0) == ($psize$1 | 0)) {
     $T$0$lcssa = $73;
     break L205;
    } else {
     $K575$0270 = $K575$0270 << 1;
     $T$0269 = $73;
    }
   }
   if ($arrayidx591 >>> 0 < (HEAP32[1006] | 0) >>> 0) _abort(); else {
    HEAP32[$arrayidx591 >> 2] = $p$0;
    HEAP32[$p$0 + 24 >> 2] = $T$0269;
    HEAP32[$p$0 + 12 >> 2] = $p$0;
    HEAP32[$p$0 + 8 >> 2] = $p$0;
    break L199;
   }
  } while (0);
  $fd609 = $T$0$lcssa + 8 | 0;
  $76 = HEAP32[$fd609 >> 2] | 0;
  $77 = HEAP32[1006] | 0;
  if ($T$0$lcssa >>> 0 < $77 >>> 0) _abort();
  if ($76 >>> 0 < $77 >>> 0) _abort(); else {
   HEAP32[$76 + 12 >> 2] = $p$0;
   HEAP32[$fd609 >> 2] = $p$0;
   HEAP32[$p$0 + 8 >> 2] = $76;
   HEAP32[$p$0 + 12 >> 2] = $T$0$lcssa;
   HEAP32[$p$0 + 24 >> 2] = 0;
   break;
  }
 } while (0);
 $dec = (HEAP32[1010] | 0) + -1 | 0;
 HEAP32[1010] = $dec;
 if (($dec | 0) == 0) $sp$0$in$i = 4464 | 0; else {
  STACKTOP = sp;
  return;
 }
 while (1) {
  $sp$0$i = HEAP32[$sp$0$in$i >> 2] | 0;
  if (($sp$0$i | 0) == 0) break; else $sp$0$in$i = $sp$0$i + 8 | 0;
 }
 HEAP32[1010] = -1;
 STACKTOP = sp;
 return;
}
function _memcpy(dest, src, num) {
 dest = dest | 0;
 src = src | 0;
 num = num | 0;
 var ret = 0;
 if ((num | 0) >= 4096) return _emscripten_memcpy_big(dest | 0, src | 0, num | 0) | 0;
 ret = dest | 0;
 if ((dest & 3) == (src & 3)) {
  while (dest & 3) {
   if ((num | 0) == 0) return ret | 0;
   HEAP8[dest >> 0] = HEAP8[src >> 0] | 0;
   dest = dest + 1 | 0;
   src = src + 1 | 0;
   num = num - 1 | 0;
  }
  while ((num | 0) >= 4) {
   HEAP32[dest >> 2] = HEAP32[src >> 2];
   dest = dest + 4 | 0;
   src = src + 4 | 0;
   num = num - 4 | 0;
  }
 }
 while ((num | 0) > 0) {
  HEAP8[dest >> 0] = HEAP8[src >> 0] | 0;
  dest = dest + 1 | 0;
  src = src + 1 | 0;
  num = num - 1 | 0;
 }
 return ret | 0;
}
function runPostSets() {}
function _memset(ptr, value, num) {
 ptr = ptr | 0;
 value = value | 0;
 num = num | 0;
 var stop = 0, value4 = 0, stop4 = 0, unaligned = 0;
 stop = ptr + num | 0;
 if ((num | 0) >= 20) {
  value = value & 255;
  unaligned = ptr & 3;
  value4 = value | value << 8 | value << 16 | value << 24;
  stop4 = stop & ~3;
  if (unaligned) {
   unaligned = ptr + 4 - unaligned | 0;
   while ((ptr | 0) < (unaligned | 0)) {
    HEAP8[ptr >> 0] = value;
    ptr = ptr + 1 | 0;
   }
  }
  while ((ptr | 0) < (stop4 | 0)) {
   HEAP32[ptr >> 2] = value4;
   ptr = ptr + 4 | 0;
  }
 }
 while ((ptr | 0) < (stop | 0)) {
  HEAP8[ptr >> 0] = value;
  ptr = ptr + 1 | 0;
 }
 return ptr - num | 0;
}
function _main() {
 var $ax4$09$i = 0, $j$010$i = 0, $sumx4$011$i = SIMD_float32x4(0, 0, 0, 0), sp = 0;
 sp = STACKTOP;
 $ax4$09$i = 8;
 $j$010$i = 0;
 $sumx4$011$i = SIMD_float32x4_splat(Math_fround(0));
 while (1) {
  $sumx4$011$i = SIMD_float32x4_add($sumx4$011$i, SIMD_float32x4_load(buffer, $ax4$09$i));
  $j$010$i = $j$010$i + 4 | 0;
  if (!($j$010$i >>> 0 < 1e3)) break; else $ax4$09$i = $ax4$09$i + 16 | 0;
 }
 STACKTOP = sp;
 return ~~((+$sumx4$011$i.w + (+$sumx4$011$i.z + (+$sumx4$011$i.x + +$sumx4$011$i.y))) / 1.0e3) | 0;
}
function copyTempDouble(ptr) {
 ptr = ptr | 0;
 HEAP8[tempDoublePtr >> 0] = HEAP8[ptr >> 0];
 HEAP8[tempDoublePtr + 1 >> 0] = HEAP8[ptr + 1 >> 0];
 HEAP8[tempDoublePtr + 2 >> 0] = HEAP8[ptr + 2 >> 0];
 HEAP8[tempDoublePtr + 3 >> 0] = HEAP8[ptr + 3 >> 0];
 HEAP8[tempDoublePtr + 4 >> 0] = HEAP8[ptr + 4 >> 0];
 HEAP8[tempDoublePtr + 5 >> 0] = HEAP8[ptr + 5 >> 0];
 HEAP8[tempDoublePtr + 6 >> 0] = HEAP8[ptr + 6 >> 0];
 HEAP8[tempDoublePtr + 7 >> 0] = HEAP8[ptr + 7 >> 0];
}
function copyTempFloat(ptr) {
 ptr = ptr | 0;
 HEAP8[tempDoublePtr >> 0] = HEAP8[ptr >> 0];
 HEAP8[tempDoublePtr + 1 >> 0] = HEAP8[ptr + 1 >> 0];
 HEAP8[tempDoublePtr + 2 >> 0] = HEAP8[ptr + 2 >> 0];
 HEAP8[tempDoublePtr + 3 >> 0] = HEAP8[ptr + 3 >> 0];
}
function stackAlloc(size) {
 size = size | 0;
 var ret = 0;
 ret = STACKTOP;
 STACKTOP = STACKTOP + size | 0;
 STACKTOP = STACKTOP + 15 & -16;
 return ret | 0;
}
function setThrew(threw, value) {
 threw = threw | 0;
 value = value | 0;
 if ((__THREW__ | 0) == 0) {
  __THREW__ = threw;
  threwValue = value;
 }
}
function _strlen(ptr) {
 ptr = ptr | 0;
 var curr = 0;
 curr = ptr;
 while (HEAP8[curr >> 0] | 0) curr = curr + 1 | 0;
 return curr - ptr | 0;
}
function setTempRet0(value) {
 value = value | 0;
 tempRet0 = value;
}
function stackRestore(top) {
 top = top | 0;
 STACKTOP = top;
}
function getTempRet0() {
 return tempRet0 | 0;
}
function stackSave() {
 return STACKTOP | 0;
}

// EMSCRIPTEN_END_FUNCS
  

    return { _strlen: _strlen, _free: _free, _main: _main, _memset: _memset, _malloc: _malloc, _memcpy: _memcpy, runPostSets: runPostSets, stackAlloc: stackAlloc, stackSave: stackSave, stackRestore: stackRestore, setThrew: setThrew, setTempRet0: setTempRet0, getTempRet0: getTempRet0 };
  })
  // EMSCRIPTEN_END_ASM
  ({ "Math": Math, "Int8Array": Int8Array, "Int16Array": Int16Array, "Int32Array": Int32Array, "Uint8Array": Uint8Array, "Uint16Array": Uint16Array, "Uint32Array": Uint32Array, "Float32Array": Float32Array, "Float64Array": Float64Array, "SIMD": SIMD }, { "abort": abort, "assert": assert, "min": Math_min, "_fflush": _fflush, "_abort": _abort, "___setErrNo": ___setErrNo, "_sbrk": _sbrk, "_time": _time, "_emscripten_memcpy_big": _emscripten_memcpy_big, "_sysconf": _sysconf, "___errno_location": ___errno_location, "STACKTOP": STACKTOP, "STACK_MAX": STACK_MAX, "tempDoublePtr": tempDoublePtr, "ABORT": ABORT, "NaN": NaN, "Infinity": Infinity }, buffer);
  var _strlen = Module["_strlen"] = asm["_strlen"];
var _free = Module["_free"] = asm["_free"];
var _main = Module["_main"] = asm["_main"];
var _memset = Module["_memset"] = asm["_memset"];
var _malloc = Module["_malloc"] = asm["_malloc"];
var _memcpy = Module["_memcpy"] = asm["_memcpy"];
var runPostSets = Module["runPostSets"] = asm["runPostSets"];
  
  Runtime.stackAlloc = asm['stackAlloc'];
  Runtime.stackSave = asm['stackSave'];
  Runtime.stackRestore = asm['stackRestore'];
  Runtime.setTempRet0 = asm['setTempRet0'];
  Runtime.getTempRet0 = asm['getTempRet0'];
  

// Warning: printing of i64 values may be slightly rounded! No deep i64 math used, so precise i64 code not included
var i64Math = null;

// === Auto-generated postamble setup entry stuff ===

if (memoryInitializer) {
  if (typeof Module['locateFile'] === 'function') {
    memoryInitializer = Module['locateFile'](memoryInitializer);
  } else if (Module['memoryInitializerPrefixURL']) {
    memoryInitializer = Module['memoryInitializerPrefixURL'] + memoryInitializer;
  }
  if (ENVIRONMENT_IS_NODE || ENVIRONMENT_IS_SHELL) {
    var data = Module['readBinary'](memoryInitializer);
    HEAPU8.set(data, STATIC_BASE);
  } else {
    addRunDependency('memory initializer');
    Browser.asyncLoad(memoryInitializer, function(data) {
      HEAPU8.set(data, STATIC_BASE);
      removeRunDependency('memory initializer');
    }, function(data) {
      throw 'could not load memory initializer ' + memoryInitializer;
    });
  }
}

function ExitStatus(status) {
  this.name = "ExitStatus";
  this.message = "Program terminated with exit(" + status + ")";
  this.status = status;
};
ExitStatus.prototype = new Error();
ExitStatus.prototype.constructor = ExitStatus;

var initialStackTop;
var preloadStartTime = null;
var calledMain = false;

dependenciesFulfilled = function runCaller() {
  // If run has never been called, and we should call run (INVOKE_RUN is true, and Module.noInitialRun is not false)
  if (!Module['calledRun'] && shouldRunNow) run();
  if (!Module['calledRun']) dependenciesFulfilled = runCaller; // try this again later, after new deps are fulfilled
}

Module['callMain'] = Module.callMain = function callMain(args) {
  assert(runDependencies == 0, 'cannot call main when async dependencies remain! (listen on __ATMAIN__)');
  assert(__ATPRERUN__.length == 0, 'cannot call main when preRun functions remain to be called');

  args = args || [];

  ensureInitRuntime();

  var argc = args.length+1;
  function pad() {
    for (var i = 0; i < 4-1; i++) {
      argv.push(0);
    }
  }
  var argv = [allocate(intArrayFromString(Module['thisProgram']), 'i8', ALLOC_NORMAL) ];
  pad();
  for (var i = 0; i < argc-1; i = i + 1) {
    argv.push(allocate(intArrayFromString(args[i]), 'i8', ALLOC_NORMAL));
    pad();
  }
  argv.push(0);
  argv = allocate(argv, 'i32', ALLOC_NORMAL);

  initialStackTop = STACKTOP;

  try {

    var ret = Module['_main'](argc, argv, 0);


    // if we're not running an evented main loop, it's time to exit
    exit(ret);
  }
  catch(e) {
    if (e instanceof ExitStatus) {
      // exit() throws this once it's done to make sure execution
      // has been stopped completely
      return;
    } else if (e == 'SimulateInfiniteLoop') {
      // running an evented main loop, don't immediately exit
      Module['noExitRuntime'] = true;
      return;
    } else {
      if (e && typeof e === 'object' && e.stack) Module.printErr('exception thrown: ' + [e, e.stack]);
      throw e;
    }
  } finally {
    calledMain = true;
  }
}




function run(args) {
  args = args || Module['arguments'];

  if (preloadStartTime === null) preloadStartTime = Date.now();

  if (runDependencies > 0) {
    Module.printErr('run() called, but dependencies remain, so not running');
    return;
  }

  preRun();

  if (runDependencies > 0) return; // a preRun added a dependency, run will be called later
  if (Module['calledRun']) return; // run may have just been called through dependencies being fulfilled just in this very frame

  function doRun() {
    if (Module['calledRun']) return; // run may have just been called while the async setStatus time below was happening
    Module['calledRun'] = true;

    if (ABORT) return; 

    ensureInitRuntime();

    preMain();

    if (ENVIRONMENT_IS_WEB && preloadStartTime !== null) {
      Module.printErr('pre-main prep time: ' + (Date.now() - preloadStartTime) + ' ms');
    }

    if (Module['_main'] && shouldRunNow) {
      Module['callMain'](args);
    }

    postRun();
  }

  if (Module['setStatus']) {
    Module['setStatus']('Running...');
    setTimeout(function() {
      setTimeout(function() {
        Module['setStatus']('');
      }, 1);
      doRun();
    }, 1);
  } else {
    doRun();
  }
}
Module['run'] = Module.run = run;

function exit(status) {
  if (Module['noExitRuntime']) {
    return;
  }

  ABORT = true;
  EXITSTATUS = status;
  STACKTOP = initialStackTop;

  // exit the runtime
  exitRuntime();

  if (ENVIRONMENT_IS_NODE) {
    // Work around a node.js bug where stdout buffer is not flushed at process exit:
    // Instead of process.exit() directly, wait for stdout flush event.
    // See https://github.com/joyent/node/issues/1669 and https://github.com/kripken/emscripten/issues/2582
    // Workaround is based on https://github.com/RReverser/acorn/commit/50ab143cecc9ed71a2d66f78b4aec3bb2e9844f6
    process['stdout']['once']('drain', function () {
      process['exit'](status);
    });
    console.log(' '); // Make sure to print something to force the drain event to occur, in case the stdout buffer was empty.
    // Work around another node bug where sometimes 'drain' is never fired - make another effort
    // to emit the exit status, after a significant delay (if node hasn't fired drain by then, give up)
    setTimeout(function() {
      process['exit'](status);
    }, 500);
  } else
  if (ENVIRONMENT_IS_SHELL && typeof quit === 'function') {
    quit(status);
  }
  // if we reach here, we must throw an exception to halt the current execution
  throw new ExitStatus(status);
}
Module['exit'] = Module.exit = exit;

function abort(text) {
  if (text) {
    Module.print(text);
    Module.printErr(text);
  }

  ABORT = true;
  EXITSTATUS = 1;

  var extra = '\nIf this abort() is unexpected, build with -s ASSERTIONS=1 which can give more information.';

  throw 'abort() at ' + stackTrace() + extra;
}
Module['abort'] = Module.abort = abort;

// {{PRE_RUN_ADDITIONS}}

if (Module['preInit']) {
  if (typeof Module['preInit'] == 'function') Module['preInit'] = [Module['preInit']];
  while (Module['preInit'].length > 0) {
    Module['preInit'].pop()();
  }
}

// shouldRunNow refers to calling main(), not run().
var shouldRunNow = true;
if (Module['noInitialRun']) {
  shouldRunNow = false;
}


run();

// {{POST_RUN_ADDITIONS}}






// {{MODULE_ADDITIONS}}






