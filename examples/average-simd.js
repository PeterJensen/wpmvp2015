function average(data) {
  var sum = SIMD.float32x4.splat(0.0);
  for (var i = 0, l = data.length; i < l; i = i+4) {
    sum = SIMD.float32x4.add(
      sum, SIMD.float32x4.load(data, i));
  }
  var total = sum.x + sum.y + sum.z + sum.w;
  return total/data.length;
}

var a = new Float32Array(100000);

function init(data) {
  for (var i = 0, l = data.length; i < l; i++) {
    data[i] = 1.0;
  }
}

function main() {
  init(a);
  var s = 0.0;
  for (var j = 0; j < 50000; j++) {
    s += average(a);
  }
  print("s = " + s);
}

main();
