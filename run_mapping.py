import h5reader
import map1d

y = h5reader.hydrof(model='s2.8si_128_hires_new2')
map1d.map2d(y, y.xzl(), y.xzr(), y.yzl(), y.yzr(), y.zzl(), y.zzr(), 0, y.xzn()[-1], -y.xzn()[-1], y.xzn()[-1], 100, 1, 1, 1e-3)
