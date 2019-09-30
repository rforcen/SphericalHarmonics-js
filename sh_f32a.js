//
// Spherical Harmonics class in Float32Array
//

const {
  sin,
  cos,
  pow,
  max,
  // sqrt,
  isinf,
  PI
} = Math


const size_of_location = 11, // coord(3)+normal(3)+color(3)+calc_texture(2)
  loc_offset = [3, 3, 3, 2],
  loc_names = ['coord', 'normal', 'color', 'calc_texture'],
  TWOPI = PI * 2

class SphericalHarmonics {
  constructor(resolution = 256, colorMap = 7, code_index = 0) {
    this.colorMap = colorMap
    this.resolution = resolution
    this.res2 = this.resolution * this.resolution
    this.n_vertex = this.res2
    this.m = []

    // mesh sotrage
    const nc = this.res2 * 3,
      nt = this.res2 * 2
    this.mesh = {
      coords: new Float32Array(nc),
      colors: new Float32Array(nc),
      normals: new Float32Array(nc),
      textures: new Float32Array(nt)
    }
    this.max_val = -1e10

    this.read_code(code_index)
    this.calc_mesh()
  }

  read_code(code_index) {
    this.code = sh_codes[code_index % sh_codes.length]
    this.calc_coeffs()
  }

  calc_coeffs() { // in 'm'
    var c = this.code.toString()
    c = '0'.repeat(8 - c.length) + c
    for (let i = 0; i < 8; i++) this.m[i] = parseInt(c.charAt(i))
  }

  powint(x, y) { // x ^ int, y is m[] so int in 1..8 range
    switch (y) {
      case 0:
        return 1
      case 1:
        return x
      case 2:
        return x * x
      case 3:
        return x * x * x
      case 4:
        return x * x * x * x
      case 5:
        return x * x * x * x * x
      case 6:
        return x * x * x * x * x * x
      case 7:
        return x * x * x * x * x * x * x
      case 8:
        return x * x * x * x * x * x * x * x
      default:
        for (let i = 1; i < y; i++) x *= x
        return x
    }
  }

  powfilt(x, y) { // general filtered power
    if (y == 0) return 1
    p = pow(x, y)
    return (isinf(p)) ? 0 : p
  }

  calc_coord(theta, phi) {
    const m = this.m
    var r = this.powint(sin(m[0] * phi), m[1])
    r += this.powint(cos(m[2] * phi), m[3])
    r += this.powint(sin(m[4] * theta), m[5])
    r += this.powint(cos(m[6] * theta), m[7])

    return [r * sin(phi) * cos(theta), r * cos(phi), r * sin(phi) * sin(theta)]
  }

  calc_texture(u, v) {
    return [u, v]
  }

  calc_normals(p0, p1, p2) {
    return normalize(cross(sub(p1, p2), sub(p1, p0)))
  }

  calc_color(v, vmin, vmax, type) {
    let dv = 0,
      vmid = 0
    let c = [1.0, 1.0, 1.0]
    let c1 = [],
      c2 = [],
      c3 = []
    let ratio = 0

    if (vmax < vmin) {
      dv = vmin
      vmin = vmax
      vmax = dv
    }
    if (vmax - vmin < 0.000001) {
      vmin -= 1
      vmax += 1
    }

    if (v < vmin) v = vmin
    if (v > vmax) v = vmax
    dv = vmax - vmin

    switch (type) {
      case 1:
        if (v < (vmin + 0.25 * dv)) {
          c[0] = 0
          c[1] = 4 * (v - vmin) / dv
          c[2] = 1
        } else if (v < (vmin + 0.5 * dv)) {
          c[0] = 0
          c[1] = 1
          c[2] = 1 + 4 * (vmin + 0.25 * dv - v) / dv
        } else if (v < (vmin + 0.75 * dv)) {
          c[0] = 4 * (v - vmin - 0.5 * dv) / dv
          c[1] = 1
          c[2] = 0
        } else {
          c[0] = 1
          c[1] = 1 + 4 * (vmin + 0.75 * dv - v) / dv
          c[2] = 0
        }
        break
      case 2:
        c[0] = (v - vmin) / dv
        c[1] = 0
        c[2] = (vmax - v) / dv
        break
      case 3:
        c[0] = (v - vmin) / dv
        c[2] = c[0]
        c[1] = c[0]
        break
      case 4:
        if (v < (vmin + dv / 6.0)) {
          c[0] = 1
          c[1] = 6 * (v - vmin) / dv
          c[2] = 0
        } else if (v < (vmin + 2.0 * dv / 6.0)) {
          c[0] = 1 + 6 * (vmin + dv / 6.0 - v) / dv
          c[1] = 1
          c[2] = 0
        } else if (v < (vmin + 3.0 * dv / 6.0)) {
          c[0] = 0
          c[1] = 1
          c[2] = 6 * (v - vmin - 2.0 * dv / 6.0) / dv
        } else if (v < (vmin + 4.0 * dv / 6.0)) {
          c[0] = 0
          c[1] = 1 + 6 * (vmin + 3.0 * dv / 6.0 - v) / dv
          c[2] = 1
        } else if (v < (vmin + 5.0 * dv / 6.0)) {
          c[0] = 6 * (v - vmin - 4.0 * dv / 6.0) / dv
          c[1] = 0
          c[2] = 1
        } else {
          c[0] = 1
          c[1] = 0
          c[2] = 1 + 6 * (vmin + 5.0 * dv / 6.0 - v) / dv
        }
        break
      case 5:
        c[0] = (v - vmin) / (vmax - vmin)
        c[1] = 1
        c[2] = 0
        break
      case 6:
        c[0] = (v - vmin) / (vmax - vmin)
        c[1] = (vmax - v) / (vmax - vmin)
        c[2] = c[0]
        break
      case 7:
        if (v < (vmin + 0.25 * dv)) {
          c[0] = 0
          c[1] = 4 * (v - vmin) / dv
          c[2] = 1 - c[1]
        } else if (v < (vmin + 0.5 * dv)) {
          c[0] = 4 * (v - vmin - 0.25 * dv) / dv
          c[1] = 1 - c[0]
          c[2] = 0
        } else if (v < (vmin + 0.75 * dv)) {
          c[1] = 4 * (v - vmin - 0.5 * dv) / dv
          c[0] = 1 - c[1]
          c[2] = 0
        } else {
          c[0] = 0
          c[2] = 4 * (v - vmin - 0.75 * dv) / dv
          c[1] = 1 - c[2]
        }
        break
      case 8:
        if (v < (vmin + 0.5 * dv)) {
          c[0] = 2 * (v - vmin) / dv
          c[1] = c[0]
          c[2] = c[0]
        } else {
          c[0] = 1 - 2 * (v - vmin - 0.5 * dv) / dv
          c[1] = c[0]
          c[2] = c[0]
        }
        break
      case 9:
        if (v < (vmin + dv / 3)) {
          c[2] = 3 * (v - vmin) / dv
          c[1] = 0
          c[0] = 1 - c[2]
        } else if (v < (vmin + 2 * dv / 3)) {
          c[0] = 0
          c[1] = 3 * (v - vmin - dv / 3) / dv
          c[2] = 1
        } else {
          c[0] = 3 * (v - vmin - 2 * dv / 3) / dv
          c[1] = 1 - c[0]
          c[2] = 1
        }
        break
      case 10:
        if (v < (vmin + 0.2 * dv)) {
          c[0] = 0
          c[1] = 5 * (v - vmin) / dv
          c[2] = 1
        } else if (v < (vmin + 0.4 * dv)) {
          c[0] = 0
          c[1] = 1
          c[2] = 1 + 5 * (vmin + 0.2 * dv - v) / dv
        } else if (v < (vmin + 0.6 * dv)) {
          c[0] = 5 * (v - vmin - 0.4 * dv) / dv
          c[1] = 1
          c[2] = 0
        } else if (v < (vmin + 0.8 * dv)) {
          c[0] = 1
          c[1] = 1 - 5 * (v - vmin - 0.6 * dv) / dv
          c[2] = 0
        } else {
          c[0] = 1
          c[1] = 5 * (v - vmin - 0.8 * dv) / dv
          c[2] = 5 * (v - vmin - 0.8 * dv) / dv
        }
        break
      case 11:
        c1.r = 200 / 255.0
        c1.g = 60 / 255.0
        c1.b = 0 / 255.0
        c2.r = 250 / 255.0
        c2.g = 160 / 255.0
        c2.b = 110 / 255.0
        c[0] = (c2.r - c1.r) * (v - vmin) / dv + c1.r
        c[1] = (c2.g - c1.g) * (v - vmin) / dv + c1.g
        c[2] = (c2.b - c1.b) * (v - vmin) / dv + c1.b
        break
      case 12:
        c1.r = 55 / 255.0
        c1.g = 55 / 255.0
        c1.b = 45 / 255.0
        /* c2.r = 200 / 255.0 c2.g =  60 / 255.0 c2.b =   0 / 255.0 */
        c2.r = 235 / 255.0
        c2.g = 90 / 255.0
        c2.b = 30 / 255.0
        c3.r = 250 / 255.0
        c3.g = 160 / 255.0
        c3.b = 110 / 255.0
        ratio = 0.4
        vmid = vmin + ratio * dv
        if (v < vmid) {
          c[0] = (c2.r - c1.r) * (v - vmin) / (ratio * dv) + c1.r
          c[1] = (c2.g - c1.g) * (v - vmin) / (ratio * dv) + c1.g
          c[2] = (c2.b - c1.b) * (v - vmin) / (ratio * dv) + c1.b
        } else {
          c[0] = (c3.r - c2.r) * (v - vmid) / ((1 - ratio) * dv) + c2.r
          c[1] = (c3.g - c2.g) * (v - vmid) / ((1 - ratio) * dv) + c2.g
          c[2] = (c3.b - c2.b) * (v - vmid) / ((1 - ratio) * dv) + c2.b
        }
        break
      case 13:
        c1.r = 0 / 255.0
        c1.g = 255 / 255.0
        c1.b = 0 / 255.0
        c2.r = 255 / 255.0
        c2.g = 150 / 255.0
        c2.b = 0 / 255.0
        c3.r = 255 / 255.0
        c3.g = 250 / 255.0
        c3.b = 240 / 255.0
        ratio = 0.3
        vmid = vmin + ratio * dv
        if (v < vmid) {
          c[0] = (c2.r - c1.r) * (v - vmin) / (ratio * dv) + c1.r
          c[1] = (c2.g - c1.g) * (v - vmin) / (ratio * dv) + c1.g
          c[2] = (c2.b - c1.b) * (v - vmin) / (ratio * dv) + c1.b
        } else {
          c[0] = (c3.r - c2.r) * (v - vmid) / ((1 - ratio) * dv) + c2.r
          c[1] = (c3.g - c2.g) * (v - vmid) / ((1 - ratio) * dv) + c2.g
          c[2] = (c3.b - c2.b) * (v - vmid) / ((1 - ratio) * dv) + c2.b
        }
        break
      case 14:
        c[0] = 1
        c[1] = 1 - (v - vmin) / dv
        c[2] = 0
        break
      case 15:
        if (v < (vmin + 0.25 * dv)) {
          c[0] = 0
          c[1] = 4 * (v - vmin) / dv
          c[2] = 1
        } else if (v < (vmin + 0.5 * dv)) {
          c[0] = 0
          c[1] = 1
          c[2] = 1 - 4 * (v - vmin - 0.25 * dv) / dv
        } else if (v < (vmin + 0.75 * dv)) {
          c[0] = 4 * (v - vmin - 0.5 * dv) / dv
          c[1] = 1
          c[2] = 0
        } else {
          c[0] = 1
          c[1] = 1
          c[2] = 4 * (v - vmin - 0.75 * dv) / dv
        }
        break
      case 16:
        if (v < (vmin + 0.5 * dv)) {
          c[0] = 0.0
          c[1] = 2 * (v - vmin) / dv
          c[2] = 1 - 2 * (v - vmin) / dv
        } else {
          c[0] = 2 * (v - vmin - 0.5 * dv) / dv
          c[1] = 1 - 2 * (v - vmin - 0.5 * dv) / dv
          c[2] = 0.0
        }
        break
      case 17:
        if (v < (vmin + 0.5 * dv)) {
          c[0] = 1.0
          c[1] = 1 - 2 * (v - vmin) / dv
          c[2] = 2 * (v - vmin) / dv
        } else {
          c[0] = 1 - 2 * (v - vmin - 0.5 * dv) / dv
          c[1] = 2 * (v - vmin - 0.5 * dv) / dv
          c[2] = 1.0
        }
        break
      case 18:
        c[0] = 0
        c[1] = (v - vmin) / (vmax - vmin)
        c[2] = 1
        break
      case 19:
        c[0] = (v - vmin) / (vmax - vmin)
        c[1] = c[0]
        c[2] = 1
        break
      case 20:
        c1.r = 0 / 255.0
        c1.g = 160 / 255.0
        c1.b = 0 / 255.0
        c2.r = 180 / 255.0
        c2.g = 220 / 255.0
        c2.b = 0 / 255.0
        c3.r = 250 / 255.0
        c3.g = 220 / 255.0
        c3.b = 170 / 255.0
        ratio = 0.3
        vmid = vmin + ratio * dv
        if (v < vmid) {
          c[0] = (c2.r - c1.r) * (v - vmin) / (ratio * dv) + c1.r
          c[1] = (c2.g - c1.g) * (v - vmin) / (ratio * dv) + c1.g
          c[2] = (c2.b - c1.b) * (v - vmin) / (ratio * dv) + c1.b
        } else {
          c[0] = (c3.r - c2.r) * (v - vmid) / ((1 - ratio) * dv) + c2.r
          c[1] = (c3.g - c2.g) * (v - vmid) / ((1 - ratio) * dv) + c2.g
          c[2] = (c3.b - c2.b) * (v - vmid) / ((1 - ratio) * dv) + c2.b
        }
        break
      case 21:
        c1.r = 255 / 255.0
        c1.g = 255 / 255.0
        c1.b = 200 / 255.0
        c2.r = 150 / 255.0
        c2.g = 150 / 255.0
        c2.b = 255 / 255.0
        c[0] = (c2.r - c1.r) * (v - vmin) / dv + c1.r
        c[1] = (c2.g - c1.g) * (v - vmin) / dv + c1.g
        c[2] = (c2.b - c1.b) * (v - vmin) / dv + c1.b
        break
      case 22:
        c[0] = 1 - (v - vmin) / dv
        c[1] = 1 - (v - vmin) / dv
        c[2] = (v - vmin) / dv
        break
      case 23:
        if (v < (vmin + 0.5 * dv)) {
          c[0] = 1
          c[1] = 2 * (v - vmin) / dv
          c[2] = c[1]
        } else {
          c[0] = 1 - 2 * (v - vmin - 0.5 * dv) / dv
          c[1] = c[0]
          c[2] = 1
        }
        break
      case 24:
        if (v < (vmin + 0.5 * dv)) {
          c[0] = 2 * (v - vmin) / dv
          c[1] = c[0]
          c[2] = 1 - c[0]
        } else {
          c[0] = 1
          c[1] = 1 - 2 * (v - vmin - 0.5 * dv) / dv
          c[2] = 0
        }
        break
      case 25:
        if (v < (vmin + dv / 3)) {
          c[0] = 0
          c[1] = 3 * (v - vmin) / dv
          c[2] = 1
        } else if (v < (vmin + 2 * dv / 3)) {
          c[0] = 3 * (v - vmin - dv / 3) / dv
          c[1] = 1 - c[0]
          c[2] = 1
        } else {
          c[0] = 1
          c[1] = 0
          c[2] = 1 - 3 * (v - vmin - 2 * dv / 3) / dv
        }
        break
    }
    return c
  }

  calc_location(i, j) { // return cnct item

    const du = TWOPI / this.resolution, // Theta
      dv = PI / this.resolution, // Phi
      dx = 1. / this.resolution,
      u = du * i,
      v = dv * j,
      du10 = du / 10.,
      dv10 = dv / 10,
      idx = i * dx,
      jdx = j * dx

    const coords = this.calc_coord(u, v)

    // update max_val
    this.max_val = max(this.max_val, maxv(coords))

    return { // coord, norm ,col, calc_texture
      coords: coords,
      normals: this.calc_normals(coords, this.calc_coord(u + du10, v), this.calc_coord(u, v + dv10)),
      colors: this.calc_color(u, 0, TWOPI, this.colorMap),
      textures: this.calc_texture(idx, jdx)
    }
  }

  calc_mesh() { // calc cnct
    const res = this.resolution

    for (let i = 0, ins_loc = 0, ins_txt = 0; i < res; i++)
      for (let j = 0; j < res; j++, ins_loc += 3, ins_txt += 2) {

        let loc = this.calc_location(j, i)

        this.mesh.coords.set(loc.coords, ins_loc)
        this.mesh.normals.set(loc.normals, ins_loc)
        this.mesh.colors.set(loc.colors, ins_loc)
        this.mesh.textures.set(loc.textures, ins_txt)
      }

    this.scale_coords(1. / this.max_val)
  }

  scale_coords(scale) {
    this.mesh.coords.forEach((_, i, v) => {
      v[i] *= scale
    })
    return this.mesh.coords
  }
}