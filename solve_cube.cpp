// g++ solve_cube.cpp -O3 -std=c++20 -o solve_cube

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <iomanip>
#include <span>
#include <iostream>
#include <fstream>

using uint128 = __uint128_t;

// https://stackoverflow.com/a/61109975/362904
template <class T>
T reverse_bits(T n) {
    short bits = sizeof(n) * 8; 
    T mask = ~T(0); // equivalent to uint32_t mask = 0b11111111111111111111111111111111;
    
    while (bits >>= 1) {
        mask ^= mask << (bits); // will convert mask to 0b00000000000000001111111111111111;
        n = (n & ~mask) >> bits | (n & mask) << bits; // divide and conquer
    }

    return n;
}

template <class T>
T find_first_zero_bit(T value) {
  for (auto i = 0; i < sizeof(T) * 8; ++i) {
    auto bit = T(1) << i;
    if (!(value & bit)) {
      return bit;
    }
  }
  return 0;
}

uint128 gen_piece(int x, int y, int z, int dir, int rot, int dis) {
  uint128 res = 0;
  for (auto i = 0; i < 4; ++i) {
    int sx = x;
    int sy = y;
    int sz = z;
    
    int* p = std::array{&sx, &sy, &sz}[dir / 2];
    int d = (dir % 2) ? 1 : -1;
    *p += d * i;

    if (sx < 0 || sy < 0 || sz < 0 ||
        sx >= 5 || sy >= 5 || sz >= 5) 
    {
      return 0;
    }

    res |= uint128(1) << (sz * 5 * 5 + sy * 5 + sx);

    if (i == 1 + dis) {
      p = std::array{std::array{&sy, &sz}, std::array{&sx, &sz}, std::array{&sx, &sy}}[dir / 2][rot / 2];
      d = (rot % 2) ? 1 : -1;
      *p += d;

      if (sx < 0 || sy < 0 || sz < 0 ||
          sx >= 5 || sy >= 5 || sz >= 5) 
      {
        return 0;
      }

      res |= uint128(1) << (sz * 5 * 5 + sy * 5 + sx);
    }
  }

  return res;
}

bool find_solution(const std::span<uint128> pieces, std::span<short> output) {
  int o = 0;
  uint128 cube = 0;
  uint128 bit = find_first_zero_bit(cube);

  for (int i = 0; i < pieces.size(); ++i) {
    for (; i < pieces.size(); ++i) {
      if (!(pieces[i] & bit)) {
        // skip:
        //continue

        // more big brain skip:
        if (pieces[i] & (bit - 1)) {
          continue;
        }
        break;
      }
      auto& piece = pieces[i];
      if (cube & piece) {
        // skip
        continue;
      }

      output[o++] = i;
      if (o == output.size()) {
        return true;
      }

      cube |= pieces[i];
      bit = find_first_zero_bit(cube);
    }
    if (o > 0) {
      // pop 1 from stack
      i = output[--o];
      cube &= ~pieces[i];
      bit = find_first_zero_bit(cube);
    }
  }
  return false;
}

void export_stl(int r, uint128 piece) {
  std::ofstream cube("solved_cube_" + std::to_string(r) + ".stl", std::ios::out | std::ios::trunc);
  cube << "solid piece_" << r << std::endl;
  static std::array triangles{
    -1.0f,-1.0f,-1.0f, // triangle 1 : begin
    -1.0f,-1.0f, 1.0f,
    -1.0f, 1.0f, 1.0f, // triangle 1 : end
    1.0f, 1.0f,-1.0f, // triangle 2 : begin
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f,-1.0f, // triangle 2 : end
    1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f,-1.0f,
    1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    -1.0f,-1.0f, 1.0f,
    1.0f,-1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f,-1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f, 1.0f,-1.0f,
    -1.0f, 1.0f,-1.0f,
    1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f, 1.0f,
    1.0f,-1.0f, 1.0f
  };
  static std::array normals{
    0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f
  };

  for (auto j = 0; j < triangles.size(); j += 3*3) {
    cube << "facet normal " 
      << normals[j/3 + 0] << " "
      << normals[j/3 + 1] << " "
      << normals[j/3 + 2] << std::endl;
    cube << "outer loop" << std::endl;
    for (auto x = 0; x < 5; ++x)
    for (auto y = 0; y < 5; ++y)
    for (auto z = 0; z < 5; ++z) {
      if (piece & (uint128(1) << (z * 5 * 5 + y * 5 + x))) {
        for (auto k = j; k < j + 3*3; k += 3) {
          cube << "vertex "
            << triangles[k + 0] / 2 + x << " "
            << triangles[k + 1] / 2 + y << " "
            << triangles[k + 2] / 2 + z << std::endl;
        }
      }
    }
    cube << "endloop" << std::endl;
    cube << "endfacet" << std::endl;
  }
  cube << "endsolid piece_" << r << std::endl;
  cube.close();
}

int main() {
  std::cout << "GEN PIECES" << std::endl;
  std::vector<uint128> pieces;

  for (auto x = 0; x < 5; ++x)
  for (auto y = 0; y < 5; ++y)
  for (auto z = 0; z < 5; ++z)
  for (auto d = 0; d < 6; ++d)
  for (auto r = 0; r < 4; ++r)
  for (auto t = 0; t < 2; ++t) {
    auto res = gen_piece(x, y, z, d, r, t);
    if (res == 0) {
      continue;
    }
    if (std::any_of(pieces.begin(), pieces.end(), [res](auto& x){ return x == res; })) {
      continue;
    }
    pieces.push_back(res);
  }
  std::cout << pieces.size() << " PIECES" << std::endl;
  std::cout << "SORT PIECES" << std::endl;

  // sepcial sort !
  std::for_each(pieces.begin(), pieces.end(), [](auto& x) {
    x = reverse_bits(x);
  });
  std::sort(pieces.begin(), pieces.end());
  std::reverse(pieces.begin(), pieces.end());
  std::for_each(pieces.begin(), pieces.end(), [](auto& x) {
    x = reverse_bits(x);
  });

  std::cout << "SEARCH SOLUTION" << std::endl;
  std::vector<short> res;
  res.resize(25);
  find_solution(pieces, res);

  // output
  for (auto r : res) {
    auto piece = pieces[r];
    std::cout << std::setw(3) << int(r) << ":";
    for (auto i = 0; i < 5*5*5; ++i) {
      if (piece & (uint128(1) << i)) {
        std::cout << "???";
      } else {
        std::cout << " ";
      }
    }
    std::cout << std::endl;

    // multiple STLs ..
    // Blender doesn't import "solids" as seperated meshes
    export_stl(r, piece);
  }
  
  std::cout << "DONE" << std::endl;
}