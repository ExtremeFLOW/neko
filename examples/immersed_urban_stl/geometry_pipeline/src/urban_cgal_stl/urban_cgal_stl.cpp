#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

struct FaceInfo {
  int nesting_level = -1;
  bool in_domain() const { return nesting_level % 2 == 1; }
};

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_2<K>;
using Fbb = CGAL::Triangulation_face_base_with_info_2<FaceInfo, K>;
using Fb = CGAL::Constrained_triangulation_face_base_2<K, Fbb>;
using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TDS, CGAL::Exact_predicates_tag>;
using Point2 = CDT::Point;
using FaceHandle = CDT::Face_handle;
using Edge = CDT::Edge;

struct Point3 {
  double x;
  double y;
  double z;
};

struct Triangle {
  Point3 a;
  Point3 b;
  Point3 c;
};

struct Building {
  double height;
  std::vector<Point2> ring;
};

struct Input {
  double xmin;
  double ymin;
  double xmax;
  double ymax;
  int nx;
  int ny;
  std::vector<double> xs;
  std::vector<double> ys;
  std::vector<double> z;
  std::vector<Building> buildings;
};

static double to_double(const Point2& p, int axis) {
  return CGAL::to_double(axis == 0 ? p.x() : p.y());
}

static void expect_tag(std::istream& in, const std::string& expected) {
  std::string got;
  in >> got;
  if (got != expected) {
    throw std::runtime_error("Expected tag '" + expected + "', got '" + got + "'");
  }
}

static Input read_input(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Cannot open input " + path);
  }

  Input input{};
  expect_tag(in, "bbox");
  in >> input.xmin >> input.ymin >> input.xmax >> input.ymax;
  expect_tag(in, "grid");
  in >> input.nx >> input.ny;

  expect_tag(in, "xs");
  input.xs.resize(input.nx);
  for (double& x : input.xs) {
    in >> x;
  }
  expect_tag(in, "ys");
  input.ys.resize(input.ny);
  for (double& y : input.ys) {
    in >> y;
  }
  expect_tag(in, "z");
  input.z.resize(static_cast<size_t>(input.nx) * static_cast<size_t>(input.ny));
  for (double& z : input.z) {
    in >> z;
  }

  int building_count = 0;
  expect_tag(in, "buildings");
  in >> building_count;
  input.buildings.reserve(building_count);
  for (int b = 0; b < building_count; ++b) {
    Building building{};
    int n = 0;
    expect_tag(in, "height");
    in >> building.height;
    expect_tag(in, "ring");
    in >> n;
    building.ring.reserve(n);
    for (int i = 0; i < n; ++i) {
      double x = 0.0;
      double y = 0.0;
      in >> x >> y;
      building.ring.emplace_back(x, y);
    }
    if (building.ring.size() >= 3 && building.height > 0.0) {
      input.buildings.push_back(std::move(building));
    }
  }

  return input;
}

static double sample_z(const Input& input, double x, double y) {
  auto xi = std::upper_bound(input.xs.begin(), input.xs.end(), x);
  auto yi = std::upper_bound(input.ys.begin(), input.ys.end(), y);
  int i = static_cast<int>(std::distance(input.xs.begin(), xi)) - 1;
  int j = static_cast<int>(std::distance(input.ys.begin(), yi)) - 1;
  i = std::max(0, std::min(i, input.nx - 2));
  j = std::max(0, std::min(j, input.ny - 2));
  double x0 = input.xs[i];
  double x1 = input.xs[i + 1];
  double y0 = input.ys[j];
  double y1 = input.ys[j + 1];
  double tx = x1 == x0 ? 0.0 : (x - x0) / (x1 - x0);
  double ty = y1 == y0 ? 0.0 : (y - y0) / (y1 - y0);
  auto at = [&](int jj, int ii) {
    return input.z[static_cast<size_t>(jj) * static_cast<size_t>(input.nx) + static_cast<size_t>(ii)];
  };
  double z00 = at(j, i);
  double z10 = at(j, i + 1);
  double z01 = at(j + 1, i);
  double z11 = at(j + 1, i + 1);
  return (1.0 - tx) * (1.0 - ty) * z00 + tx * (1.0 - ty) * z10 +
         (1.0 - tx) * ty * z01 + tx * ty * z11;
}

static double signed_area(const std::vector<Point2>& ring) {
  double area = 0.0;
  for (size_t i = 0; i < ring.size(); ++i) {
    const Point2& a = ring[i];
    const Point2& b = ring[(i + 1) % ring.size()];
    area += to_double(a, 0) * to_double(b, 1) - to_double(b, 0) * to_double(a, 1);
  }
  return 0.5 * area;
}

static void insert_closed_constraint(CDT& cdt, const std::vector<Point2>& ring) {
  std::vector<CDT::Vertex_handle> handles;
  handles.reserve(ring.size());
  for (const Point2& p : ring) {
    handles.push_back(cdt.insert(p));
  }
  for (size_t i = 0; i < handles.size(); ++i) {
    cdt.insert_constraint(handles[i], handles[(i + 1) % handles.size()]);
  }
}

static void mark_domains(CDT& ct, FaceHandle start, int index, std::list<Edge>& border) {
  if (start->info().nesting_level != -1) {
    return;
  }
  std::list<FaceHandle> queue;
  queue.push_back(start);
  while (!queue.empty()) {
    FaceHandle fh = queue.front();
    queue.pop_front();
    if (fh->info().nesting_level != -1) {
      continue;
    }
    fh->info().nesting_level = index;
    for (int i = 0; i < 3; ++i) {
      Edge edge(fh, i);
      FaceHandle neighbor = fh->neighbor(i);
      if (neighbor->info().nesting_level == -1) {
        if (ct.is_constrained(edge)) {
          border.push_back(edge);
        } else {
          queue.push_back(neighbor);
        }
      }
    }
  }
}

static void mark_domains(CDT& cdt) {
  for (auto fit = cdt.all_faces_begin(); fit != cdt.all_faces_end(); ++fit) {
    fit->info().nesting_level = -1;
  }
  std::list<Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border);
  while (!border.empty()) {
    Edge edge = border.front();
    border.pop_front();
    FaceHandle neighbor = edge.first->neighbor(edge.second);
    if (neighbor->info().nesting_level == -1) {
      mark_domains(cdt, neighbor, edge.first->info().nesting_level + 1, border);
    }
  }
}

static CDT build_cdt(const Input& input) {
  CDT cdt;
  std::vector<Point2> outer = {
      Point2(input.xmin, input.ymin),
      Point2(input.xmax, input.ymin),
      Point2(input.xmax, input.ymax),
      Point2(input.xmin, input.ymax),
  };
  insert_closed_constraint(cdt, outer);

  for (const Building& building : input.buildings) {
    insert_closed_constraint(cdt, building.ring);
  }

  for (double y : input.ys) {
    for (double x : input.xs) {
      cdt.insert(Point2(x, y));
    }
  }

  mark_domains(cdt);
  return cdt;
}

static Point3 p3(const Input& input, const Point2& p, double z) {
  return {to_double(p, 0) - input.xmin, to_double(p, 1) - input.ymin, z};
}

static double orientation_xy(const Point3& a, const Point3& b, const Point3& c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

static void add_up_triangle(std::vector<Triangle>& tris, Point3 a, Point3 b, Point3 c) {
  if (orientation_xy(a, b, c) < 0.0) {
    std::swap(b, c);
  }
  tris.push_back({a, b, c});
}

static void add_down_triangle(std::vector<Triangle>& tris, Point3 a, Point3 b, Point3 c) {
  if (orientation_xy(a, b, c) > 0.0) {
    std::swap(b, c);
  }
  tris.push_back({a, b, c});
}

static void add_quad(std::vector<Triangle>& tris, Point3 a, Point3 b, Point3 c, Point3 d) {
  tris.push_back({a, b, c});
  tris.push_back({a, c, d});
}

static std::vector<Triangle> terrain_triangles(const Input& input, const CDT& cdt) {
  std::vector<Triangle> tris;
  for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
    if (!fit->info().in_domain()) {
      continue;
    }
    Point2 pa = fit->vertex(0)->point();
    Point2 pb = fit->vertex(1)->point();
    Point2 pc = fit->vertex(2)->point();
    Point3 a = p3(input, pa, sample_z(input, to_double(pa, 0), to_double(pa, 1)));
    Point3 b = p3(input, pb, sample_z(input, to_double(pb, 0), to_double(pb, 1)));
    Point3 c = p3(input, pc, sample_z(input, to_double(pc, 0), to_double(pc, 1)));
    add_up_triangle(tris, a, b, c);
  }
  return tris;
}

static std::vector<std::array<int, 3>> ear_clip(const std::vector<Point2>& ring) {
  std::vector<Point2> pts = ring;
  if (signed_area(pts) < 0.0) {
    std::reverse(pts.begin(), pts.end());
  }
  std::vector<int> remaining(pts.size());
  for (size_t i = 0; i < pts.size(); ++i) {
    remaining[i] = static_cast<int>(i);
  }
  auto point_in_tri = [&](const Point2& p, const Point2& a, const Point2& b, const Point2& c) {
    double px = to_double(p, 0), py = to_double(p, 1);
    double ax = to_double(a, 0), ay = to_double(a, 1);
    double bx = to_double(b, 0), by = to_double(b, 1);
    double cx = to_double(c, 0), cy = to_double(c, 1);
    double v0x = cx - ax, v0y = cy - ay;
    double v1x = bx - ax, v1y = by - ay;
    double v2x = px - ax, v2y = py - ay;
    double denom = v0x * v1y - v1x * v0y;
    if (std::abs(denom) < 1e-12) {
      return false;
    }
    double u = (v2x * v1y - v1x * v2y) / denom;
    double v = (v0x * v2y - v2x * v0y) / denom;
    return u >= -1e-12 && v >= -1e-12 && u + v <= 1.0 + 1e-12;
  };

  std::vector<std::array<int, 3>> out;
  int guard = 0;
  while (remaining.size() > 3 && guard < static_cast<int>(pts.size() * pts.size())) {
    ++guard;
    bool clipped = false;
    for (size_t pos = 0; pos < remaining.size(); ++pos) {
      int prev = remaining[(pos + remaining.size() - 1) % remaining.size()];
      int curr = remaining[pos];
      int next = remaining[(pos + 1) % remaining.size()];
      double cross = (to_double(pts[curr], 0) - to_double(pts[prev], 0)) *
                         (to_double(pts[next], 1) - to_double(pts[prev], 1)) -
                     (to_double(pts[curr], 1) - to_double(pts[prev], 1)) *
                         (to_double(pts[next], 0) - to_double(pts[prev], 0));
      if (cross <= 1e-10) {
        continue;
      }
      bool contains = false;
      for (int idx : remaining) {
        if (idx == prev || idx == curr || idx == next) {
          continue;
        }
        if (point_in_tri(pts[idx], pts[prev], pts[curr], pts[next])) {
          contains = true;
          break;
        }
      }
      if (contains) {
        continue;
      }
      out.push_back({prev, curr, next});
      remaining.erase(remaining.begin() + static_cast<long>(pos));
      clipped = true;
      break;
    }
    if (!clipped) {
      out.clear();
      for (size_t i = 1; i + 1 < pts.size(); ++i) {
        out.push_back({0, static_cast<int>(i), static_cast<int>(i + 1)});
      }
      return out;
    }
  }
  if (remaining.size() == 3) {
    out.push_back({remaining[0], remaining[1], remaining[2]});
  }
  return out;
}

static std::vector<Triangle> building_surface_triangles(const Input& input, bool closed_bottoms) {
  std::vector<Triangle> tris;
  for (Building building : input.buildings) {
    if (signed_area(building.ring) < 0.0) {
      std::reverse(building.ring.begin(), building.ring.end());
    }
    std::vector<double> base_z;
    base_z.reserve(building.ring.size());
    for (const Point2& p : building.ring) {
      base_z.push_back(sample_z(input, to_double(p, 0), to_double(p, 1)));
    }

    auto roof_tris = ear_clip(building.ring);
    for (const auto& tri : roof_tris) {
      Point3 a = p3(input, building.ring[tri[0]], base_z[tri[0]] + building.height);
      Point3 b = p3(input, building.ring[tri[1]], base_z[tri[1]] + building.height);
      Point3 c = p3(input, building.ring[tri[2]], base_z[tri[2]] + building.height);
      add_up_triangle(tris, a, b, c);
      if (closed_bottoms) {
        Point3 ba = p3(input, building.ring[tri[0]], base_z[tri[0]]);
        Point3 bb = p3(input, building.ring[tri[1]], base_z[tri[1]]);
        Point3 bc = p3(input, building.ring[tri[2]], base_z[tri[2]]);
        add_down_triangle(tris, ba, bb, bc);
      }
    }

    for (size_t i = 0; i < building.ring.size(); ++i) {
      size_t j = (i + 1) % building.ring.size();
      Point3 b0 = p3(input, building.ring[i], base_z[i]);
      Point3 b1 = p3(input, building.ring[j], base_z[j]);
      Point3 t1 = p3(input, building.ring[j], base_z[j] + building.height);
      Point3 t0 = p3(input, building.ring[i], base_z[i] + building.height);
      add_quad(tris, b0, b1, t1, t0);
    }
  }
  return tris;
}

static bool point_in_ring(const std::vector<Point2>& ring, double x, double y) {
  bool inside = false;
  for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
    double xi = to_double(ring[i], 0);
    double yi = to_double(ring[i], 1);
    double xj = to_double(ring[j], 0);
    double yj = to_double(ring[j], 1);
    bool intersects = ((yi > y) != (yj > y)) &&
                      (x < (xj - xi) * (y - yi) / ((yj - yi) == 0.0 ? 1e-30 : (yj - yi)) + xi);
    if (intersects) {
      inside = !inside;
    }
  }
  return inside;
}

static double building_height_at(const Input& input, double x, double y) {
  double height = 0.0;
  for (const Building& building : input.buildings) {
    if (point_in_ring(building.ring, x, y)) {
      height = std::max(height, building.height);
    }
  }
  return height;
}

static double face_building_height(const Input& input, FaceHandle face) {
  Point2 pa = face->vertex(0)->point();
  Point2 pb = face->vertex(1)->point();
  Point2 pc = face->vertex(2)->point();
  double cx = (to_double(pa, 0) + to_double(pb, 0) + to_double(pc, 0)) / 3.0;
  double cy = (to_double(pa, 1) + to_double(pb, 1) + to_double(pc, 1)) / 3.0;
  return building_height_at(input, cx, cy);
}

static std::vector<Triangle> cdt_building_triangles(const Input& input, const CDT& cdt, bool closed_bottoms) {
  std::vector<Triangle> tris;

  for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
    if (fit->info().in_domain() || fit->info().nesting_level <= 1) {
      continue;
    }
    double height = face_building_height(input, fit);
    if (height <= 0.0) {
      continue;
    }
    Point2 pa = fit->vertex(0)->point();
    Point2 pb = fit->vertex(1)->point();
    Point2 pc = fit->vertex(2)->point();
    Point3 a = p3(input, pa, sample_z(input, to_double(pa, 0), to_double(pa, 1)) + height);
    Point3 b = p3(input, pb, sample_z(input, to_double(pb, 0), to_double(pb, 1)) + height);
    Point3 c = p3(input, pc, sample_z(input, to_double(pc, 0), to_double(pc, 1)) + height);
    add_up_triangle(tris, a, b, c);
    if (closed_bottoms) {
      Point3 ba = p3(input, pa, sample_z(input, to_double(pa, 0), to_double(pa, 1)));
      Point3 bb = p3(input, pb, sample_z(input, to_double(pb, 0), to_double(pb, 1)));
      Point3 bc = p3(input, pc, sample_z(input, to_double(pc, 0), to_double(pc, 1)));
      add_down_triangle(tris, ba, bb, bc);
    }
  }

  for (auto eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit) {
    FaceHandle f = eit->first;
    int i = eit->second;
    FaceHandle n = f->neighbor(i);
    bool f_terrain = f->info().in_domain();
    bool n_terrain = !cdt.is_infinite(n) && n->info().in_domain();
    bool f_building = !f_terrain && f->info().nesting_level > 1;
    bool n_building = !n_terrain && !cdt.is_infinite(n) && n->info().nesting_level > 1;
    if (!((f_terrain && n_building) || (n_terrain && f_building))) {
      continue;
    }

    Point2 pa = f->vertex(cdt.cw(i))->point();
    Point2 pb = f->vertex(cdt.ccw(i))->point();
    double height = face_building_height(input, f_building ? f : n);
    if (height <= 0.0) {
      continue;
    }
    double za = sample_z(input, to_double(pa, 0), to_double(pa, 1));
    double zb = sample_z(input, to_double(pb, 0), to_double(pb, 1));
    Point3 a = p3(input, pa, za);
    Point3 b = p3(input, pb, zb);
    Point3 bt = p3(input, pb, zb + height);
    Point3 at = p3(input, pa, za + height);
    add_quad(tris, a, b, bt, at);
  }

  return tris;
}

static std::vector<Triangle> top_cap_triangles(const Input& input, const CDT& cdt, double top_z) {
  std::vector<Triangle> tris;
  for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
    int domain = fit->info().nesting_level;
    if (domain <= 0) {
      continue;
    }
    Point3 a = p3(input, fit->vertex(0)->point(), top_z);
    Point3 b = p3(input, fit->vertex(1)->point(), top_z);
    Point3 c = p3(input, fit->vertex(2)->point(), top_z);
    add_down_triangle(tris, a, b, c);
  }
  return tris;
}

static double area2_3d(const Triangle& t) {
  double ux = t.b.x - t.a.x, uy = t.b.y - t.a.y, uz = t.b.z - t.a.z;
  double vx = t.c.x - t.a.x, vy = t.c.y - t.a.y, vz = t.c.z - t.a.z;
  double nx = uy * vz - uz * vy;
  double ny = uz * vx - ux * vz;
  double nz = ux * vy - uy * vx;
  return nx * nx + ny * ny + nz * nz;
}

static void remove_degenerate(std::vector<Triangle>& tris) {
  tris.erase(
      std::remove_if(
          tris.begin(),
          tris.end(),
          [](const Triangle& t) {
            return area2_3d(t) < 1e-18;
          }),
      tris.end());
}

using EdgeKey = std::tuple<long long, long long, long long, long long, long long, long long>;

static long long q(double v) {
  return std::llround(v * 1000000.0);
}

static EdgeKey edge_key(Point3 a, Point3 b) {
  std::array<long long, 3> qa = {q(a.x), q(a.y), q(a.z)};
  std::array<long long, 3> qb = {q(b.x), q(b.y), q(b.z)};
  if (qb < qa) {
    std::swap(qa, qb);
  }
  return {qa[0], qa[1], qa[2], qb[0], qb[1], qb[2]};
}

static std::vector<Triangle> close_boundaries_to_top(const std::vector<Triangle>& surface, double top_z) {
  std::map<EdgeKey, int> counts;
  std::map<EdgeKey, std::pair<Point3, Point3>> oriented;
  for (const Triangle& t : surface) {
    std::array<std::pair<Point3, Point3>, 3> edges = {{{t.a, t.b}, {t.b, t.c}, {t.c, t.a}}};
    for (const auto& edge : edges) {
      EdgeKey key = edge_key(edge.first, edge.second);
      counts[key] += 1;
      oriented[key] = edge;
    }
  }

  std::vector<Triangle> tris;
  for (const auto& item : counts) {
    if (item.second != 1) {
      continue;
    }
    auto edge = oriented[item.first];
    Point3 a = edge.first;
    Point3 b = edge.second;
    Point3 bt{b.x, b.y, top_z};
    Point3 at{a.x, a.y, top_z};
    add_quad(tris, a, b, bt, at);
  }
  return tris;
}

static Point3 normal(const Triangle& t) {
  double ux = t.b.x - t.a.x, uy = t.b.y - t.a.y, uz = t.b.z - t.a.z;
  double vx = t.c.x - t.a.x, vy = t.c.y - t.a.y, vz = t.c.z - t.a.z;
  Point3 n{uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx};
  double len = std::sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
  if (len == 0.0) {
    return {0.0, 0.0, 0.0};
  }
  return {n.x / len, n.y / len, n.z / len};
}

static void write_float(std::ofstream& out, float v) {
  out.write(reinterpret_cast<const char*>(&v), sizeof(float));
}

static void write_binary_stl(const std::string& path, const std::vector<Triangle>& tris) {
  std::ofstream out(path, std::ios::binary);
  if (!out) {
    throw std::runtime_error("Cannot write " + path);
  }
  std::string header = "urban_cgal_stl";
  header.resize(80, ' ');
  out.write(header.data(), 80);
  uint32_t count = static_cast<uint32_t>(tris.size());
  out.write(reinterpret_cast<const char*>(&count), sizeof(uint32_t));
  for (const Triangle& t : tris) {
    Point3 n = normal(t);
    for (double v : {n.x, n.y, n.z, t.a.x, t.a.y, t.a.z, t.b.x, t.b.y, t.b.z, t.c.x, t.c.y, t.c.z}) {
      write_float(out, static_cast<float>(v));
    }
    uint16_t attr = 0;
    out.write(reinterpret_cast<const char*>(&attr), sizeof(uint16_t));
  }
}

int main(int argc, char** argv) {
  if (argc != 7) {
    std::cerr << "usage: urban_cgal_stl input.txt terrain.stl buildings.stl combined.stl boxed.stl top_z\n";
    return 2;
  }

  try {
    Input input = read_input(argv[1]);
    double top_z = std::stod(argv[6]);
    CDT cdt = build_cdt(input);

    std::vector<Triangle> terrain = terrain_triangles(input, cdt);
    std::vector<Triangle> buildings_closed = cdt_building_triangles(input, cdt, true);
    std::vector<Triangle> building_skin = cdt_building_triangles(input, cdt, false);

    std::vector<Triangle> combined = terrain;
    combined.insert(combined.end(), building_skin.begin(), building_skin.end());

    std::vector<Triangle> boxed = combined;
    std::vector<Triangle> sides = close_boundaries_to_top(combined, top_z);
    std::vector<Triangle> cap = top_cap_triangles(input, cdt, top_z);
    boxed.insert(boxed.end(), sides.begin(), sides.end());
    boxed.insert(boxed.end(), cap.begin(), cap.end());

    remove_degenerate(terrain);
    remove_degenerate(buildings_closed);
    remove_degenerate(combined);
    remove_degenerate(boxed);

    write_binary_stl(argv[2], terrain);
    write_binary_stl(argv[3], buildings_closed);
    write_binary_stl(argv[4], combined);
    write_binary_stl(argv[5], boxed);

    std::cout << "{"
              << "\"terrain_facets\":" << terrain.size() << ","
              << "\"buildings_facets\":" << buildings_closed.size() << ","
              << "\"combined_facets\":" << combined.size() << ","
              << "\"boxed_facets\":" << boxed.size() << ","
              << "\"cdt_vertices\":" << cdt.number_of_vertices() << ","
              << "\"cdt_faces\":" << cdt.number_of_faces() << ","
              << "\"building_regions\":" << input.buildings.size()
              << "}\n";
  } catch (const std::exception& e) {
    std::cerr << "urban_cgal_stl: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
