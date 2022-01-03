#include <bits/stdc++.h>
#include <random>
#ifdef _MSC_VER
#include <ppl.h>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#else
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#endif

/** compro_io **/

/* tuple */
// out
namespace aux {
    template<typename T, unsigned N, unsigned L>
    struct tp {
        static void output(std::ostream& os, const T& v) {
            os << std::get<N>(v) << ", ";
            tp<T, N + 1, L>::output(os, v);
        }
    };
    template<typename T, unsigned N>
    struct tp<T, N, N> {
        static void output(std::ostream& os, const T& v) { os << std::get<N>(v); }
    };
}
template<typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& t) {
    os << '[';
    aux::tp<std::tuple<Ts...>, 0, sizeof...(Ts) - 1>::output(os, t);
    return os << ']';
}

template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x);

/* pair */
// out
template<class S, class T>
std::ostream& operator<<(std::ostream& os, const std::pair<S, T>& p) {
    return os << "[" << p.first << ", " << p.second << "]";
}
// in
template<class S, class T>
std::istream& operator>>(std::istream& is, const std::pair<S, T>& p) {
    return is >> p.first >> p.second;
}

/* container */
// out
template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x) {
    bool f = true;
    os << "[";
    for (auto& y : x) {
        os << (f ? "" : ", ") << y;
        f = false;
    }
    return os << "]";
}
// in
template <
    class T,
    class = decltype(std::begin(std::declval<T&>())),
    class = typename std::enable_if<!std::is_same<T, std::string>::value>::type
>
std::istream& operator>>(std::istream& is, T& a) {
    for (auto& x : a) is >> x;
    return is;
}

/* struct */
template<typename T>
auto operator<<(std::ostream& out, const T& t) -> decltype(out << t.stringify()) {
    out << t.stringify();
    return out;
}

/* setup */
struct IOSetup {
    IOSetup(bool f) {
        if (f) { std::cin.tie(nullptr); std::ios::sync_with_stdio(false); }
        std::cout << std::fixed << std::setprecision(15);
    }
} iosetup(true);

/** string formatter **/
template<typename... Ts>
std::string format(const std::string& f, Ts... t) {
    size_t l = std::snprintf(nullptr, 0, f.c_str(), t...);
    std::vector<char> b(l + 1);
    std::snprintf(&b[0], l + 1, f.c_str(), t...);
    return std::string(&b[0], &b[0] + l);
}

template<typename T>
std::string stringify(const T& x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

/* dump */
#define ENABLE_DUMP
#ifdef ENABLE_DUMP
#define DUMPOUT std::cerr
std::ostringstream DUMPBUF;
#define dump(...) do{DUMPBUF<<"  ";DUMPBUF<<#__VA_ARGS__<<" :[DUMP - "<<__LINE__<<":"<<__FUNCTION__<<"]"<<std::endl;DUMPBUF<<"    ";dump_func(__VA_ARGS__);DUMPOUT<<DUMPBUF.str();DUMPBUF.str("");DUMPBUF.clear();}while(0);
void dump_func() { DUMPBUF << std::endl; }
template <class Head, class... Tail> void dump_func(Head&& head, Tail&&... tail) { DUMPBUF << head; if (sizeof...(Tail) == 0) { DUMPBUF << " "; } else { DUMPBUF << ", "; } dump_func(std::move(tail)...); }
#else
#define dump(...) void(0);
#endif

/* timer */
class Timer {
    double t = 0, paused = 0, tmp;
public:
    Timer() { reset(); }
    static double time() {
#ifdef _MSC_VER
        return __rdtsc() / 3.0e9;
#else
        unsigned long long a, d;
        __asm__ volatile("rdtsc"
            : "=a"(a), "=d"(d));
        return (d << 32 | a) / 3.0e9;
#endif
    }
    void reset() { t = time(); }
    void pause() { tmp = time(); }
    void restart() { paused += time() - tmp; }
    double elapsed_ms() { return (time() - t - paused) * 1000.0; }
} timer;

/* rand */
struct Xorshift {
    uint64_t x = 88172645463325252LL;
    void set_seed(unsigned seed, int rep = 100) { x = uint64_t((seed + 1) * 10007); for (int i = 0; i < rep; i++) next_int(); }
    unsigned next_int() { x = x ^ (x << 7); return x = x ^ (x >> 9); }
    unsigned next_int(unsigned mod) { x = x ^ (x << 7); x = x ^ (x >> 9); return x % mod; }
    unsigned next_int(unsigned l, unsigned r) { x = x ^ (x << 7); x = x ^ (x >> 9); return x % (r - l + 1) + l; } // inclusive
    double next_double() { return double(next_int()) / UINT_MAX; }
} rnd;

/* shuffle */
template<typename T>
void shuffle_vector(std::vector<T>& v, Xorshift& rnd) {
    int n = v.size();
    for (int i = n - 1; i >= 1; i--) {
        int r = rnd.next_int(i);
        std::swap(v[i], v[r]);
    }
}

/* split */
std::vector<std::string> split(std::string str, const std::string& delim) {
    for (char& c : str) if (delim.find(c) != std::string::npos) c = ' ';
    std::istringstream iss(str);
    std::vector<std::string> parsed;
    std::string buf;
    while (iss >> buf) parsed.push_back(buf);
    return parsed;
}

template<typename A, size_t N, typename T> inline void Fill(A(&array)[N], const T& val) {
    std::fill((T*)array, (T*)(array + N), val);
}

template<typename T> bool chmax(T& a, const T& b) { if (a < b) { a = b; return true; } return false; }
template<typename T> bool chmin(T& a, const T& b) { if (a > b) { a = b; return true; } return false; }



using pii = std::pair<int, int>;

constexpr int H = 50;
constexpr int W = 50;
constexpr int K = 2500; // num turns
constexpr char d2c[] = "RULD";
int c2d[256];
constexpr int dr[] = { 0, -1, 0, 1 };
constexpr int dc[] = { 1, 0, -1, 0 };

void init() {
    c2d['R'] = 0; c2d['U'] = 1; c2d['L'] = 2; c2d['D'] = 3;
}

struct Food {
    int id, r, c, value, decay;
    Food(int id, int r, int c, int value, int decay) : id(id), r(r), c(c), value(value), decay(decay) {}
    std::string stringify() const {
        return format("Food [id=%d, r=%d, c=%d, value=%d, decay=%d]", id, r, c, value, decay);
    }
};
std::ostream& operator<<(std::ostream& o, Food* f) {
    o << (f ? f->stringify() : "null");
    return o;
}

struct TestCase {

    std::vector<std::string> board;
    int N;

    std::vector<Food> foods;
    std::vector<std::vector<Food*>> food_map;

    // last: start, otherwise: food
    // 始点が末尾なのは food の id 管理上都合がよいため
    std::vector<pii> points;

    std::vector<std::vector<int>> dist; // (u, v) の最短距離
    std::vector<std::vector<std::string>> route; // (u, v) の最短経路 (RULD)

    TestCase(std::istream& in) {
        { int buf; in >> buf >> buf >> buf; }
        int sr, sc;
        in >> sr >> sc;
        sr--; sc--;
        board.resize(H);
        in >> board;
        food_map.resize(H, std::vector<Food*>(W, nullptr));
        in >> N;
        for (int i = 0; i < N; i++) {
            int fr, fc, F, D;
            in >> fr >> fc >> F >> D;
            foods.emplace_back(i, fr - 1, fc - 1, F, D);
            points.emplace_back(fr - 1, fc - 1);
        }
        points.emplace_back(sr, sc);
        for (int i = 0; i < N; i++) {
            auto [id, fr, fc, F, D] = foods[i];
            food_map[fr][fc] = &foods[i];
        }
        calc_shortest_paths();
    }

    void calc_shortest_paths() {
        static constexpr int inf = INT_MAX / 2;
        int V = points.size();
        dist.resize(V, std::vector<int>(V, inf));
        route.resize(V, std::vector<std::string>(V));

        int dist_map[H][W];
        char prev_map[H][W];

        for (int u = 0; u < V; u++) {
            // init
            for (int i = 0; i < H; i++) {
                for (int j = 0; j < W; j++) {
                    dist_map[i][j] = inf;
                    prev_map[i][j] = '$';
                }
            }
            // bfs
            auto [sr, sc] = points[u];
            std::queue<pii> qu;
            qu.emplace(sr, sc);
            dist_map[sr][sc] = 0;
            while (!qu.empty()) {
                auto [r, c] = qu.front(); qu.pop();
                for (int d = 0; d < 4; d++) {
                    int nr = r + dr[d], nc = c + dc[d];
                    if (board[nr][nc] == '#' || dist_map[nr][nc] != inf) continue;
                    dist_map[nr][nc] = dist_map[r][c] + 1;
                    prev_map[nr][nc] = d2c[d];
                    qu.emplace(nr, nc);
                }
            }
            // route
            for (int v = 0; v < V; v++) {
                if (u == v) continue;
                // v -> u
                auto [r, c] = points[v];
                auto ch = prev_map[r][c];
                std::string path;
                while (ch != '$') {
                    path += ch;
                    int d = (c2d[ch] + 2) & 3; // 逆方向
                    r += dr[d]; c += dc[d];
                    ch = prev_map[r][c];
                }
                std::reverse(path.begin(), path.end());
                route[u][v] = path;
                dist[u][v] = path.size();
            }
        }
    }

    struct Result {
        int score;
        int highest_score;
        int turn_to_truncate;
    };

    Result evaluate(const std::string& cmds) const {
        Result res = { 0, 0, -1 };
        int turn = -1;
        auto [r, c] = points.back();
        std::vector<bool> used(N, false);
        for (char cmd : cmds) {
            turn++;
            if (cmd == '-') continue;
            int nr = r + dr[c2d[cmd]], nc = c + dc[c2d[cmd]];
            if (board[nr][nc] == '#') continue;
            if (food_map[nr][nc] && !used[food_map[nr][nc]->id]) {
                const auto& food = foods[food_map[nr][nc]->id];
                res.score += food.value - food.decay * turn;
                used[food.id] = true;
                if (res.highest_score < res.score) {
                    res.highest_score = res.score;
                    res.turn_to_truncate = turn;
                }
            }
            r = nr; c = nc;
            if (turn == K - 1) break;
        }
        return res;
    }

    Result evaluate(const std::vector<int>& ids) const {
        Result res = { 0, 0, -1 };
        int turn = -1;
        int prev_id = ids[0]; // start point
        auto [r, c] = points[prev_id];
        std::vector<bool> used(N, false);
        for (int i = 1; i < ids.size(); i++) {
            int id = ids[i];
            for (char cmd : route[prev_id][id]) {
                turn++;
                if (cmd == '-') continue;
                int nr = r + dr[c2d[cmd]], nc = c + dc[c2d[cmd]];
                if (board[nr][nc] == '#') continue;
                if (food_map[nr][nc] && !used[food_map[nr][nc]->id]) {
                    const auto& food = foods[food_map[nr][nc]->id];
                    res.score += food.value - food.decay * turn;
                    used[food.id] = true;
                    if (res.highest_score < res.score) {
                        res.highest_score = res.score;
                        res.turn_to_truncate = turn;
                    }
                }
                r = nr; c = nc;
                if (turn == K - 1) break;
            }
            if (turn == K - 1) break;
            prev_id = id;
        }
        return res;
    }

    Result evaluate(const std::vector<int>& ids, const std::vector<bool>& skip_flag) const {
        Result res = { 0, 0, -1 };
        int turn = -1;
        int prev_id = ids[0]; // start point
        auto [r, c] = points[prev_id];
        std::vector<bool> used(N, false);
        for (int i = 1; i < ids.size(); i++) {
            if (skip_flag[i]) continue;
            int id = ids[i];
            for (char cmd : route[prev_id][id]) {
                turn++;
                if (cmd == '-') continue;
                int nr = r + dr[c2d[cmd]], nc = c + dc[c2d[cmd]];
                if (board[nr][nc] == '#') continue;
                if (food_map[nr][nc] && !used[food_map[nr][nc]->id]) {
                    const auto& food = foods[food_map[nr][nc]->id];
                    res.score += food.value - food.decay * turn;
                    used[food.id] = true;
                    if (res.highest_score < res.score) {
                        res.highest_score = res.score;
                        res.turn_to_truncate = turn;
                    }
                }
                r = nr; c = nc;
                if (turn == K - 1) break;
            }
            if (turn == K - 1) break;
            prev_id = id;
        }
        return res;
    }

    std::string cvt(const std::vector<int>& ids, const std::vector<bool>& skip_flag) const {
        std::string cmds;
        int prev_id = ids[0];
        for (int i = 1; i < ids.size(); i++) {
            if (skip_flag[i]) continue;
            int id = ids[i];
            cmds += route[prev_id][id];
            if (cmds.size() >= K) break;
            prev_id = id;
        }
        if (cmds.size() < K) {
            cmds += std::string(K - cmds.size(), '-');
        }
        else {
            cmds = cmds.substr(0, K);
        }
        return cmds;
    }

};

struct DistanceTSP {

    Xorshift rnd;

    int V;
    std::vector<std::vector<int>> dist;

    int cost;
    std::vector<int> path;

    DistanceTSP(const TestCase& tc) : V(tc.points.size()), dist(tc.dist) {
        cost = 0;
        path.push_back(tc.points.size() - 1);
        for (int u = 0; u < tc.foods.size(); u++) {
            cost += dist[path.back()][u];
            path.push_back(u);
        }
    }

    int two_opt_diff(int idx1, int idx2) const {
        if (idx2 == V - 1) {
            // 終点自由端
            int p1 = path[idx1], p2 = path[idx1 + 1], p3 = path[idx2];
            return dist[p1][p3] - dist[p1][p2];
        }
        else {
            int p1 = path[idx1], p2 = path[idx1 + 1], p3 = path[idx2], p4 = path[idx2 + 1];
            return dist[p1][p3] + dist[p2][p4] - dist[p1][p2] - dist[p3][p4];
        }
    }

    double get_temp(double start_temp, double end_temp, int loop, int num_loop) {
        return end_temp + (start_temp - end_temp) * (num_loop - loop) / num_loop;
    }

    void two_opt_annealing(int num_loop) {
        int V = path.size();
        int min_cost = cost;
        std::vector<int> best_path;
        for (int loop = 0; loop < num_loop; loop++) {
            // "パスの"添字 (点の id ではない)
            int idx1 = rnd.next_int(V - 1);
            int idx2 = rnd.next_int(V);
            if (idx1 == idx2) idx2++;
            if (idx1 > idx2) std::swap(idx1, idx2);
            int diff = two_opt_diff(idx1, idx2);
            double temp = get_temp(1.0, 0.01, loop, num_loop);
            double prob = std::exp(-diff / temp);
            if (prob > rnd.next_double()) {
                // accept
                std::reverse(path.begin() + idx1 + 1, path.begin() + idx2 + 1);
                cost += diff;
                if (cost < min_cost) {
                    min_cost = cost;
                    best_path = path;
                }
            }
            if (!(loop & 0x7FFFFF)) {
                dump(loop, min_cost, cost);
            }
        }
        dump(num_loop, min_cost, cost);
        cost = min_cost;
        path = best_path;
    }

};

int main() {

#ifdef _MSC_VER
    std::ifstream ifs("C:\\dev\\heuristic\\tasks\\RCO2017qualB\\tester\\in\\6.in");
    std::ofstream ofs("C:\\dev\\heuristic\\tasks\\RCO2017qualB\\tester\\out\\6.out");
    std::istream& in = ifs;
    std::ostream& out = ofs;
#else
    std::istream& in = std::cin;
    std::ostream& out = std::cout;
#endif

    init();

    const auto tc = TestCase(in);

    dump(timer.elapsed_ms());

    DistanceTSP dtsp(tc);
    dtsp.two_opt_annealing(100000000);

    dump(timer.elapsed_ms());

    auto path = dtsp.path;
    // 価値の低い餌を取る必要はない
    std::vector<bool> skip_flag(path.size(), false);
    auto best_res = tc.evaluate(path);

    dump(best_res.highest_score);

    while (timer.elapsed_ms() < 9900) {
        auto inner_res = best_res;
        int skip_id = -1;
        for (int id = 1; id < path.size(); id++) {
            if (skip_flag[id]) continue;
            skip_flag[id] = true;
            auto res = tc.evaluate(path, skip_flag);
            if (inner_res.highest_score < res.highest_score) {
                inner_res = res;
                skip_id = id;
            }
            skip_flag[id] = false;
        }
        if (skip_id == -1) break;
        best_res = inner_res;
        skip_flag[skip_id] = true;
    }

    dump(best_res.highest_score);

    auto ans = tc.cvt(path, skip_flag);
    for (int i = best_res.turn_to_truncate + 1; i < K; i++) ans[i] = '-';

    out << ans << std::endl;

    dump(timer.elapsed_ms());

    return 0;
}