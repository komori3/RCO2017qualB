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

    int sr, sc;
    std::vector<std::string> board;
    int N;
    std::vector<Food> foods;
    std::vector<std::vector<Food*>> food_map;

    TestCase(std::istream& in) {
        { int buf; in >> buf >> buf >> buf; }
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
        }
        for (int i = 0; i < N; i++) {
            auto [id, fr, fc, F, D] = foods[i];
            food_map[fr][fc] = &foods[i];
        }
    }

    struct Result {
        int score;
        int best_score;
        int best_turn;
    };

    Result evaluate(const std::string& cmds) const {
        Result res = {0, 0, -1};
        assert(cmds.size() == K);
        int turn = -1;
        int r = sr, c = sc;
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
                if (res.best_score < res.score) {
                    res.best_score = res.score;
                    res.best_turn = turn;
                }
            }
            r = nr; c = nc;
        }
        return res;
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

    std::string ans;
    std::string dirs = "RULD";
    for (int i = 0; i < K; i++) ans += dirs[rnd.next_int(4)];

    auto [score, best_score, best_turn] = tc.evaluate(ans);
    dump(score, best_score, best_turn);

    // ans の best_turn 文字目までは採用
    // -> best_turn + 1 文字目以降は全て '-'
    for (int i = best_turn + 1; i < K; i++) ans[i] = '-';

    out << ans << std::endl;

    return 0;
}