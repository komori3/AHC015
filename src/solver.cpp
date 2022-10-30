#define _CRT_NONSTDC_NO_WARNINGS
#define _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING
#include <bits/stdc++.h>
#include <random>
#include <unordered_set>
#include <array>
#include <optional>
#ifdef _MSC_VER
#include <opencv2/core.hpp>
#include <opencv2/core/utils/logger.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <conio.h>
#include <ppl.h>
#include <filesystem>
#include <intrin.h>
//#include <boost/multiprecision/cpp_int.hpp>
int __builtin_clz(unsigned int n)
{
    unsigned long index;
    _BitScanReverse(&index, n);
    return 31 - index;
}
int __builtin_ctz(unsigned int n)
{
    unsigned long index;
    _BitScanForward(&index, n);
    return index;
}
namespace std {
    inline int __lg(int __n) { return sizeof(int) * 8 - 1 - __builtin_clz(__n); }
}
//using __uint128_t = boost::multiprecision::uint128_t;
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
std::istream& operator>>(std::istream& is, std::pair<S, T>& p) {
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

std::ostream& operator<<(std::ostream& os, const std::vector<bool>& v) {
    std::string s(v.size(), ' ');
    for (int i = 0; i < (int)v.size(); i++) s[i] = v[i] + '0';
    os << s;
    return os;
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
#define DUMPOUT std::cerr
std::ostringstream DUMPBUF;
#define dump(...) do{DUMPBUF<<"  ";DUMPBUF<<#__VA_ARGS__<<" :[DUMP - "<<__LINE__<<":"<<__FUNCTION__<<"]"<<std::endl;DUMPBUF<<"    ";dump_func(__VA_ARGS__);DUMPOUT<<DUMPBUF.str();DUMPBUF.str("");DUMPBUF.clear();}while(0);
void dump_func() { DUMPBUF << std::endl; }
template <class Head, class... Tail> void dump_func(Head&& head, Tail&&... tail) { DUMPBUF << head; if (sizeof...(Tail) == 0) { DUMPBUF << " "; } else { DUMPBUF << ", "; } dump_func(std::move(tail)...); }

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
    double elapsed_ms() const { return (time() - t - paused) * 1000.0; }
};

/* rand */
struct Xorshift {
    static constexpr uint64_t M = INT_MAX;
    static constexpr double e = 1.0 / M;
    uint64_t x = 88172645463325252LL;
    Xorshift() {}
    Xorshift(uint64_t seed) { reseed(seed); }
    inline void reseed(uint64_t seed) { x = 0x498b3bc5 ^ seed; for (int i = 0; i < 20; i++) next(); }
    inline uint64_t next() { x = x ^ (x << 7); return x = x ^ (x >> 9); }
    inline int next_int() { return next() & M; }
    inline int next_int(int mod) { return next() % mod; }
    inline int next_int(int l, int r) { return l + next_int(r - l + 1); }
    inline double next_double() { return next_int() * e; }
};

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

template<typename T, typename ...Args> auto make_vector(T x, int arg, Args ...args) { if constexpr (sizeof...(args) == 0)return std::vector<T>(arg, x); else return std::vector(arg, make_vector<T>(x, args...)); }
template<typename T> bool chmax(T& a, const T& b) { if (a < b) { a = b; return true; } return false; }
template<typename T> bool chmin(T& a, const T& b) { if (a > b) { a = b; return true; } return false; }

using ll = long long;
using ld = double;
//using ld = boost::multiprecision::cpp_bin_float_quad;
using pii = std::pair<int, int>;
using pll = std::pair<ll, ll>;

using std::cin, std::cout, std::cerr, std::endl, std::string, std::vector, std::array;



constexpr int T = 100;
constexpr int N = 10;
using Board = int[N][N];

constexpr int di[] = { 0, -1, 0, 1 };
constexpr int dj[] = { 1, 0, -1, 0 };
const string d2c = "RFLB";
char c2d[256];

void init() {
    c2d['R'] = 0;
    c2d['F'] = 1;
    c2d['L'] = 2;
    c2d['B'] = 3;
}

struct State {

    int fs[T] = {};
    int den = 0;
    Board board = {};
    int t = 0;

    State() {}

    void init(std::istream& in) {
        in >> fs;
        int ctr[4] = {};
        for (int f : fs) ctr[f]++;
        den = ctr[1] * ctr[1] + ctr[2] * ctr[2] + ctr[3] * ctr[3];
    }

    void load(int p) {
        int q = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (!board[i][j]) {
                    q++;
                    if (q == p) {
                        board[i][j] = fs[t];
                        t++;
                        return;
                    }
                }
            }
        }
        assert(false);
    }

    void load_random(Xorshift& rnd) {
        vector<int> perm(T);
        std::iota(perm.begin(), perm.end(), 0);
        shuffle_vector(perm, rnd);
        for (int p : perm) {
            int i = p / N, j = p % N;
            if (!board[i][j]) {
                board[i][j] = fs[t];
                t++;
                return;
            }
        }
    }

    void apply_move(Board& board, char c) const {
        if (c == 'L') {
            for (int i = 0; i < N; i++) {
                int k = 0;
                for (int j = 0; j < N; j++) {
                    if (board[i][j]) {
                        board[i][k] = board[i][j];
                        if (k != j) board[i][j] = 0;
                        k++;
                    }
                }
            }
        }
        else if (c == 'R') {
            for (int i = 0; i < N; i++) {
                int k = N - 1;
                for (int j = N - 1; j >= 0; j--) {
                    if (board[i][j]) {
                        board[i][k] = board[i][j];
                        if (k != j) board[i][j] = 0;
                        k--;
                    }
                }
            }
        }
        else if (c == 'F') {
            for (int j = 0; j < N; j++) {
                int k = 0;
                for (int i = 0; i < N; i++) {
                    if (board[i][j]) {
                        board[k][j] = board[i][j];
                        if (k != i) board[i][j] = 0;
                        k++;
                    }
                }
            }
        }
        else {
            for (int j = 0; j < N; j++) {
                int k = N - 1;
                for (int i = N - 1; i >= 0; i--) {
                    if (board[i][j]) {
                        board[k][j] = board[i][j];
                        if (k != i) board[i][j] = 0;
                        k--;
                    }
                }
            }
        }
    }

    void query(std::ostream& out, char c) {
        out << c << endl;
        apply_move(board, c);
    }

    void query_greedy(std::ostream& out) {
        int best_score = -1;
        char best_dir = ' ';
        for (char c : d2c) {
            Board cboard;
            std::memcpy(cboard, board, sizeof(int) * N * N);
            apply_move(cboard, c);
            if (chmax(best_score, compute_score(cboard))) {
                best_dir = c;
            }
        }
        query(out, best_dir);
    }

    int simulate(Xorshift& rnd, const string& ops) const {
        State cstate(*this);
        for (char c : ops) {
            cstate.load_random(rnd);
            cstate.apply_move(cstate.board, c);
            if (cstate.t == T) break;
        }
        return compute_score(cstate.board);
    }

    void query_simulate(std::ostream& out, Xorshift& rnd, int len, double duration) {
        Timer timer;
        double start_time = timer.elapsed_ms(), end_time = start_time + duration;

        int dir_ctr[4] = {};
        int dir_sum[4] = {};

        string ops(len, ' ');
        int loop = 0;
        while (timer.elapsed_ms() < end_time) {
            for (int i = 0; i < len; i++) ops[i] = d2c[rnd.next_int(4)];
            int dir = c2d[ops[0]];
            dir_ctr[dir]++;
            dir_sum[dir] += simulate(rnd, ops);
            loop++;
        }
        
        int best_dir = -1;
        double best_avg = -1;
        for (int d = 0; d < 4; d++) {
            if (chmax(best_avg, (double)dir_sum[d] / dir_ctr[d])) {
                best_dir = d;
            }
        }
        query(out, d2c[best_dir]);
        dump(t, loop, compute_score(board));
    }

    int compute_score(const Board& board) const {
        int s2 = 0;
        bool used[N][N] = {};
        for (int si = 0; si < N; si++) {
            for (int sj = 0; sj < N; sj++) {
                if (!board[si][sj] || used[si][sj]) continue;
                int s = 0, v = board[si][sj];
                std::queue<pii> qu;
                qu.emplace(si, sj);
                used[si][sj] = true;
                s++;
                while (!qu.empty()) {
                    auto [i, j] = qu.front(); qu.pop();
                    for (int d = 0; d < 4; d++) {
                        int ni = i + di[d], nj = j + dj[d];
                        if (ni < 0 || ni >= N || nj < 0 || nj >= N || used[ni][nj] || board[ni][nj] != v) continue;
                        qu.emplace(ni, nj);
                        used[ni][nj] = true;
                        s++;
                    }
                }
                s2 += s * s;
            }
        }
        return (int)round(1e6 * s2 / den);
    }

};

int solve(std::istream& in, std::ostream& out) {
    Xorshift rnd;

    State state;
    state.init(in);
    for (int t = 0; t < T; t++) {
        int p;
        in >> p;
        state.load(p);
        state.query_simulate(out, rnd, 10, 180);
    }
    return state.compute_score(state.board);
}

#ifdef _MSC_VER
// マルチテストケース実行
void batch_test(int seed_begin = 0, int num_seed = 100, int step = 1) {

    constexpr int batch_size = 5;
    const int block_size = batch_size * step;
    int seed_end = seed_begin + num_seed;

    vector<int> scores(num_seed);

    concurrency::critical_section mtx;
    for (int batch_begin = seed_begin; batch_begin < seed_end; batch_begin += block_size) {
        int batch_end = std::min(batch_begin + block_size, seed_end);
        vector<int> seeds;
        for (int seed = batch_begin; seed < batch_end; seed += step) seeds.push_back(seed);
        concurrency::parallel_for_each(seeds.begin(), seeds.end(), [&mtx, &scores](int seed) {
            std::ifstream ifs(format("tools_win/in/%04d.txt", seed));
            std::istream& in = ifs;
            std::ofstream ofs(format("tools_win/out/%04d.txt", seed));
            std::ostream& out = ofs;

            int score = solve(in, out);

            ifs.close();
            ofs.close();

            {
                mtx.lock();
                scores[seed] = score;
                cerr << format("seed=%d, score=%d\n", seed, scores[seed]);
                mtx.unlock();
            }
        });
    }

    dump(std::accumulate(scores.begin(), scores.end(), 0) * 2);
}
#endif



int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

#ifdef _MSC_VER
    std::ifstream ifs(R"(tools_win\in\0000.txt)");
    std::ofstream ofs(R"(tools_win\out\0000.txt)");
    std::istream& in = ifs;
    std::ostream& out = ofs;
#else
    std::istream& in = cin;
    std::ostream& out = cout;
#endif

    init();

#if 1
    int score = solve(in, out);
    dump(score);
#else
    batch_test();
#endif

    return 0;
}