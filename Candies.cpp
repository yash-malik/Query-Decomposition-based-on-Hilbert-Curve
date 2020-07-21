// Submission to https://codingcompetitions.withgoogle.com/kickstart/round/000000000019ff43/0000000000337b4d
// based on Query Decomposition using Hilbert Curve in MO's Algorithm


// #pragma comment(linker, "/stack:200000000")
// #pragma GCC optimize("O3")
// #pragma GCC optimize("Ofast")
// #pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
// #pragma GCC optimize("unroll-loops")

#include<bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#ifdef on_linux
#include <sys/resource.h>
#endif

using namespace std;
using namespace __gnu_pbds;

#ifdef dobby_is_a_free_elf
#define TRACE
#endif
#ifdef TRACE
#define trace(...) __f(#__VA_ARGS__, __VA_ARGS__)
template <typename Arg1>
void __f(const char* name, Arg1&& arg1) {
	cout << name << " : " << arg1 << std::endl;
}
template <typename Arg1, typename... Args>
void __f(const char* names, Arg1&& arg1, Args&&... args) {
	const char* comma = strchr(names + 1, ','); cout.write(names, comma - names) << " : " << arg1 << " | "; __f(comma + 1, args...);
}
#else
#define trace(...) 0;
#endif

#define sz(x)                   (int)((x).size())
#define all(x)                  (x).begin(), (x).end()
#define mem0(x)                 memset(x, 0, sizeof (x))
#define mem1(x)                 memset(x, -1, sizeof (x))
#define sleep(x)                std::this_thread::sleep_for(std::chrono::milliseconds(x))
#define DECIMAL(n)              std::cout << std::fixed << std::setprecision(n);
#define rep(i,a,b)              for(int32_t i = (a); i < (b); ++i)
#define between(i,x,y)          ((i) >= (x) && (i) <= (y))
// #define clamp(i,x,y)            (((i) < (x)) ? (x) : ((y) < (i)) ? (y) : (i)); assert((x) <= (y));
#define sqr(a)                  ((a) * (a))
#define pii                     pair<int32_t, int32_t>
#define pll                     pair<long long, long long>
#define mp                      make_pair
#define fi                      first
#define sc                      second
#define pb                      push_back
#define ppb                     pop_back
#define pf                      push_front
#define ppf                     pop_front
#define lbd                     lower_bound
#define ubd                     upper_bound
#define int                     long long
#define ll                      long long

template<typename T, typename U> istream& operator>>(istream& in, pair<T, U> &a) {in >> a.first >> a.second; return in;}
template<typename T, typename U> ostream& operator<<(ostream& out, pair<T, U> a) {out << a.first << " " << a.second; return out;}
template<typename T, typename U> static inline void remax(T &x, U y) {if (x < y) {x = y;}}
template<typename T, typename U> static inline void remin(T &x, U y) {if (y < x) {x = y;}}

static mt19937_64 gen(chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count());

inline long long toint(const string &s) {stringstream ss; ss << s; long long x; ss >> x; return x;}
inline string tostring(long long number) {stringstream ss; ss << number; return ss.str();}
inline string tobin(long long x) {return bitset<63>(x).to_string();}

void alllower(string &u) {transform(u.begin(), u.end(), u.begin(), ::tolower); return;}
void allupper(string &u) {transform(u.begin(), u.end(), u.begin(), ::toupper); return;}

template <typename T>
using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
// find_by_order(k)  returns iterator to kth element starting from 0;
// order_of_key(k) returns count of elements strictly smaller than k;
// erase, insert same as normal set

const long double PI = 3.14159265358979323846264338;
const long double eps = 1e-10;
const long long fftmod = 998244353;
const long long MOD = 1000000007;
const long long INF = 1e18;
const int32_t N = 1e5 + 10;
const int32_t NN = 1e6 + 10;

inline int64_t hilbertOrder(int x, int y, int pow, int rotate) {
	if (pow == 0) {
		return 0;
	}
	int hpow = 1 << (pow - 1);
	int seg = (x < hpow) ? ((y < hpow) ? 0 : 3) : ((y < hpow) ? 1 : 2);
	seg = (seg + rotate) & 3;
	const int rotateDelta[4] = {3, 0, 0, 1};
	int nx = x & (x ^ hpow), ny = y & (y ^ hpow);
	int nrot = (rotate + rotateDelta[seg]) & 3;
	int64_t subSquareSize = int64_t(1) << (2 * pow - 2);
	int64_t ans = seg * subSquareSize;
	int64_t add = hilbertOrder(nx, ny, pow - 1, nrot);
	ans += (seg == 1 || seg == 2) ? add : (subSquareSize - add - 1);
	return ans;
}

const int block = cbrtl(1.0 * N * N);
struct Query {
	int t, l, r, idx;
	pair<int, ll> ord;

	Query(int a, int b, int c) {
		t = a, l = b, r = c;
		calcOrder();
	}
	inline void calcOrder() {
		int tblock  = t / block;
		int lblock = ((tblock) % 4 <= 1) ? -(l / block) : +(l / block);
		int rblock = ((lblock) & 1) ? -r : +r;
		ord = {tblock, hilbertOrder(l, r, 18, 0)};
	}
};

inline bool operator<(const Query &a, const Query &b) {
	return a.ord < b.ord;
}


inline void dobbysolver(int testcase)
{
	int n, q;
	cin >> n >> q;
	int a[n];
	for (int i = 0; i < n; ++i)
		cin >> a[i];

	vector<Query> qry;
	vector<pii> upd;
	int tt = -1;				// update query ctr
	for (int i = 0; i < q; ++i) {
		char type;
		cin >> type;
		if (type == 'Q') {
			int l, r;
			cin >> l >> r;
			l--, r--;
			qry.emplace_back(tt, l, r);
		} else {
			int pos, val;
			cin >> pos >> val;
			pos--;
			tt++;
			upd.pb({pos, val});
		}
	}

	sort(all(qry));

	int l = 1, r = 0, t = -1;
	int ans = 0, oddans = 0, evenans = 0, oddsum = 0, evensum = 0;

	for (Query Q : qry) {
		while (t < Q.t) {
			t++;
			int pos = upd[t].fi;
			int nval = upd[t].sc;
			upd[t].sc = a[pos];
			if (between(pos, l, r)) {
				int c = (pos - l + 1);
				if (c & 1) {
					oddans -= c * a[pos];
					oddans += c * nval;
					oddsum -= a[pos];
					oddsum += nval;
				} else {
					evenans += c * a[pos];
					evenans -= c * nval;
					evensum -= a[pos];
					evensum += nval;
				}
			}
			a[pos] = nval;
		}
		while (t > Q.t) {
			int pos = upd[t].fi;
			int nval = upd[t].sc;
			upd[t].sc = a[pos];
			if (between(pos, l, r)) {
				int c = (pos - l + 1);
				if (c & 1) {
					oddans -= c * a[pos];
					oddans += c * nval;
					oddsum -= a[pos];
					oddsum += nval;
				} else {
					evenans += c * a[pos];
					evenans -= c * nval;
					evensum -= a[pos];
					evensum += nval;
				}
			}
			a[pos] = nval;
			t--;
		}
		while (l > Q.l) {
			l--;
			swap(oddans, evenans);
			oddans *= -1;
			evenans *= -1;
			swap(oddsum, evensum);
			oddans += oddsum;
			evenans -= evensum;

			oddans += a[l];
			oddsum += a[l];
		}
		while (r < Q.r) {
			r++;
			int c = (r - l + 1);
			if (c & 1) {
				oddans += c * a[r];
				oddsum += a[r];
			} else {
				evenans -= c * a[r];
				evensum += a[r];
			}
		}
		while (l < Q.l) {
			swap(oddans, evenans);
			oddans *= -1;
			evenans *= -1;
			swap(oddsum, evensum);
			oddans -= oddsum;
			evenans += evensum;

			evensum -= a[l];
			l++;
		}
		while (r > Q.r) {
			int c = (r - l + 1);
			if (c & 1) {
				oddans -= c * a[r];
				oddsum -= a[r];
			} else {
				evenans += c * a[r];
				evensum -= a[r];
			}
			r--;
		}
		ans += (oddans + evenans);
	}
	cout << ans;
	return;
}

int32_t main()
{
#ifdef on_linux
	rlimit cpu_time {.rlim_cur = 2, .rlim_max = RLIM_INFINITY};
	setrlimit(RLIMIT_CPU, &cpu_time);
#endif
	ios_base::sync_with_stdio(false);
	cin.tie(nullptr);
	cout.tie(nullptr);
	int TESTS = 1;
	cin >> TESTS;
	for (int i = 1; i <= TESTS; ++i) {
		cout << "Case #" << i << ": ";
		dobbysolver(i);
		if (i != TESTS)	cout << '\n';
	}
	return 0;
}
