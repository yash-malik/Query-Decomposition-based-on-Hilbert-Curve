// Submission to https://codingcompetitions.withgoogle.com/kickstart/round/000000000019ff43/0000000000337b4d
// based on Query Decomposition using Hilbert Curve in MO's Algorithm


#include<bits/stdc++.h>

using namespace std;

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

const int32_t N = 1e5 + 10;
const int block = cbrtl(1.0 * N * N);

struct Query {
	int t, l, r, idx;
	pair<int, int64_t> ord;

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

inline bool operator<(const Query & a, const Query & b) {
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
	vector<pair<int, int>> upd;
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
			upd.push_back({pos, val});
		}
	}

	sort(qry.begin(), qry.end());

	int l = 1, r = 0, t = -1;
	int64_t ans = 0, oddans = 0, evenans = 0, oddsum = 0, evensum = 0;

	for (Query Q : qry) {
		while (t < Q.t) {
			t++;
			int pos = upd[t].first;
			int nval = upd[t].second;
			upd[t].second = a[pos];
			if (l <= pos && pos <= r) {
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
			int pos = upd[t].first;
			int nval = upd[t].second;
			upd[t].second = a[pos];
			if (l <= pos && pos <= r) {
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

int32_t main() {
	ios_base::sync_with_stdio(false);
	cin.tie(nullptr); cout.tie(nullptr);
	int TESTS = 1;
	cin >> TESTS;
	for (int i = 1; i <= TESTS; ++i) {
		cout << "Case #" << i << ": ";
		dobbysolver(i);
		if (i != TESTS)	cout << '\n';
	}
	return 0;
}
