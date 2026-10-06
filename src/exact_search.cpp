#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// Exact interval DP. N <= 2047 and 1 <= k <= 9 keep even the
// total cost of a chain below 2^121, safely inside signed __int128.
using I = __int128_t;

std::string dec(I x) {
    if (x == 0) return "0";
    bool negative = x < 0;
    if (negative) x = -x;
    std::string s;
    while (x) { s += char('0' + x % 10); x /= 10; }
    if (negative) s += '-';
    std::reverse(s.begin(), s.end());
    return s;
}

int main(int argc, char** argv) {
    // argv[1]: maximum interval size; argv[2]: output JSONL path.
    // One run computes every j <= N for k = 1,...,9.
    if (argc > 3) {
        std::cerr << "Usage: exact_search [N] [output.jsonl]\n";
        return 1;
    }
    try {
    const int N = argc > 1 ? std::stoi(argv[1]) : 2047;
    const std::string path = argc > 2 ? argv[2] : "exact_results.jsonl";
    if (N < 1 || N > 2047) { std::cerr << "N must be in [1, 2047]\n"; return 1; }
    std::ofstream out(path);
    if (!out) { std::cerr << "Cannot create output\n"; return 1; }
    const int stride = N + 2;
    const size_t cells = size_t(stride) * stride;
    auto id = [stride](int l, int r) { return size_t(l) * stride + r; };
    const I inf = I(1) << 126;
    int binom[10][10] = {};
    for (int a = 0; a <= 9; ++a) {
        binom[a][0] = binom[a][a] = 1;
        for (int b = 1; b < a; ++b) binom[a][b] = binom[a-1][b-1] + binom[a-1][b];
    }
    for (int k = 1; k <= 9; ++k) {
        std::vector<I> cost(cells), baseline(cells), powers(N+2, 1);
        std::vector<int> counts(cells);
        std::vector<uint16_t> root(cells);
        for (int q = 1; q <= N; ++q)
            for (int r = 0; r < k; ++r) powers[q] *= q;
        // Empty intervals have cost zero. Process shorter intervals first.
        for (int length = 1; length <= N; ++length) {
            for (int l = 1; l + length - 1 <= N; ++l) {
                const int r = l + length - 1;
                const size_t pos = id(l, r);
                I best = inf;
                int best_count = 0, best_q = 0;
                for (int q = l; q <= r; ++q) {
                    // Every target in this interval pays q^k for the current query.
                    I candidate = length * powers[q] + cost[id(l,q-1)] + cost[id(q+1,r)];
                    // Ties: fewer total queries, then the smallest root.
                    // Ascending q order resolves the last tie automatically.
                    if (candidate <= best) {
                        int count = length + counts[id(l,q-1)] + counts[id(q+1,r)];
                        if (candidate < best || count < best_count) {
                            best = candidate; best_count = count; best_q = q;
                        }
                    }
                }
                cost[pos] = best; counts[pos] = best_count; root[pos] = best_q;
                // Baseline: the lower midpoint for even-sized intervals.
                const int mid = (l+r)/2;
                baseline[pos] = length*powers[mid] + baseline[id(l,mid-1)] + baseline[id(mid+1,r)];
                assert(cost[pos] <= baseline[pos]);
            }
        }
        // Recover moments of the selected optimal tree for each prefix.
        for (int j = 1; j <= N; ++j) {
            std::vector<I> moments(k+1);
            std::vector<std::pair<int,int>> pending{{1,j}};
            while (!pending.empty()) {
                auto [l,r] = pending.back(); pending.pop_back();
                if (l > r) continue;
                const int q = root[id(l,r)], length = r-l+1;
                I power = 1;
                for (int s = 0; s <= k; ++s) {
                    moments[s] += length * power;
                    power *= j-q;
                }
                pending.emplace_back(l,q-1);
                pending.emplace_back(q+1,r);
            }
            I check = 0;
            for (int s = 0; s <= k; ++s) {
                I power = 1;
                for (int r = 0; r < k-s; ++r) power *= j;
                I term = binom[k][s] * moments[s] * power;
                check += (s%2 ? -term : term);
            }
            assert(check == cost[id(1,j)]);
            assert(moments[0] == counts[id(1,j)]);
            out << "{\"j\":" << j << ",\"k\":" << k
                << ",\"optimal_sum\":\"" << dec(cost[id(1,j)])
                << "\",\"binary_sum\":\"" << dec(baseline[id(1,j)])
                << "\",\"root\":" << root[id(1,j)] << ",\"moments\":[";
            for (int s = 0; s <= k; ++s) {
                if (s) out << ',';
                out << '"' << dec(moments[s]) << '"';
            }
            out << "]}\n";
        }
        out.flush();
        std::cerr << "Completed k=" << k << ", j=1.." << N << '\n';
    }
    if (!out) {
        std::cerr << "Output write failed\n";
        return 1;
    }
    return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
