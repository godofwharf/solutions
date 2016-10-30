/*
Estimating the vertex connectivity of a directed graph by finding maximum flow using Dinic's algorithm

Solution Complexity: O(|V|^2 E)
*/

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <cctype>
#include <cstdio>
#include <string>
#include <vector>
#include <cmath>
#include <queue>
#include <stack>
#include <set>
#include <map>
using namespace std;

const int MAXN = 400;
int cap[MAXN][MAXN];
int deg[MAXN];
int adj[MAXN][MAXN];
int q[MAXN];
int prev[MAXN];
int us[MAXN];
int vs[MAXN];
const int inf = 9999999;

int dinic(int n, int s, int t) {
    int flow = 0;
    while (true) {
        for (int z = 0; z < n; z++)
            prev[z] = -1;
        int qf = 0, qb = 0;
        prev[q[qb++] = s] = -2;
        while (qb > qf && prev[t] == -1)
            for (int u = q[qf++], i = 0, v; i < deg[u]; i++)
                if (prev[v = adj[u][i]] == -1 && cap[u][v])
                    prev[q[qb++] = v] = u;
        if (prev[t] < 0)
            break;
        for (int z = 0; z < n; z++)
            if (cap[z][t] && prev[z] != -1) {
                int bot = cap[z][t];
                for (int v = z, u = prev[v]; u >= 0; v = u, u = prev[v])
                    bot = min(bot, cap[u][v]);
                if (!bot)
                    continue;
                cap[z][t] -= bot, cap[t][z] += bot;
                for (int v = z, u = prev[v]; u >= 0; v = u, u = prev[v])
                    cap[u][v] -= bot, cap[v][u] += bot;
                flow += bot;
            }
    }

    return flow;
}

int main() {

    int t ;
    cin >> t;
    while (t--) {
        int n, m;

        scanf("%d %d",&n, &m);

        for(int i = 0; i < m; i++) {
            scanf("%d %d", &us[i], &vs[i]);
        }
        int ans = inf;
        for(int s = 0; s < n; s++)
            for(int t = s + 1; t < n; t++)
            {
                for(int i = 0; i < 2 * n; i++) {
                    deg[i] = 0;
                    for(int j = 0; j < 2 * n; j++)
                        cap[i][j] = adj[i][j] = 0;
                }
                for(int i = 0; i < m; i++)
                    cap[us[i] + n][vs[i]] = cap[vs[i] + n][us[i]] = inf;
                for(int i = 0; i < n; i++)
                    cap[i][i + n] = 1;
                cap[s][s + n] = cap[t][t + n] = inf;
                for(int i = 0; i < 2 * n; i++)
                    for(int j = 0; j < 2 * n; j++)
                        if (cap[i][j])
                            adj[i][deg[i]++] = j;
                int cur = dinic(2 * n, s, t + n);
                ans = min(ans, cur);
            }
        if(ans == inf)
            ans = n;
        printf("%d\n", ans);
    }
    return 0;
}
