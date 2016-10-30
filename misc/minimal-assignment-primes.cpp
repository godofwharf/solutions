/* Solution to a problem which I conceived myself.
Walt recently discovered his passion for Mathematics. He is particularly interested in prime numbers. 
He has become so obsessed that he wants every batch to contain exactly p pounds of crystal meth, where p is a prime. 
Jesse has already cooked N batches. All the N batches contain some integer amount of pounds of crystal meth. 
In one move, Walt can transfer 1 pound from one batch to another batch. He wants to figure out the minimum number of moves 
required so as to make all the batches contain prime amount of pounds of crystal meth. He knows that it is not possible always. 
Will you help him?

Input:
The first line contains the number of testcases T.
1 <= T <= 100

Each testcase contains an integer N - no. of batches
1 <= N <= 50

The next line contains N space separated integers which correspond to amounts of meth ni in pounds in the i-th batch.
1 <= ni <= 200

The sum of amounts in all the batches won’t exceed 200.

Output:
Print -1 if it is impossible. 
Otherwise print the minimum number of moves necessary to make all batches contain prime amount of pounds of crystal meth

Sample Input:
2
3
12 19 11
3
4 5 4

Sample Output:
18
-1

Solution Complexity: O(SlgS(lglgS) + L * N^3 + S * lg(S/lgS)) 
S - Sum of the amounts in all batches
L - Number of ways to paritition S into a set of distinct primes

#include <iostream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdio>
#define MAXS 500
#define MAXN 50
using namespace std;

const int inf=9999999;
bool isprime[MAXS+1];
vector<vector<int> > comb;
vector<int> primes;

typedef vector<double> VD;
typedef vector<VD> VVD;
typedef vector<int> VI;

void sieve()
{
    memset(isprime, true, sizeof(isprime));
    isprime[0] = false;
    isprime[1] = false;
    for(int i = 2; i * i <= MAXS; i++)
    {
        if(isprime[i])
        {
            for(int j = 2 * i; j <= MAXS; j += i)
                isprime[j]=false;
        }
    }
    for(int i = 0; i <= MAXS; i++)
    {
        if(isprime[i] && i != 2)
            primes.push_back(i);
    }
}

void findCombination(int e, int S, int N, vector<int> cur)
{
    if(N == 0)
    {
        comb.push_back(cur);
        return;
    }
    int a, b;
    int avg = (int) ceil((float) S / N);
    a = lower_bound(primes.begin(), primes.begin() + e, avg) - primes.begin();// returns index of first element in range whose value is >=avg

    int sum = 0;
    for(int j = 0; j < N - 1; j++)
        sum += primes[j];
    int x = S - sum;

    b = upper_bound(primes.begin(), primes.begin() + e, x) - primes.begin() - 1;// returns index of first element in range whose value is > x

    for(int i = a; i <= b; i++)
    {
        cur.push_back(primes[i]);
        findCombination(i, S - primes[i], N - 1, cur);
        cur.pop_back();
    }
}

double MinCostMatching(const VVD &cost, VI &Lmate, VI &Rmate) {
    int n = int(cost.size());

    // construct dual feasible solution
    VD u(n);
    VD v(n);
    for (int i = 0; i < n; i++) {
        u[i] = cost[i][0];
        for (int j = 1; j < n; j++) u[i] = min(u[i], cost[i][j]);
    }
    for (int j = 0; j < n; j++) {
        v[j] = cost[0][j] - u[0];
        for (int i = 1; i < n; i++) v[j] = min(v[j], cost[i][j] - u[i]);
    }

    // construct primal solution satisfying complementary slackness
    Lmate = VI(n, -1);
    Rmate = VI(n, -1);
    int mated = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (Rmate[j] != -1) continue;
            if (fabs(cost[i][j] - u[i] - v[j]) < 1e-10) {
                Lmate[i] = j;
                Rmate[j] = i;
                mated++;
                break;
            }
        }
    }

    VD dist(n);
    VI dad(n);
    VI seen(n);

    // repeat until primal solution is feasible
    while (mated < n) {

        // find an unmatched left node
        int s = 0;
        while (Lmate[s] != -1) s++;

        // initialize Dijkstra
        fill(dad.begin(), dad.end(), -1);
        fill(seen.begin(), seen.end(), 0);
        for (int k = 0; k < n; k++)
            dist[k] = cost[s][k] - u[s] - v[k];

        int j = 0;
        while (true) {
            // find closest
            j = -1;
            for (int k = 0; k < n; k++) {
                if (seen[k]) continue;
                if (j == -1 || dist[k] < dist[j]) j = k;
            }
            seen[j] = 1;

            // termination condition
            if (Rmate[j] == -1) break;

            // relax neighbors
            const int i = Rmate[j];
            for (int k = 0; k < n; k++) {
                if (seen[k]) continue;
                const double new_dist = dist[j] + cost[i][k] - u[i] - v[k];
                if (dist[k] > new_dist) {
                    dist[k] = new_dist;
                    dad[k] = j;
                }
            }
        }

        // update dual variables
        for (int k = 0; k < n; k++) {
            if (k == j || !seen[k]) continue;
            const int i = Rmate[k];
            v[k] += dist[k] - dist[j];
            u[i] -= dist[k] - dist[j];
        }
        u[s] += dist[j];

        // augment along path
        while (dad[j] >= 0) {
            const int d = dad[j];
            Rmate[j] = Rmate[d];
            Lmate[Rmate[j]] = j;
            j = d;
        }
        Rmate[j] = s;
        Lmate[s] = j;

        mated++;
    }

    double value = 0;
    for (int i = 0; i < n; i++)
        value += cost[i][Lmate[i]];

    return value;
}

int main()
{
    int t;
    cin >> t;
    sieve();
    while(t--)
    {
        int S = 0, N, temp;
        cin >> N;
        vector<int> amounts;
        for(int i = 0; i < N; i++)
        {
            scanf("%d",&temp);
            amounts.push_back(temp);
            S += temp;
        }

        int two = 0;
        vector<int> cur;
        if(S % 2 == 0 && N % 2 == 1)
        {
            two = 1;
            cur.push_back(2);
        }
        findCombination(primes.size(), S - two * 2, N - two, cur);
        int ans = inf;
        for(int i = 0; i < comb.size(); i++) {
            vector<vector<double> > cost;
            for(int j = 0; j < comb[0].size(); j++) {
                vector<double> v;
                for(int k = 0; k < amounts.size(); k++) v.push_back(abs(amounts[k] - comb[i][j]));
                cost.push_back(v);
            }
            vector<int> lpart;
            vector<int> rpart;
            lpart = comb[i];
            rpart = amounts;
            ans = min(ans, (int) MinCostMatching(cost, lpart, rpart));
        }
        if(ans == inf) ans = -1;
        cout << ans << endl;
        comb.clear();
    }
    return 0;
}
