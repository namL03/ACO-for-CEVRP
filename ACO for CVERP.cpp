#include<bits/stdc++.h>
#include <iostream>
#include <set>
#include <map>
#include <cmath>
#include <random>
#include <iterator>
using namespace std;
using namespace std::chrono;
double optimal_value, energy_consumption;
double p_ = 0.98, alpha = 1, beta = 2, pr = 0.05;
int num_vehicles, num_customers, num_stations, capacity, energy_capacity, num_nodes;
double d[1010][1010], p[1010][1010];
int demand[1010];
vector<int>best_sol;
double best_fitness = 1e18;
struct Node
{
	int x;
	int y;
} node[1010];
//int getRandomNumber(int n) {
//    std::random_device rd;  
//    std::mt19937 gen(rd()); 
//    std::uniform_int_distribution<> distr(1, n); 
//    return distr(gen); 
//}
int getRandomNumber(int n) {
    // Get the current time as a seed value
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    // Initialize the random number generator with the seed
    std::mt19937 gen(seed);

    // Define the distribution range from 1 to n
    std::uniform_int_distribution<> distr(1, n);

    // Generate and return a random number
    return distr(gen);
}
double get_distance(int x, int y)
{
	int xd = node[x].x - node[y].x;
	int yd = node[x].y - node[y].y;
	return sqrt(xd * xd + yd * yd);
}
double fitness_eval(vector<int>v)
{
	double ans = 0;
	for(int i = 0; i < v.size() - 1; i++)
	{
		ans += d[v[i]][v[i + 1]];
	}
	return ans;
}
void update_pheromone(int i, int j, double k)
{
	if(best_sol.size() != 0)
	{
		double p_max = 1.0 / ((1 - p_) * best_fitness); 
		double p_min = p_max * (1 - pow(pr, 1.0 / num_customers)) / 
						((num_customers / 2.0 - 1) * pow(pr, 1.0 / num_customers));
		if(k > p_max) k = p_max;
		if(k < p_min) k = p_min;
	}
	p[i][j] = k;
}
int roulette_wheel_selection(int u, const std::set<int>& s) {
    std::map<int, double> probabilities;
    double sum = 0.0;

    for (int v : s) {
        if (d[u][v] == 0) continue; // Handling division by zero
        probabilities[v] = pow(p[u][v], alpha) / pow(d[u][v], beta);
        sum += probabilities[v];
    }
    
    if (sum == 0.0) return -1; // Handle zero sum case

    for (auto& elem : probabilities) {
        elem.second /= sum;
    }

    // Enhanced seed using time and random_device
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() + std::random_device()();
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(0, 1);
    double random_number = dis(gen);
    double cumulative_probability = 0.0;

    for (auto& elem : probabilities) {
        cumulative_probability += elem.second;
        if (random_number <= cumulative_probability) {
            //std::cout << elem.first << std::endl;
            return elem.first;
        }
    }

    return -1; // Fallback in case of rounding errors
}
vector<int>r, new_r, stations;
set<int>s;
long long V[1010];
int P[1010];
double sum_dis[2010];
int check[2010];
vector<int> route_construction()
{
	r.clear();
	new_r.clear();
	for(int i = 1; i <= num_customers; i++)
	{
		s.insert(i);
	}
//	for(int i:s)
//	{
//		cout << i << '\n';
//	}
	int st = getRandomNumber(num_customers);
	//cout << st << '\n';
	r.push_back(st);
	s.erase(st);
	while(!s.empty())
	{
		int u = r[r.size() - 1];
		int v = roulette_wheel_selection(u, s);
		r.push_back(v);
		s.erase(v);
	}
//	for(int i:r)
//	{
//		cout << i << " ";
//	}
//	cout << '\n';
	int flag = 0;
	for(int i = 0; i < r.size(); i++)
	{
		if(r[i] == 1) flag = 1;
		if(flag) new_r.push_back(r[i]);
	}
	for(int i = 0; i < r.size(); i++)
	{
		if(r[i] == 1) break;
		new_r.push_back(r[i]);
	}
//	for(int i = 0; i < new_r.size(); i++)
//	{
//		cout << new_r[i] << " ";
//	}
//	cout << '\n';
	V[0] = 0;
	for(int i = 1; i < num_customers; i++)
	{
		V[i] = 1e18;
	}
	for(int i = 0; i < num_customers; i++)
	{
		P[i] = 0;
	}
	for(int i = 1; i < num_customers; i++)
	{
		int j = i;
		int tc = 0;
		double td = 0;
		while(j < num_customers && tc <= capacity)
		{
			tc += demand[new_r[j]];
			if(i == j) td = d[1][new_r[j]] + d[new_r[j]][1];
			else td += d[new_r[j - 1]][new_r[j]] + d[new_r[j]][1] - d[new_r[j - 1]][1];
			if(tc <= capacity)
			{
				if(V[i - 1] + td < V[j])
				{
					V[j] = V[i - 1] + td;
					P[j] = i - 1;
				}
				j++;
			}
		}
	}
	vector<int>route;
	int las = num_customers - 1;
	while(las != 0)
	{
		int st = P[las];
		route.push_back(1);
		for(int i = las; i > st; i--)
		{
			route.push_back(new_r[i]);
		}
		las = st;
	}
	route.push_back(1);
//	for(int i = 0; i < route.size(); i++)
//	{
//		cout << route[i] << " ";
//	}
//	cout << '\n';
//	double current_c = 0;
//	cout << capacity << '\n';
//	for(int i = 1; i < route.size(); i++)
//	{
//		current_c += demand[route[i]];
//		cout << i << " " << current_c << '\n';
//		if(route[i] == 1) current_c = 0;
//	}
//	cout << '\n';
	double D = energy_capacity / energy_consumption;
	double ld = D;
	vector<int>final_route;
	for(int i = 0; i < route.size() - 1; i++)
	{
		if(route[i] == 1)
		{
			ld = D;
		}
		final_route.push_back(route[i]);
		double best_dis = 1e18;
		int chosen;
		for(int j = num_customers + 1; j <= num_nodes; j++)
		{
			if(d[route[i]][j] < ld)
			{
				if(d[route[i]][j] + d[j][route[i + 1]] < best_dis)
				{
					best_dis = d[route[i]][j] + d[j][route[i + 1]];
					chosen = j;
				}
			}
		}
		final_route.push_back(chosen);
		ld = D - d[chosen][route[i + 1]];
	}
	final_route.push_back(1);
	route.clear();
	double current_dis = 0, to_be_dis = 0;
	int check_point = 0;
	check[0] = 0;
	for(int i = 1; i < final_route.size(); i++)
	{
		check[i] = 0;
		if(check_point == 0) current_dis += d[final_route[i - 1]][final_route[i]];
		else to_be_dis += d[final_route[i - 1]][final_route[i]];
		if(final_route[i] > num_customers)
		{
			if(check_point == 0) check_point = i;
			else
			{
				if(current_dis + to_be_dis 
				- d[final_route[check_point - 1]][final_route[check_point]]
				- d[final_route[check_point]][final_route[check_point + 1]]
				+ d[final_route[check_point - 1]][final_route[check_point + 1]] < D)
				{
					check[check_point] = 1;
					current_dis = current_dis + to_be_dis
					- d[final_route[check_point - 1]][final_route[check_point]]
					- d[final_route[check_point]][final_route[check_point + 1]]
					+ d[final_route[check_point - 1]][final_route[check_point + 1]];
					to_be_dis = 0;
				}
				else
				{
					current_dis = to_be_dis;
					to_be_dis = 0;
				}
				check_point = i;
			}
		}
		if(final_route[i] == 1)
		{
			if(check_point != 0)
			{
				if(current_dis + to_be_dis 
				- d[final_route[check_point - 1]][final_route[check_point]]
				- d[final_route[check_point]][final_route[check_point + 1]]
				+ d[final_route[check_point - 1]][final_route[check_point + 1]] < D)
				{
					check[check_point] = 1;
				}
			}
			current_dis = 0;
			to_be_dis = 0;
			check_point = 0;
		}
	}
	for(int i = 0; i < final_route.size(); i++)
	{
		if(check[i] == 0)
		{
			route.push_back(final_route[i]);
		}
	}
	//cout << d[5][2] + d[2][30] << '\n'; 
//	for(int i : final_route)
//	{
//		cout << i << " ";
//	}
//	cout << '\n';
//	for(int i : route)
//	{
//		cout << i << " ";
//	}
//	cout << '\n';
//	double current_d = 0, current_c = 0;
//	cout << energy_capacity << '\n';
//	cout << capacity << '\n';
//	for(int i = 1; i < route.size(); i++)
//	{
//		current_d += d[route[i - 1]][route[i]] * energy_consumption;
//		current_c += demand[route[i]];
//		cout << i << " " << current_d << " " << current_c << '\n';
//		if(route[i] == 1 || route[i] > num_customers) current_d = 0;
//		if(route[i] == 1) current_c = 0;
//
//	}
	return route;
}

int main()
{
	freopen("X-n351-k40.evrp", "r", stdin);
	freopen("X-n351-k40.out", "w", stdout);
	ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
	string s;
	for(int i = 1; i <= 3; i++)
	{
		getline(cin, s);
	}
	cin >> s >> optimal_value;
	cin >> s >> num_vehicles;
	cin >> s >> num_customers;
	cin >> s >> num_stations;
	cin >> s >> capacity;
	cin >> s >> energy_capacity;
	cin >> s >> energy_consumption;
	getline(cin, s);
	getline(cin, s);
	getline(cin, s);
	//cout << s;
	num_nodes = num_customers + num_stations;
//	cout << optimal_value << '\n';
//	cout << num_vehicles << '\n';
//	cout << num_customers << '\n';
//	cout << num_stations << '\n';
//	cout << capacity << '\n';
//	cout << energy_capacity << '\n';
//	cout << energy_consumption << '\n';
	for(int i = 1; i <= num_nodes; i++)
	{
		int id, x, y;
		cin >> id >> x >> y;
		node[id].x = x;
		node[id].y = y;
		//cout << id << " " << node[id].x << " " << node[id].y << '\n'; 
	}
	getline(cin, s);
	getline(cin, s);
	for(int i = 1; i <= num_customers; i++)
	{
		int id, x;
		cin >> id >> x;
		demand[id] = x;
		//cout << id << " " << demand[id] << '\n';
	}
	vector<double> all_sol;
	double time_limit;
	if(num_nodes < 100) time_limit = num_nodes * 1.0 / 100;
	else if(num_nodes < 1000) time_limit = num_nodes * 2.0 / 100;
	else time_limit = num_nodes * 3.0 / 100;
	for(int i = 1; i <= num_nodes; i++)
	{
		for(int j = 1; j <= num_nodes; j++)
		{
			d[i][j] = get_distance(i, j);
			//cout << i << " " << j << " " << d[i][j] << '\n';
		}
	}
	for(int num_run = 1; num_run <= 20; num_run++)
	{
		auto start = high_resolution_clock::now();
		best_sol.clear();
		best_fitness = 1e18;
		for(int i = 1; i <= num_customers; i++)
		{
			for(int j = 1; j <= num_customers; j++)
			{
				p[i][j] = 1e-4;
				//cout << i << " " << j << " " << p[i][j] << '\n';
			}
		}
		while (duration_cast<seconds>(high_resolution_clock::now() - start).count() < time_limit * 60)
		{
			vector<int>iter_best;
			double iter_fit = 1e18;
			for(int ant = 1; ant <= num_customers; ant++)
			{
				vector<int>res = route_construction();
				double fit = fitness_eval(res);
				if(fit < iter_fit)
				{
					iter_fit = fit;
					iter_best.clear();
					for(int i:res)
					{
						iter_best.push_back(i);
					}
				}
			}
			if(iter_fit < best_fitness)
			{
				best_fitness = iter_fit;
				best_sol.clear();
				for(int i:iter_best)
				{
					best_sol.push_back(i);
				}
			}
			for(int i = 1; i <= num_customers; i++)
			{
				for(int j = 1; j <= num_customers; j++)
				{
					double new_k = p[i][j] * p_;
					update_pheromone(i, j, new_k);
				}
			}
			int cur = 1;
			for(int i = 1; i < iter_best.size(); i++)
			{
				if(iter_best[i] > num_customers || iter_best[i] == 1) continue;
				double new_k = p[cur][iter_best[i]] + 1 / iter_fit;
				update_pheromone(cur, iter_best[i], new_k);
				cur = iter_best[i]; 
			}
			double new_k = p[cur][1] + 1 / iter_fit;
			update_pheromone(cur, 1, new_k);
		}
		all_sol.push_back(best_fitness);
		cout << best_fitness << '\n';
//		for(int i:best_sol)
//		{
//			cout << i << " ";
//		}
//		cout << '\n';
//		cout << '\n' << capacity << '\n' << energy_capacity << '\n';
//		double current_d = 0;
//		double current_c = 0;
//		for(int i = 1; i < best_sol.size(); i++)
//		{
//			current_d += d[best_sol[i - 1]][best_sol[i]] * energy_consumption;
//			current_c += demand[best_sol[i]];
//			cout << i << " " << current_d << " " << current_c << '\n';
//			if(best_sol[i] == 1 || best_sol[i] > num_customers) current_d = 0;
//			if(best_sol[i] == 1) current_c = 0;
//		}
	}
	double min_sol = 1e18, mean_sol = 0, std_sol = 0;
	for(double sol: all_sol)
	{
		if(sol < min_sol)
		{
			min_sol = sol;
		}
		mean_sol += sol;
	}
	mean_sol /= 20;
	for(double sol: all_sol)
	{
		std_sol += (sol - mean_sol) * (sol - mean_sol);
	}
	std_sol = sqrt(std_sol);
	cout << "Min: " << min_sol << '\n';
	cout << "Mean: " << mean_sol << '\n';
	cout << "Std: " << std_sol << '\n';
	
}
