// Dijkstra's Algorthim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <queue> 
#include <list>
#include <cstdlib>
#include <ctime>
#include <fstream>
using namespace std;


# define INF 0x3f3f3f3f 


typedef pair<int, int> iPair;

// Source of GeeksForGeeks
//https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-using-priority_queue-stl/


class Graph
{
    int V;    // No. of vertices 

    // In a weighted graph, we need to store vertex 
    // and weight pair for every edge 
    list< pair<int, int> >* adj;

public:
    Graph(int V);  // Constructor 

    // function to add an edge to graph 
    void addEdge(int u, int v, int w);

    // prints shortest path from s 
    int shortestPath(int s, int d);
};

// Allocates memory for adjacency list 
Graph::Graph(int V)
{
    this->V = V;
    adj = new list<iPair>[V];
}

void Graph::addEdge(int u, int v, int w)
{
    adj[u].push_back(make_pair(v, w));
    adj[v].push_back(make_pair(u, w));
}

// Prints shortest paths from src to all other vertices 
int Graph::shortestPath(int src, int d)
{
    // Create a priority queue to store vertices that 
    // are being preprocessed. This is weird syntax in C++. 
    // Refer below link for details of this syntax 
    // https://www.geeksforgeeks.org/implement-min-heap-using-stl/ 
    priority_queue< iPair, vector <iPair>, greater<iPair> > pq;

    // Create a vector for distances and initialize all 
    // distances as infinite (INF) 
    vector<int> dist(V, INF);

    // Insert source itself in priority queue and initialize 
    // its distance as 0. 
    pq.push(make_pair(0, src));
    dist[src] = 0;

    /* Looping till priority queue becomes empty (or all
      distances are not finalized) */
    while (!pq.empty())
    {
        // The first vertex in pair is the minimum distance 
        // vertex, extract it from priority queue. 
        // vertex label is stored in second of pair (it 
        // has to be done this way to keep the vertices 
        // sorted distance (distance must be first item 
        // in pair) 
        int u = pq.top().second;
        pq.pop();

        // 'i' is used to get all adjacent vertices of a vertex 
        list< pair<int, int> >::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
        {
            // Get vertex label and weight of current adjacent 
            // of u. 
            int v = (*i).first;
            int weight = (*i).second;

            //  If there is shorted path to v through u. 
            if (dist[v] > dist[u] + weight)
            {
                // Updating distance of v 
                dist[v] = dist[u] + weight;
                pq.push(make_pair(dist[v], v));
            }
        }
    }

    for (int i = 0; i < V; ++i) {
        if (i == d)
            return dist[i];
    }
}



class Point
{
    int x;
    int y;
public:
    Point(int _x, int _y)
    {
        x = _x;
        y = _y;
    }
    int getX() const { return x; }
    int getY() const { return y; }
};

class myComparator
{
public:
    int operator() (const Point& p1, const Point& p2)
    {
        return p1.getX() > p2.getX();
    }
};


// My implementation of the Dijkstra Algorthim
// It calaculates the shortest distance between two given nodes in a graph
// C : Graph as vector based Source - destination - weight implementation
// S : Source node
// d : Destination node

int find_shortest_path(vector <tuple<int, int, int>> C, int s, int d) {
    priority_queue <Point, vector<Point>, myComparator > pq;
    vector <int> k;
    pq.push(Point(0, s));
    int curr_path_weight = 0;
    int curr_node = s;
    while (s != d) { // while destination not reached
        vector <tuple<int, int, int>> curr;

        // The for loop pushes the current neighbours 
        for (int i = 0; i < C.size(); i++) {
            if (get<0>(C[i]) == s)
                curr.push_back(make_tuple(get<0>(C[i]), get<1>(C[i]), get<2>(C[i])));
        }
        // Pushes the neighbours as undirected graph
        // The od one was a directed graph
        for (int i = 0; i < C.size(); i++) {
            if (get<1>(C[i]) == s)
                curr.push_back(make_tuple(get<1>(C[i]), get<0>(C[i]), get<2>(C[i])));
        }
        int least = 0;

        // it checks for least distance   
        for (int i = 0; i < curr.size(); i++) {
            if (get<2>(curr[i]) < get<2>(curr[least])) {
                least = i;
            }
        }

        // Push the current destinations and distances 
        for (int i = 0; i < curr.size(); i++) {
            pq.push(Point(get<2>(curr[i]) + curr_path_weight, get<1>(curr[i])));
        }
        // pop the  new expected source distance
        s = pq.top().getY();
        // pop the new expected source least distance
        int temp_path_weight = pq.top().getX();
        // remove the current node
        pq.pop();
        //if there is no neghibours
        if (curr.size() == 0) {
            curr_path_weight = temp_path_weight;
            continue;
        }
        //compare between two distances "Current and observed"
        if (get<2>(curr[least]) + curr_path_weight < temp_path_weight || temp_path_weight == 0) {
            curr_path_weight = get<2>(curr[least]);
            s = get<1>(curr[least]);

        }
        else {
            curr_path_weight = temp_path_weight;
        }

    }
    //If reached return the current weight which is the smallest
    if (s == d)
        return curr_path_weight;
    else
        return 0;
}
int main()
{

    ofstream R("Output.txt");
    for (int k = 0; k < 50; k++) {
        srand(time(0));
        int n = 20;
        Graph g(n);
        vector <tuple<int, int, int>> v;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                int w = (rand() % 10 + 1 + (k * 5) + 10 + k + i + j) % 7 + 1;
                int prob = rand() % 10;
                //cout << "i = " << i <<", j = " << j << ", w = "<< w << endl;
                if (w > 0) {
                    v.push_back(make_tuple(i, j, w));
                    g.addEdge(i, j, w);
                }
            }
        }
        int s = 0;
        int d = 10;


        // Stress Testing 
        if (g.shortestPath(s, d) != find_shortest_path(v, s, d)) {
            R << "solution differ\n";
            R << "internet solution " << g.shortestPath(s, d) << ", my solution" << find_shortest_path(v, s, d) << endl;
            int l = 0;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    R << "i = " << i << " , j = " << j << " , w = " << get<2>(v[l]) << endl;
                    l++;
                }
            }
        }
        else {
            R << "solution ok" << endl;
        }

    }



}







