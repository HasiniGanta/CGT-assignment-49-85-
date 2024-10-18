#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>
#include <time.h>

// Function prototypes
void findFundamentalCutsetsAndCircuits(int adjMatrix[][100], int weights[][100], int parent[], int n);
int minKey(int key[], bool mstSet[], int n);
void DFS(int adjMatrix[][100], int u, bool visited[], int n);
bool isConnected(int adjMatrix[][100], int n);
void articulationPointsUtil(int adjMatrix[][100], int u, bool visited[], int disc[], int low[], int parent[], bool ap[], int n);
void findArticulationPoints(int adjMatrix[][100], int n);
int findEdgeConnectivity(int adjMatrix[][100], int n);
void checkKConnectivity(int adjMatrix[][100], int n);
int edgeConnectivityUtil(int adjMatrix[][100], int u, int v, int n);

// Structure for a heap node
typedef struct {
    int degree;
    int vertex;
} HeapNode;

// Function to swap two heap nodes
void swap(HeapNode* a, HeapNode* b) {
    HeapNode temp = *a;
    *a = *b;
    *b = temp;
}

// Function to heapify a subtree rooted with node i
void heapify(HeapNode heap[], int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < n && heap[left].degree > heap[largest].degree)
        largest = left;

    if (right < n && heap[right].degree > heap[largest].degree)
        largest = right;

    if (largest != i) {
        swap(&heap[i], &heap[largest]);
        heapify(heap, n, largest);
    }
}

// Function to build a max heap
void buildHeap(HeapNode heap[], int n) {
    for (int i = n / 2 - 1; i >= 0; i--) {
        heapify(heap, n, i);
    }
}

// Function to delete the max element from the heap
HeapNode delMax(HeapNode heap[], int* n) {
    HeapNode maxNode = heap[0];
    heap[0] = heap[*n - 1];
    (*n)--;
    heapify(heap, *n, 0);
    return maxNode;
}

// Function to insert a new element into the heap
void insertHeap(HeapNode heap[], int* n, int degree, int vertex) {
    heap[*n].degree = degree;
    heap[*n].vertex = vertex;
    (*n)++;
    for (int i = *n / 2 - 1; i >= 0; i--) {
        heapify(heap, *n, i);
    }
}

// Function to merge two heaps
void mergeHeaps(HeapNode heap[], int* n, HeapNode h1[], int* n1) {
    for (int i = 0; i < *n1; i++) {
        insertHeap(heap, n, h1[i].degree, h1[i].vertex);
    }
}

// Function to generate a random weight
int generateRandomWeight() {
    return rand() % 10 + 1;  // Random weight between 1 and 10
}

// Function to implement Bellman-Ford algorithm
void bellmanFord(int adjMatrix[][100], int weights[][100], int n, int src) {
    int dist[n];

    // Step 1: Initialize distances from src to all other vertices as INFINITE
    for (int i = 0; i < n; i++) {
        dist[i] = INT_MAX;
    }
    dist[src] = 0;

    // Step 2: Relax all edges |V| - 1 times
    for (int i = 0; i < n - 1; i++) {
        for (int u = 0; u < n; u++) {
            for (int v = 0; v < n; v++) {
                if (adjMatrix[u][v] == 1 && dist[u] != INT_MAX && dist[u] + weights[u][v] < dist[v]) {
                    dist[v] = dist[u] + weights[u][v];
                }
            }
        }
    }

    // Step 3: Check for negative-weight cycles
    for (int u = 0; u < n; u++) {
        for (int v = 0; v < n; v++) {
            if (adjMatrix[u][v] == 1 && dist[u] != INT_MAX && dist[u] + weights[u][v] < dist[v]) {
                printf("Graph contains negative weight cycle\n");
                return;
            }
        }
    }

    // Print the shortest distances
    printf("Vertex Distance from Source %d:\n", src);
    for (int i = 0; i < n; i++) {
        if (dist[i] == INT_MAX) {
            printf("Vertex %d: INF\n", i);
        } else {
            printf("Vertex %d: %d\n", i, dist[i]);
        }
    }
}

// Function to find the vertex with the minimum key value
int minKey(int key[], bool mstSet[], int n) {
    int min = INT_MAX, minIndex;

    for (int v = 0; v < n; v++) {
        if (mstSet[v] == false && key[v] < min) {
            min = key[v];
            minIndex = v;
        }
    }

    return minIndex;
}

// Function to implement Prim's algorithm for MST and print the MST path
void primMST(int adjMatrix[][100], int weights[][100], int n) {
    int parent[n]; // Array to store constructed MST
    int key[n];    // Key values to pick minimum weight edges
    bool mstSet[n]; // To keep track of vertices included in the MST

    // Initialize all keys as INFINITE and mstSet[] as false
    for (int i = 0; i < n; i++) {
        key[i] = INT_MAX;
        mstSet[i] = false;
    }

    key[0] = 0;   // Start from the first vertex
    parent[0] = -1; // Root of the MST

    // The MST will have n-1 edges
    for (int count = 0; count < n - 1; count++) {
        // Pick the minimum key vertex not yet included in the MST
        int u = minKey(key, mstSet, n);

        // Add the picked vertex to the MST set
        mstSet[u] = true;

        // Update the key values and parent of the adjacent vertices
        for (int v = 0; v < n; v++) {
            if (adjMatrix[u][v] && mstSet[v] == false && weights[u][v] < key[v]) {
                parent[v] = u;
                key[v] = weights[u][v];
            }
        }
    }

    // Print the MST path (edges and their weights)
    printf("MST Path (Edges and Weights):\n");
    for (int i = 1; i < n; i++) {
        printf("%d - %d: %d\n", parent[i], i, weights[i][parent[i]]);
    }

    // Find fundamental cutsets and circuits
    findFundamentalCutsetsAndCircuits(adjMatrix, weights, parent, n);
}

// Function to find fundamental cutsets and circuits
void findFundamentalCutsetsAndCircuits(int adjMatrix[][100], int weights[][100], int parent[], int n) {
    printf("\nFundamental Cutsets and Circuits:\n");

    // Iterate through each non-MST edge
    for (int u = 0; u < n; u++) {
        for (int v = 0; v < n; v++) {
            if (adjMatrix[u][v] && u != parent[v] && v != parent[u]) {
                // This is a non-MST edge, it forms a fundamental circuit when added to the MST

                printf("Fundamental Circuit for edge (%d, %d):\n", u, v);

                // Ensure the parent is valid (i.e., not -1) before printing the cycle
                if (parent[u] != -1 && parent[v] != -1) {
                    printf("Cycle: %d -> %d -> %d\n", u, v, parent[v]);
                } else {
                    printf("Cycle: %d -> %d (No further parent)\n", u, v);
                }

                // Print fundamental cutset
                printf("Fundamental Cutset for edge (%d, %d): { %d, %d }\n", u, v, u, v);
            }
        }
    }
}

// DFS to check if the graph is connected
void DFS(int adjMatrix[][100], int u, bool visited[], int n) {
    visited[u] = true;
    for (int v = 0; v < n; v++) {
        if (adjMatrix[u][v] && !visited[v]) {
            DFS(adjMatrix, v, visited, n);
        }
    }
}

// Check if the graph is connected
bool isConnected(int adjMatrix[][100], int n) {
    bool visited[n];
    for (int i = 0; i < n; i++) {
        visited[i] = false;
    }

    DFS(adjMatrix, 0, visited, n);

    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            return false;
        }
    }
    return true;
}

// Utility function to find edge connectivity using max-flow like technique
int edgeConnectivityUtil(int adjMatrix[][100], int u, int v, int n) {
    int cutset = 0;
    bool visited[n];
    for (int i = 0; i < n; i++) {
        visited[i] = false;
    }

    DFS(adjMatrix, u, visited, n);
    if (!visited[v]) {
        return 1;  // Direct cut found
    }

    return cutset;
}

// Find the edge connectivity λ of the graph
int findEdgeConnectivity(int adjMatrix[][100], int n) {
    int minCut = INT_MAX;

    for (int u = 0; u < n; u++) {
        for (int v = u + 1; v < n; v++) {
            if (adjMatrix[u][v] == 1) {
                int cut = edgeConnectivityUtil(adjMatrix, u, v, n);
                if (cut < minCut) {
                    minCut = cut;
                }
            }
        }
    }

    return minCut;
}

// Utility function for articulation points (find vertex connectivity κ)
void articulationPointsUtil(int adjMatrix[][100], int u, bool visited[], int disc[], int low[], int parent[], bool ap[], int n) {
    static int time = 0;
    int children = 0;

    visited[u] = true;
    disc[u] = low[u] = ++time;

    for (int v = 0; v < n; v++) {
        if (adjMatrix[u][v]) {
            if (!visited[v]) {
                children++;
                parent[v] = u;
                articulationPointsUtil(adjMatrix, v, visited, disc, low, parent, ap, n);

                low[u] = (low[u] < low[v]) ? low[u] : low[v];

                if (parent[u] == -1 && children > 1) {
                    ap[u] = true;
                }

                if (parent[u] != -1 && low[v] >= disc[u]) {
                    ap[u] = true;
                }
            } else if (v != parent[u]) {
                low[u] = (low[u] < disc[v]) ? low[u] : disc[v];
            }
        }
    }
}

// Find articulation points (vertex connectivity κ)
void findArticulationPoints(int adjMatrix[][100], int n) {
    bool visited[n];
    int disc[n];
    int low[n];
    int parent[n];
    bool ap[n];

    for (int i = 0; i < n; i++) {
        visited[i] = false;
        parent[i] = -1;
        ap[i] = false;
    }

    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            articulationPointsUtil(adjMatrix, i, visited, disc, low, parent, ap, n);
        }
    }

    int vertexConnectivity = 0;
    printf("\nArticulation points (Critical Vertices):\n");
    for (int i = 0; i < n; i++) {
        if (ap[i]) {
            printf("Vertex %d\n", i);
            vertexConnectivity++;
        }
    }

    if (vertexConnectivity == 0) {
        printf("The graph is fully connected (No articulation points).\n");
    }

    printf("Vertex Connectivity (κ) = %d\n", vertexConnectivity);
}

// Check the K-connectivity of the graph
void checkKConnectivity(int adjMatrix[][100], int n) {
    int edgeConnectivity = findEdgeConnectivity(adjMatrix, n);
    findArticulationPoints(adjMatrix, n);
    printf("\nEdge Connectivity (λ) = %d\n", edgeConnectivity);

    int K = edgeConnectivity;  // Typically K-connectivity is the same as edge connectivity
    printf("The graph is %d-connected (K-connectivity).\n", K);
}

// Main function to implement the algorithm
int main() {
    srand(time(0));  // Seed for random number generation

    int n;
    printf("Enter the number of vertices: ");
    scanf("%d", &n);

    int degrees[n];
    HeapNode heap[n];
    HeapNode h1[n];
    int adjMatrix[100][100]; // Assuming maximum 100 vertices for simplicity
    int weights[100][100];   // Weights for the edges
    int heapSize = 0;
    int h1Size = 0;

    printf("Enter the degree sequence: ");
    for (int i = 0; i < n; i++) {
        scanf("%d", &degrees[i]);
        heap[i].degree = degrees[i];
        heap[i].vertex = i;
        heapSize++;
    }

    // Initialize the adjacency matrix with zeros and assign random weights
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            adjMatrix[i][j] = 0;
            weights[i][j] = 0;
        }
    }

    // Build the initial heap
    buildHeap(heap, heapSize);

    while (heapSize > 0) {
        HeapNode maxNode = delMax(heap, &heapSize);
        int d = maxNode.degree;
        int u = maxNode.vertex;

        h1Size = 0;
        for (int i = 0; i < d; i++) {
            if (heapSize == 0) {
                printf("Invalid degree sequence\n");
                return 1;
            }
            HeapNode maxNeighbor = delMax(heap, &heapSize);
            int r = maxNeighbor.degree;
            int v = maxNeighbor.vertex;

            adjMatrix[u][v] = adjMatrix[v][u] = 1; // Create an edge
            int weight = generateRandomWeight();   // Assign random weight
            weights[u][v] = weights[v][u] = weight;

            if (r > 1) {
                h1[h1Size].degree = r - 1;
                h1[h1Size].vertex = v;
                h1Size++;
            }
        }

        mergeHeaps(heap, &heapSize, h1, &h1Size);
    }

    // Print the adjacency matrix and weights
    printf("Adjacency Matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", adjMatrix[i][j]);
        }
        printf("\n");
    }

    printf("\nWeights Matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", weights[i][j]);
        }
        printf("\n");
    }

    // Perform Bellman-Ford algorithm from a given source vertex
    int src;
    printf("\nEnter the source vertex for Bellman-Ford: ");
    scanf("%d", &src);
    bellmanFord(adjMatrix, weights, n, src);

    // Perform Prim's algorithm to find the MST and fundamental cutsets/circuits
    primMST(adjMatrix, weights, n);

    // Find the edge connectivity, vertex connectivity, and K-connectivity
    checkKConnectivity(adjMatrix, n);

    return 0;
}