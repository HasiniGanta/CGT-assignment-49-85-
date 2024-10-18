#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>
#include <time.h>

void findFundamentalCutsetsAndCircuits(int adjMatrix[][100], int weights[][100], int parent[], int n);
int minKey(int key[], bool mstSet[], int n);
void DFS(int adjMatrix[][100], int u, bool visited[], int n);

typedef struct {
    int degree;
    int vertex;
} HeapNode;

void swap(HeapNode* a, HeapNode* b) {
    HeapNode temp = *a;
    *a = *b;
    *b = temp;
}

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

// to build max heap
void buildHeap(HeapNode heap[], int n) {
    for (int i = n / 2 - 1; i >= 0; i--) {
        heapify(heap, n, i);
    }
}

//delete the max element from the heap
HeapNode delMax(HeapNode heap[], int* n) {
    HeapNode maxNode = heap[0];
    heap[0] = heap[*n - 1];
    (*n)--;
    heapify(heap, *n, 0);
    return maxNode;
}

void insertHeap(HeapNode heap[], int* n, int degree, int vertex) {
    heap[*n].degree = degree;
    heap[*n].vertex = vertex;
    (*n)++;
    for (int i = *n / 2 - 1; i >= 0; i--) {
        heapify(heap, *n, i);
    }
}

void mergeHeaps(HeapNode heap[], int* n, HeapNode h1[], int* n1) {
    for (int i = 0; i < *n1; i++) {
        insertHeap(heap, n, h1[i].degree, h1[i].vertex);
    }
}

int generateRandomWeight() {
    return rand() % 10 + 1;  // Random weight between 1 and 10
}

//Bellman-Ford algorithm
void bellmanFord(int adjMatrix[][100], int weights[][100], int n, int src) {
    int dist[n];

    // Step 1: initialize distances from src to all other vertices as infinite
    for (int i = 0; i < n; i++) {
        dist[i] = INT_MAX;
    }
    dist[src] = 0;

    // Step 2: iterate all edges |V| - 1 times
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
    printf("Vertex Distance from Source %d:\n", src);
    for (int i = 0; i < n; i++) {
        if (dist[i] == INT_MAX) {
            printf("Vertex %d: INF\n", i);
        } else {
            printf("Vertex %d: %d\n", i, dist[i]);
        }
    }
}

//finding the vertex with the minimum key value
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
    int parent[n]; //array to store constructed MST
    int key[n];    //to pick minimum weight edges
    bool mstSet[n]; //track of vertices included in the MST

    for (int i = 0; i < n; i++) {
        key[i] = INT_MAX;
        mstSet[i] = false;
    }

    key[0] = 0;   
    parent[0] = -1; //root
    for (int count = 0; count < n - 1; count++) {
        // Pick the minimum key vertex not included in the MST
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

    // fundamental cutsets and circuits
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
                printf("Fundamental Cutset for edge (%d, %d): { %d, %d }\n", u, v, u, v);
            }
        }
    }
}
void sortDesc(int* degrees, int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - 1 - i; j++) {
            if (degrees[j] < degrees[j + 1]) {
                int temp = degrees[j];
                degrees[j] = degrees[j + 1];
                degrees[j + 1] = temp;
            }
        }
    }
}

//minimum degree in the array
int findMinDegree(int* degrees, int n) {
    int minDegree = degrees[0];
    for (int i = 1; i < n; i++) {
        if (degrees[i] < minDegree) {
            minDegree = degrees[i];
        }
    }
    return minDegree;
}
// Havel-Hakimi algorithm
int havelHakimi(int* degrees, int n) {
    int* temp = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        temp[i] = degrees[i];
    }
    while (1) {
        sortDesc(temp, n);
        if (temp[0] == 0) {
            free(temp);
            return 1; // Sequence is graphical
        }
        int d = temp[0];
        temp[0] = 0;
        if (d >= n || d < 0) {
            free(temp);
            return 0; // Sequence is not graphical
        }
        for (int i = 1; i <= d; i++) {
            temp[i]--;
            if (temp[i] < 0) {
                free(temp);
                return 0; // Sequence is not graphical
            }
        }
    }
}
// K-connectivity
int kConnectivity(int* degrees, int n) {
    int minDegree = findMinDegree(degrees, n);
    return minDegree - 1; // k-connectivity is minimum degree - 1
}


int main() {
    srand(time(0));  // Seed for random number generation

    int n;
    printf("Enter the number of vertices: ");
    scanf("%d", &n);

    int degrees[n];
    HeapNode heap[n];
    HeapNode h1[n];
    int adjMatrix[100][100]; 
    int weights[100][100];   
    int heapSize = 0;
    int h1Size = 0;

    printf("Enter the degree sequence:\n");
    for (int i = 0; i < n; i++) {
        scanf("%d", &degrees[i]);
        heap[i].degree = degrees[i];
        heap[i].vertex = i;
        heapSize++;
    }

    // Initialize the adjacency matrix
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
    // Perform Bellman-Ford algorithm 
    int src;
    printf("\nEnter the source vertex for Bellman-Ford: ");
    scanf("%d", &src);
    bellmanFord(adjMatrix, weights, n, src);

    // Perform Prim's algorithm 
    primMST(adjMatrix, weights, n);
    if (havelHakimi(degrees, n)) {
        printf("The degree sequence is graphical.\n");
        int edgeConnectivity = findMinDegree(degrees, n);
        int vertexConnectivity = findMinDegree(degrees, n);
        int kConnectivityValue = kConnectivity(degrees, n);

        printf("Edge Connectivity: %d\n", edgeConnectivity);
        printf("Vertex Connectivity: %d\n", vertexConnectivity);
        printf("k-Connectivity: %d\n", kConnectivityValue);
    } else {
        printf("The degree sequence is not graphical.\n");
    }


    return 0;
}