#include <assert.h>
#include <iostream>
#include <omp.h>
#define DIMEN (2)
#define THREADS (8)
#define PSO_C1 (2)
#define PSO_C2 (2)
#define PSO_MAX_W (0.9)
#define PSO_MIN_W (0.4)
#define random(from, to) (rand() % ((to) - (from)) + (from))
using namespace std;

/**
 * Type definitions
 */
typedef int Valtype;
typedef struct Coord
{
  Valtype ref[DIMEN];
} Coord;
typedef struct Node
{
  double veloc[DIMEN];
  Coord coord;
  Coord best_coord;
  Valtype best_value;
} Node;

/**
 * Funtion definitions
 */
Valtype f(Coord coord);
double rand_double();
bool compare(Valtype a, Valtype b);
Node find_global_best(Node *nodes, unsigned int nodes_amount);
Valtype pso_next_coord(Node node, unsigned int dimen);
double pso_next_veloc(Node node, unsigned int dimen);
Valtype random_coord(Node node, unsigned int dimen);
double random_veloc(Node node, unsigned int dimen);
void assign(Node *nodes,
            unsigned int nodes_amount,
            bool calc_best,
            double (*next_veloc)(Node node, unsigned int dimen),
            Valtype (*next_coord)(Node node, unsigned int dimen));

/**
 * Global variables
 */
double w;
Node global_best;
const Coord INIT_COORD_RANGE[2] = {{-10, -10}, {10, 10}};

int main()
{
  unsigned int times, nodes_amount;
  cout << "Execution times: ";
  cin >> times;
  cout << "Number of nodes: ";
  cin >> nodes_amount;
  Node nodes[nodes_amount];

  // initialization
  assign(nodes, nodes_amount, false, random_veloc, random_coord);
  for (unsigned int i = 0; i < nodes_amount; i++)
    nodes[i].best_value = f(nodes[i].coord), nodes[i].best_coord = nodes[i].coord;
  global_best = find_global_best(nodes, nodes_amount);

  // pso
  for (unsigned int gen = 0; gen < times; gen++)
  {
    w = ((PSO_MAX_W - PSO_MIN_W) * (times - gen) / times + PSO_MIN_W);
    assign(nodes, nodes_amount, true, pso_next_veloc, pso_next_coord);
    global_best = find_global_best(nodes, nodes_amount);
  }

  // print
  cout << "Result: (" << global_best.best_coord.ref[0];
  for (unsigned int i = 1; i < DIMEN; i++)
    cout << ", " << global_best.best_coord.ref[i];
  cout << ").";

  return 0;
}

Valtype f(Coord coord)
{
  Valtype ret = 0;
  for (unsigned int i = 0; i < DIMEN; i++)
    ret += coord.ref[i] * coord.ref[i];
  return ret;
}

double rand_double()
{
  return (double)rand() / RAND_MAX;
};

bool compare(Valtype a, Valtype b)
{
  return a < b;
}

Node find_global_best(Node *nodes, unsigned int nodes_amount)
{
  assert(nodes_amount > 0);
  Node best = nodes[0];
  for (unsigned int i = 1; i < nodes_amount; i++)
    best = compare(best.best_value, nodes[i].best_value) ? best : nodes[i];
  return best;
}

void assign(Node *nodes,
            unsigned int nodes_amount,
            bool calc_best,
            double (*next_veloc)(Node node, unsigned int dimen),
            Valtype (*next_coord)(Node node, unsigned int dimen))
{
#pragma omp parallel for num_threads(THREADS)
  for (unsigned int i = 0; i < nodes_amount; i++)
  {
    for (unsigned int j = 0; j < DIMEN; j++)
    {
      nodes[i].veloc[j] = next_veloc(nodes[i], j);
      nodes[i].coord.ref[j] = next_coord(nodes[i], j);
    }
    if (calc_best)
    {
      Valtype value = f(nodes[i].coord);
      if (compare(value, nodes[i].best_value))
        nodes[i].best_value = value, nodes[i].best_coord = nodes[i].coord;
    }
  }
}

Valtype random_coord(Node _node, unsigned int dimen)
{
  return random(INIT_COORD_RANGE[0].ref[dimen], INIT_COORD_RANGE[1].ref[dimen]);
}

double random_veloc(Node _node, unsigned int dimen)
{
  return rand_double();
};

Valtype pso_next_coord(Node node, unsigned int dimen)
{
  return node.coord.ref[dimen] + node.veloc[dimen];
};

double pso_next_veloc(Node node, unsigned int dimen)
{
  return w * node.veloc[dimen] +
         rand_double() * PSO_C1 * (node.best_coord.ref[dimen] - node.coord.ref[dimen]) +
         rand_double() * PSO_C2 * (global_best.best_coord.ref[dimen] - node.coord.ref[dimen]);
};
