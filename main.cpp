#include <assert.h>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "main.hpp"
#define rand_double() ((double)rand() / RAND_MAX)

int main()
{
  std::cout << "Execution times: ";
  std::cin >> times;
  std::cout << "Number of nodes: ";
  std::cin >> nodes_amount;
  Node nodes[nodes_amount];
  time_t start_ts = time(0);

  init_nodes(nodes);
  for (unsigned int gen = 0; gen < times; gen++)
    next_gen(nodes, (PSO_MAX_W - PSO_MIN_W) * (times - gen) / times + PSO_MIN_W);

  // print
  std::cout << "Cost: " << time(0) - start_ts << std::endl
            << "Value: " << global_best.best_value << std::endl
            << "Coordinate: (" << global_best.best_coord.ref[0];
  for (unsigned int i = 1; i < DIMEN; i++)
    std::cout << ", " << global_best.best_coord.ref[i];
  std::cout << ")" << std::endl;

  return 0;
}

bool compare(Valtype a, Valtype b)
{
  return a < b;
}

Valtype coord_in_range(Valtype coord, unsigned int dimen)
{
  return std::max(std::min(coord, COORD_RANGE[1].ref[dimen]), COORD_RANGE[0].ref[dimen]);
}

Valtype f(Coord coord)
{
  Valtype ret = 0.0;
  for (unsigned int i = 0; i < DIMEN; i++)
    ret += coord.ref[i] * coord.ref[i] - cos(M_PI * 2.0 * coord.ref[i]) + 10.0;
  return ret;
}

void init_nodes(Node *nodes)
{
  for (unsigned int i = 0; i < nodes_amount; i++)
  {
    for (unsigned int j = 0; j < DIMEN; j++)
      nodes[i].veloc[j] = random_veloc(nodes[i], j),
      nodes[i].coord.ref[j] = random_coord(nodes[i], j);
    nodes[i].best_value = f(nodes[i].best_coord = nodes[i].coord);
    if (!i || compare(nodes[i].best_value, global_best.best_value))
      global_best = nodes[i];
  }
}

void next_gen(Node *nodes, double w)
{
  Valtype value;
#pragma omp parallel for num_threads(THREADS) private(value)
  for (unsigned int i = 0; i < nodes_amount; i++)
  {
    for (unsigned int j = 0; j < DIMEN; j++)
      nodes[i].veloc[j] = pso_next_veloc(nodes[i], j, w),
      nodes[i].coord.ref[j] = pso_next_coord(nodes[i], j, w);
    if (compare(value = f(nodes[i].coord), nodes[i].best_value))
    {
      nodes[i].best_value = value, nodes[i].best_coord = nodes[i].coord;
      if (compare(value, global_best.best_value))
#pragma omp critical
        global_best = nodes[i];
    }
  }
}

Valtype pso_next_coord(Node node, unsigned int dimen, double w)
{
  return node.coord.ref[dimen] + node.veloc[dimen];
};

double pso_next_veloc(Node node, unsigned int dimen, double w)
{
  return coord_in_range(
      w * node.veloc[dimen] +
          rand_double() * PSO_C1 * (node.best_coord.ref[dimen] - node.coord.ref[dimen]) +
          rand_double() * PSO_C2 * (global_best.best_coord.ref[dimen] - node.coord.ref[dimen]),
      dimen);
};

Valtype random_coord(Node _node, unsigned int dimen)
{
  Valtype range = COORD_RANGE[1].ref[dimen] - COORD_RANGE[0].ref[dimen];
  return rand_double() * range + COORD_RANGE[0].ref[dimen];
}

double random_veloc(Node _node, unsigned int dimen)
{
  return rand_double();
};
