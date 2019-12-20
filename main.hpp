#ifndef PSO_OPENMP
#define PSO_OPENMP 1

#define THREADS (8)
#define DIMEN (10)
#define PSO_C1 (2)
#define PSO_C2 (2)
#define PSO_MAX_W (0.9)
#define PSO_MIN_W (0.4)

/**
 * Type definitions
 */
typedef double Valtype;
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
bool compare(Valtype a, Valtype b);
Valtype coord_in_range(Valtype coord, unsigned int dimen);
Valtype f(Coord coord);
void init_nodes(Node *nodes);
void next_gen(Node *nodes, double w);
Valtype pso_next_coord(Node node, unsigned int dimen, double w);
double pso_next_veloc(Node node, unsigned int dimen, double w);
Valtype random_coord(Node node, unsigned int dimen);
double random_veloc(Node node, unsigned int dimen);

/**
 * Global variables
 */
Node global_best;
unsigned int times, nodes_amount;
const Coord COORD_RANGE[2] = {{-5.12, -5.12, -5.12, -5.12, -5.12, -5.12, -5.12, -5.12, -5.12, -5.12},
                              {5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12}};

#endif /* PSO_OPENMP */
