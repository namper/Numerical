#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

/**
 * Define the constants here
 */
double initial_x, final_x, initial_t, final_t;
int partition_x, partition_t;

double h = (final_x - initial_x)/partition_x;
double tau = (final_t - initial_t)/partition_t;

double gamma = tau/(h*h);

double boundary_1(double t);
double boundary_2(double t);
double initial(double x);

double source(double x, double t);

double exact_function(double x, double t);

vector<vector<double>> createGrid(int x_partitions, int t_partitions);
void fillTheGridWithInitialConditions(vector<vector<double>> & grid);

double F_source(vector<vector<double>> & grid , int index_x, int index_t);
double alpha(int index);
double beta(int index, int index_t, vector<vector<double>> & grid);

void approximateTheNextLine(vector<vector<double>> & grid, int index_t);

int main(int argc, char const *argv[])
{
    vector<vector<double>> y = createGrid(partition_x, partition_t);
    fillTheGridWithInitialConditions(y);
    for(int k = 1; k <= partition_t; k++)
    {
        approximateTheNextLine(y);
    }
    return 0;
}

/**
 * Codes of helper functions
 */

vector<vector<double>> createGrid(int x_partitions, int t_partitions)
{
    /**
     * Creates (x_partitions+1)X(t_partitions+1) grid full of zeros
     */

    vector<double> zeros;
    for(int i = 0; i <= t_partitions; i++)
    {
        zeros.push_back(0);
    }
    vector<vector<double>> matrix;
    for(int i = 0; i <= x_partitions; i++)
    {
        matrix.push_back(zeros);
    }
    return matrix;
}

void fillTheGridWithInitialConditions(vector<vector<double>> & grid)
{
    /**
     * Fills the grid with initial conditions defined above
     */

    for(int k = 0; k <= partition_x; k++)
    {
        grid[k][0] = initial(initial_x + h*k); 
    }
    for(int k = 0; k < partition_t; k++)
    {
        grid[0][k] = boundary_1(initial_t + tau*k);
        grid[partition_x][k] = boundary_2(initial_t + tau*k);
    }
}

double F_source(vector<vector<double>> & grid , int index_x, int index_t)
{
    return -grid[index_x][index_t - 1] - tau*source(initial_x + index_x*h, initial_t + index_t*tau);
}

double alpha(int index)
{
    if (index == 1)
        return 0;
    else
        return gamma / (1 + 2 * gamma - gamma * alpha(index - 1)); 
}

double beta(int index, int index_t, vector<vector<double>> & grid)
{
    if (index == 1)
        return boundary_1(initial_t + index_t*tau);
    else
        return (gamma* beta(index - 1, index_t) - F_source(grid, index - 1, index_t))/(1 + 2*gamma - gamma*alpha(index - 1));
}

void approximateTheNextLine(vector<vector<double>> & grid, int index_t)
{
    for (int i = partition_x - 1; i >= 1; i--)
        grid[i][index_t] = alpha(i + 1)*grid[i+1][index_t] + beta(i+1, index_t, grid);
}
