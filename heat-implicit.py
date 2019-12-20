import numpy as np 

initial_x = 0
final_x = 1

initial_t = 0
final_t = 1

x_partition = 5 
t_partition = 15

h = (final_x - initial_x)/x_partition
tau = (final_t - initial_t)/t_partition

gamma = tau/(h**2)
def boundary_1(t):
    return 2 + t

def boundary_2(t):
    return 2 + t

def initial_cond(x):
    return 2

def source (x, t):
    return x**2 + t

def create_grid():
    y = []
    for i in range(x_partition + 1):
        y.append(np.zeros(shape=t_partition + 1))
    return y

def approximate_solution():
    solutions = create_grid()
    for i in range( t_partition + 1 ):
        solutions[0][i] = boundary_1(initial_t + tau*i)
        solutions[x_partition][i] = boundary_2(initial_t + tau*i)

    for i in range( x_partition + 1):
        solutions[i][0] = initial_cond(initial_x + h*i)

    def F(index_i, index_j):
        return -solutions[index_i][index_j - 1] - tau*source(initial_x + h*index_i, initial_t + tau*index_j)

    def alpha(index):
        if index == 1:
            return 0
        return gamma/(1 + 2*gamma - gamma*alpha(index - 1)) 
    
    def beta(index, index_j):
        if index == 1:
            return boundary_1(index_j - 1)
        return (gamma*beta(index - 1, index_j) - F(index - 1, index_j))/(1 + 2*gamma - gamma*alpha(index - 1))

    for index_j in range(1, t_partition + 1):
        for index_i in range(x_partition - 1, 0, -1):
            solutions[index_i][index_j] = alpha(index_i + 1)*solutions[index_i + 1][index_j] + beta(index_i+1,index_j)

    return solutions
if __name__ == "__main__":
    matrix = approximate_solution()
    for i in range(x_partition+1):
        print(matrix[i])