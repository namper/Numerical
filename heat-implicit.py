import numpy as np


class ApproximatedSolution:
    initial_x = 0
    final_x = 1

    initial_t = 0
    final_t = 1

    x_partition = 5
    t_partition = 15

    h = (final_x - initial_x) / x_partition
    tau = (final_t - initial_t) / t_partition

    gamma = tau / (h ** 2)

    def __init__(self):
        self.grid = np.zeros(shape=(self.x_partition + 1, self.t_partition + 1))
        self.approximate_solution()

    @staticmethod
    def boundary_1(t):
        return 2 + t

    @staticmethod
    def boundary_2(t):
        return 2 + t

    @staticmethod
    def initial_cond(x):
        return 2

    @staticmethod
    def source(x, t):
        return x ** 2 + t

    def f(self, index_i, index_j):
        return -self.grid[index_i][index_j - 1] - self.tau * \
               self.source(self.initial_x + self.h * index_i,
                           self.initial_t + self.tau * index_j)

    def alpha(self, index):
        if index == 1:
            return 0
        return self.gamma / (1 + 2 * self.gamma - self.gamma * self.alpha(index - 1))

    def beta(self, index, index_j):
        if index == 1:
            return self.boundary_1(index_j - 1)
        return (self.gamma * self.beta(index - 1, index_j) - self.f(index - 1, index_j)) / (
                1 + 2 * self.gamma - self.gamma * self.alpha(index - 1))

    def initialize_boundary_conditions(self):
        self.grid.T[0] = np.array([self.initial_cond(self.initial_x + self.h * i) for i in range(self.x_partition + 1)])
        print(self.x_partition)
        self.grid[self.x_partition] = np.array(
            [self.boundary_2(self.initial_t + self.tau * i) for i in range(self.t_partition + 1)]
        )
        self.grid[0] = np.array([self.boundary_1(self.initial_t + self.tau * i) for i in range(self.t_partition + 1)])

    def approximate_solution(self):
        self.initialize_boundary_conditions()
        for j in range(1, self.t_partition + 1):
            for i in range(self.x_partition - 1, 0, -1):
                self.grid[i][j] = self.alpha(i + 1) * self.grid[i + 1][j] + self.beta(i + 1, j)


if __name__ == "__main__":
    approximated_matrix = ApproximatedSolution().grid
    print(approximated_matrix)
