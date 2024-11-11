import random
import math
import time

class Delta:
    def __init__(self, x: int = 0, y: int = 0):
        self.x = x
        self.y = y

    def __str__(self):
        return '(%d, %d)' % (self.x, self.y)

def compute(sequence, selected) -> float:
    n = len(sequence)
    total_distance = 0.0
    for idx in range(n):
        start_pt = selected[sequence[idx]]
        end_pt = selected[sequence[(idx + 1) % n]]
        total_distance += math.hypot(start_pt.x - end_pt.x, start_pt.y - end_pt.y)
    return total_distance

def genetic_algorithm(collection, population_size=100, generations=500):
    total_clusters = len(collection)

    # 初始化种群
    population = []
    for _ in range(population_size):
        individual = []
        for cluster in collection:
            idx = random.randint(0, len(cluster) - 1)
            individual.append(idx)
        population.append(individual)

    best_distance = float('inf')
    best_individual = None
    best_sequence = None

    for generation in range(generations):
        # 计算适应度
        fitness = []
        for individual in population:
            selected_points = [cluster[individual[i]] for i, cluster in enumerate(collection)]
            # 生成访问顺序
            sequence = list(range(total_clusters))
            random.shuffle(sequence)
            total_distance = compute(sequence, selected_points)
            fitness.append((total_distance, individual, sequence))

        # 更新最佳解
        fitness.sort(key=lambda x: x[0])
        if fitness[0][0] < best_distance:
            best_distance = fitness[0][0]
            best_individual = fitness[0][1][:]
            best_sequence = fitness[0][2][:]

        # 选择操作（锦标赛选择）
        new_population = []
        for _ in range(population_size):
            i1, i2 = random.sample(range(population_size), 2)
            winner = fitness[i1] if fitness[i1][0] < fitness[i2][0] else fitness[i2]
            new_population.append(winner[1])

        # 交叉操作（部分匹配交叉 PMX）
        for i in range(0, population_size - 1, 2):
            parent1 = new_population[i]
            parent2 = new_population[i + 1]
            if random.random() < 0.8:  # 80%的概率进行交叉
                point = random.randint(1, total_clusters - 2)
                child1 = parent1[:point] + parent2[point:]
                child2 = parent2[:point] + parent1[point:]
                new_population[i] = child1
                new_population[i + 1] = child2

        # 变异操作
        for individual in new_population:
            if random.random() < 0.1:  # 10%的概率进行变异
                idx = random.randint(0, total_clusters - 1)
                individual[idx] = random.randint(0, len(collection[idx]) - 1)

        population = new_population

    # 最后使用2-opt优化路径
    selected_points = [collection[i][best_individual[i]] for i in range(total_clusters)]
    best_sequence, best_distance = two_opt(best_sequence, selected_points)

    output = '[' + ', '.join([str(selected_points[v]) for v in best_sequence]) + ']'
    print(output)
    print(f'Best total distance: {best_distance}')

def two_opt(sequence, selected):
    n = len(sequence)
    best_sequence = sequence[:]
    best_distance = compute(sequence, selected)
    improved = True
    while improved:
        improved = False
        for i in range(n - 1):
            for j in range(i+2, n):
                if j - i == 1: continue
                new_sequence = best_sequence[:]
                new_sequence[i+1:j+1] = reversed(best_sequence[i+1:j+1])
                new_distance = compute(new_sequence, selected)
                if new_distance < best_distance - 1e-6:
                    best_sequence = new_sequence
                    best_distance = new_distance
                    improved = True
        if not improved:
            break
    return best_sequence, best_distance

if __name__ == "__main__":
    random.seed(time.time())
    collection = []
    for row in [eval(w) for w in input().split("@")]:
        tmp = []
        for (x, y) in row:
            tmp.append(Delta(x, y))
        collection.append(tmp)
    genetic_algorithm(collection, population_size=200, generations=1000)
