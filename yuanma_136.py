import random as alpha
import math as beta
import time as gamma


class SAParam(object):

    def __init__(self, start: float, end: float, reduction: float):
        self.start = start
        self.end = end
        self.reduction = reduction


class Delta(object):
    def __init__(self, x: int = 0, y: int = 0):
        self.x = x
        self.y = y

    def __str__(self):
        return '(%d, %d)' % (self.x, self.y)

    def __hash__(self):
        return hash((self.x, self.y))

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y


def compute(sequence, selected) -> float:
    n = len(sequence)
    total_distance = 0.0
    for idx in range(n):
        start_pt = selected[sequence[idx]]
        end_pt = selected[sequence[(idx + 1) % n]]
        total_distance += beta.hypot(start_pt.x - end_pt.x, start_pt.y - end_pt.y)
    return total_distance

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
                if new_distance < best_distance:
                    best_sequence = new_sequence
                    best_distance = new_distance
                    improved = True
        if not improved:
            break
    return best_sequence, best_distance

def SA(config: SAParam) -> None:
    collection = []
    for row in [eval(w) for w in input().split("@")]:
        tmp = []
        for (x, y) in row:
            tmp.append(Delta(x, y))
        collection.append(tmp)

    total = len(collection)

    alpha.seed(gamma.time())
    min_dist = beta.inf
    best_points = []
    best_order = []

    temperature = config.start
    epsilon = config.end
    reduction = config.reduction

    # 初始化选定的点
    current_selected = []
    for i in range(total):
        idx = alpha.randint(0, len(collection[i]) - 1)
        current_selected.append(collection[i][idx])

    # 初始化访问顺序
    current_sequence = list(range(total))
    alpha.shuffle(current_sequence)

    current_total_distance = compute(current_sequence, current_selected)

    min_dist = current_total_distance
    best_points = current_selected[:]
    best_order = current_sequence[:]

    while temperature > epsilon:
        # 生成邻居解
        new_selected = current_selected[:]
        new_sequence = current_sequence[:]

        if alpha.random() < 0.5:
            # 更改选定的点
            cluster_idx = alpha.randint(0, total -1)
            new_idx = alpha.randint(0, len(collection[cluster_idx]) -1)
            new_selected[cluster_idx] = collection[cluster_idx][new_idx]
        else:
            # 更改访问顺序
            idx1, idx2 = alpha.sample(range(total), 2)
            new_sequence[idx1], new_sequence[idx2] = new_sequence[idx2], new_sequence[idx1]

        # 确保选定的点是唯一的
        if len(set(new_selected)) != total:
            temperature *= reduction
            continue

        new_total_distance = compute(new_sequence, new_selected)
        delta = new_total_distance - current_total_distance

        if delta < 0:
            # 接受新解
            current_selected = new_selected
            current_sequence = new_sequence
            current_total_distance = new_total_distance

            if new_total_distance < min_dist:
                min_dist = new_total_distance
                best_points = new_selected[:]
                best_order = new_sequence[:]
        else:
            # 以一定概率接受更差的解
            prob = beta.exp(-delta / temperature)
            if alpha.random() < prob:
                current_selected = new_selected
                current_sequence = new_sequence
                current_total_distance = new_total_distance

        temperature *= reduction

    # 使用2-opt算法优化最佳路径
    best_order, min_dist = two_opt(best_order, best_points)

    output = '[' + ', '.join([str(best_points[v]) for v in best_order]) + ']'
    print(output)

    print(f'Best total distance: {min_dist}')


if __name__ == "__main__":
    SA(SAParam(
        start=1000,
        end=0.01,
        reduction=0.9998
    ))
