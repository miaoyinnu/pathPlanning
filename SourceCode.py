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


def SA(config: SAParam) -> None:
    collection = []
    for row in [eval(w) for w in input().split("@")]:
        collection.append([Delta(x, y) for x, y in row])

    total = len(collection)

    alpha.seed(gamma.time())
    min_dist = beta.inf
    best_points = []
    best_order = []

    temperature = config.start
    epsilon = config.end
    reduction = config.reduction

    def make_matrix(selected):
        ret = [[beta.hypot(selected[i].x - selected[j].x, selected[i].y - selected[j].y) for j in range(total)] for i in range(total)]
        return ret

    while temperature > epsilon:
        temperature *= reduction

        selected = []
        for i in range(total):
            idx = alpha.randint(0, len(collection[i]) - 1)
            selected.append(collection[i][idx])

        if len(set(selected)) != total:
            continue

        distance_matrix = make_matrix(selected)

        # 使用启发式方法改进序列生成
        sequence = [0]
        visited = [False] * total
        visited[0] = True

        for _ in range(1, total):
            current = sequence[-1]
            next_idx = min((j for j in range(total) if not visited[j]), key=lambda x: distance_matrix[current][x])
            sequence.append(next_idx)
            visited[next_idx] = True

        # 2-opt 优化
        for _ in range(5):  # 限制优化次数
            for i in range(1, total - 2):
                for j in range(i + 1, total - 1):
                    if j - i == 1:
                        continue
                    new_seq = sequence[:]
                    new_seq[i:j] = sequence[j - 1:i - 1:-1]
                    if compute(new_seq, selected) < compute(sequence, selected):
                        sequence = new_seq

        total_distance = compute(sequence, selected)

        if total_distance < min_dist:
            min_dist = total_distance
            best_points = selected
            best_order = sequence

    output = '[' + ', '.join([str(best_points[v]) for v in best_order]) + ']'
    print(output)
    print(f'Best total distance: {min_dist}')



if __name__ == "__main__":
    SA(SAParam(
        start=100,
        end=0.01,
        reduction=0.9998
    ))