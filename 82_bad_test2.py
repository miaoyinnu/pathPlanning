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

def three_opt(sequence, selected):
    n = len(sequence)
    best_sequence = sequence[:]
    best_distance = compute(sequence, selected)
    improved = True

    while improved:
        improved = False
        for i in range(n - 2):
            for j in range(i + 1, n - 1):
                for k in range(j + 1, n):
                    new_sequences = [
                        best_sequence[:i+1] + best_sequence[j+1:k+1] + best_sequence[i+1:j+1] + best_sequence[k+1:],
                        best_sequence[:i+1] + best_sequence[j+1:k+1][::-1] + best_sequence[i+1:j+1] + best_sequence[k+1:],
                        best_sequence[:i+1] + best_sequence[j+1:k+1] + best_sequence[i+1:j+1][::-1] + best_sequence[k+1:],
                        best_sequence[:i+1] + best_sequence[i+1:j+1][::-1] + best_sequence[j+1:k+1] + best_sequence[k+1:],
                    ]
                    for new_seq in new_sequences:
                        new_distance = compute(new_seq, selected)
                        if new_distance < best_distance - 1e-6:
                            best_sequence = new_seq
                            best_distance = new_distance
                            improved = True
                            break
                if improved:
                    break
            if improved:
                break
    return best_sequence, best_distance

def optimize_point_selection(selected, sequence, collection):
    improved = True
    while improved:
        improved = False
        for cluster_idx in range(len(collection)):
            current_point = selected[cluster_idx]
            best_point = current_point
            best_distance = compute(sequence, selected)
            for candidate_point in collection[cluster_idx]:
                if candidate_point == current_point or candidate_point in selected:
                    continue
                selected_candidate = selected[:]
                selected_candidate[cluster_idx] = candidate_point
                if len(set(selected_candidate)) != len(selected_candidate):
                    continue
                total_distance = compute(sequence, selected_candidate)
                if total_distance < best_distance - 1e-6:
                    best_point = candidate_point
                    best_distance = total_distance
                    improved = True
            selected[cluster_idx] = best_point
    return selected

def SA(config: SAParam) -> None:
    collection = []
    str1 = "[(41, 31), (42, 31), (43, 31), (44, 31), (45, 31), (45, 31), (45, 32), (45, 33), (45, 34), (45, 35), (45, 36), (45, 36), (44, 36), (43, 36), (42, 36), (41, 36), (41, 36), (41, 35), (41, 34), (41, 33), (41, 32), (41, 31)]@[(48, 36), (49, 36), (50, 36), (51, 36), (52, 36), (53, 36), (53, 36), (53, 37), (53, 38), (53, 39), (53, 40), (53, 41), (53, 42), (53, 42), (52, 42), (51, 42), (50, 42), (49, 42), (48, 42), (48, 42), (48, 41), (48, 40), (48, 39), (48, 38), (48, 37), (48, 36)]@[(35, 15), (36, 15), (37, 15), (38, 15), (38, 15), (38, 16), (38, 17), (38, 18), (38, 19), (38, 20), (38, 20), (37, 20), (36, 20), (35, 20), (35, 20), (35, 19), (35, 18), (35, 17), (35, 16), (35, 15)]@[(17, 35), (18, 35), (19, 35), (20, 35), (21, 35), (22, 35), (22, 35), (22, 36), (22, 37), (22, 38), (22, 39), (22, 39), (21, 39), (20, 39), (19, 39), (18, 39), (17, 39), (17, 39), (17, 38), (17, 37), (17, 36), (17, 35)]@[(8, 36), (9, 36), (10, 36), (11, 36), (12, 36), (12, 36), (12, 37), (12, 38), (12, 39), (12, 40), (12, 41), (12, 41), (11, 41), (10, 41), (9, 41), (8, 41), (8, 41), (8, 40), (8, 39), (8, 38), (8, 37), (8, 36)]@[(9, 49), (10, 49), (11, 49), (12, 49), (12, 49), (12, 50), (12, 51), (12, 52), (12, 52), (11, 52), (10, 52), (9, 52), (9, 52), (9, 51), (9, 50), (9, 49)]@[(10, 8), (11, 8), (12, 8), (13, 8), (14, 8), (14, 8), (14, 9), (14, 10), (14, 11), (14, 12), (14, 12), (13, 12), (12, 12), (11, 12), (10, 12), (10, 12), (10, 11), (10, 10), (10, 9), (10, 8)]@[(52, 21), (51, 23), (49, 25), (47, 26), (44, 25), (42, 23), (42, 21), (42, 18), (44, 16), (47, 16), (49, 16), (51, 18)]@[(49, 10), (48, 12), (46, 13), (44, 13), (42, 13), (41, 11), (41, 8), (42, 6), (44, 6), (46, 6), (48, 7)]@[(27, 18), (26, 20), (23, 22), (20, 22), (17, 20), (17, 18), (17, 15), (20, 13), (23, 13), (26, 15)]@[(37, 48), (38, 48), (39, 48), (40, 48), (40, 48), (40, 49), (40, 50), (40, 51), (40, 51), (39, 51), (38, 51), (37, 51), (37, 51), (37, 50), (37, 49), (37, 48)]@[(21, 55), (22, 55), (23, 55), (24, 55), (24, 55), (24, 56), (24, 57), (24, 58), (24, 58), (23, 58), (22, 58), (21, 58), (21, 58), (21, 57), (21, 56), (21, 55)]@[(5, 5), (6, 5), (7, 5), (8, 5), (8, 6), (8, 7), (8, 8), (7, 8), (6, 8), (5, 8), (5, 7), (5, 6), (5, 5)]@[(20, 20), (21, 20), (22, 20), (23, 20), (23, 21), (23, 22), (23, 23), (22, 23), (21, 23), (20, 23), (20, 22), (20, 21), (20, 20)]@[(30, 10), (31, 10), (32, 10), (33, 10), (33, 11), (33, 12), (33, 13), (32, 13), (31, 13), (30, 13), (30, 12), (30, 11), (30, 10)]@[(12, 35), (13, 35), (14, 35), (15, 35), (15, 36), (15, 37), (15, 38), (14, 38), (13, 38), (12, 38), (12, 37), (12, 36), (12, 35)]@[(40, 45), (41, 45), (42, 45), (43, 45), (43, 46), (43, 47), (43, 48), (42, 48), (41, 48), (40, 48), (40, 47), (40, 46), (40, 45)]@[(15, 25), (16, 25), (17, 25), (18, 25), (18, 26), (18, 27), (18, 28), (17, 28), (16, 28), (15, 28), (15, 27), (15, 26), (15, 25)]@[(50, 50), (51, 50), (52, 50), (53, 50), (53, 51), (53, 52), (53, 53), (52, 53), (51, 53), (50, 53), (50, 52), (50, 51), (50, 50)]@[(35, 15), (36, 15), (37, 15), (38, 15), (38, 16), (38, 17), (38, 18), (37, 18), (36, 18), (35, 18), (35, 17), (35, 16), (35, 15)]@[(25, 25), (26, 25), (27, 25), (28, 25), (28, 26), (28, 27), (28, 28), (27, 28), (26, 28), (25, 28), (25, 27), (25, 26), (25, 25)]"
    str2 = "[(41, 31), (42, 31), (43, 31), (44, 31), (45, 31), (45, 31), (45, 32), (45, 33), (45, 34), (45, 35), (45, 36), (45, 36), (44, 36), (43, 36), (42, 36), (41, 36), (41, 36), (41, 35), (41, 34), (41, 33), (41, 32), (41, 31)]@[(48, 36), (49, 36), (50, 36), (51, 36), (52, 36), (53, 36), (53, 36), (53, 37), (53, 38), (53, 39), (53, 40), (53, 41), (53, 42), (53, 42), (52, 42), (51, 42), (50, 42), (49, 42), (48, 42), (48, 42), (48, 41), (48, 40), (48, 39), (48, 38), (48, 37), (48, 36)]@[(35, 15), (36, 15), (37, 15), (38, 15), (38, 15), (38, 16), (38, 17), (38, 18), (38, 19), (38, 20), (38, 20), (37, 20), (36, 20), (35, 20), (35, 20), (35, 19), (35, 18), (35, 17), (35, 16), (35, 15)]@[(17, 35), (18, 35), (19, 35), (20, 35), (21, 35), (22, 35), (22, 35), (22, 36), (22, 37), (22, 38), (22, 39), (22, 39), (21, 39), (20, 39), (19, 39), (18, 39), (17, 39), (17, 39), (17, 38), (17, 37), (17, 36), (17, 35)]@[(8, 36), (9, 36), (10, 36), (11, 36), (12, 36), (12, 36), (12, 37), (12, 38), (12, 39), (12, 40), (12, 41), (12, 41), (11, 41), (10, 41), (9, 41), (8, 41), (8, 41), (8, 40), (8, 39), (8, 38), (8, 37), (8, 36)]@[(9, 49), (10, 49), (11, 49), (12, 49), (12, 49), (12, 50), (12, 51), (12, 52), (12, 52), (11, 52), (10, 52), (9, 52), (9, 52), (9, 51), (9, 50), (9, 49)]@[(10, 8), (11, 8), (12, 8), (13, 8), (14, 8), (14, 8), (14, 9), (14, 10), (14, 11), (14, 12), (14, 12), (13, 12), (12, 12), (11, 12), (10, 12), (10, 12), (10, 11), (10, 10), (10, 9), (10, 8)]@[(52, 21), (51, 23), (49, 25), (47, 26), (44, 25), (42, 23), (42, 21), (42, 18), (44, 16), (47, 16), (49, 16), (51, 18)]@[(49, 10), (48, 12), (46, 13), (44, 13), (42, 13), (41, 11), (41, 8), (42, 6), (44, 6), (46, 6), (48, 7)]@[(27, 18), (26, 20), (23, 22), (20, 22), (17, 20), (17, 18), (17, 15), (20, 13), (23, 13), (26, 15)]@[(37, 48), (38, 48), (39, 48), (40, 48), (40, 48), (40, 49), (40, 50), (40, 51), (40, 51), (39, 51), (38, 51), (37, 51), (37, 51), (37, 50), (37, 49), (37, 48)]@[(21, 55), (22, 55), (23, 55), (24, 55), (24, 55), (24, 56), (24, 57), (24, 58), (24, 58), (23, 58), (22, 58), (21, 58), (21, 58), (21, 57), (21, 56), (21, 55)]"

    for row in [eval(w) for w in str2.split("@")]:
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

    # Initialize selected points
    current_selected = []
    for i in range(total):
        idx = alpha.randint(0, len(collection[i]) - 1)
        current_selected.append(collection[i][idx])

    # Initialize sequence
    current_sequence = list(range(total))
    alpha.shuffle(current_sequence)

    current_total_distance = compute(current_sequence, current_selected)

    min_dist = current_total_distance
    best_points = current_selected[:]
    best_order = current_sequence[:]

    iteration = 0
    max_iterations = int(1e6)

    accepted = 0
    attempted = 0

    while temperature > epsilon and iteration < max_iterations:
        iteration += 1
        # Generate neighbor solution
        new_selected = current_selected[:]
        new_sequence = current_sequence[:]

        action = alpha.random()
        if action < 0.2:
            # Change the selected point in one cluster
            cluster_idx = alpha.randint(0, total - 1)
            new_idx = alpha.randint(0, len(collection[cluster_idx]) - 1)
            new_selected[cluster_idx] = collection[cluster_idx][new_idx]
        elif action < 0.4:
            # Swap two points in the sequence
            idx1, idx2 = alpha.sample(range(total), 2)
            new_sequence[idx1], new_sequence[idx2] = new_sequence[idx2], new_sequence[idx1]
        elif action < 0.6:
            # Reverse a subsequence
            idx1, idx2 = sorted(alpha.sample(range(total), 2))
            new_sequence[idx1:idx2+1] = reversed(new_sequence[idx1:idx2+1])
        elif action < 0.8:
            # Change multiple selected points
            num_changes = alpha.randint(2, 5)
            clusters_to_change = alpha.sample(range(total), num_changes)
            for cluster_idx in clusters_to_change:
                new_idx = alpha.randint(0, len(collection[cluster_idx]) - 1)
                new_selected[cluster_idx] = collection[cluster_idx][new_idx]
        else:
            # Swap selected points between two clusters
            cluster_idx1, cluster_idx2 = alpha.sample(range(total), 2)
            new_selected[cluster_idx1], new_selected[cluster_idx2] = new_selected[cluster_idx2], new_selected[cluster_idx1]

        # Ensure selected points are unique
        if len(set(new_selected)) != total:
            continue

        new_total_distance = compute(new_sequence, new_selected)
        delta = new_total_distance - current_total_distance

        attempted += 1
        if delta < 0 or alpha.random() < beta.exp(-delta / temperature):
            accepted += 1
            # Accept new solution
            current_selected = new_selected
            current_sequence = new_sequence
            current_total_distance = new_total_distance

            if new_total_distance < min_dist - 1e-6:
                min_dist = new_total_distance
                best_points = new_selected[:]
                best_order = new_sequence[:]
        # Adaptive cooling
        if iteration % 1000 == 0:
            acceptance_rate = accepted / attempted if attempted > 0 else 0
            if acceptance_rate > 0.95:
                temperature *= 0.5  # Increase cooling if acceptance rate is too high
            elif acceptance_rate < 0.05:
                temperature *= 1.05  # Decrease cooling if acceptance rate is too low
            accepted = 0
            attempted = 0

        temperature *= reduction

    # Optimize point selection
    best_points = optimize_point_selection(best_points, best_order, collection)

    # Use 3-opt algorithm to optimize the best tour
    best_order, min_dist = three_opt(best_order, best_points)

    output = '[' + ', '.join([str(best_points[v]) for v in best_order]) + ']'
    print(output)
    print(f'Best total distance: {min_dist}')

if __name__ == "__main__":
    SA(SAParam(
        start=1e5,    # Increased initial temperature
        end=1e-6,     # Lower end temperature
        reduction=0.999995  # Slower cooling schedule
    ))
