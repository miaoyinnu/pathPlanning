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
    max_iterations = int(1e5)

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
