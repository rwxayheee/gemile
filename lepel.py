import networkx as nx

def find_longest_valid_path(graph):
    longest_path = []

    # Helper function to perform DFS
    def dfs(current_node, current_path, visited_edges, node_visits):
        nonlocal longest_path

        # Check if the path is valid (even length) and if it's the longest valid path so far
        if (len(current_path) == len(node_visits) - 1 and  # Path length
            len(current_path) - 1 > len(longest_path) and  # Longer than previous longest
            all(v <= (2 if node == current_path[0] else 1) for node, v in node_visits.items()) and  # Start node can be visited twice
            all(graph.degree[node] == 1 if node_visits[node] == 2 else True for node in current_path)  # Degree 1 nodes as start/end
           ):
            # Update longest path, making sure to copy current_path
            longest_path = list(current_path)
            print(f"New Longest Valid Path Found: {longest_path}")

        # Traverse neighbors
        for neighbor in sorted(graph.neighbors(current_node)):  # Sorted for consistent traversal
            edge = frozenset([current_node, neighbor])  # Treat edge as unordered set

            if edge not in visited_edges:
                # Visit edge
                visited_edges.add(edge)
                current_path.append(neighbor)
                node_visits[neighbor] += 1

                # Recursive DFS
                dfs(neighbor, current_path, visited_edges, node_visits)

                # Backtrack
                visited_edges.remove(edge)
                node_visits[neighbor] -= 1
                current_path.pop()

        # Special condition: check if returning to the start node with even length
        if (current_node == current_path[0] and
            len(current_path) % 2 == 0 and
            len(current_path) - 1 > len(longest_path)):
            longest_path = list(current_path)
            print(f"Loop Back to Start with Even Length: {longest_path}")

    # Initialize node visit counts
    node_visits = {node: 0 for node in graph.nodes}

    # Start DFS from all nodes (since start/end nodes can be visited differently)
    for start_node in graph.nodes:
        node_visits[start_node] += 1
        dfs(start_node, [start_node], set(), node_visits)
        node_visits[start_node] -= 1

    return longest_path


def find_and_sort_isolated_subgraphs(graph):
    # Find all connected components (isolated subgraphs) in the graph
    connected_components = list(nx.connected_components(graph))
    
    # Each component is a set of nodes representing an isolated subgraph
    isolated_subgraphs = [graph.subgraph(component).copy() for component in connected_components]
    
    # Sort isolated subgraphs by the number of nodes from high to low
    sorted_subgraphs = sorted(isolated_subgraphs, key=lambda subgraph: len(subgraph.nodes()), reverse=True)
    
    return sorted_subgraphs


# Create the graph and edges
G = nx.Graph()
edges = [
    (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1),
    (1, 7), (7, 8), (8, 9), (9, 10), (10, 2),
    (10, 11), (11, 12), (12, 13), (13, 3), 
    (13, 14)
]
G.add_edges_from(edges)

# Find the longest valid path
corrected_longest_path = find_longest_valid_path(G)
print("Final Longest Valid Path:", corrected_longest_path)
