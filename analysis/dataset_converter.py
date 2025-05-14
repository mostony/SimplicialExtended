def convert(vertices_path, edges_path, output_path):
    with open(vertices_path, 'r') as f:
        lengths = [int(line.strip()) for line in f]

    with open(edges_path, 'r') as f:
        vertices = [int(line.strip()) for line in f]
    assert sum(lengths) == len(vertices)

    with open(output_path, 'w') as f:
       cur_index = 0
       for size in lengths:
           hyperedge = vertices[cur_index:cur_index+size]
           cur_index += size
           f.write(','.join(map(str, hyperedge)) + '\n')
       

if __name__ == '__main__':
    datasets = ["email-Enron", "email-Eu"]
    for d in datasets:
        input_vertices_path = f"data/{d}/{d}-nverts.txt"
        input_edges_path = f"data/{d}/{d}-simplices.txt"
        output_path = f"data/{d}/hyperedges-{d}.txt"
        convert(input_vertices_path, input_edges_path, output_path)

