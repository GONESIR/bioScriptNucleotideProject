def DeBruijn(patterns):
    graph = {}
    for i in range(len(patterns)):
        prefix = patterns[i][:-1]
        suffix = patterns[i][1:]

        if prefix not in graph:
            graph[prefix] = [suffix]
        else:
            graph[prefix].append(suffix)   
    return graph

patterns = []
with open('rosalind_ba3e.txt', 'r') as myFile:
    for line in myFile:
        patterns.append(line.strip())

deBruijnGraph = DeBruijn(patterns)

with open("output.txt", 'w') as writer:
    for key, values in sorted(deBruijnGraph.items()):
        result = ",".join(values)
        writer.write(f"{key} -> {result}\n")