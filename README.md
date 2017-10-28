# go louvain

This package implements the louvain algorithm in go.

# usage

```bash
./louvain_runner -i [input file]
```

An example is shown below.
```bash
git clone https://github.com/ken57/go_louvain.git

cd go_louvain
go build src/louvain_runner.go
./louvain_runner -i src/louvain/resource/karate.txt
```

and it will be obtained.
```
Modularity Q: 0.418803
Nodes to communities.
nodeId: 10 communityId: 0
nodeId: 34 communityId: 1
nodeId: 14 communityId: 0
nodeId: 15 communityId: 1
...
```

Louvain method can hierachical structured community detection.
CommunityIds of each layer are obtained by -l option.
```
./louvain_runner -l -i src/louvain/resource/karate.txt
Best Modularity Q: 0.418803
[NodeId] [CommunityId in each layer]
10 [0 0 0]
34 [1 1 1]
14 [0 0 0]
```

# A format of input file

The input file is in the following format.
```
[source],[dest]
```
An example is shown below.
```
1,2
1,3
1,4
1,5
1,6
1,7
...
```

The input file is Interpreted as an undirected graph.
# Test

```
go test ./src/louvain
```
