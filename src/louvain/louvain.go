package louvain

type Community struct {
	inWeight    WeightType //内部权重2倍加上自环边1倍
	totalWeight WeightType //等于inWeight+连接外部节点的权重
}

type Level struct {
	graph         Graph
	communities   []Community //社区的社区信息
	inCommunities []int       //节点所属社区
}

type Louvain struct {
	level   []Level
	current *Level
	m2      WeightType
	nodeNum []int
}

func NewLouvain(graph Graph) *Louvain {
	louvain := Louvain{}
	louvain.level = make([]Level, 1)
	louvain.current = &louvain.level[0]
	louvain.current.graph = graph
	louvain.current.communities = make([]Community, louvain.current.graph.GetNodeSize())
	louvain.current.inCommunities = make([]int, louvain.current.graph.GetNodeSize())
	louvain.m2 = WeightType(0.0)
	for nodeId := 0; nodeId < louvain.current.graph.GetNodeSize(); nodeId++ {
		louvain.current.inCommunities[nodeId] = nodeId
		neigh := WeightType(0.0)
		for _, edge := range louvain.current.graph.GetIncidentEdges(nodeId) {
			neigh += edge.weight
		}
		self := WeightType(louvain.current.graph.GetSelfWeight(nodeId))
		louvain.current.communities[nodeId] = Community{self, neigh + self}
		louvain.m2 += (neigh + 2*self)
	}
	return &louvain
}

func (this *Louvain) BestModularity() WeightType {
	return this.Modularity(len(this.level) - 1)
}

func (this *Louvain) Modularity(level int) WeightType {
	return this.level[level].modularity(this.m2)
}

func (this *Level) modularity(m2 WeightType) WeightType {
	q := WeightType(0.0)
	for _, comm := range this.communities {
		q += comm.inWeight/m2 - (comm.totalWeight/m2)*(comm.totalWeight/m2)
	}
	return q
}

// 合并社区(一次层级内社区合并，要移动的节点的所属分片标签inCommunity变为移动到的分片的分片编号)
func (this *Louvain) merge() bool {
	improved := false

	q := make([]int, this.current.graph.GetNodeSize())
	mark := make([]bool, this.current.graph.GetNodeSize())
	for nodeId := 0; nodeId < this.current.graph.GetNodeSize(); nodeId++ {
		q[nodeId] = nodeId
	}

	for len(q) != 0 { //q为判断队列，对所有节点进行操作，加入到deltaQ最大的社区
		nodeId := q[0]
		q = q[1:] // pop_front，弹出一个节点进行操作
		mark[nodeId] = true

		neighWeights := map[int]WeightType{}

		self := WeightType(this.current.graph.GetSelfWeight(nodeId))
		totalWeight := self

		neighWeightsKeys := make([]int, 0, len(neighWeights))
		for _, edge := range this.current.graph.GetIncidentEdges(nodeId) {
			destCommId := this.current.inCommunities[edge.destId]
			if _, exists := neighWeights[destCommId]; !exists {
				neighWeightsKeys = append(neighWeightsKeys, destCommId)
			}
			neighWeights[destCommId] += edge.weight
			totalWeight += edge.weight
		}

		prevCommunity := this.current.inCommunities[nodeId]
		prevNeighWeight := WeightType(neighWeights[prevCommunity])
		this.remove(nodeId, prevCommunity, 2*prevNeighWeight+self, totalWeight)

		maxInc := WeightType(0.0) //max_detaQ
		bestCommunity := prevCommunity
		bestNeighWeight := WeightType(prevNeighWeight) //指与当前节点相连的社区的权重
		for _, community := range neighWeightsKeys {
			weight := neighWeights[community]
			inc := WeightType(weight - this.current.communities[community].totalWeight*totalWeight/this.m2) //detaQ
			if inc > maxInc {
				maxInc = inc
				bestCommunity = community
				bestNeighWeight = weight
			}
		}

		this.insert(nodeId, bestCommunity, 2*bestNeighWeight+self, totalWeight)

		if bestCommunity != prevCommunity {
			improved = true
			for _, edge := range this.current.graph.GetIncidentEdges(nodeId) { //移动节点的邻居需要重新判断deltaQ，因此加入判断队列q
				if mark[edge.destId] {
					q = append(q, edge.destId)
					mark[edge.destId] = false
				}
			}
		}
	}

	return improved
}

func (this *Louvain) insert(nodeId int, community int, inWeight WeightType, totalWeight WeightType) {
	this.current.inCommunities[nodeId] = community
	this.current.communities[community].inWeight += inWeight
	this.current.communities[community].totalWeight += totalWeight
}

func (this *Louvain) remove(nodeId int, community int, inWeight WeightType, totalWeight WeightType) {
	this.current.inCommunities[nodeId] = -1
	this.current.communities[community].inWeight -= inWeight
	this.current.communities[community].totalWeight -= totalWeight
}

func (this *Louvain) rebuild() {

	renumbers := map[int]int{}
	num := 0 //新社区的数量
	//循环前移动节点的所属分片标签inCommunity为移动到的分片的分片编号，这个循环的目的是为了给新的社区重新编号
	for nodeId, inCommunity := range this.current.inCommunities {
		if commId, exists := renumbers[inCommunity]; !exists {
			//fmt.Println("1-nodeId", nodeId, "inCommunity:", inCommunity, " num:", num)
			renumbers[inCommunity] = num
			this.current.inCommunities[nodeId] = num
			num++
		} else {
			//fmt.Println("2-nodeId", nodeId, "inCommunity:", inCommunity, " num:", num)
			this.current.inCommunities[nodeId] = commId
		}
	}

	//移动旧编号社区的社区信息到新编号社区
	newCommunities := make([]Community, num)
	for nodeId := 0; nodeId < len(this.current.communities); nodeId++ {
		//fmt.Println("nodeId:", nodeId, " inCommunity:", this.current.inCommunities[nodeId])
		if comm, exists := renumbers[nodeId]; exists {
			newCommunities[comm] = this.current.communities[nodeId]
			//fmt.Println(("comm:"), comm)
		}
	}

	//二维切片，表示新社区所含的节点ID
	communityNodes := make([][]int, len(newCommunities))
	for nodeId := 0; nodeId < this.current.graph.GetNodeSize(); nodeId++ {
		communityNodes[this.current.inCommunities[nodeId]] = append(communityNodes[this.current.inCommunities[nodeId]], nodeId)
	}

	//一个新社区作为一个点，构造新层级的图
	newGraph := Graph{make(Edges, len(newCommunities)), make([]WeightType, len(newCommunities))}

	//新图中添加边
	for commId := 0; commId < newGraph.GetNodeSize(); commId++ {
		newEdges := map[int]WeightType{}
		selfWeight := WeightType(0.0)
		for _, nodeId := range communityNodes[commId] {
			edges := this.current.graph.GetIncidentEdges(nodeId)
			for _, edge := range edges {
				newEdges[this.current.inCommunities[edge.destId]] += edge.weight
			}
			selfWeight += this.current.graph.GetSelfWeight(nodeId)
		}

		newGraph.AddSelfEdge(commId, newEdges[commId]+selfWeight)
		for nodeId, weight := range newEdges {
			if nodeId != commId {
				newGraph.AddDirectedEdge(commId, nodeId, weight)
			}
		}
	}

	newInCommunities := make([]int, newGraph.GetNodeSize())
	for nodeId := 0; nodeId < len(newInCommunities); nodeId++ {
		newInCommunities[nodeId] = nodeId
	}

	this.level = append(this.level, Level{newGraph, newCommunities, newInCommunities})
	this.current = &this.level[len(this.level)-1]
}

func (this *Louvain) GetLevel(n int) *Level {
	return &this.level[n]
}

func (this *Louvain) Compute() {
	for this.merge() { //不断的合并重构直到没有节点移动，没有增益
		this.rebuild()
	}
}

func (this *Louvain) GetPertition(level int) ([]int, []int) {
	nodeToCommunity := make([]int, this.level[0].graph.GetNodeSize())
	nodeNum := make([]int, this.GetCommunitiesNum())
	for nodeId, _ := range nodeToCommunity {
		commId := nodeId
		for l := 0; l != level; l++ {
			commId = this.level[l].inCommunities[commId]
		}
		nodeToCommunity[nodeId] = commId
		nodeNum[commId] = nodeNum[commId] + 1
	}
	return nodeToCommunity, nodeNum
}

func (this *Louvain) GetBestPertition() ([]int, []int) {
	return this.GetPertition(len(this.level))
}

func (this *Louvain) GetNodeToCommunityInEachLevel(nodeId int) []int {
	nodeToCommunityInEachLevel := make([]int, len(this.level))
	commId := nodeId
	for l := 0; l != len(this.level); l++ {
		commId = this.level[l].inCommunities[commId]
		nodeToCommunityInEachLevel[l] = commId
	}
	return nodeToCommunityInEachLevel
}

// 输出最后的社区数量
func (this *Louvain) GetCommunitiesNum() int {
	return len(this.current.inCommunities)
}

// 输出最后的社区ID和拥有的节点数量
func (this *Louvain) GetCommunitiesNodeNum() []int {
	communities_num := this.GetCommunitiesNum()
	communities := make([]int, communities_num)
	for nodeId := 0; nodeId < this.current.graph.GetNodeSize(); nodeId++ {
		commId := this.current.inCommunities[nodeId]
		communities[commId] = communities[commId] + 1
	}

	return communities
}

// 返回社区的m2(无自环边的情况下等于边权重的2倍)
func (this *Louvain) GetM2() WeightType {
	return this.m2
}
