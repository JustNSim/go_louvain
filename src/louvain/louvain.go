package louvain

import (
	"fmt"
	"math"
)

type Community struct {
	inWeight    WeightType //内部权重2倍加上自环边1倍
	totalWeight WeightType //等于inWeight+连接外部节点的权重
	inTXnum     int        //社区内部交易数量
	exTXnum     int        //社区外部交易数量
	inShard     int        //所属分片
	nodeNum     int        //社区拥有的账户数量
}

type Level struct {
	graph         Graph
	communities   []Community //社区的社区信息
	inCommunities []int       //节点所属社区
}

type Shard struct {
	nodeNumInShard int
	performance    float32
	processTime    float32
	inTXnum        int
	exTXnum        int
}
type Louvain struct {
	level    []Level
	current  *Level
	m2       WeightType
	shard    []Shard
	shardNum int
	lambda   float32
}

func NewLouvain(graph Graph, shardPerformance []float32, random bool, lambda float32, shardNum int) *Louvain {
	louvain := Louvain{}
	louvain.level = make([]Level, 1)
	louvain.current = &louvain.level[0]
	louvain.current.graph = graph
	louvain.current.communities = make([]Community, louvain.current.graph.GetNodeSize())
	louvain.current.inCommunities = make([]int, louvain.current.graph.GetNodeSize())
	louvain.m2 = WeightType(0.0)
	louvain.shard = make([]Shard, len(shardPerformance))
	louvain.shardNum = shardNum
	louvain.lambda = lambda
	for i := 0; i < len(shardPerformance); i++ {
		louvain.shard[i].performance = shardPerformance[i]
	}
	var round int = 0
	for nodeId := 0; nodeId < louvain.current.graph.GetNodeSize(); nodeId++ {
		louvain.current.inCommunities[nodeId] = nodeId
		neigh := WeightType(0.0)
		for _, edge := range louvain.current.graph.GetIncidentEdges(nodeId) {
			neigh += edge.weight
		}
		self := WeightType(louvain.current.graph.GetSelfWeight(nodeId))
		if !random { //先不分配分片，模块度社区划分后再分配分片
			louvain.current.communities[nodeId] = Community{self, neigh + self, 0, int(neigh), -1, 1}
		} else { //随机均匀分配到各分片
			if len(louvain.shard) == 0 {
				fmt.Println("shard number is 0,err!")
			} else {
				louvain.current.communities[nodeId] = Community{self, neigh + self, 0, int(neigh), round, 1}
				round = (round + 1) % len(louvain.shard)
			}
		}
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
func (this *Louvain) Merge(isModularity bool, nodeToShard []int) (bool, []int) {
	improved := false
	if !isModularity {
		this.current = &this.level[0]
		//打印this.current.inCommunities
		// for i := 0; i < len(this.current.inCommunities); i++ {
		// 	fmt.Printf("nodeId: %d inCommunity: %d  Community.InShard: %d\n", i, this.current.inCommunities[i], this.current.communities[this.current.inCommunities[i]].inShard)
		// }
	}

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

		//neighWeightsKeys存邻居节点所在社区
		neighWeightsKeys := make([]int, 0, len(neighWeights))
		otherShardNeighList := make([]int, 0, len(neighWeights))
		otherShardNeighShardIndex := make([]int, 0, len(neighWeights))
		for _, edge := range this.current.graph.GetIncidentEdges(nodeId) {
			destCommId := this.current.inCommunities[edge.destId]
			if _, exists := neighWeights[destCommId]; !exists {
				neighWeightsKeys = append(neighWeightsKeys, destCommId)
			}
			neighWeights[destCommId] += edge.weight
			totalWeight += edge.weight

			if !isModularity {
				neighShard := nodeToShard[edge.destId]
				myShard := nodeToShard[nodeId]
				if myShard != neighShard {
					otherShardNeighList = append(otherShardNeighList, destCommId)
					otherShardNeighShardIndex = append(otherShardNeighShardIndex, neighShard)
				}
			}

		}

		prevCommunity := this.current.inCommunities[nodeId]
		prevNeighWeight := WeightType(neighWeights[prevCommunity])
		// if !isModularity {
		// 	// fmt.Println("nodeId:", nodeId)
		// 	// fmt.Println("prevCommunity:", prevCommunity)
		// 	// fmt.Println("len of communities:", len(this.current.communities))
		// 	//打印this.current.inCommunities
		// 	// for i := 0; i < len(this.current.inCommunities); i++ {
		// 	// 	fmt.Printf("nodeId: %d inCommunity: %d \n", i, this.current.inCommunities[i])
		// 	// }
		// 	// fmt.Println("shard index:", this.current.communities[prevCommunity].inShard)
		// 	this.shard[this.current.communities[prevCommunity].inShard].nodeNumInShard -= 1
		// }
		this.remove(nodeId, prevCommunity, 2*prevNeighWeight+self, totalWeight)
		this.removeTX(nodeId, prevCommunity, int(prevNeighWeight), int(totalWeight-prevNeighWeight-self))

		maxInc := WeightType(0.0) //max_detaQ
		bestCommunity := prevCommunity
		bestNeighWeight := WeightType(prevNeighWeight) //指与当前节点相连的社区的权重

		if isModularity {
			for _, community := range neighWeightsKeys {
				weight := neighWeights[community]
				inc := WeightType(weight - this.current.communities[community].totalWeight*totalWeight/this.m2) //detaQ
				if inc > maxInc {
					maxInc = inc
					bestCommunity = community
					bestNeighWeight = weight
				}
			}
		} else { //以处理时间为依据判断是否移动节点，
			myShard := nodeToShard[nodeId]
			bestShard := myShard
			myShardProcessTime := this.shard[myShard].processTime
			myShardInTxNum := this.shard[myShard].inTXnum
			myShardExTxNum := this.shard[myShard].exTXnum
			bestShardProcesstime := this.shard[myShard].processTime
			bestShardInTxNum := this.shard[myShard].inTXnum
			bestShardExTxNum := this.shard[myShard].exTXnum
			for index, community := range otherShardNeighList { //此时因为是第一层，所以community就是邻居节点编号
				//判断不在一个分片的邻居节点
				neighShard := otherShardNeighShardIndex[index]
				if myShard != neighShard {
					//计算节点移动到邻居所在社区后，两个分片的处理时间的最大值是否下降，
					nodeToMyShardTXnum := this.CalCommShardTXnum(community, myShard, &this.level[0])
					nodeToOtherShardTXnum := int(totalWeight-self) - nodeToMyShardTXnum
					sourceShardProcessTime := this.shard[myShard].processTime
					destShardProcessTime := this.shard[neighShard].processTime
					var biggerProcessTime float32 = 0.0
					if sourceShardProcessTime > destShardProcessTime {
						biggerProcessTime = sourceShardProcessTime
					} else {
						biggerProcessTime = destShardProcessTime
					}
					newSourceExTXnum := this.shard[myShard].exTXnum - nodeToOtherShardTXnum + nodeToMyShardTXnum
					newSourceInTXnum := this.shard[myShard].inTXnum - nodeToMyShardTXnum
					newSourceProcessTime := (float32(newSourceInTXnum) + this.lambda*float32(newSourceExTXnum)) / this.shard[myShard].performance
					//计算节点移动到邻居所在社区后，目标分片的处理时间
					nodeToDestShardTXnum := this.CalCommShardTXnum(nodeId, neighShard, &this.level[0])
					nodeToNoDestShardTXnum := int(totalWeight-self) - nodeToDestShardTXnum
					newDestExTXnum := this.shard[neighShard].exTXnum + nodeToNoDestShardTXnum - nodeToDestShardTXnum
					newDestInTXnum := this.shard[neighShard].inTXnum + nodeToDestShardTXnum
					newDestProcessTime := (float32(newDestInTXnum) + this.lambda*float32(newDestExTXnum)) / this.shard[neighShard].performance
					var newBiggerProcessTime float32 = 0.0
					if newSourceProcessTime > newDestProcessTime {
						newBiggerProcessTime = newSourceProcessTime
					} else {
						newBiggerProcessTime = newDestProcessTime
					}
					//如果移动节点后，两个分片的处理时间的最大值下降，则移动节点到邻居所在社区
					if newBiggerProcessTime < biggerProcessTime {
						//移动节点到邻居所在社区
						bestCommunity = community
						bestShard = neighShard
						bestNeighWeight = neighWeights[community]
						myShardProcessTime = newSourceProcessTime
						myShardInTxNum = newSourceInTXnum
						myShardExTxNum = newSourceExTXnum
						bestShardInTxNum = newDestInTXnum
						bestShardExTxNum = newDestExTXnum
						bestShardProcesstime = newDestProcessTime
					}
				}
			}
			if bestCommunity != prevCommunity { //节点移动后更新分片的处理时间、节点所属分片
				this.shard[bestShard].inTXnum = bestShardInTxNum
				this.shard[bestShard].exTXnum = bestShardExTxNum
				this.shard[bestShard].processTime = bestShardProcesstime
				this.shard[bestShard].nodeNumInShard += 1
				this.shard[myShard].inTXnum = myShardInTxNum
				this.shard[myShard].exTXnum = myShardExTxNum
				this.shard[myShard].processTime = myShardProcessTime
				this.shard[myShard].nodeNumInShard -= 1
				nodeToShard[nodeId] = bestShard
			}
		}

		this.insert(nodeId, bestCommunity, 2*bestNeighWeight+self, totalWeight)
		this.insertTX(nodeId, bestCommunity, int(bestNeighWeight), int(totalWeight-bestNeighWeight-self))

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
	if !isModularity {
		this.current = &this.level[len(this.level)-1]
	}
	return improved, nodeToShard
}

func (this *Louvain) insert(nodeId int, community int, inWeight WeightType, totalWeight WeightType) {
	//this.current.communities[nodeId].inShard = this.current.communities[community].inShard
	this.current.inCommunities[nodeId] = community
	this.current.communities[community].inWeight += inWeight
	this.current.communities[community].totalWeight += totalWeight
}

func (this *Louvain) insertTX(nodeId int, community int, toThisCommTXnum int, toOtherCommTXnum int) {
	this.current.communities[community].inTXnum += toThisCommTXnum
	this.current.communities[community].exTXnum = this.current.communities[community].exTXnum + toOtherCommTXnum - toThisCommTXnum
}

func (this *Louvain) remove(nodeId int, community int, inWeight WeightType, totalWeight WeightType) {
	this.current.inCommunities[nodeId] = -1
	this.current.communities[community].inWeight -= inWeight
	this.current.communities[community].totalWeight -= totalWeight
}

func (this *Louvain) removeTX(nodeId int, community int, toThisCommTXnum int, toOtherCommTXnum int) {
	this.current.communities[community].inTXnum -= toThisCommTXnum
	this.current.communities[community].exTXnum = this.current.communities[community].exTXnum - toOtherCommTXnum + toThisCommTXnum
}

func (this *Louvain) rebuild() { //添加 节点的inShard更新，好像复制时已经复制了

	renumbers := map[int]int{}
	num := 0 //新社区的数量
	//循环前移动节点的所属分片标签inCommunity为移动到的社区的社区编号，这个循环的目的是为了给新的社区重新编号
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

func (this *Louvain) Compute(isModularity bool) {
	nodeToShard := make([]int, 1)
	improved, _ := this.Merge(isModularity, nodeToShard) //在模块度划分时不会用上
	for improved {                                       //不断的合并重构直到没有节点移动，没有增益
		this.rebuild()
		improved, _ = this.Merge(isModularity, nodeToShard)
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
		this.current.communities[commId].nodeNum += 1
		nodeNum[commId] = nodeNum[commId] + 1
		//this.level[0].communities[commId].inShard = this.current.communities[commId].inShard
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

// 返回社区的m2(无自环边的情况下等于边权重的2倍)
func (this *Louvain) GetM2() WeightType {
	return this.m2
}

// 返回最后的社区数量
func (this *Louvain) GetCommunitiesNum() int {
	return len(this.current.inCommunities)
}

// 返回最后的社区ID和拥有的节点数量,用getbestpartition
func (this *Louvain) GetCommunitiesNodeNum() ([]int, []int) {
	_, commNodeNum := this.GetBestPertition()
	communities_num := len(commNodeNum)
	//对数组communities_num中各标签对应的值进行排序，得到的标签顺序存在descRank中
	descRank := make([]int, communities_num)
	descNodeNum := make([]int, communities_num)
	sourceIndex := make([]int, communities_num)
	for i := 0; i < communities_num; i++ {
		sourceIndex[i] = i
	}
	for i := 0; i < communities_num; i++ {
		max := commNodeNum[i]
		maxIndex := sourceIndex[i]
		Index := i
		for j := i + 1; j < communities_num; j++ {
			if commNodeNum[j] > max {
				max = commNodeNum[j]
				maxIndex = sourceIndex[j]
				Index = j
			}
		}
		descRank[i] = maxIndex
		descNodeNum[i] = max
		commNodeNum[Index] = commNodeNum[i]
		sourceIndex[Index] = sourceIndex[i]
	}

	// 打印出commNodeNum和descRank
	//fmt.Println(descNodeNum)
	//fmt.Println(descRank)
	return descNodeNum, descRank
}

// 在分片内只有一个社区时进行计算分片的处理时间
func (this *Louvain) CalShardProcessTime(lambda float32) (int, int) {
	var maxProcessTimeShard int = -1
	var maxProcessTime float32 = 0.0
	var minProcessTimeShard int = -1
	var minProcessTime float32 = math.MaxFloat32
	for _, comm := range this.current.communities {
		if comm.inShard == -1 {
			continue
		}
		this.shard[comm.inShard].processTime = (float32(comm.inTXnum) + lambda*float32(comm.exTXnum)) / this.shard[comm.inShard].performance
		if this.shard[comm.inShard].processTime > maxProcessTime {
			maxProcessTime = this.shard[comm.inShard].processTime
			maxProcessTimeShard = comm.inShard
			// if minProcessTime == 0.0 {
			// 	minProcessTime = maxProcessTime
			// 	minProcessTimeShard = comm.inShard
			// }
		}
		if this.shard[comm.inShard].processTime < minProcessTime {
			minProcessTime = this.shard[comm.inShard].processTime
			minProcessTimeShard = comm.inShard
		}
	}
	return maxProcessTimeShard, minProcessTimeShard
}

// 寻找处理时间最短的分片
func (this *Louvain) FindMinProcessTimeShard() int {
	var minProcessTimeShard int = -1
	var minProcessTime float32 = math.MaxFloat32
	for index, shard := range this.shard {
		if shard.processTime < minProcessTime {
			minProcessTime = shard.processTime
			minProcessTimeShard = index
		}
	}
	return minProcessTimeShard
}

// 计算某个社区与某个分片的交易数量
func (this *Louvain) CalCommShardTXnum(community int, shard int, levelIndex *Level) int {
	var commToShardTXnum int = 0
	for _, edge := range levelIndex.graph.GetIncidentEdges(community) {
		if levelIndex.communities[edge.destId].inShard == shard {
			commToShardTXnum += int(edge.weight)
		}
	}
	return commToShardTXnum
}

// 分配社区到分片。将节点数最多的社区分配到前几个分片，如果社区数多于分片数，后续社区加到处理时间最短的分片
func (this *Louvain) CommToShard(lambda float32) {
	//返回最后一层（current）的社区按拥有的节点数降序排列，descRank是社区原有的标签，descNodeNum是社区拥有的节点数
	descNodeNum, descRank := this.GetCommunitiesNodeNum()
	//打印descNodeNum和descRank
	//fmt.Println("descNodeNum", descNodeNum)
	//fmt.Println("descRank", descRank)
	if len(descRank) < len(this.shard) {
		fmt.Println("社区划分后的社区数小于分片数！会导致最后几个分片分配到的账户为空")
		//将社区分到前几个分片
		fmt.Println("len(descRank):", len(descRank))
		for i := 0; i < len(descRank); i++ {
			this.current.communities[descRank[i]].inShard = i
			fmt.Println("descRank[i]:", descRank[i], " i:", i)
			//更新分片的节点数、内部交易、外部交易、处理时间
			this.shard[i].nodeNumInShard = descNodeNum[i]
			this.shard[i].inTXnum = this.current.communities[descRank[i]].inTXnum
			this.shard[i].exTXnum = this.current.communities[descRank[i]].exTXnum
			this.shard[i].processTime = (float32(this.shard[i].inTXnum) + lambda*float32(this.shard[i].exTXnum)) / this.shard[i].performance
		}
		for i := len(descRank); i < len(this.shard); i++ {
			this.shard[i].processTime = 0.0
			this.shard[i].nodeNumInShard = 0
			this.shard[i].inTXnum = 0
			this.shard[i].exTXnum = 0
		}

	} else { //分片的节点数和之后打印的不一样，以及为什么第三步merge没有节点移动
		for i := 0; i < len(this.shard); i++ {
			this.current.communities[descRank[i]].inShard = i
			this.shard[i].nodeNumInShard = descNodeNum[i]
			//fmt.Println("commID:", descRank[i], "shardID i:", i, "commNodeNum:", descNodeNum[i])
			this.shard[i].inTXnum = this.current.communities[descRank[i]].inTXnum
			this.shard[i].exTXnum = this.current.communities[descRank[i]].exTXnum
			this.shard[i].processTime = (float32(this.shard[i].inTXnum) + lambda*float32(this.shard[i].exTXnum)) / this.shard[i].performance
		}
		//_, minProcessTimeShard := this.CalShardProcessTime(lambda)
		minProcessTimeShard := this.FindMinProcessTimeShard()
		//打印descRank
		// for i := 0; i < len(descRank); i++ {
		// 	fmt.Println("i:", i, "     descRank[i]:", descRank[i])
		// }
		for i := len(this.shard); i < len(descRank); i++ {
			commID := descRank[i]
			//加入到处理时间最短的分片
			this.current.communities[commID].inShard = minProcessTimeShard
			//计算该社区与要加入的分片的交易数量
			commToShardTXnum := this.CalCommShardTXnum(commID, minProcessTimeShard, this.current)
			//更新分片的交易数量
			//fmt.Println("commID:", commID, " minProcessTimeShard:", minProcessTimeShard)
			this.shard[minProcessTimeShard].inTXnum += commToShardTXnum + this.current.communities[commID].inTXnum
			this.shard[minProcessTimeShard].exTXnum += this.current.communities[commID].exTXnum - commToShardTXnum
			//更新分片的节点数量、处理时间
			this.shard[minProcessTimeShard].nodeNumInShard += descNodeNum[i]
			//fmt.Println("communityID=descRank[i]:", descRank[i], "shardID :", minProcessTimeShard, "commNodeNum:", descNodeNum[i])
			this.shard[minProcessTimeShard].processTime = (float32(this.shard[minProcessTimeShard].inTXnum) +
				lambda*float32(this.shard[minProcessTimeShard].exTXnum)) / this.shard[minProcessTimeShard].performance
			//更新最短处理时间的分片
			minProcessTimeShard = this.FindMinProcessTimeShard()

		}
	}
}

// 获取louvain的lambda
func (this *Louvain) GetLambda() float32 {
	return this.lambda
}

// 打印分片的处理时间和节点数量
func (this *Louvain) PrintShard() {
	for i := 0; i < len(this.shard); i++ {
		fmt.Printf("shard: %d nodeNumInShard: %d processTime: %f \n", i, this.shard[i].nodeNumInShard, this.shard[i].processTime)
	}
}

// 根据节点-分片映射表得到各分片包含的节点列表
func (this *Louvain) GetShardNodeList(nodeToShard []int) [][]int {
	shardNodeList := make([][]int, len(this.shard))
	for i := 0; i < len(nodeToShard); i++ {
		shardId := nodeToShard[i]
		shardNodeList[shardId] = append(shardNodeList[shardId], i)
	}
	return shardNodeList
}

// 获取louvain的lambda
func (this *Louvain) GetLevelNum() int {
	return len(this.level)
}

// 将节点在levelIndex(0开始计数）层所在社区所属的分片标签传递到第一层
func (this *Louvain) UpdateNodeShardIndex() ([]int, []int) {
	levelIndex := len(this.level) - 1
	nodeToCommunity := make([]int, this.level[0].graph.GetNodeSize())
	nodeToShard := make([]int, this.level[0].graph.GetNodeSize())
	//nodeNum := make([]int, this.GetCommunitiesNum())
	for nodeId, _ := range nodeToCommunity {
		commId := this.level[0].inCommunities[nodeId]
		firstLevelCommId := commId
		for l := 1; l < len(this.level); l++ {
			commId = this.level[l].inCommunities[commId]
		}
		// nodeToCommunity[nodeId] = commId
		// this.current.communities[commId].nodeNum += 1
		// nodeNum[commId] = nodeNum[commId] + 1
		this.level[0].communities[firstLevelCommId].inShard = this.level[levelIndex].communities[commId].inShard
		nodeToShard[nodeId] = this.level[levelIndex].communities[commId].inShard
	}
	return nodeToCommunity, nodeToShard
}
