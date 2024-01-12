package main

import (
	"flag"
	"fmt"
	"time"

	"./louvain"
)

func main() {
	var inputFilename = flag.String("i", "", "Input graph filename for loouvain. [Required]")
	//var showCommunityIdOfEachLayer = flag.Bool("l", false, "Show community id of each layer as result. Default setting is false.")
	flag.Parse()

	if *inputFilename == "" {
		fmt.Print("Input filename [-i] is required.")
		return
	}

	startTime1 := time.Now()
	graphReader := louvain.NewGraphReader()
	graph := graphReader.Load(*inputFilename)

	startTime2 := time.Now()
	//设置每个社区性能值
	var shardNum = 4
	shardPerformance := []float32{250, 250, 200, 100}
	var israndom bool = false
	var lambda float32 = 2.0

	louvain := louvain.NewLouvain(graph, shardPerformance, israndom, lambda, shardNum)
	var isModularity bool = true
	louvain.Compute(isModularity) //社区划分
	//打印社区划分结果
	fmt.Printf("Number of nodes: %d\n", graph.GetNodeSize())
	fmt.Printf("Number of communities: %d\n", louvain.GetCommunitiesNum())
	//打印每个节点所属的社区,GetBestPertition()中梳理了节点最后对应的社区。更新了节点在第一层的社区所属的分片
	// nodeToCommunity, nodeNum := louvain.GetBestPertition()
	// for nodeId, commId := range nodeToCommunity {
	// 	fmt.Printf("nodeId: %s communityId: %d \n", graphReader.GetNodeLabel(nodeId), commId)
	// }
	// //打印每个社区的节点数
	// for commId, nodeNum := range nodeNum {
	// 	fmt.Printf("commId: %d nodeNum: %d \n", commId, nodeNum)
	// }

	//将社区分配到分片
	fmt.Println("---------start the second stage comm to shard")
	louvain.CommToShard(louvain.GetLambda())
	_, nodeToShard := louvain.UpdateNodeShardIndex()
	louvain.PrintShard()
	//按处理时间判断节点是否移动
	isModularity = false
	result := true
	for i := 1; result; i++ {
		fmt.Println("start the third stage node Merge: ", i, " times")
		result, nodeToShard = louvain.Merge(isModularity, nodeToShard)
		if i > 20 {
			break
		}
	}

	louvain.PrintShard()
	//louvain.UpdateNodeShardIndex(0)
	shardNodeList := louvain.GetShardNodeList(nodeToShard)
	//打印二维数组shardNodeList
	for i := 0; i < len(shardNodeList); i++ {
		//fmt.Printf("shard[%d] %d nodes: %v\n", i, len(shardNodeList[i]), shardNodeList[i])
		fmt.Printf("shard[%d] %d nodes\n", i, len(shardNodeList[i]))
	}

	// if *showCommunityIdOfEachLayer == false {
	// 	fmt.Printf("Nodes to communities.\n")

	// 	//打印每个节点所属的社区
	// 	// nodeToCommunity, nodeNum := louvain.GetBestPertition()
	// 	// for nodeId, commId := range nodeToCommunity {
	// 	// 	fmt.Printf("nodeId: %s communityId: %d \n", graphReader.GetNodeLabel(nodeId), commId)
	// 	// }
	// 	//打印每个社区的节点数
	// 	_, nodeNum := louvain.GetBestPertition()
	// 	for commId, nodeNum := range nodeNum {
	// 		fmt.Printf("commId: %d nodeNum: %d \n", commId, nodeNum)
	// 	}
	// } else {
	// 	fmt.Println("[NodeId] [CommunityId in each layer]")
	// 	for nodeId := 0; nodeId != graph.GetNodeSize(); nodeId++ {
	// 		fmt.Print(graphReader.GetNodeLabel(nodeId) + " ")
	// 		fmt.Println(louvain.GetNodeToCommunityInEachLevel(nodeId))
	// 	}
	// }

	endTime := time.Now()
	// 计算时间差
	duration1 := startTime2.Sub(startTime1)
	duration2 := endTime.Sub(startTime2)
	duration3 := endTime.Sub(startTime1)

	// 打印时间差
	fmt.Printf("graph load 耗时: %v\n", duration1)
	fmt.Printf("社区划分耗时: %v\n", duration2)
	fmt.Printf("总耗时: %v\n", duration3)

}
