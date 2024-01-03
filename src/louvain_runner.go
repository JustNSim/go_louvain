package main

import (
	"flag"
	"fmt"
	"time"

	"./louvain"
)

func main() {
	var inputFilename = flag.String("i", "", "Input graph filename for loouvain. [Required]")
	var showCommunityIdOfEachLayer = flag.Bool("l", false, "Show community id of each layer as result. Default setting is false.")
	flag.Parse()

	if *inputFilename == "" {
		fmt.Print("Input filename [-i] is required.")
		return
	}

	startTime1 := time.Now()
	graphReader := louvain.NewGraphReader()
	graph := graphReader.Load(*inputFilename)

	startTime2 := time.Now()
	louvain := louvain.NewLouvain(graph)
	louvain.Compute()

	fmt.Printf("Best Modularity Q: %f\n", louvain.BestModularity())
	fmt.Printf("Number of nodes: %d\n", graph.GetNodeSize())
	fmt.Printf("Number of communities: %d\n", louvain.GetCommunitiesNum())

	if *showCommunityIdOfEachLayer == false {
		fmt.Printf("Nodes to communities.\n")

		//打印每个节点所属的社区
		//nodeToCommunity, nodeNum := louvain.GetBestPertition()
		// for nodeId, commId := range nodeToCommunity {
		// 	fmt.Printf("nodeId: %s communityId: %d \n", graphReader.GetNodeLabel(nodeId), commId)
		// }

		//louvain.GetBestPertition()

		//打印每个社区的节点数
		_, nodeNum := louvain.GetBestPertition()
		for commId, nodeNum := range nodeNum {
			fmt.Printf("commId: %d nodeNum: %d \n", commId, nodeNum)
		}
	} else {
		fmt.Println("[NodeId] [CommunityId in each layer]")
		for nodeId := 0; nodeId != graph.GetNodeSize(); nodeId++ {
			fmt.Print(graphReader.GetNodeLabel(nodeId) + " ")
			fmt.Println(louvain.GetNodeToCommunityInEachLevel(nodeId))
		}
	}

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
