47a48,79
> /// コンストラクタ(ファイルロード用)
> BCMOctree::BCMOctree(RootGrid* rootGrid, const std::vector<Pedigree>& pedigrees)
>  : rootGrid(rootGrid), divider(0), ordering(PEDIGREELIST)
> {
> 	using namespace std;
> 	int nRoot = rootGrid->getSize();
> 	rootNodes = new Node*[nRoot];
> 	for(int i = 0; i < nRoot; i++){ rootNodes[i] = new Node(i); }
> 	
> 	leafNodeArray.clear();
> 	leafNodeArray.reserve(pedigrees.size());
> 
> 	// Build Octree form Pedigree List
> 	for(vector<Pedigree>::const_iterator ped = pedigrees.begin(); ped != pedigrees.end(); ++ped) {
> 		unsigned int rootId = ped->getRootID();
> 		Node* node = rootNodes[rootId];
> 		for(int l = 1; l <= ped->getLevel(); l++) {
> 			if(node->isLeafNode()) {
> 				node->makeChildNodes();
> 				for(int cid = 0; cid < 8; cid++) {
> 					node->getChild(cid)->setActive(false);
> 				}
> 			}
> 			int cid = ped->getChildId(l);
> 			node = node->getChild(cid);
> 		}
> 		node->setActive(true);
> 		node->setBlockID(leafNodeArray.size());
> 		leafNodeArray.push_back(node);
> 	}
> 
> }
52c84,85
<   for (int id = 0; id < rootGrid->getSize(); id++) deleteNode(rootNodes[id]);
---
>   //for (int id = 0; id < rootGrid->getSize(); id++) deleteNode(rootNodes[id]);
>   for (int id = 0; id < rootGrid->getSize(); id++) delete rootNodes[id];
384c417,419
<     if (!neighbor) continue;
---
>     if (!neighbor){
> 		continue;
> 	}
