
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <algorithm>
//#include "tree.hh"

// get_genus(x)
// are_same_genus(x,y)
// is_ancestor(x,y)
// get_lca(x,y)

class TaxonomyTree{
	public:
		TaxonomyTree(){}
		
		TaxonomyTree(std::string taxonomy_tree_prefix){

			std::ifstream names_file(taxonomy_tree_prefix + "/names.dmp");
			std::ifstream nodes_file(taxonomy_tree_prefix + "/nodes.dmp");

			
			// --- read names file to map ---
			for( std::string line; getline( names_file, line ); )
			{
				std::vector<std::string> tokens;
				std::istringstream iss(line);
				std::string token;
				while(std::getline(iss, token, '\t')){   // but we can specify a different one
					tokens.push_back(token);
				}
				if(!(tokens[6] == "scientific name" || tokens[6] == "genbank common name")){
					continue;
				}
				//std::cout << tokens[0] << std::endl;
				//std::cout << tokens[2] << std::endl;
				//std::cout << tokens[6] << std::endl;
				names_map[tokens[0]] = std::make_pair(tokens[2],tokens[6]);
			}
			// --- read nodes file to map ---
			for( std::string line; getline( nodes_file, line ); )
			{
				std::vector<std::string> tokens;
				std::istringstream iss(line);
				std::string token;
				while(std::getline(iss, token, '\t')){   // but we can specify a different one
					tokens.push_back(token);
				}
				nodes_map[tokens[0]] = std::make_pair(tokens[2],tokens[4]);
			}
		}

	std::string get_parent(std::string node){
		return std::get<0>(nodes_map[node]);
	}
	std::string get_genus(std::string node){
		while(std::get<0>(nodes_map[node]) != "1"){
			if(std::get<1>(nodes_map[node]) == "genus"){
				return std::get<0>(nodes_map[node]);
			} else{
				node = std::get<0>(nodes_map[node]); 
			}
		}
		return "-1";
	}
	// are_same_genus(x,y)
	bool are_same_genus(std::string node_x, std::string node_y){
		if(get_genus(node_x) == get_genus(node_y) && get_genus(node_x) != "-1"){
			return true;
		}
		return false;
	}
	bool is_ancestor(std::string node_x, std::string node_y){ // whether x is ancestor of y
		if(node_x == node_y){
			return true;
		} // would makes sense to keep though. But commenting now for debug (counting) purposes
		while(std::get<0>(nodes_map[node_y]) != "1"){
			if(std::get<0>(nodes_map[node_y]) == node_x){
				return true;
			} else{
				node_y = std::get<0>(nodes_map[node_y]); 
			}
		}
		return false;
	}
	std::string get_lca(std::string node_x, std::string node_y){
		std::vector<std::string> node_x_ancestors;
		while(true){
			node_x_ancestors.push_back(node_x);
			if(node_x == "1"){
				break;
			}
			node_x = std::get<0>(nodes_map[node_x]);
		}
		while(true){
			if(std::find(node_x_ancestors.begin(), node_x_ancestors.end(), node_y) != node_x_ancestors.end()){
				return node_y;
			}
			if(node_y == "1"){
				break;
			}
			node_y = std::get<0>(nodes_map[node_y]);
		}
	}
	std::string get_node_id_from_name(std::string name){
		for (auto it = names_map.begin(); it != names_map.end(); it++){
			if(std::get<0>(it->second) == name){
				return it->first;
			}
		}
	}
	private:
		std::map<std::string,std::tuple<std::string,std::string>> names_map; // id -> name, name type
		std::map<std::string,std::tuple<std::string,std::string>> nodes_map; // id -> parent id, taxonomic rank
};
/*
int main(){
	std::map<std::string,std::tuple<std::string,std::string>> names_map; // id -> name, name type
	std::map<std::string,std::tuple<std::string,std::string>> nodes_map; // id -> parent id, taxonomic rank
	std::string taxonomy_tree_prefix = "/projects/btl_scratch/tgoktas/miBF-LINKS-project/tests/databases/refseq_with_annotations/taxonomy";

	TaxonomyTree t(taxonomy_tree_prefix);
	//std::cout << "here 1" << std::endl;
	//std::cout << (t.get_parent("20")) << std::endl;
	//std::cout << (t.get_genus("33")) << std::endl;
	std::cout << (t.are_same_genus("1262452","1241582")) << std::endl;
	//std::cout << (t.is_ancestor("31","34")) << std::endl;
	std::cout << (t.get_lca("1262452","1241582")) << std::endl;
	//std::cout << (t.get_node_id_from_name("Buchnera aphidicola")) << std::endl;

}
*/