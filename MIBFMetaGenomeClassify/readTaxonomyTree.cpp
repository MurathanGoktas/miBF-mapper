
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

int main(){
	std::map<std::string,std::tuple<std::string,std::string>> names_map; // id -> name, name type
	std::map<std::string,std::tuple<std::string,std::string>> nodes_map; // id -> name, name type
	std::string taxonomy_tree_prefix = "/projects/btl_scratch/tgoktas/miBF-LINKS-project/tests/databases/refseq_with_annotations/taxonomy";

	std::ifstream names_file(taxonomy_tree_prefix + "/names.dmp");
	std::ifstream nodes_file(taxonomy_tree_prefix + "/nodes.dmp");

/*
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
*/
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