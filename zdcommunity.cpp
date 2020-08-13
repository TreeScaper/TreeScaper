#include "zdcommunity.h"


template<class T>
void print_comm_array(Matrix<T> &arr, int n, File &output, bool arr_is_covariance, double highfreq, double lowfreq, int &covariance_freeid_size, int &covariance_nonfree_id_size, int *covariance_freeid, int *covariance_nonfree_id)
{
	covariance_freeid_size = 0;
	covariance_nonfree_id_size = 0;

	if (arr.get_row() != 0 && arr.get_col() != 0)
	{
		if (arr_is_covariance)
		{
			for (int i = 0; i < n; i++)
			{
				bool crit;
				crit = (arr(i, i) <= lowfreq * (1 - lowfreq) || arr(i, i) <= highfreq * (1 - highfreq));

				if (crit)   //if((*arr)[i][i] != 0)
				{
					covariance_freeid[covariance_freeid_size] = i;
					covariance_freeid_size++;
				}
				else
				{
					covariance_nonfree_id[covariance_nonfree_id_size] = i;
					covariance_nonfree_id_size++;
				}
			}
			output << ">" << endl;
			for (int i = 0; i < covariance_nonfree_id_size; i++)
			{
				output << i << " " << i << endl;
				//                output << covariance_nonfree_id[i] << " " << covariance_nonfree_id[i] << endl;
			}
			output << ">" << endl;
			output << "1 1" << endl;
			output << ">" << endl;
			for (int i = 0; i < covariance_nonfree_id_size; i++)
			{
				for (int j = 0; j < covariance_nonfree_id_size; j++)
				{
					output << j << " " << i << " " << arr(covariance_nonfree_id[j], covariance_nonfree_id[i]) << " " << "1" << endl;
					//                    output << covariance_nonfree_id[j]  << " " << covariance_nonfree_id[i] << " " << (*arr)[covariance_nonfree_id[j]][covariance_nonfree_id[i]] << " " << "1" << endl;
				}
			}
		}
		else
		{
			output << ">" << endl;
			for (int i = 0; i < n; i++)
			{
				output << i << " " << i << endl;
			}
			output << ">" << endl;
			output << "1 1" << endl;
			output << ">" << endl;

			for (int i = 0; i < n; i++)
				covariance_nonfree_id[i] = i;
			covariance_nonfree_id_size = n;

			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					output << j << " " << i << " " << arr(i, j) << " " << "1" << endl;
		}
	}
	output.close();
}

int read_conf(char* filename, int* &conf, int* &sign)
{
	//ifstream finput;
	//finput.open(filename, fstream::in | fstream::binary);
	//if (finput.fail())
	//{
	//	cout << "Could not read config file: " << filename << "\n\n";
	//	exit(-1);
	//}

	File file_Conf(filename);
	if (!file_Conf.is_open())
	{
		cout << "Could not read config file: " << filename << "\n\n";
		exit(-1);
	}
	int pos_header = file_Conf.end_header();
	file_Conf.close();
	ifstream finput;
	finput.open(filename, fstream::in | fstream::binary);
	finput.seekg(pos_header, ios::beg);

	int nb_layers;
	//read number of layers
	finput.read((char *)&nb_layers, sizeof(int));
	conf = (int*)malloc((long)nb_layers * sizeof(int));
	sign = (int*)malloc((long)nb_layers * sizeof(int));

	// read conf
	finput.read((char *)conf, nb_layers * sizeof(int));
	if (finput.fail() || finput.eof())
	{
		cout << "Error!\n\n";
		exit(-1);
		//throw new exception();
	}
	finput.read((char *)sign, nb_layers * sizeof(int));

	return nb_layers;
}

void create_resolution(double lp, double ln, int nb_layers, int* sign, double* &lambda)
{
	for (int i = 0; i < nb_layers; i++)
	{
		if (sign[i] == POSITIVE)
			lambda[i] = lp;
		else if (sign[i] == NEGATIVE)
			lambda[i] = ln;
		else
			cout << "Error: Layer not recognized!\n\n";
	}
}

bool community_detection_automatically(Matrix<double> &mat, map<String, String> &paras) {
	string highfreq = (char*)paras["-hf"];
	string lowfreq = (char*)paras["-lf"];
	int modelType = 0;
	if (paras["-cm"] == (String) "CNM")
		modelType = 3;
	else
		if (paras["-cm"] == (String) "CPM")
			modelType = 4;
		else
			if (paras["-cm"] == (String) "ERNM")
				modelType = 2;
			else
				if (paras["-cm"] == (String) "NNM")
					modelType = 1;
	int size = atoi((char*)paras["-size"]);
	bool label_flag = (paras["-ft"] == String("Cova")); 
	// label_flag is true when labels are bipartition and false when labels are trees.
	


	srand(time(NULL));
	int covariance_freeid_size = 0;
	int *covariance_freeid = new int[size];
	int covariance_nonfree_id_size = 0;
	int *covariance_nonfree_id = new int[size];

	Header_info info;
	File file_src(paras["-f"]);
	file_src >> info;
	String adj_src_file = info["source"];
	info.insert("created", time_stamp());
	info.insert("size", paras["-size"]);
	info.insert("source", paras["-f"]);

	String temp_file = paras["-path"]; 
	temp_file += "CDtemp_temp.txt";
	File file_Comm_temp(temp_file);
	file_Comm_temp.clean();


	double highfrequence = atof(highfreq.c_str());
	double lowfrequence = atof(lowfreq.c_str());
	if (label_flag && (highfrequence > 1.0 || highfrequence < 0.0
		|| lowfrequence > 1.0 || lowfrequence < 0.0 || (highfrequence - lowfrequence) <= 0.0))
	{
		cout << "Warning: The high and low frequencies must be between 0 and 1!\n\n";
		return false;
	}
	print_comm_array(mat, size, file_Comm_temp, label_flag, highfrequence, lowfrequence, covariance_freeid_size, covariance_nonfree_id_size, covariance_freeid, covariance_nonfree_id);
	info.insert("active_node_size", to_string(covariance_nonfree_id_size));

	char *infile = strdup((char*)temp_file);

	String outname_Graph = paras["-path"];
	outname_Graph += "CDtemp.bin";
	char *outfile = (char*)outname_Graph;

	String outname_Node = paras["-path"];
	outname_Node += "CDtemp_node_map.txt";
	char *node_map_file = (char*)outname_Node;
	String outname_Conf = paras["-path"];
	outname_Conf += "CDtemp.conf";
	char *conf_file = (char*)outname_Conf;
	int is_weighted = 1;
	int is_directed = 1;
	int is_single_slice = 0;
	double interslice_weight = 1.0;

	int* conf = NULL;
	int* sign = NULL;
	double* lambda = NULL;

	Slicer s(infile, (double)interslice_weight, is_directed, is_weighted, is_single_slice);
	Graph* g = s.get_graph();
	g->display_binary(outfile);
	delete g;
	s.display_node_mapping(node_map_file);
	s.display_conf(conf_file, modelType);

	int layers = read_conf(conf_file, conf, sign);
	lambda = new double[layers];


	//----------------fixed lambda neg, find two lambda_+ ----------
	//----------------such that the numbers of community  ----------
	//----------------are minimum or maximum.             ----------

	double lambda_neg = 0;
	double lambda_pos_min = -1, lambda_pos_max = 1;
	Community * community;
	map<double, Community *> LamCommunities;
	map<double, Community *>::iterator it, it2;
	map<double, double> mods;
	int times = 0;
	int numNodes = 0;

	cout << "Searching minimum of lambda_pos";

	while (times <= 20)
	{
		cout << ".";
		create_resolution(lambda_pos_min, lambda_neg, layers, sign, lambda);
		community = new Community(outfile, conf, sign, lambda);
		int stochastic = 0;
		GreedyLouvain::iterate_randomly = stochastic;
		GreedyLouvain::detect_communities(community);
		if (community->nb_comm == 1)
		{
			numNodes = community->g->nb_nodes;
			mods[lambda_pos_min] = community->modularity();
			LamCommunities[lambda_pos_min] = community;
			break;
		}
		else
		{
			delete community;
			community = NULL;
			lambda_pos_min *= 2;
		}
		times++;
	}
	if (times == 11)
	{
		cout << "Error: Cannot find lambda pos min!\n\n";
		return false;
	}
	cout << " Done!\n";

	Matrix<double> mat_adj;
	if (paras["-dm"] == (String) "URF") {
		File file_Adj(adj_src_file);
		if (!file_Adj.is_open()) {
			cout << "Warning! Unable to open adjacency matrix file: /" << (char*)adj_src_file << "/ for determining lambda_pos_max.\n";
			cout << "This may cause failure in CD of adjacency matrix from unweighted RF distance due to massive repeated distance.\n";
		}
		else {
			mat_adj.resize(size, size);
			file_Adj.end_header();
			for (int i = 0; i < size; i++)
			{
				for (int j = 0; j <= i; j++)
				{
					file_Adj >> mat_adj.matrix[i][j];
					mat_adj.matrix[j][i] = mat_adj.matrix[i][j];
				}
			}
			file_Adj.close();
		}
	}
	times = 0;

	cout << "Searching maximum of lambda_pos";
	while (times <= 20)
	{
		cout << ".";
		create_resolution(lambda_pos_max, lambda_neg, layers, sign, lambda);
		community = new Community(outfile, conf, sign, lambda);
		int stochastic = 0;
		GreedyLouvain::iterate_randomly = stochastic;
		GreedyLouvain::detect_communities(community);

		if (paras["-dm"] == (String) "URF" && !mat_adj.is_empty())//Affinity matrix from Unweighted RF distance is read.
		{
			
			bool allsametopo = true;
			double samevalue = mat_adj(0, 0);
			for (int i = 0; i < community->nb_comm; i++)
			{
				int repidx = -1;
				for (int j = 0; j < covariance_nonfree_id_size; j++)
				{
					if (community->n2c[j] == i)
					{
						if (repidx == -1)
						{
							repidx = j;
						}
						else
						{
							if (fabs(mat_adj(repidx, j) - samevalue) >= 1e-10)
							{
								allsametopo = false;
							}
						}
					}
				}
			}
			if (allsametopo)
			{
				mods[lambda_pos_max] = community->modularity();
				LamCommunities[lambda_pos_max] = community;
				break;
			}
		}
		if (community->nb_comm == covariance_nonfree_id_size)
		{
			mods[lambda_pos_max] = community->modularity();
			LamCommunities[lambda_pos_max] = community;
			break;
		}
		else
		{
			delete community;
			community = NULL;
			lambda_pos_max *= 2;
		}
		times++;
	}
	if (times == 11)
	{
		cout << "Error: Cannot find lambda pos max!\n\n";
		return false;
	}
	cout << " Done!\n";

	cout << "lambda neg: " << lambda_neg << endl;
	cout << "lambda pos min: " << lambda_pos_min << ", lambda pos max: " << lambda_pos_max << endl;

	queue<pair<double, double> > qq;
	pair<double, double> p, pp;
	vector<pair<double, double> > plateausLb(0);
	vector<pair<double, double> > plateausUb(0);
	bool SameAsFirst, SameAsSecond;
	double lambda_pos, plateaubound = 0;
	int atleastsearchnum = 0;
	p = make_pair(lambda_pos_min, lambda_pos_max);
	qq.push(p);

	// width first search
	cout << "\nThe testing values of lambda pos are:" << endl;
	cout << lambda_pos_min << ", " << lambda_pos_max << ", ";
	while (!qq.empty())
	{
		pp = qq.front();
		qq.pop();
		lambda_pos = (pp.first + pp.second) / 2;
		cout << lambda_pos << ", ";
		std::cout.flush();
		create_resolution(lambda_pos, lambda_neg, layers, sign, lambda);
		community = new Community(outfile, conf, sign, lambda);
		int stochastic = 0;
		GreedyLouvain::iterate_randomly = stochastic;
		GreedyLouvain::detect_communities(community);

		SameAsFirst = false;
		if (lambda_pos_min == pp.first)
		{
			if (community->nb_comm == 1)
			{
				SameAsFirst = true;
				lambda_pos_min = lambda_pos;
			}
		}
		else
		{
			it = LamCommunities.find(pp.first);
			if (it != LamCommunities.end())
			{
				if (Community::IsSameCommunity(community, it->second))
					SameAsFirst = true;
			}
		}
		SameAsSecond = false;
		if (lambda_pos_max == pp.second)
		{
			if (community->nb_comm == covariance_nonfree_id_size)
			{
				SameAsSecond = true;
				lambda_pos_max = lambda_pos;
			}
		}
		else
		{
			it = LamCommunities.find(pp.second);
			if (it != LamCommunities.end())
			{
				if (Community::IsSameCommunity(community, it->second))
					SameAsSecond = true;
			}
		}
		mods[lambda_pos] = community->modularity();
		LamCommunities[lambda_pos] = community;

		// added new interval for testing
		if (!SameAsFirst && ((lambda_pos - pp.first) > plateaubound || plateaubound == 0 || LamCommunities.size() + qq.size() < atleastsearchnum))
		{
			qq.push(make_pair(pp.first, lambda_pos));
		}
		if (!SameAsSecond && ((pp.second - lambda_pos) > plateaubound || plateaubound == 0 || LamCommunities.size() + qq.size() < atleastsearchnum))
		{
			qq.push(make_pair(lambda_pos, pp.second));
		}

		// find plateaus
		if (SameAsFirst && community->nb_comm != 1)
		{
			if (plateaubound == 0)
				plateaubound = (lambda_pos - pp.first) / 2;
			int i;
			for (i = 0; i < plateausLb.size(); i++)
			{
				if (plateausLb[i].second == pp.first)
				{
					plateausLb[i].second = lambda_pos;
					break;
				}
			}
			if (i == plateausLb.size())
			{
				plateausLb.resize(plateausLb.size() + 1);
				plateausLb[plateausLb.size() - 1] = make_pair(pp.first, lambda_pos);
			}
		}
		if (SameAsSecond && community->nb_comm != covariance_nonfree_id_size)
		{
			if (plateaubound == 0)
				plateaubound = (pp.second - lambda_pos) / 2;
			int i;
			for (i = 0; i < plateausLb.size(); i++)
			{
				if (plateausLb[i].first == pp.second)
				{
					plateausLb[i].first = lambda_pos;
					break;
				}
			}
			if (i == plateausLb.size())
			{
				plateausLb.resize(plateausLb.size() + 1);
				plateausLb[plateausLb.size() - 1] = make_pair(lambda_pos, pp.second);
			}
		}
	}
	cout << endl;

	cout << "\nThe found plateaus are:" << endl;
	for (int i = 0; i < plateausLb.size(); i++)
		cout << "[" << plateausLb[i].first << ", " << plateausLb[i].second
		<< "], length: " << plateausLb[i].second - plateausLb[i].first << ", number of communities: " <<
		LamCommunities[plateausLb[i].first]->nb_comm << endl;

	cout << "\nDetailed check lambda pos are (if necessary):" << endl;

	plateausUb.resize(plateausLb.size());
	for (int i = 0; i < plateausLb.size(); i++)
		plateausUb[i] = make_pair(plateausLb[i].first, plateausLb[i].second);

	double extstart, extend, max_length = 0;
	vector<int> totestidix(0);

	do {
		totestidix.resize(0);

		for (int i = 0; i < plateausUb.size(); i++)
		{
			extstart = plateausLb[i].first - (plateausLb[i].second - plateausLb[i].first);
			extend = plateausLb[i].second + (plateausLb[i].second - plateausLb[i].first);
			for (it = LamCommunities.begin(); it != LamCommunities.end(); it++)
			{
				if (it->first > extstart && it->first < plateausLb[i].first)
					extstart = it->first;
				if (it->first < extend && it->first > plateausLb[i].second)
					extend = it->first;
			}
			plateausUb[i].first = extstart;
			plateausUb[i].second = extend;

			if (plateausLb[i].second - plateausLb[i].first > max_length)
			{
				max_length = plateausLb[i].second - plateausLb[i].first;
			}
		}

		for (int i = 0; i < plateausLb.size(); i++)
		{
			if (((plateausUb[i].second - plateausUb[i].first) - max_length) >= 0.000001)
			{
				totestidix.resize(totestidix.size() + 1);
				totestidix[totestidix.size() - 1] = i;
			}
		}
		if (totestidix.size() <= 1)
			break;
		else
		{
			for (int i = 0; i < totestidix.size(); i++)
			{
				int idx = totestidix[i];

				lambda_pos = (plateausUb[idx].first + plateausLb[idx].first) / 2;
				if (LamCommunities.find(lambda_pos) == LamCommunities.end())
				{
					cout << lambda_pos << ", ";
					std::cout.flush();
					create_resolution(lambda_pos, lambda_neg, layers, sign, lambda);
					community = new Community(outfile, conf, sign, lambda);
					int stochastic = 0;
					GreedyLouvain::iterate_randomly = stochastic;
					GreedyLouvain::detect_communities(community);
					mods[lambda_pos] = community->modularity();
					LamCommunities[lambda_pos] = community;
				}
				it = LamCommunities.find(lambda_pos);
				it2 = LamCommunities.find(plateausLb[idx].first);
				if (Community::IsSameCommunity(it2->second, it->second))
					plateausLb[idx].first = lambda_pos;
				else
					plateausUb[idx].first = lambda_pos;
				lambda_pos = (plateausLb[idx].second + plateausUb[idx].second) / 2;
				if (LamCommunities.find(lambda_pos) == LamCommunities.end())
				{
					cout << lambda_pos << ", ";
					std::cout.flush();
					create_resolution(lambda_pos, lambda_neg, layers, sign, lambda);
					community = new Community(outfile, conf, sign, lambda);
					int stochastic = 0;
					GreedyLouvain::iterate_randomly = stochastic;
					GreedyLouvain::detect_communities(community);
					mods[lambda_pos] = community->modularity();
					LamCommunities[lambda_pos] = community;
				}
				it = LamCommunities.find(lambda_pos);
				it2 = LamCommunities.find(plateausLb[idx].second);
				if (Community::IsSameCommunity(it2->second, it->second))
					plateausLb[idx].second = lambda_pos;
				else
					plateausUb[idx].second = lambda_pos;
			}
			cout << endl;
		}
	} while (true);

	cout << "\nThe found plateaus are:" << endl;
	for (int i = 0; i < plateausLb.size(); i++)
	{
		cout << "Transition lower bound: [" << plateausUb[i].first << ", " << plateausLb[i].first << "]" << endl;
		cout << "Transition upper bound: [" << plateausLb[i].second << ", " << plateausUb[i].second << "]" << endl;
		cout << "Updated plateau: [" << plateausLb[i].first << ", " << plateausLb[i].second << "], length: " <<
			plateausLb[i].second - plateausLb[i].first << ", number of communities: " <<
			LamCommunities[plateausLb[i].first]->nb_comm << endl << endl;
	}

	
	

	int com_info_col = LamCommunities.size() + 1;
	double **com_info = new double *[covariance_nonfree_id_size + 4];
	for (int i = 0; i < covariance_nonfree_id_size + 4; i++)
		com_info[i] = new double[com_info_col];

	com_info[0][0] = numNodes;

	com_info[1][0] = lambda_neg;
	com_info[2][0] = 0;
	com_info[3][0] = 0;
	for (int i = 0; i < covariance_nonfree_id_size; i++)
		com_info[i + 4][0] = covariance_nonfree_id[i];

	// compute labels of communities
	int k = 0;
	com_info[0][1] = 0;
	it = LamCommunities.begin();
	community = it->second;

	com_info[1][1] = it->first;
	com_info[2][1] = it->second->nb_comm;
	com_info[3][1] = mods[it->first];//---it->second->modularity();
	for (int i = 4; i < covariance_nonfree_id_size + 4; i++)
		com_info[i][1] = it->second->n2c[i - 4];

	it++;
	for (; it != LamCommunities.end(); it++)
	{
		k++;
		if (Community::IsSameCommunity(community, it->second)) // problem
			com_info[0][k + 1] = com_info[0][k];
		else
			com_info[0][k + 1] = com_info[0][k] + 1;
		community = it->second;
		com_info[1][k + 1] = it->first;
		com_info[2][k + 1] = it->second->nb_comm;
		com_info[3][k + 1] = mods[it->first];
		for (int i = 4; i < covariance_nonfree_id_size + 4; i++)
			com_info[i][k + 1] = it->second->n2c[i - 4];
	}

	//plateaus
	String outname_pla("CD_Plateaus");
	outname_pla.make_stdname(paras);
	info.insert("output_type", "Plateaus of Community detection result");
	
	info.insert("node_type", paras["-node"]);

	if (modelType == 3)
		info.insert("CD_model","Configuration Null Model");
	else if (modelType == 4)
		info.insert("CD_model", "Constant Potts Model");
	else if (modelType == 2)
		info.insert("CD_model", "Erdos-Renyi Null Model");
	else if (modelType == 1)
		info.insert("CD_model", "No Null Model");
	info.insert("tuning", "automatically");

	File file_Pla(outname_pla);
	file_Pla.clean();
	file_Pla << info;

	file_Pla << "lengthLB:\t";
	for (int i = 0; i < plateausLb.size(); i++)
		file_Pla << plateausLb[i].second - plateausLb[i].first << "\t";
	file_Pla << "\n";
	file_Pla << "lengthUB:\t";
	for (int i = 0; i < plateausLb.size(); i++)
		file_Pla << plateausUb[i].second - plateausUb[i].first << "\t";
	file_Pla << "\n";
	file_Pla << "startUB:\t";
	for (int i = 0; i < plateausLb.size(); i++)
		file_Pla << plateausUb[i].first << "\t";
	file_Pla << "\n";
	file_Pla << "startLB:\t";
	for (int i = 0; i < plateausLb.size(); i++)
		file_Pla << plateausLb[i].first << "\t";
	file_Pla << "\n";
	file_Pla << "endLB:\t";
	for (int i = 0; i < plateausLb.size(); i++)
		file_Pla << plateausLb[i].second << "\t";
	file_Pla << "\n";
	file_Pla << "endUB:\t";
	for (int i = 0; i < plateausLb.size(); i++)
		file_Pla << plateausUb[i].second << "\t";
	file_Pla << "\n";
	if (label_flag)
		file_Pla << "Bipartition index" << "\t" << "Community Index" << "\n";
	else
		file_Pla << "Tree index" << "\t" << "Community Index" << "\n";
	for (int i = 0; i < covariance_nonfree_id_size; i++)
	{
		file_Pla << covariance_nonfree_id[i] << "\t";
		for (int j = 0; j < plateausLb.size(); j++)
		{
			file_Pla << LamCommunities[plateausLb[j].first]->n2c[i] << "\t";
		}
		file_Pla << "\n";
	}

	// output communities information
	String outname_CD("Community");
	outname_CD.make_stdname(paras);

	info.insert("output_type", "Community detection result");
	File file_CD(outname_CD);
	file_CD.clean();
	file_CD << info;

	for (int i = 0; i < covariance_nonfree_id_size + 4; i++)
	{
		if (i == 0)
			file_CD << "Same community as previous or not (first number is number of " << paras["-node"] << ")\n";
		else if (i == 1)
			file_CD << "Value of lambda: " << "\n";
		else if (i == 2)
			file_CD << "Number of communities: " << "\n";
		else if (i == 3)
			file_CD << "Value of modularity: " << "\n";
		else if (i == 4)
			file_CD << "Community index (first column is " << paras["-node"] << " index):\n";
		for (int j = 0; j < com_info_col; j++)
			file_CD << com_info[i][j] << "\t";
		file_CD << "\n";
	}

	//delete temporary file
	free(infile);

	free(conf);
	free(sign);
	free(covariance_freeid);
	free(covariance_nonfree_id);
	if (lambda != NULL)
		delete[] lambda;

	for (it = LamCommunities.begin(); it != LamCommunities.end(); it++)
		delete it->second;

	const char *tempfile0 = (char*)temp_file;
	remove(tempfile0);
	const char *tempfile1 = outfile;
	remove(tempfile1);
	const char *tempfile2 = node_map_file;
	remove(tempfile2);
	const char *tempfile3 = conf_file;
	remove(tempfile3);

	cout << "Output community results to file: " << outname_CD << endl;
	cout << "and " << outname_pla << "\n\n";
	return true;
}

bool community_detection_manually(Matrix<double> &mat, map<String, String> &paras)
{
	string highfreq = (char*)paras["-hf"];
	string lowfreq = (char*)paras["-lf"];
	int modelType = 0;
	if (paras["-cm"] == (String) "CNM")
		modelType = 3;
	else
		if (paras["-cm"] == (String) "CPM")
			modelType = 4;
		else
			if (paras["-cm"] == (String) "ERNM")
				modelType = 2;
			else
				if (paras["-cm"] == (String) "NNM")
					modelType = 1;
	int size = atoi((char*)paras["-size"]);
	bool label_flag = (paras["-node"] == String("Bipartition"));
	Array<double> lambda_pos_list;
	double lpiv = atof((char *)paras["-lpiv"]);
	double lp = atof((char *)paras["-lp"]);
	double lps = atof((char *)paras["-lps"]);
	double lpe = atof((char *)paras["-lpe"]);

	if (0 == lpiv)
	{
		lambda_pos_list.resize(1);
		lambda_pos_list[0] = lp;
	}
	else
	{
		int size = (int)((lpe - lps) / lpiv + 1);
		lambda_pos_list.resize(size);
		for (int i = 0; i < size; i++)
			lambda_pos_list[i] = lps + i * lpiv;
	}

	Array<double> lambda_neg_list;
	double lniv = atof((char *)paras["-lniv"]);
	double ln = atof((char *)paras["-ln"]);
	double lns = atof((char *)paras["-lns"]);
	double lne = atof((char *)paras["-lne"]);

	if (0 == lniv)
	{
		lambda_neg_list.resize(1);
		lambda_neg_list[0] = ln;
	}
	else
	{
		int size = (int)((lne - lns) / lniv + 1);
		lambda_neg_list.resize(size);
		for (int i = 0; i < size; i++)
			lambda_neg_list[i] = lns + i * lniv;
	}

	srand(time(NULL));
	int covariance_freeid_size = 0;
	int *covariance_freeid = new int[size];
	int covariance_nonfree_id_size = 0;
	int *covariance_nonfree_id = new int[size];


	Header_info info;
	File file_src(paras["-f"]);
	file_src >> info;
	info.insert("created", time_stamp());
	info.insert("size", paras["-size"]);
	info.insert("source", paras["-f"]);

	String temp_file = paras["-path"];
	temp_file += "CDtemp_temp.txt";
	File file_Comm_temp(temp_file);
	file_Comm_temp.clean();

	double highfrequence = atof(highfreq.c_str());
	double lowfrequence = atof(lowfreq.c_str());
	if (label_flag && (highfrequence > 1.0 || highfrequence < 0.0
		|| lowfrequence > 1.0 || lowfrequence < 0.0 || (highfrequence - lowfrequence) <= 0.0))
	{
		cout << "Warning: The high and low frequencies must be between 0 and 1!\n\n";
		return false;
	}

	print_comm_array(mat, size, file_Comm_temp, label_flag, highfrequence, lowfrequence, covariance_freeid_size, covariance_nonfree_id_size, covariance_freeid, covariance_nonfree_id);
	info.insert("active_node_size", to_string(covariance_nonfree_id_size));

	char *infile = strdup((char*)temp_file);

	String outname_Graph = paras["-path"];
	outname_Graph += "CDtemp.bin";
	char *outfile = (char*)outname_Graph;

	String outname_Node = paras["-path"];
	outname_Node += "CDtemp_node_map.txt";
	char *node_map_file = (char*)outname_Node;
	String outname_Conf = paras["-path"];
	outname_Conf += "CDtemp.conf";
	char *conf_file = (char*)outname_Conf;

	int is_weighted = 1;
	int is_directed = 1;
	int is_single_slice = 0;
	double interslice_weight = 1.0;

	int* conf = NULL;
	int* sign = NULL;
	double* lambda = NULL;

	Slicer s(infile, (double)interslice_weight, is_directed, is_weighted, is_single_slice);
	Graph* g = s.get_graph();
	g->display_binary(outfile);
	delete g;
	s.display_node_mapping(node_map_file);
	s.display_conf(conf_file, modelType);

	int layers = read_conf(conf_file, conf, sign);
	lambda = new double[layers];
	double lambda_pos;// = atof(lambda_pos_list.c_str());
	double lambda_neg;// = atof(lambda_neg_list.c_str());

	info.insert("node_type", paras["-node"]);

	if (modelType == 3)
		info.insert("CD_model", "Configuration Null Model");
	else if (modelType == 4)
		info.insert("CD_model", "Constant Potts Model");
	else if (modelType == 2)
		info.insert("CD_model", "Erdos-Renyi Null Model");
	else if (modelType == 1)
		info.insert("CD_model", "No Null Model");
	info.insert("tuning", "manually");
	String lnlist = to_string(lambda_neg_list[0]);
	for (int i=1;i<lambda_neg_list.get_length();i++){
		lnlist += ",";
		lnlist += to_string(lambda_neg_list[i]);
	}
	info.insert("lambda_neg", lnlist);
	String lplist = to_string(lambda_pos_list[0]);
	for (int i = 1; i<lambda_pos_list.get_length(); i++) {
		lplist += ",";
		lplist += to_string(lambda_pos_list[i]);
	}
	info.insert("lambda_pos", lplist);

	String outname_CD("Community");
	outname_CD.make_stdname(paras);
	info.insert("output_type", "Community detection result");
	File file_CD(outname_CD);
	file_CD.clean();
	file_CD << info;
	if (!file_CD.is_open())
	{
		cout << "Unable to open the file: /" << outname_CD << "/\n\n";
		return false;
	}
	int lambdasize = 0;

	if (lambda_pos_list.get_length() == 1)
		lambdasize = lambda_neg_list.get_length();
	else if (lambda_neg_list.get_length() == 1)
		lambdasize = lambda_pos_list.get_length();
	else
	{
		cout << "Error: The length of the array of lambda negative or the length of the array lambda positive should be one!\n\n";
		return false;
	};
	int com_info_col = lambdasize + 1;

	double **com_info = new double *[covariance_nonfree_id_size + 4];
	for (int i = 0; i < covariance_nonfree_id_size + 4; i++)
		com_info[i] = new double[com_info_col];

	Community** communities = new Community *[lambdasize];

	for (int k = 0; k < lambdasize; k++)
	{
		if (lambda_pos_list.get_length() == 1)
		{
			lambda_pos = lambda_pos_list[0];
			lambda_neg = lambda_neg_list[k];
		}
		else
		{
			lambda_pos = lambda_pos_list[k];
			lambda_neg = lambda_neg_list[0];
		}
		if (modelType != 1)
		{
			cout << "Using lambda+ " << lambda_pos << ", lambda- " << lambda_neg << endl;
		}
		create_resolution(lambda_pos, lambda_neg, layers, sign, lambda);
		communities[k] = new Community(outfile, conf, sign, lambda);
		int stochastic = 0;
		GreedyLouvain::iterate_randomly = stochastic;
		if (stochastic)
			cout << "Using random node order" << endl;
		cout << "Network has "
			<< communities[k]->g->nb_nodes << " nodes, "
			<< communities[k]->g->nb_links << " links, " << endl;

		GreedyLouvain::detect_communities(communities[k]);
		double mod = communities[k]->modularity();

		cout << "Value of modularity:" << mod << endl;
		cout << "Number of communities: " << communities[k]->nb_comm << endl;
		//co->display_comm2node();

		if (covariance_freeid != NULL && covariance_freeid_size != 0)
		{
			for (int i = 0; i < communities[k]->nb_comm; i++)
			{
				std::cout << "Community " << i + 1 << " includes nodes: ";
				for (int j = 0; j < communities[k]->size; j++)
				{
					if (communities[k]->n2c[j] == i)
						std::cout << covariance_nonfree_id[j] << ",";
				}
				std::cout << "\n";
			}
			std::cout << "Free node index:" << endl;
			for (int i = 0; i < covariance_freeid_size; i++)
				std::cout << covariance_freeid[i] << ", ";
			std::cout << "\n";
		}
		else
		{
			for (int i = 0; i < communities[k]->nb_comm; i++)
			{
				std::cout << "Community " << i + 1 << " includes nodes: ";
				for (int j = 0; j < communities[k]->size; j++)
				{
					if (communities[k]->n2c[j] == i)
						std::cout << j << ",";
				}
				std::cout << "\n";
			}
		}

		if (k == 0)
		{
			com_info[0][0] = communities[k]->g->nb_nodes;

			if (lambda_pos_list.get_length() == 1)
				com_info[1][0] = lambda_pos_list[0];
			else
				com_info[1][0] = lambda_neg_list[0];
			com_info[2][0] = 0;
			com_info[3][0] = 0;
			for (int i = 0; i < covariance_nonfree_id_size; i++)
				com_info[i + 4][0] = covariance_nonfree_id[i];
		}

		if (k == 0)
			com_info[0][1] = 0;
		else
		{
			if (Community::IsSameCommunity(communities[k - 1], communities[k]))
				com_info[0][k + 1] = com_info[0][k];
			else
				com_info[0][k + 1] = com_info[0][k] + 1;
		}
		if (lambda_pos_list.get_length() == 1)
			com_info[1][k + 1] = lambda_neg_list[k];
		else
			com_info[1][k + 1] = lambda_pos_list[k];

		com_info[2][k + 1] = communities[k]->nb_comm;
		com_info[3][k + 1] = mod;

		for (int i = 4; i < covariance_nonfree_id_size + 4; i++)
			com_info[i][k + 1] = communities[k]->n2c[i - 4];
	}

	for (int i = 0; i < covariance_nonfree_id_size + 4; i++)
	{
		if (i == 0)
		{
			file_CD << "Same community as previous or not (first number is number of "<< paras["-node"] <<")\n";
		}
		else if (i == 1)
			file_CD << "Value of lambda: " << "\n";
		else if (i == 2)
			file_CD << "Number of communities: " << "\n";
		else if (i == 3)
			file_CD << "Value of modularity: " << "\n";
		else if (i == 4)
		{
			file_CD << "Community index (first column is " << paras["-node"] << " index): \n";
		}
		for (int j = 0; j < lambdasize + 1; j++)
			file_CD << com_info[i][j] << "\t";
		file_CD << "\n";
	}

	//delete temporary file
	free(infile);

	free(conf);
	free(sign);
	free(covariance_freeid);
	free(covariance_nonfree_id);
	if (lambda != NULL)
		delete[] lambda;


	const char *tempfile0 = (char*)temp_file;
	remove(tempfile0);
	const char *tempfile1 = outfile;
	remove(tempfile1);
	const char *tempfile2 = node_map_file;
	remove(tempfile2);
	const char *tempfile3 = conf_file;
	remove(tempfile3);

	cout << "Output community results to file: " << outname_CD << endl;
	return true;
}