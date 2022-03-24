#include "zdcommunity.hpp"

void check_self_covariance(SpecMat::LowerTri<PRECISION>& adjacency, double hf, double lf, int& covariance_freeid_size, int* covariance_freeid, int& covariance_nonfree_id_size, int* covariance_nonfree_id)
{
	int n = adjacency.dimension();
	adjacency.form_row_ptr();
	assert(n != 0);

	for (int i = 0; i < n; i++)
	{
		bool crit;
		crit = (adjacency[i][i] <= lf * (1 - lf) || adjacency[i][i] <= hf * (1 - hf));

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
}

Graph *SymmetricOneSliceAdjaceny2Graph(SpecMat::LowerTri<PRECISION> &adjacency, int *node_id, int node_size, double weight_threshold)
{
	adjacency.form_row_ptr();

	int pos_node_size_per_node;
	int *pos_node_id_per_node;
	double *pos_node_weight_per_node;
	pos_node_id_per_node = new int[node_size];
	pos_node_weight_per_node = new double[node_size];

	int neg_node_size_per_node;
	int *neg_node_id_per_node;
	double *neg_node_weight_per_node;
	neg_node_id_per_node = new int[node_size];
	neg_node_weight_per_node = new double[node_size];

	double *r_i;
	int total_degree = 0; // nb_links
	int acc_degree = 0;

	Array<int> degrees(0, 10 * node_size);
	Array<int> links(0, 10 * node_size);
	Array<double> weights(0, 10 * node_size);
	double w_ij = 0;

	for(int i = 0; i < node_size; i++)
	{
		pos_node_size_per_node = 0;
		neg_node_size_per_node = 0;
		// r_i = adjacency[node_id[i]];
		for(int j = 0; j < node_size; j++)
		{
			w_ij = node_id[j] <= node_id[i] ? adjacency[node_id[i]][node_id[j]] : adjacency[node_id[j]][node_id[i]];
			if (w_ij < -weight_threshold) // Negative weights
			{
				neg_node_id_per_node[neg_node_size_per_node] = node_id[j];
				neg_node_weight_per_node[neg_node_size_per_node++] = -w_ij;
				total_degree++;	
			}
			else if (w_ij > weight_threshold) // Positive weights
			{
				pos_node_id_per_node[pos_node_size_per_node] = node_id[j];
				pos_node_weight_per_node[pos_node_size_per_node++] = w_ij;
				total_degree++;	
			}
		}
		acc_degree += pos_node_size_per_node;
		degrees.push(acc_degree); // pos_in
		acc_degree += pos_node_size_per_node;
		degrees.push(acc_degree); // pos_out;
		acc_degree += neg_node_size_per_node;
		degrees.push(acc_degree); // neg_in
		acc_degree += neg_node_size_per_node;
		degrees.push(acc_degree); // neg_out
		degrees.push(acc_degree); // inter_in;
		degrees.push(acc_degree); // inter_out;

		for (int k = 0; k < 2; k++) // one for pos_in, one for pos_out
		{
			for(int j = 0; j < pos_node_size_per_node; j++)
			{
				links.push(pos_node_id_per_node[j]);
				weights.push(pos_node_weight_per_node[j]);
			}
		}

		for (int k = 0; k < 2; k++) // one for neg_in, one for neg_out
		{
			for(int j = 0; j < neg_node_size_per_node; j++)
			{
				links.push(neg_node_id_per_node[j]);
				weights.push(neg_node_weight_per_node[j]);
			}
		}
	}


	degrees.align();
	links.align();
	weights.align();

	auto degrees_ptr = degrees.pop_vec();
	auto links_ptr = links.pop_vec();
	auto weights_ptr = weights.pop_vec();

	Graph* g = new Graph(node_size, 3, total_degree, degrees_ptr, links_ptr, weights_ptr);
	return g;
}

// int read_conf(char* filename, int* &conf, int* &sign)
// {
// 	//ifstream finput;
// 	//finput.open(filename, fstream::in | fstream::binary);
// 	//if (finput.fail())
// 	//{
// 	//	cout << "Could not read config file: " << filename << "\n\n";
// 	//	exit(-1);
// 	//}

// 	std::ifstream file_Conf;
// 	file_Conf.open(filename);
// 	if (!file_Conf.is_open())
// 	{
// 		cout << "Could not read config file: " << filename << "\n\n";
// 		exit(-1);
// 	}
// 	int pos_header = end_header(file_Conf);
// 	file_Conf.close();
// 	ifstream finput;
// 	finput.open(filename, fstream::in | fstream::binary);
// 	finput.seekg(pos_header, ios::beg);

// 	int nb_layers;
// 	//read number of layers
// 	finput.read((char *)&nb_layers, sizeof(int));
// 	conf = (int*)malloc((long)nb_layers * sizeof(int));
// 	sign = (int*)malloc((long)nb_layers * sizeof(int));

// 	// read conf
// 	finput.read((char *)conf, nb_layers * sizeof(int));
// 	if (finput.fail() || finput.eof())
// 	{
// 		cout << "Error!\n\n";
// 		exit(-1);
// 		//throw new exception();
// 	}
// 	finput.read((char *)sign, nb_layers * sizeof(int));

// 	return nb_layers;
// }

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

bool community_detection_automatically(SpecMat::LowerTri<PRECISION> &adj, map<String, String> &paras) {
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
	bool label_flag = (paras["-ft"] == (String) "Cova"); 
	// label_flag is true when labels are bipartition and false when labels are trees.
	


	srand(time(NULL));
	int covariance_freeid_size = 0;
	int *covariance_freeid = NULL;
	int covariance_nonfree_id_size = 0;
	int *covariance_nonfree_id = NULL;
	covariance_freeid = new int[adj.dimension()];
	covariance_nonfree_id = new int[adj.dimension()];

	String temp = paras["-post"];
	paras["-post"] = "";
	String info_item[4] = { "created","output_type","size","source" };
	String info_content[4] = { time_stamp(),"Comminty detection temp file", paras["-size"], paras["-f"] };
	Header_info info(info_item, info_content, 3);

	// String temp_file = make_stdname("Community_temp", paras);
	// paras["-post"] = temp;
	// File file_Comm_temp(temp_file);
	// file_Comm_temp.clean();


	double highfrequence = atof(highfreq.c_str());
	double lowfrequence = atof(lowfreq.c_str());
	if (label_flag && (highfrequence > 1.0 || highfrequence < 0.0
		|| lowfrequence > 1.0 || lowfrequence < 0.0 || (highfrequence - lowfrequence) <= 0.0))
	{
		cout << "Warning: The high and low frequencies must be between 0 and 1!\n\n";
		return false;
	}
	// print_comm_array(mat, size, file_Comm_temp, label_flag, highfrequence, lowfrequence);
	if (label_flag)
		check_self_covariance(adj, highfrequence, lowfrequence, covariance_freeid_size, covariance_freeid, covariance_nonfree_id_size, covariance_nonfree_id);
	else
	{
		covariance_nonfree_id_size = adj.dimension();
		for(int i = 0; i < adj.dimension(); i++)
			covariance_nonfree_id[i] = i;
	}

	Graph *g = SymmetricOneSliceAdjaceny2Graph(adj, covariance_nonfree_id, covariance_nonfree_id_size, 0.0);
	

	int is_weighted = 1;
	int is_directed = 1;
	int is_single_slice = 0;
	double interslice_weight = 1.0;

	int* conf = new int[3];
	conf[0] = modelType;
	conf[1] = modelType;
	conf[2] = 1;
	int* sign = new int[3];
	sign[0] = 1;
	sign[1] = -1;
	sign[2] = 1;

	double* lambda = NULL;


	int layers = 3;
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

	Graph *gc = nullptr;

	while (times <= 20)
	{
		create_resolution(lambda_pos_min, lambda_neg, layers, sign, lambda);
		if (gc != nullptr)
			delete gc;
		gc = new Graph(*g);
		community = new Community(gc, conf, sign, lambda);
		// community = new Community(outfile, conf, sign, lambda);
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

	times = 0;
	while (times <= 20)
	{
		create_resolution(lambda_pos_max, lambda_neg, layers, sign, lambda);
		if (gc != nullptr)
			delete gc;
		gc = new Graph(*g);
		community = new Community(gc, conf, sign, lambda);
		// community = new Community(outfile, conf, sign, lambda);
		int stochastic = 0;
		GreedyLouvain::iterate_randomly = stochastic;
		GreedyLouvain::detect_communities(community);

		if (paras["-dm"] == (String) "URF")//Affinity matrix from Unweighted RF distance is read.
		{
			bool allsametopo = true;
			// double samevalue = mat(0, 0);
			double samevalue = adj[0][0];
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
							// if (fabs(mat(repidx, j) - samevalue) >= 1e-10)
							double check = (repidx >= j) ? adj[repidx][j] - samevalue : adj[j][repidx] - samevalue;
							if (fabs(check) >= 1e-10)
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
		if (gc != nullptr)
			delete gc;
		gc = new Graph(*g);
		community = new Community(gc, conf, sign, lambda);
		// community = new Community(outfile, conf, sign, lambda);
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
					if (gc != nullptr)
						delete gc;
					gc = new Graph(*g);
					community = new Community(gc, conf, sign, lambda);
					// community = new Community(outfile, conf, sign, lambda);
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
					if (gc != nullptr)
						delete gc;
					gc = new Graph(*g);
					community = new Community(gc, conf, sign, lambda);
					// community = new Community(outfile, conf, sign, lambda);
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
	String outname_Pla("CD_Plateaus");
	outname_Pla.make_stdname(paras);
	info.insert("output_type", "Plateaus of Community detection result");
	
	if (label_flag)
		info.insert("label_type", "Bipartition");
	else
		info.insert("label_type", "Tree");
	info.insert("label_feature", "todo");

	if (modelType == 3)
		info.insert("CD_model","Configuration Null Model");
	else if (modelType == 4)
		info.insert("CD_model", "Constant Potts Model");
	else if (modelType == 2)
		info.insert("CD_model", "Erdos-Renyi Null Model");
	else if (modelType == 1)
		info.insert("CD_model", "No Null Model");
	info.insert("tuning", "automatica");

	std::ofstream file_Pla;
	file_Pla.open((char *)outname_Pla);
	if (!file_Pla.is_open())
		throw(1);
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
	file_Pla.close();

	// output communities information
	String outname_CD("Community");
	outname_CD.make_stdname(paras);

	info.insert("output_type", "Community detection result");
	std::ofstream file_CD;
	file_CD.open((char *)outname_CD);
	file_CD << info;

	for (int i = 0; i < covariance_nonfree_id_size + 4; i++)
	{
		if (i == 0)
		{
			if (label_flag)
				file_CD << "Same community as previous or not (first number is number of bipartitions)" << "\n";
			else
				file_CD << "Same community as previous or not (first number is number of trees)" << "\n";
		}
		else if (i == 1)
			file_CD << "Value of lambda: " << "\n";
		else if (i == 2)
			file_CD << "Number of communities: " << "\n";
		else if (i == 3)
			file_CD << "Value of modularity: " << "\n";
		else if (i == 4)
		{
			if (label_flag)
				file_CD << "Community index (first column is bipartition index): " << "\n";
			else
				file_CD << "Community index (first column is tree index): " << "\n";
		}
		for (int j = 0; j < com_info_col; j++)
			file_CD << com_info[i][j] << "\t";
		file_CD << "\n";
	}
	file_CD.close();

	//delete temporary file
	// free(infile);
	// free(outfile);
	// free(node_map_file);
	// free(conf_file);

	free(conf);
	free(sign);
	if (lambda != NULL)
		delete[] lambda;

	for (it = LamCommunities.begin(); it != LamCommunities.end(); it++)
		delete it->second;

	cout << "Output community results to file: " << outname_CD << endl;
	cout << "and " << outname_Pla << "\n\n";
	return true;
}
