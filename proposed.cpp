#include "class.h"

//データ読み込み
void set_data(std::vector<Data> &z);
void set_data_place(std::vector<Data> &z);

void insert(int node, Data z);

void get_leaf_node_id(Node node);
value_type dist_nodes(Node node1, Node node2);
int share_parent_node(Node node1, Node node2);
value_type dist_data_to_node(Data data, Node node2);
void get_all_node(int node);
int find_node(int node, Data z);

vector<Data> objects;
vector<Data> objects2;
map<int, Node> set_of_nodes;
vector<int> child_node_id;
vector<int> child_node_id2;

int T = 0;
int node_id = 2;

value_type calc_spatial_score(Data z1, Data z2, double dist_max);
value_type calc_textual_score(Data z1, Data z2);

int main(int argc, char *argv[])
{
	value_type dist_max;																										//位置スコア計算用
	std::chrono::system_clock::time_point start, end, mesure_time_start, mesure_time_end, mesure_time_start1, mesure_time_end1; //実行時間計測用　変数

	//データの初期化
	for (int k = 0; k < num; k++)
	{
		objects.push_back({
			k,
			0,
			0,
		});

		objects2.push_back({
			k,
			0,
			0,
		});
	}

	//データ読み込み
	set_data(objects2);	//Twitter
	//set_data_place(objects2);	//Place

	Range whole_region;

	whole_region.x_max = 50.0;
	whole_region.x_min = 25.0;
	whole_region.y_max = -60.0;
	whole_region.y_min = -125.0;

	dist_max = sqrt(pow(abs(whole_region.x_max - whole_region.x_min), 2.0) + pow(abs(whole_region.y_max - whole_region.y_min), 2.0));

	Range region = {0, 0, 0, 0};

	Rtree rtree;
	Node root_node;
	Node init_node;

	init_node.MBR = whole_region;
	init_node.isleaf = true;
	init_node.id = 1;
	init_node.depth = 1;
	set_of_nodes[init_node.id] = init_node;

	root_node.MBR = whole_region;
	root_node.isleaf = false;
	root_node.child_node[1] = init_node.id;
	root_node.id = 0;
	root_node.parent_node = -1;
	root_node.depth = 0;
	set_of_nodes[root_node.id] = root_node;
	set_of_nodes[init_node.id].parent_node = root_node.id;

	rtree.root_node = &set_of_nodes[root_node.id];

	//パラメータ出力
	cout << "proposed" << endl;
	cout << "TOP_K = " << TOP_K << endl;
	cout << "num = " << num << endl;
	cout << "NODE_MAX = " << NODE_MAX << endl;
	cout << "alpha = " << alpha << endl;
	cout << "l = " << L << endl;
	cout << "preprocessing" << endl;

	//vector<string> distinct_key;

	mt19937 rnd(1);
	multiset<AnswerPair> answer;
	value_type tau;
	map<string, int> key_frecency;
	multiset<Key_tf> key_tf;
	int current_node_id;

	int rand_k = 1;
	int rand1, rand2;
	value_type spacial_score, textual_score, score;
	AnswerPair pair1;
	int num_of_calclation = 0;
	int num_of_calclation1 = 0;
	int num_of_calclation2 = 0;
	int num_of_calclation3 = 0;
	int num_of_calclation4 = 0;

	int num_share_key = 0;
	int max_num_share_key = 0;
	float jaccard_max = 0;
	int jaccard_cnt = 0;

	value_type time;

	//Rtree作成時間の計測開始
	start = std::chrono::system_clock::now();
	mesure_time_start = std::chrono::system_clock::now();

	for (int i = 0; i < num; i++)
	{
		rand1 = rnd() % num;

		objects[i] = objects2[rand1];
	}

	//Rtree作成
	for (int i = 0; i < num; i++)
	{
		insert(root_node.id, objects[i]);
	}

	//Rtree作成時間の出力
	mesure_time_end = std::chrono::system_clock::now();
	time = std::chrono::duration_cast<std::chrono::microseconds>(mesure_time_end - mesure_time_start).count();
	cout << "insert time(rtree)= " << time / 1000000 << "sec" << endl;

	//「ノード情報の保存」の時間計測開始
	mesure_time_start = std::chrono::system_clock::now();

	//葉ノード　データのキーワード集合の和集合の保存
	for (auto itr = set_of_nodes.begin(); itr != set_of_nodes.end(); itr++)
	{
		if (itr->second.isleaf == true)
		{
			for (int i = 0; i < itr->second.data_id.size(); i++)
			{
				for (int k = 0; k < objects[itr->second.data_id[i] - 1].key.size(); k++)
				{
					set_of_nodes[itr->first].key.push_back(objects[itr->second.data_id[i] - 1].key[k]);
				}
			}
			std::sort(set_of_nodes[itr->first].key.begin(), set_of_nodes[itr->first].key.end());
			set_of_nodes[itr->first].key.erase(std::unique(set_of_nodes[itr->first].key.begin(), set_of_nodes[itr->first].key.end()), set_of_nodes[itr->first].key.end());
		}
	}

	//葉ノード　取りうるJaccard類似度の最大値の保存
	for (auto itr = set_of_nodes.begin(); itr != set_of_nodes.end(); itr++)
	{
		jaccard_max = 0;
		if (itr->second.isleaf == true)
		{
			for (int i = 0; i < itr->second.data_id.size(); i++)
			{
				itr->second.share_key.push_back(0);
			}
			for (int i = 0; i < itr->second.data_id.size(); i++)
			{
				for (int j = 0; j < itr->second.data_id.size(); j++)
				{
					if (i != j)
					{
						num_share_key = 0;
						for (int k = 0; k < objects[itr->second.data_id[i] - 1].key.size(); k++)
						{
							for (int l = 0; l < objects[itr->second.data_id[j] - 1].key.size(); l++)
							{
								if (objects[itr->second.data_id[i] - 1].key[k] == objects[itr->second.data_id[j] - 1].key[l])
								{
									num_share_key++;
								}
							}
						}
						set_of_nodes[itr->first].share_key[i] = num_share_key;
						set_of_nodes[itr->first].share_key[j] = num_share_key;
						float jaccard = 2 * num_share_key;
						jaccard /= objects[itr->second.data_id[i] - 1].key.size() + objects[itr->second.data_id[j] - 1].key.size() - num_share_key; //Jaccard
						//jaccard /= sqrt(objects[itr->second.data_id[i] - 1].key.size() * objects[itr->second.data_id[j] - 1].key.size());	//Cosine
						//jaccard /= (objects[itr->second.data_id[i] - 1].key.size() + objects[itr->second.data_id[j] - 1].key.size());	//Dice
						if (jaccard > jaccard_max)
						{
							jaccard_max = jaccard;
						}
					}
				}
			}

			max_num_share_key = 0;
			for (int i = 0; i < set_of_nodes[itr->first].share_key.size(); i++)
			{
				if (set_of_nodes[itr->first].share_key[i] > max_num_share_key)
				{
					max_num_share_key = set_of_nodes[itr->first].share_key[i];
				}
			}
			set_of_nodes[itr->first].max_share_key = max_num_share_key;
			set_of_nodes[itr->first].jaccard_max = jaccard_max;

			max_num_share_key = 0;
			for (int i = 0; i < itr->second.data_id.size(); i++)
			{
				if (objects[itr->second.data_id[i] - 1].key.size() > max_num_share_key)
				{
					max_num_share_key = objects[itr->second.data_id[i] - 1].key.size();
				}
			}
			set_of_nodes[itr->first].max_num_key = max_num_share_key;
		}
	}

	//中間ノード　取りうるJaccard類似度の最大値の保存
	for (auto itr = set_of_nodes.begin(); itr != set_of_nodes.end(); itr++)
	{
		if (itr->second.isleaf == false)
		{
			child_node_id.clear();
			max_num_share_key = 0;
			jaccard_max = 0;
			get_leaf_node_id(itr->second);
			for (int i = 0; i < child_node_id.size(); i++)
			{
				if (set_of_nodes[child_node_id[i]].max_share_key > max_num_share_key)
				{
					max_num_share_key = set_of_nodes[child_node_id[i]].max_share_key;
				}
				if (set_of_nodes[child_node_id[i]].jaccard_max > jaccard_max)
				{
					jaccard_max = set_of_nodes[child_node_id[i]].jaccard_max;
				}
			}
			set_of_nodes[itr->first].max_share_key = max_num_share_key;
			set_of_nodes[itr->first].jaccard_max = jaccard_max;

			if (jaccard_max != 1)
			{
				double text_max = 0;
				for (int i = 0; i < child_node_id.size(); i++)
				{
					for (int j = 0; j < set_of_nodes[child_node_id[i]].data_id.size(); j++)
					{
						int pivot = floor((1 - jaccard_max) * objects[set_of_nodes[child_node_id[i]].data_id[j] - 1].key.size()) + 1;
						for (int k = 0; k < pivot; k++)
						{
							for (int l = 0; l < set_of_nodes[child_node_id[i]].data_id.size(); l++)
							{
								for (int p = 0; p < objects[set_of_nodes[child_node_id[i]].data_id[l] - 1].key.size(); p++)
								{
									if (objects[set_of_nodes[child_node_id[i]].data_id[j] - 1].key[k] == objects[set_of_nodes[child_node_id[i]].data_id[l] - 1].key[p])
									{
										textual_score = calc_textual_score(objects[set_of_nodes[child_node_id[i]].data_id[j] - 1], objects[set_of_nodes[child_node_id[i]].data_id[l] - 1]);
										if (textual_score > text_max)
										{
											text_max = textual_score;
										}
									}
								}
							}
						}
					}
				}
				set_of_nodes[itr->first].jaccard_max = text_max;
			}

			max_num_share_key = 0;
			for (int i = 0; i < child_node_id.size(); i++)
			{
				if (set_of_nodes[child_node_id[i]].max_num_key > max_num_share_key)
				{
					max_num_share_key = set_of_nodes[child_node_id[i]].max_num_key;
				}
			}
			set_of_nodes[itr->first].max_num_key = max_num_share_key;
		}
	}

	//「ノード情報の保存」の時間出力
	mesure_time_end = std::chrono::system_clock::now();
	time = std::chrono::duration_cast<std::chrono::microseconds>(mesure_time_end - mesure_time_start).count();
	cout << "preprocessing time = " << time / 1000000 << "sec" << endl;

	//閾値の初期値取得　計測開始
	cout << "calc init knn" << endl;
	mesure_time_start = std::chrono::system_clock::now();

	//ランダムなkペアから閾値の初期値取得
	/*while (rand_k <= TOP_K)
	{

		rand1 = rnd() % num;
		rand2 = rnd() % num;

		spacial_score = calc_spatial_score(objects[rand1], objects[rand2], dist_max);
		textual_score = calc_textual_score(objects[rand1], objects[rand2]);
		score = alpha * spacial_score + (1 - alpha) * textual_score;

		pair1 = {objects[rand1].id, objects[rand2].id, score, spacial_score, textual_score};
		answer.insert(pair1);
		rand_k++;
	}*/

	auto itr_ans = answer.begin();
	tau = itr_ans->score;

	//init knn,threshold
	int number_of_key = 0;
	set<string> type_of_key;
	value_type max_score = 0;
	int max_score_node;
	vector<NodeScore> node_order;
	multiset<NodeScore> node_order2;
	NodeScore node_score;
	for (auto itr = set_of_nodes.begin(); itr != set_of_nodes.end(); itr++)
	{
		if (itr->second.isleaf == true)
		{
			number_of_key = 0;
			type_of_key.clear();
			spacial_score = 1 - (sqrt(pow(abs(itr->second.MBR.x_max - itr->second.MBR.x_min), 2.0) + pow(abs(itr->second.MBR.y_max - itr->second.MBR.y_min), 2.0)) / dist_max);
			for (int i = 0; i < itr->second.data_id.size(); i++)
			{
				for (auto itr2 = objects[itr->second.data_id[i] - 1].key.begin(); itr2 != objects[itr->second.data_id[i] - 1].key.end(); itr2++)
				{
					number_of_key++;
					type_of_key.insert(*itr2);
				}
			}
			textual_score = itr->second.jaccard_max;
			score = alpha * spacial_score + (1 - alpha) * textual_score;
			if (score > max_score)
			{
				max_score = score;
				max_score_node = itr->second.id;
			}
			node_score.node_id = itr->second.id;
			node_score.score = score;
			node_score.flag = false;
			node_order2.insert(node_score);
		}
	}
	auto itr_order1 = node_order2.end();
	for (int j = 0; j < L; j++)
	{
		itr_order1--;
		for (int i = 0; i < set_of_nodes[itr_order1->node_id].data_id.size(); i++)
		{
			for (int j = 0; j < set_of_nodes[itr_order1->node_id].data_id.size(); j++)
			{
				if (i != j)
				{
					spacial_score = calc_spatial_score(objects[set_of_nodes[itr_order1->node_id].data_id[i] - 1], objects[set_of_nodes[itr_order1->node_id].data_id[j] - 1], dist_max);
					textual_score = calc_textual_score(objects[set_of_nodes[itr_order1->node_id].data_id[i] - 1], objects[set_of_nodes[itr_order1->node_id].data_id[j] - 1]);
					score = alpha * spacial_score + (1 - alpha) * textual_score;

					if (answer.size() < TOP_K)
					{
						pair1 = {objects[set_of_nodes[itr_order1->node_id].data_id[i] - 1].id, objects[set_of_nodes[itr_order1->node_id].data_id[j] - 1].id, score, spacial_score, textual_score};
						answer.insert(pair1);
					}
					else
					{
						pair1 = {objects[set_of_nodes[itr_order1->node_id].data_id[i] - 1].id, objects[set_of_nodes[itr_order1->node_id].data_id[j] - 1].id, score, spacial_score, textual_score};
						answer.insert(pair1);
						auto itr_ans1 = answer.begin();
						answer.erase(itr_ans1);
						itr_ans1 = answer.begin();
						tau = itr_ans1->score;
					}
				}
			}
		}
	}

	//閾値の初期値、閾値の初期値取得時間　出力
	mesure_time_end = std::chrono::system_clock::now();
	time = std::chrono::duration_cast<std::chrono::microseconds>(mesure_time_end - mesure_time_start).count();
	cout << "tau = " << fixed << setprecision(6) << tau << endl;
	cout << "init knn calclation time = " << fixed << setprecision(6) << time / 1000000 << "sec" << endl;

	//Top-k spatio-textual similarity joinの実行時間　計測開始
	cout << "exe join" << endl;
	Node current_node;
	bool current_node_flag;
	value_type node_dist;
	vector<NodePair> node_pair_list;
	NodePair node_pair;
	int max_num_key;
	value_type time_node_to_node, time_node_to_data, time_calc_join, time_thread;
	int count_thread;

	node_order.clear();
	for (auto itr = set_of_nodes.begin(); itr != set_of_nodes.end(); itr++)
	{
		NodeScore temp_node = {itr->second.id, 1, false};
		if (itr->second.isleaf == true)
		{
			node_order.push_back(temp_node);
		}
	}
	/*for (auto itr = node_order2.begin(); itr != node_order2.end(); itr++)
	{
		node_order.push_back(*itr);
	}

	for (int j = 0; j < L; j++)
	{
		auto itr_join = node_order.end();
		itr_join--;
		itr_join->flag = true;
	}*/

	auto itr_ans1 = answer.begin();

	while (node_order.size() != 0)
	{
		mesure_time_start = std::chrono::system_clock::now();
		auto itr_order = node_order.end();
		itr_order--;
		current_node = set_of_nodes[itr_order->node_id];
		current_node_flag = itr_order->flag;
		node_order.erase(itr_order);
		node_pair_list.clear();
		count_thread = 0;

		for (int i = 0; i < node_order.size(); i++)
		{
			//mesure_time_start1 = std::chrono::system_clock::now();
			int max_num_key1;
			NodePair node_pair1;
			double node_dist1 = dist_nodes(current_node, set_of_nodes[node_order[i].node_id]);
			double spacial_score1 = 1 - (node_dist1 / dist_max);

			int share_parent_node_id = share_parent_node(current_node, set_of_nodes[node_order[i].node_id]);
			if (current_node.max_num_key < set_of_nodes[node_order[i].node_id].max_num_key)
			{
				max_num_key1 = set_of_nodes[node_order[i].node_id].max_num_key;
			}
			if (current_node.max_num_key > set_of_nodes[node_order[i].node_id].max_num_key)
			{
				max_num_key1 = current_node.max_num_key;
			}
			if (max_num_key1 == 0)
			{
				max_num_key1 = 1;
			}
			double textual_score1 = set_of_nodes[share_parent_node_id].jaccard_max;
			if (textual_score1 > 1)
			{
				textual_score1 = 1;
			}
			//cout << textual_score << endl;
			double score1 = alpha * spacial_score1 + (1 - alpha) * textual_score1;
			//num_of_calclation1++;
			if (score1 > tau)
			{
				node_pair1.node_id1 = current_node.id;
				node_pair1.node_id2 = set_of_nodes[node_order[i].node_id].id;
				node_pair1.score = score1;
			}
		}

		mesure_time_end = std::chrono::system_clock::now();
		time = std::chrono::duration_cast<std::chrono::microseconds>(mesure_time_end - mesure_time_start).count();

		time_node_to_node = time_node_to_node + time;

		if (current_node_flag == false)
		{
			node_pair.node_id1 = current_node.id;
			node_pair.node_id2 = current_node.id;
			node_pair.score = 1;
			node_pair_list.push_back(node_pair);
		}

		for (auto itr = node_pair_list.begin(); itr != node_pair_list.end(); itr++)
		{
			if (itr->score > tau)
			{
				if (itr->node_id1 != itr->node_id2)
				{
					for (int i = 0; i < set_of_nodes[itr->node_id1].data_id.size(); i++)
					{
						mesure_time_start = std::chrono::system_clock::now();

						node_dist = dist_data_to_node(objects[set_of_nodes[itr->node_id1].data_id[i] - 1], set_of_nodes[itr->node_id2]);
						spacial_score = 1 - (node_dist / dist_max);
						for (int j = 0; j < set_of_nodes[itr->node_id2].key.size(); j++)
						{
							for (int k = 0; k < objects[set_of_nodes[itr->node_id1].data_id[i] - 1].key.size(); k++)
							{
								if (set_of_nodes[itr->node_id2].key[j] == objects[set_of_nodes[itr->node_id1].data_id[i] - 1].key[k])
								{
									textual_score++;
								}
							}
						}
						textual_score /= objects[set_of_nodes[itr->node_id1].data_id[i] - 1].key.size();
						score = alpha * spacial_score + (1 - alpha) * textual_score;
						num_of_calclation3++;
						mesure_time_end = std::chrono::system_clock::now();
						time = std::chrono::duration_cast<std::chrono::microseconds>(mesure_time_end - mesure_time_start).count();
						time_node_to_data = time_node_to_data + time;
						if (score > tau)
						{
							mesure_time_start = std::chrono::system_clock::now();
							num_of_calclation4++;
							for (int j = 0; j < set_of_nodes[itr->node_id2].data_id.size(); j++)
							{
								spacial_score = calc_spatial_score(objects[set_of_nodes[itr->node_id1].data_id[i] - 1], objects[set_of_nodes[itr->node_id2].data_id[j] - 1], dist_max);
								textual_score = calc_textual_score(objects[set_of_nodes[itr->node_id1].data_id[i] - 1], objects[set_of_nodes[itr->node_id2].data_id[j] - 1]);
								score = alpha * spacial_score + (1 - alpha) * textual_score;
								num_of_calclation++;
								if (score > tau)
								{
									if (answer.size() < TOP_K)
									{
										pair1 = {objects[set_of_nodes[itr->node_id1].data_id[i] - 1].id, objects[set_of_nodes[itr->node_id2].data_id[j] - 1].id, score, spacial_score, textual_score};
										answer.insert(pair1);
									}
									else
									{
										pair1 = {objects[set_of_nodes[itr->node_id1].data_id[i] - 1].id, objects[set_of_nodes[itr->node_id2].data_id[j] - 1].id, score, spacial_score, textual_score};
										answer.insert(pair1);
										itr_ans1 = answer.begin();
										answer.erase(itr_ans1);
										itr_ans1 = answer.begin();
										tau = itr_ans1->score;
									}
								}
							}
							mesure_time_end = std::chrono::system_clock::now();
							time = std::chrono::duration_cast<std::chrono::microseconds>(mesure_time_end - mesure_time_start).count();
							time_calc_join = time_calc_join + time;
						}
					}
				}
				else
				{
					mesure_time_start = std::chrono::system_clock::now();
					for (int i = 0; i < set_of_nodes[itr->node_id1].data_id.size(); i++)
					{
						for (int j = 0; j < set_of_nodes[itr->node_id2].data_id.size(); j++)
						{
							if (i != j)
							{
								spacial_score = calc_spatial_score(objects[set_of_nodes[itr->node_id1].data_id[i] - 1], objects[set_of_nodes[itr->node_id2].data_id[j] - 1], dist_max);
								textual_score = calc_textual_score(objects[set_of_nodes[itr->node_id1].data_id[i] - 1], objects[set_of_nodes[itr->node_id2].data_id[j] - 1]);
								score = alpha * spacial_score + (1 - alpha) * textual_score;
								num_of_calclation++;
								if (score > tau)
								{
									if (answer.size() < TOP_K)
									{
										pair1 = {objects[set_of_nodes[itr->node_id1].data_id[i] - 1].id, objects[set_of_nodes[itr->node_id2].data_id[j] - 1].id, score, spacial_score, textual_score};
										answer.insert(pair1);
									}
									else
									{
										pair1 = {objects[set_of_nodes[itr->node_id1].data_id[i] - 1].id, objects[set_of_nodes[itr->node_id2].data_id[j] - 1].id, score, spacial_score, textual_score};
										answer.insert(pair1);
										itr_ans1 = answer.begin();
										answer.erase(itr_ans1);
										itr_ans1 = answer.begin();
										tau = itr_ans1->score;
									}
								}
							}
						}
					}
					mesure_time_end = std::chrono::system_clock::now();
					time = std::chrono::duration_cast<std::chrono::microseconds>(mesure_time_end - mesure_time_start).count();
					time_calc_join = time_calc_join + time;
				}
			}
		}
	}

	end = std::chrono::system_clock::now();
	time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

	cout << "overall time = " << time << "sec" << endl;
	cout << "filtering time(node to node) = " << fixed << setprecision(6) << time_node_to_node / 1000000 << "sec" << endl;
	cout << "filtering time(node to data) = " << fixed << setprecision(6) << time_node_to_data / 1000000 << "sec" << endl;
	cout << "join calclation time = " << fixed << setprecision(6) << time_calc_join / 1000000 << "sec" << endl;
	cout << "tau = " << fixed << setprecision(6) << tau << endl;

	int num_knn = 1;
	/*for (auto itr = answer.begin(); itr != answer.end(); itr++)
	{
		cout << num_knn << " : data1 = " << itr->data1_id;
		cout << ", data2 = " << itr->data2_id;
		cout << ", score = " << itr->score << endl;
		num_knn++;
	}*/

	return 0;
}