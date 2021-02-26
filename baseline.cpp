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
    value_type dist_max;                                                                                                        //位置スコア計算用
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
    set_data(objects2); //Twitter
    //set_data_place(objects2);	//Place

    Range whole_region;

    whole_region.x_max = 55.0;
    whole_region.x_min = 20.0;
    whole_region.y_max = -45.0;
    whole_region.y_min = -130.0;

    dist_max = sqrt(pow(abs(whole_region.x_max - whole_region.x_min), 2.0) + pow(abs(whole_region.y_max - whole_region.y_min), 2.0));

    Range region = {0, 0, 0, 0};

    Rtree rtree;
    Node root_node;
    Node init_node;

    init_node.MBR = region;
    init_node.isleaf = true;
    init_node.id = 1;
    init_node.depth = 1;
    set_of_nodes[init_node.id] = init_node;

    root_node.MBR = region;
    root_node.isleaf = false;
    root_node.child_node[1] = init_node.id;
    root_node.id = 0;
    root_node.parent_node = -1;
    root_node.depth = 0;
    set_of_nodes[root_node.id] = root_node;
    set_of_nodes[init_node.id].parent_node = root_node.id;

    rtree.root_node = &set_of_nodes[root_node.id];

    start = std::chrono::system_clock::now();
    cout << "baseline" << endl;
    cout << "TOP_K = " << TOP_K << endl;
    cout << "num = " << num << endl;
    cout << "NODE_MAX = " << NODE_MAX << endl;
    cout << "alpha = " << alpha << endl;
    cout << "insert rtree" << endl;

    mesure_time_start = std::chrono::system_clock::now();

    mt19937 rnd1(1);
    int rand;

    for (int i = 0; i < num; i++)
    {
        rand = rnd1() % num;

        objects[i] = objects2[rand];
    }

    for (int i = 0; i < num; i++)
    {
        //cout << i << endl;
        insert(root_node.id, objects2[i]);
    }

    mesure_time_end = std::chrono::system_clock::now();
    value_type time = std::chrono::duration_cast<std::chrono::microseconds>(mesure_time_end - mesure_time_start).count();

    cout << "insert time(rtree)= " << fixed << setprecision(6) << time / 1000000 << "sec" << endl;

    random_device rnd;
    multiset<AnswerPair> answer;
    value_type tau;
    map<string, int> key_frecency;
    multiset<Key_tf> key_tf;
    int current_node_id;

    int rand_k = 0;
    int rand1, rand2;
    value_type spacial_score, textual_score, score;
    AnswerPair pair1;
    vector<bool> init_flag;
    vector<vector<bool>> flag_pair;
    int num_of_calclation = 0;

    for (int i = 0; i < num; i++)
    {
        init_flag.push_back(false);
    }

    for (int i = 0; i < num; i++)
    {
        flag_pair.push_back(init_flag);
    }
    for (int i = 0; i < num; i++)
    {
        flag_pair[i][i] = true;
    }

    mesure_time_start = std::chrono::system_clock::now();

    cout << "calc init knn" << endl;

    while (1)
    {
        rand1 = rnd() % 10000;
        auto itr = set_of_nodes.find(rand1);
        if (itr != set_of_nodes.end())
        {
            if (itr->second.isleaf == true)
            {
                for (int i = 0; i < itr->second.data_id.size(); i++)
                {
                    for (int j = 0; j < itr->second.data_id.size(); j++)
                    {
                        if (i > j)
                        {
                            spacial_score = calc_spatial_score(objects[itr->second.data_id[i] - 1], objects[itr->second.data_id[j] - 1], dist_max);
                            textual_score = calc_textual_score(objects[itr->second.data_id[i] - 1], objects[itr->second.data_id[j] - 1]);
                            score = alpha * spacial_score + (1 - alpha) * textual_score;

                            pair1 = {objects[itr->second.data_id[i] - 1].id, objects[itr->second.data_id[j] - 1].id, score, spacial_score, textual_score};
                            answer.insert(pair1);
                            rand_k++;
                            if (rand_k == TOP_K)
                            {
                                break;
                            }
                        }
                    }
                    if (rand_k == TOP_K)
                    {
                        break;
                    }
                }
            }
        }
        if (rand_k == TOP_K)
        {
            break;
        }
    }

    auto itr_ans = answer.begin();
    tau = itr_ans->score;

    mesure_time_end = std::chrono::system_clock::now();
    time = std::chrono::duration_cast<std::chrono::microseconds>(mesure_time_end - mesure_time_start).count();

    cout << "init knn calclation time = " << fixed << setprecision(6) << time / 1000000 << "sec" << endl;

    cout << "tau = " << fixed << setprecision(6) << tau << endl;

    cout << "exe join" << endl;
    mesure_time_start = std::chrono::system_clock::now();

    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < objects[i].key.size(); j++)
        {
            auto itr = key_frecency.find(objects[i].key[j]);
            if (itr != key_frecency.end())
            {
                itr->second++;
            }
            else
            {
                key_frecency[objects[i].key[j]] = 1;
            }
        }
    }

    for (auto itr = key_frecency.begin(); itr != key_frecency.end(); itr++)
    {
        Key_tf tf = {itr->first, itr->second};
        key_tf.insert(tf);
    }

    for (int i = 0; i < num; i++)
    {
        string key_x;
        for (int j = 0; j < objects[i].key.size() - 1; j++)
        {
            for (int k = objects[i].key.size() - 1; k > i; k--)
            {
                if (key_frecency[objects[i].key[k]] > key_frecency[objects[i].key[k - 1]])
                {
                    key_x = objects[i].key[k];
                    objects[i].key[k] = objects[i].key[k - 1];
                    objects[i].key[k - 1] = key_x;
                }
            }
        }
    }

    vector<string> pivot_term;
    map<string, map<int, vector<int>>> invert;
    map<int, vector<int>> invert_minus;
    int pivot;
    value_type node_dist;
    value_type min_node_dist;
    value_type calclated_node_dist;
    value_type lowerbound_t;
    value_type lowerbound;
    Triple triple;
    multiset<Triple> order_triple;

    for (int i = 0; i < num; i++)
    {
        //cout << i << endl;
        current_node_id = find_node(root_node.id, objects[i]);
        for (int j = 0; j < objects[i].key.size(); j++)
        {
            textual_score = (double)(objects[i].key.size() - j) / objects[i].key.size();
            spacial_score = 1;
            lowerbound = alpha * spacial_score + (1 - alpha) * textual_score;
            triple = {i + 1, objects[i].key[j], current_node_id, lowerbound};
            order_triple.insert(triple);
            // cout << "score(leaf)" << lowerbound << endl;
        }
        while (current_node_id != 0)
        {
            if (current_node_id == -1)
            {
                break;
            }
            current_node_id = set_of_nodes[current_node_id].parent_node;
            child_node_id2.clear();
            get_all_node(current_node_id);
            min_node_dist = 10000000;
            for (auto itr = child_node_id2.begin(); itr != child_node_id2.end(); itr++)
            {
                node_dist = dist_data_to_node(objects[i], set_of_nodes[*itr]);
                if (node_dist < min_node_dist)
                {
                    min_node_dist = node_dist;
                }
            }
            spacial_score = 1 - (min_node_dist / dist_max);
            for (int j = 0; j < objects[i].key.size(); j++)
            {
                textual_score = (double)(objects[i].key.size() - j) / objects[i].key.size();
                lowerbound = alpha * spacial_score + (1 - alpha) * textual_score;
                triple = {i + 1, objects[i].key[j], current_node_id, lowerbound};
                order_triple.insert(triple);
                //cout << "score" << lowerbound << endl;
            }
        }
    }

    mesure_time_end = std::chrono::system_clock::now();
    time = std::chrono::duration_cast<std::chrono::microseconds>(mesure_time_end - mesure_time_start).count();

    cout << "sort time = " << fixed << setprecision(6) << time / 1000000 << "sec" << endl;

    mesure_time_start = std::chrono::system_clock::now();

    cout << order_triple.size() << endl;
    Triple current_triple;
    bool join_flag = false;
    while (join_flag == false)
    {
        //cout << "tau = " << tau << endl;
        auto itr_triple = order_triple.end();
        itr_triple--;
        current_triple = *itr_triple;
        order_triple.erase(itr_triple);
        if (current_triple.score < tau)
        {
            join_flag = true;
        }
        for (auto itr = invert[current_triple.key][current_triple.node_id].begin(); itr != invert[current_triple.key][current_triple.node_id].end(); itr++)
        {
            if (flag_pair[current_triple.data_id - 1][*itr - 1] == false)
            {
                spacial_score = calc_spatial_score(objects[current_triple.data_id - 1], objects[*itr - 1], dist_max);
                textual_score = calc_textual_score(objects[current_triple.data_id - 1], objects[*itr - 1]);
                score = alpha * spacial_score + (1 - alpha) * textual_score;
                num_of_calclation++;
                if (score > tau)
                {
                    pair1 = {objects[current_triple.data_id - 1].id, objects[*itr - 1].id, score, spacial_score, textual_score};
                    answer.insert(pair1);
                    auto itr_ans1 = answer.begin();
                    answer.erase(itr_ans1);
                    itr_ans1 = answer.begin();
                    tau = itr_ans1->score;
                    flag_pair[current_triple.data_id - 1][*itr - 1] = true;
                    flag_pair[*itr - 1][current_triple.data_id - 1] = true;
                }
            }
        }
        invert[current_triple.key][current_triple.node_id].push_back(objects[current_triple.data_id - 1].id);
    }

    mesure_time_end = std::chrono::system_clock::now();
    time = std::chrono::duration_cast<std::chrono::microseconds>(mesure_time_end - mesure_time_start).count();

    cout << "join  = " << fixed << setprecision(6) << time / 1000000 << "sec" << endl;

    end = std::chrono::system_clock::now();
    time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    //cout << time << "sec" << fixed << setprecision(6) << endl;
    cout << "tau = " << fixed << setprecision(6) << tau << endl;
    //cout << "number_of_calclation = " << num_of_calclation << endl;

    /*int num_knn = 1;
    for (auto itr = answer.begin(); itr != answer.end(); itr++)
    {
        cout << num_knn << " : data1 = " << itr->data1_id;
        cout << ", data2 = " << itr->data2_id;
        cout << ", score = " << itr->score << endl;
        num_knn++;
    }*/

    return 0;
}