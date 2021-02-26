#include "class.h"

//データ読み込み
void set_data(std::vector<Data> &z);
void set_data_place(std::vector<Data> &z);

Range set_MBR(vector<int> data_id);
Range set_MBR(Node node);
bool collision_MBR(Range current_MBR1, Range current_MBR2);
value_type min_dist_MBR_x(Range current_MBR1, Range current_MBR2);
value_type min_dist_MBR_y(Range current_MBR1, Range current_MBR2);
void insert(int node, Data z);
void divide_leafnode(int node, Data z);
void divide_nonleafnode(int nonleafnode);

int find_node(int node, Data z);

vector<Data> objects;
map<int, Node> set_of_nodes;
vector<int> child_node_id;
int share_parent_node_id;

int node_id = 2;

int T = 0;

value_type calc_spatial_score(Data z1, Data z2, double dist_max);
value_type calc_textual_score(Data z1, Data z2);

int main(int argc, char *argv[])
{
    value_type dist_max;
    std::chrono::system_clock::time_point start, end;
    //データの初期化
    for (int k = 0; k < num; k++)
    {
        objects.push_back({
            k,
            0,
            0,
        });
    }

    //データ読み込み
    set_data(objects); //Twitter
                       //set_data_place(objects);	//Place

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
    cout << "linear" << endl;
    cout << "exe knn" << endl;

    //int a = find_node(root_node.id, objects[0]);

    random_device rnd;
    multiset<AnswerPair> answer;
    value_type tau;
    map<string, int> key_frecency;
    multiset<Key_tf> key_tf;
    int current_node_id;

    int rand_k = 1;
    int rand1, rand2;
    value_type spacial_score, textual_score, score;
    AnswerPair pair1;
    vector<bool> init_flag;
    vector<vector<bool>> flag_pair;

    start = std::chrono::system_clock::now();

    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < num; j++)
        {
            if (i != j)
            {
                spacial_score = calc_spatial_score(objects[i], objects[j], dist_max);
                textual_score = calc_textual_score(objects[i], objects[j]);
                score = alpha * spacial_score + (1 - alpha) * textual_score;
                pair1 = {objects[i].id, objects[j].id, score, spacial_score, textual_score};
                answer.insert(pair1);
                if (answer.size() > TOP_K)
                {
                    auto itr_ans1 = answer.begin();
                    answer.erase(itr_ans1);
                    itr_ans1 = answer.begin();
                    tau = itr_ans1->score;
                }
            }
        }
    }

    end = std::chrono::system_clock::now();
    value_type time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    cout << time << "sec" << endl;

    int num_knn = 1;
    for (auto itr = answer.begin(); itr != answer.end(); itr++)
    {
        cout << num_knn << " : data1 = " << itr->data1_id;
        cout << ", data2 = " << itr->data2_id;
        cout << ", score = " << itr->score << endl;
        num_knn++;
    }

    return 0;
}