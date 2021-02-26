#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <random>
#include <chrono>
#include <utility>
#include <iomanip>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>

#define MAX 2048

#define TOP_K 100
#define num 1000000
#define alpha 0.5
#define NODE_MAX 1200
#define L TOP_K

typedef double value_type;

using namespace std;

class Data
{
public:
    int id;
    value_type x;
    value_type y;
    vector<string> key;
};

class Range
{
public:
    value_type x_max;
    value_type x_min;
    value_type y_max;
    value_type y_min;
};

class AnswerPair
{
public:
    int data1_id;
    int data2_id;
    double score;
    double spacial_score;
    double textual_score;
    bool operator<(const AnswerPair &lhs) const
    {
        return score < lhs.score;
    }
};

class NodeScore
{
public:
    int node_id;
    double score;
    bool flag;
    bool operator<(const NodeScore &lhs) const
    {
        return score < lhs.score;
    }
};
class NodePair
{
public:
    int node_id1;
    int node_id2;
    value_type score;
};

class Key_tf
{
public:
    string key;
    int tf;
    bool operator<(const Key_tf &lhs) const
    {
        return tf < lhs.tf;
    }
};

class Triple
{
public:
    int data_id;
    string key;
    int node_id;
    value_type score;
    bool operator<(const Triple &lhs) const
    {
        return score < lhs.score;
    }
};

//データ グローバル変数

//ノードクラス
class Node
{
public:
    Node();
    ~Node();
    void set_MBR(Data z);
    Range set_MBR(Node node);
    Range set_MBR(Range original_MBR, Data z);
    Range set_MBR(vector<int> data_id);

    bool ismax();
    bool innerMBR(Data z);
    value_type area_MBR();
    value_type area_MBR(Range current_MBR);
    value_type increase_MBR(Data z);
    value_type min_dist_MBR_x(Range current_MBR1, Range current_MBR2);
    value_type min_dist_MBR_y(Range current_MBR1, Range current_MBR2);
    bool collision_MBR(Range current_MBR1, Range current_MBR2);

    int id;
    int depth;
    Range MBR;
    map<int, int> child_node;
    int parent_node;
    bool isleaf;
    vector<int> data_id;
    vector<int> share_key;
    vector<string> key;
    float jaccard_max;
    int max_share_key;
    int max_num_key;
};

class Rtree
{
public:
    Rtree();
    ~Rtree();
    void insert(Data z);

    Node *root_node;
};