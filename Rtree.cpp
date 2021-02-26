#include "class.h"

extern vector<Data> objects;
extern int node_id;
extern map<int, Node> set_of_nodes;
extern vector<int> child_node_id;
extern vector<int> child_node_id2;
extern int share_parent_node_id;

Node::Node() {}
Node::~Node() {}

//位置スコア計算
value_type calc_spatial_score(Data z1, Data z2, double dist_max)
{
    value_type score, dist;
    dist = sqrt(pow(abs(z1.x - z2.x), 2.0) + pow(abs(z1.y - z2.y), 2.0));
    score = 1 - (dist / dist_max);
    return score;
}

//キーワードスコア計算
value_type calc_textual_score(Data z1, Data z2)
{
    value_type score;
    set<string> set_key;
    int num_key = 0;
    set_key.clear();
    for (int i = 0; i < z1.key.size(); i++)
    {
        set_key.insert(z1.key[i]);
    }
    for (int i = 0; i < z2.key.size(); i++)
    {
        set_key.insert(z2.key[i]);
    }

    for (int i = 0; i < z1.key.size(); i++)
    {
        for (int j = 0; j < z2.key.size(); j++)
        {
            if (z1.key[i] == z2.key[j])
            {
                num_key++;
            }
        }
    }
    score = (double)num_key / (set_key.size()); //Jaccard
    //score = (double)num_key / sqrt(abs(z1.key.size() * z2.key.size()));	//Cosine
    //score = (double)2 * num_key / (z1.key.size() + z2.key.size()); //Dice
    return score;
}

void Node::set_MBR(Data z)
{
    if (MBR.x_max < z.x)
    {
        MBR.x_max = z.x;
    }
    if (MBR.x_min > z.x)
    {
        MBR.x_min = z.x;
    }
    if (MBR.y_max < z.y)
    {
        MBR.y_max = z.y;
    }
    if (MBR.y_min > z.y)
    {
        MBR.y_min = z.y;
    }
};

Range Node::set_MBR(Node node)
{
    Range init_MBR;
    auto itr1 = node.child_node.begin();
    init_MBR = set_of_nodes[itr1->second].MBR;
    for (auto itr = node.child_node.begin(); itr != node.child_node.end(); itr++)
    {
        if (init_MBR.x_max < set_of_nodes[itr->second].MBR.x_max)
        {
            init_MBR.x_max = set_of_nodes[itr->second].MBR.x_max;
        }
        if (init_MBR.x_min > set_of_nodes[itr->second].MBR.x_min)
        {
            init_MBR.x_min = set_of_nodes[itr->second].MBR.x_min;
        }
        if (init_MBR.y_max < set_of_nodes[itr->second].MBR.y_max)
        {
            init_MBR.y_max = set_of_nodes[itr->second].MBR.y_max;
        }
        if (init_MBR.y_min > set_of_nodes[itr->second].MBR.y_min)
        {
            init_MBR.y_min = set_of_nodes[itr->second].MBR.y_min;
        }
    }
    return init_MBR;
};

Range Node::set_MBR(Range original_MBR, Data z)
{
    if (original_MBR.x_max < z.x)
    {
        original_MBR.x_max = z.x;
    }
    if (original_MBR.x_min > z.x)
    {
        original_MBR.x_min = z.x;
    }
    if (original_MBR.y_max < z.y)
    {
        original_MBR.y_max = z.y;
    }
    if (original_MBR.y_min > z.y)
    {
        original_MBR.y_min = z.y;
    }
    return original_MBR;
};

Range Node::set_MBR(vector<int> data_id)
{
    Range init_MBR;
    init_MBR.x_max = objects[data_id[0]].x;
    init_MBR.x_min = objects[data_id[0]].x;
    init_MBR.y_max = objects[data_id[0]].y;
    init_MBR.y_min = objects[data_id[0]].y;
    for (int i = 1; i < data_id.size(); i++)
    {
        if (init_MBR.x_max < objects[data_id[i]].x)
        {
            init_MBR.x_max = objects[data_id[i]].x;
        }
        if (init_MBR.x_min > objects[data_id[i]].x)
        {
            init_MBR.x_min = objects[data_id[i]].x;
        }
        if (init_MBR.y_max < objects[data_id[i]].y)
        {
            init_MBR.y_max = objects[data_id[i]].y;
        }
        if (init_MBR.y_min > objects[data_id[i]].y)
        {
            init_MBR.y_min = objects[data_id[i]].y;
        }
    }
    return init_MBR;
};

bool Node::ismax()
{
    if (data_id.size() >= NODE_MAX)
    {
        return true;
    }
    else
    {
        return false;
    }
};

bool Node::innerMBR(Data z)
{
    if (MBR.x_max >= z.x && MBR.x_min <= z.x && MBR.y_max >= z.y && MBR.y_min <= z.y)
    {
        return true;
    }
    else
    {
        return false;
    }
};

value_type Node::area_MBR()
{
    value_type area;
    area = abs(MBR.x_max - MBR.x_min) * abs(MBR.y_max - MBR.y_min);
    return area;
};

value_type Node::area_MBR(Range current_MBR)
{
    value_type area;
    area = abs(current_MBR.x_max - current_MBR.x_min) * abs(current_MBR.y_max - current_MBR.y_min);
    return area;
};

value_type Node::increase_MBR(Data z)
{
    value_type original_area = area_MBR();
    Range increased_MBR = set_MBR(MBR, z);
    value_type increased_area = area_MBR(increased_MBR);
    return increased_area - original_area;
};

bool Node::collision_MBR(Range current_MBR1, Range current_MBR2)
{
    value_type center_x1, center_y1, center_x2, center_y2;
    value_type width1, height1, width2, height2;

    center_x1 = (current_MBR1.x_max + current_MBR1.x_min) / 2;
    center_y1 = (current_MBR1.y_max + current_MBR1.y_min) / 2;
    center_x2 = (current_MBR2.x_max + current_MBR2.x_min) / 2;
    center_y2 = (current_MBR2.y_max + current_MBR2.y_min) / 2;

    width1 = abs(current_MBR1.x_max - current_MBR1.x_min);
    height1 = abs(current_MBR1.y_max - current_MBR1.y_min);
    width2 = abs(current_MBR2.x_max - current_MBR2.x_min);
    height2 = abs(current_MBR2.y_max - current_MBR2.y_min);

    if (abs(center_x1 - center_x2) < width1 / 2 + width2 / 2 && abs(center_y1 - center_y2) < height1 / 2 + height2 / 2)
    {
        return true;
    }
    else
    {
        return false;
    }
};

value_type Node::min_dist_MBR_x(Range current_MBR1, Range current_MBR2)
{
    value_type min_dist;
    value_type dist1, dist2;
    if (collision_MBR(current_MBR1, current_MBR2) == true)
    {
        min_dist = 0;
        return min_dist;
    }
    else
    {
        dist1 = abs(current_MBR1.x_max - current_MBR2.x_min);
        dist2 = abs(current_MBR2.x_max - current_MBR1.x_min);
        if (dist1 < dist2)
        {
            return dist1;
        }
        else
        {
            return dist2;
        }
    }
};

value_type Node::min_dist_MBR_y(Range current_MBR1, Range current_MBR2)
{
    value_type min_dist;
    value_type dist1, dist2;
    if (collision_MBR(current_MBR1, current_MBR2) == true)
    {
        min_dist = 0;
        return min_dist;
    }
    else
    {
        dist1 = abs(current_MBR1.y_max - current_MBR2.y_min);
        dist2 = abs(current_MBR2.y_max - current_MBR1.y_min);
        if (dist1 < dist2)
        {
            return dist1;
        }
        else
        {
            return dist2;
        }
    }
};

Rtree::Rtree() {}
Rtree::~Rtree() {}

Range set_MBR(vector<int> data_id)
{
    Range init_MBR;
    init_MBR.x_max = objects[data_id[0] - 1].x;
    init_MBR.x_min = objects[data_id[0] - 1].x;
    init_MBR.y_max = objects[data_id[0] - 1].y;
    init_MBR.y_min = objects[data_id[0] - 1].y;
    for (int i = 1; i < data_id.size(); i++)
    {
        if (init_MBR.x_max < objects[data_id[i] - 1].x)
        {
            init_MBR.x_max = objects[data_id[i] - 1].x;
        }
        if (init_MBR.x_min > objects[data_id[i] - 1].x)
        {
            init_MBR.x_min = objects[data_id[i] - 1].x;
        }
        if (init_MBR.y_max < objects[data_id[i] - 1].y)
        {
            init_MBR.y_max = objects[data_id[i] - 1].y;
        }
        if (init_MBR.y_min > objects[data_id[i] - 1].y)
        {
            init_MBR.y_min = objects[data_id[i] - 1].y;
        }
    }
    return init_MBR;
};

Range set_MBR(Node node)
{
    Range init_MBR;
    auto itr1 = node.child_node.begin();
    init_MBR = set_of_nodes[itr1->second].MBR;
    for (auto itr = node.child_node.begin(); itr != node.child_node.end(); itr++)
    {
        if (init_MBR.x_max < set_of_nodes[itr->second].MBR.x_max)
        {
            init_MBR.x_max = set_of_nodes[itr->second].MBR.x_max;
        }
        if (init_MBR.x_min > set_of_nodes[itr->second].MBR.x_min)
        {
            init_MBR.x_min = set_of_nodes[itr->second].MBR.x_min;
        }
        if (init_MBR.y_max < set_of_nodes[itr->second].MBR.y_max)
        {
            init_MBR.y_max = set_of_nodes[itr->second].MBR.y_max;
        }
        if (init_MBR.y_min > set_of_nodes[itr->second].MBR.y_min)
        {
            init_MBR.y_min = set_of_nodes[itr->second].MBR.y_min;
        }
    }
    return init_MBR;
};

bool collision_MBR(Range current_MBR1, Range current_MBR2)
{
    value_type center_x1, center_y1, center_x2, center_y2;
    value_type width1, height1, width2, height2;

    center_x1 = (current_MBR1.x_max + current_MBR1.x_min) / 2;
    center_y1 = (current_MBR1.y_max + current_MBR1.y_min) / 2;
    center_x2 = (current_MBR2.x_max + current_MBR2.x_min) / 2;
    center_y2 = (current_MBR2.y_max + current_MBR2.y_min) / 2;

    width1 = abs(current_MBR1.x_max - current_MBR1.x_min);
    height1 = abs(current_MBR1.y_max - current_MBR1.y_min);
    width2 = abs(current_MBR2.x_max - current_MBR2.x_min);
    height2 = abs(current_MBR2.y_max - current_MBR2.y_min);

    if (abs(center_x1 - center_x2) < (width1 / 2) + (width2 / 2) && abs(center_y1 - center_y2) < (height1 / 2) + (height2 / 2))
    {
        return true;
    }
    else
    {
        return false;
    }
};

value_type min_dist_MBR_x(Range current_MBR1, Range current_MBR2)
{
    value_type min_dist;
    value_type dist1, dist2;
    if (collision_MBR(current_MBR1, current_MBR2) == true)
    {
        min_dist = 0;
        return min_dist;
    }
    else
    {
        dist1 = abs(current_MBR1.x_max - current_MBR2.x_min);
        dist2 = abs(current_MBR2.x_max - current_MBR1.x_min);
        if (dist1 < dist2)
        {
            return dist1;
        }
        else
        {
            return dist2;
        }
    }
};

value_type min_dist_MBR_y(Range current_MBR1, Range current_MBR2)
{
    value_type min_dist;
    value_type dist1, dist2;
    if (collision_MBR(current_MBR1, current_MBR2) == true)
    {
        min_dist = 0;
        return min_dist;
    }
    else
    {
        dist1 = abs(current_MBR1.y_max - current_MBR2.y_min);
        dist2 = abs(current_MBR2.y_max - current_MBR1.y_min);
        if (dist1 < dist2)
        {
            return dist1;
        }
        else
        {
            return dist2;
        }
    }
};

void update_child_node(int parent_node, int node)
{
    set_of_nodes[node].parent_node = parent_node;
    set_of_nodes[node].depth = set_of_nodes[node].depth + 1;
    if (set_of_nodes[node].isleaf == false)
    {
        for (auto itr = set_of_nodes[node].child_node.begin(); itr != set_of_nodes[node].child_node.end(); itr++)
        {
            update_child_node(node, itr->second);
        }
    }
}

Range calc_group_MBR(map<int, int> group_node_id, int node_id)
{
    Range init_MBR;
    auto itr1 = group_node_id.begin();
    init_MBR = set_of_nodes[itr1->second].MBR;
    for (auto itr = group_node_id.begin(); itr != group_node_id.end(); itr++)
    {
        if (init_MBR.x_max < set_of_nodes[itr->second].MBR.x_max)
        {
            init_MBR.x_max = set_of_nodes[itr->second].MBR.x_max;
        }
        if (init_MBR.x_min > set_of_nodes[itr->second].MBR.x_min)
        {
            init_MBR.x_min = set_of_nodes[itr->second].MBR.x_min;
        }
        if (init_MBR.y_max < set_of_nodes[itr->second].MBR.y_max)
        {
            init_MBR.y_max = set_of_nodes[itr->second].MBR.y_max;
        }
        if (init_MBR.y_min > set_of_nodes[itr->second].MBR.y_min)
        {
            init_MBR.y_min = set_of_nodes[itr->second].MBR.y_min;
        }
    }

    if (init_MBR.x_max < set_of_nodes[node_id].MBR.x_max)
    {
        init_MBR.x_max = set_of_nodes[node_id].MBR.x_max;
    }
    if (init_MBR.x_min > set_of_nodes[node_id].MBR.x_min)
    {
        init_MBR.x_min = set_of_nodes[node_id].MBR.x_min;
    }
    if (init_MBR.y_max < set_of_nodes[node_id].MBR.y_max)
    {
        init_MBR.y_max = set_of_nodes[node_id].MBR.y_max;
    }
    if (init_MBR.y_min > set_of_nodes[node_id].MBR.y_min)
    {
        init_MBR.y_min = set_of_nodes[node_id].MBR.y_min;
    }
    return init_MBR;
}

Range set_MBR(vector<int> data_id, int data_id1)
{
    Range init_MBR;
    init_MBR.x_max = objects[data_id[0] - 1].x;
    init_MBR.x_min = objects[data_id[0] - 1].x;
    init_MBR.y_max = objects[data_id[0] - 1].y;
    init_MBR.y_min = objects[data_id[0] - 1].y;
    for (int i = 1; i < data_id.size(); i++)
    {
        if (init_MBR.x_max < objects[data_id[i] - 1].x)
        {
            init_MBR.x_max = objects[data_id[i] - 1].x;
        }
        if (init_MBR.x_min > objects[data_id[i] - 1].x)
        {
            init_MBR.x_min = objects[data_id[i] - 1].x;
        }
        if (init_MBR.y_max < objects[data_id[i] - 1].y)
        {
            init_MBR.y_max = objects[data_id[i] - 1].y;
        }
        if (init_MBR.y_min > objects[data_id[i] - 1].y)
        {
            init_MBR.y_min = objects[data_id[i] - 1].y;
        }
    }
    if (init_MBR.x_max < objects[data_id1 - 1].x)
    {
        init_MBR.x_max = objects[data_id1 - 1].x;
    }
    if (init_MBR.x_min > objects[data_id1 - 1].x)
    {
        init_MBR.x_min = objects[data_id1 - 1].x;
    }
    if (init_MBR.y_max < objects[data_id1 - 1].y)
    {
        init_MBR.y_max = objects[data_id1 - 1].y;
    }
    if (init_MBR.y_min > objects[data_id1 - 1].y)
    {
        init_MBR.y_min = objects[data_id1 - 1].y;
    }
    return init_MBR;
};

void divide_nonleafnode(int nonleafnode)
{
    map<int, int> group1_node, group2_node;
    int max_xdist_node_id1, max_xdist_node_id2, max_ydist_node_id1, max_ydist_node_id2;
    value_type max_xdist, max_ydist;
    value_type xdist, ydist;
    value_type x_axis_dist, y_axis_dist;
    value_type norm_xdist, norm_ydist;
    x_axis_dist = abs(set_of_nodes[nonleafnode].MBR.x_max - set_of_nodes[nonleafnode].MBR.x_min);
    y_axis_dist = abs(set_of_nodes[nonleafnode].MBR.y_max - set_of_nodes[nonleafnode].MBR.y_min);

    max_xdist = 0;
    max_ydist = 0;

    value_type area1, area2, area3, all_area, calc_area, calc_area1, calc_area2;
    value_type max_area = -10000000;
    int group1_id, group2_id;
    Range calc_MBR1, calc_MBR2;
    auto itr_del1 = set_of_nodes[nonleafnode].child_node.begin();
    auto itr_del2 = set_of_nodes[nonleafnode].child_node.begin();

    all_area = (abs(set_of_nodes[nonleafnode].MBR.x_max - set_of_nodes[nonleafnode].MBR.x_min)) * (abs(set_of_nodes[nonleafnode].MBR.y_max - set_of_nodes[nonleafnode].MBR.y_min));

    for (auto itr1 = set_of_nodes[nonleafnode].child_node.begin(); itr1 != set_of_nodes[nonleafnode].child_node.end(); itr1++)
    {
        for (auto itr2 = set_of_nodes[nonleafnode].child_node.begin(); itr2 != set_of_nodes[nonleafnode].child_node.end(); itr2++)
        {
            if (itr1 != itr2)
            {
                area1 = (abs(set_of_nodes[itr1->second].MBR.x_max - set_of_nodes[itr1->second].MBR.x_min)) * (abs(set_of_nodes[itr1->second].MBR.y_max - set_of_nodes[itr1->second].MBR.y_min));
                area2 = (abs(set_of_nodes[itr2->second].MBR.x_max - set_of_nodes[itr2->second].MBR.x_min)) * (abs(set_of_nodes[itr2->second].MBR.y_max - set_of_nodes[itr2->second].MBR.y_min));
                calc_area = all_area - area1 - area2;
                if (calc_area > max_area)
                {
                    group1_id = itr1->second;
                    group2_id = itr2->second;
                    max_area = calc_area;
                    itr_del1 = itr1;
                    itr_del2 = itr2;
                }
            }
        }
    }
    group1_node[group1_id] = group1_id;
    group2_node[group2_id] = group2_id;

    if (itr_del1 != itr_del2)
    {
        set_of_nodes[nonleafnode].child_node.erase(itr_del1);
        set_of_nodes[nonleafnode].child_node.erase(itr_del2);
    }
    else
    {
        itr_del1 = set_of_nodes[nonleafnode].child_node.begin();
        itr_del2 = itr_del1++;
        set_of_nodes[nonleafnode].child_node.erase(itr_del1);
        set_of_nodes[nonleafnode].child_node.erase(itr_del2);
    }

    while (set_of_nodes[nonleafnode].child_node.size() != 0)
    {
        auto itr_child = set_of_nodes[nonleafnode].child_node.begin();
        calc_MBR1 = calc_group_MBR(group1_node, itr_child->second);
        calc_MBR2 = calc_group_MBR(group2_node, itr_child->second);
        area1 = (abs(calc_MBR1.x_max - calc_MBR1.x_min)) * (abs(calc_MBR1.y_max - calc_MBR1.y_min));
        area2 = (abs(calc_MBR2.x_max - calc_MBR2.x_min)) * (abs(calc_MBR2.y_max - calc_MBR2.y_min));
        area3 = (abs(set_of_nodes[itr_child->second].MBR.x_max - set_of_nodes[itr_child->second].MBR.x_min)) * (abs(set_of_nodes[itr_child->second].MBR.y_max - set_of_nodes[itr_child->second].MBR.y_min));
        calc_area1 = all_area - area1 - area3;
        calc_area2 = all_area - area2 - area3;
        if (calc_area1 <= calc_area2)
        {
            group1_node[itr_child->second] = itr_child->second;
        }
        if (calc_area1 > calc_area2)
        {
            group2_node[itr_child->second] = itr_child->second;
        }
        set_of_nodes[nonleafnode].child_node.erase(itr_child);
    }

    int remaining_id;
    int i = 0;

    Node group1_nodes, group2_nodes;
    if (set_of_nodes[nonleafnode].parent_node == -1)
    {
        group1_nodes.id = node_id;
        group1_nodes.child_node = group1_node;
        group1_nodes.MBR = set_MBR(group1_nodes);
        group1_nodes.isleaf = false;
        group1_nodes.parent_node = nonleafnode;
        group1_nodes.depth = set_of_nodes[nonleafnode].depth + 1;
        set_of_nodes[group1_nodes.id] = group1_nodes;
        node_id++;
        group2_nodes.id = node_id;
        group2_nodes.child_node = group2_node;
        group2_nodes.MBR = set_MBR(group2_nodes);
        group2_nodes.isleaf = false;
        group2_nodes.parent_node = nonleafnode;
        group2_nodes.depth = set_of_nodes[nonleafnode].depth + 1;
        set_of_nodes[group2_nodes.id] = group2_nodes;
        node_id++;
        for (auto itr = group1_nodes.child_node.begin(); itr != group1_nodes.child_node.end(); itr++)
        {
            update_child_node(group1_nodes.id, itr->first);
        }
        for (auto itr = group2_nodes.child_node.begin(); itr != group2_nodes.child_node.end(); itr++)
        {
            update_child_node(group2_nodes.id, itr->first);
        }
        set_of_nodes[group1_nodes.id].parent_node = nonleafnode;
        set_of_nodes[group2_nodes.id].parent_node = nonleafnode;
        set_of_nodes[nonleafnode].child_node.clear();
        set_of_nodes[nonleafnode].child_node[group1_nodes.id] = group1_nodes.id;
        set_of_nodes[nonleafnode].child_node[group2_nodes.id] = group2_nodes.id;
    }
    else
    {
        group1_nodes.id = node_id;
        group1_nodes.child_node = group1_node;
        group1_nodes.MBR = set_MBR(group1_nodes);
        group1_nodes.isleaf = false;
        group1_nodes.parent_node = set_of_nodes[nonleafnode].parent_node;
        group1_nodes.depth = set_of_nodes[set_of_nodes[nonleafnode].parent_node].depth + 1;
        set_of_nodes[group1_nodes.id] = group1_nodes;
        node_id++;
        group2_nodes.id = node_id;
        group2_nodes.child_node = group2_node;
        group2_nodes.MBR = set_MBR(group2_nodes);
        group2_nodes.isleaf = false;
        group2_nodes.parent_node = set_of_nodes[nonleafnode].parent_node;
        group2_nodes.depth = set_of_nodes[set_of_nodes[nonleafnode].parent_node].depth + 1;
        set_of_nodes[group2_nodes.id] = group2_nodes;
        node_id++;
        for (auto itr = group1_nodes.child_node.begin(); itr != group1_nodes.child_node.end(); itr++)
        {
            set_of_nodes[itr->second].parent_node = group1_nodes.id;
        }
        for (auto itr = group2_nodes.child_node.begin(); itr != group2_nodes.child_node.end(); itr++)
        {
            set_of_nodes[itr->second].parent_node = group2_nodes.id;
        }
        set_of_nodes[set_of_nodes[nonleafnode].parent_node].child_node[group1_nodes.id] = group1_nodes.id;
        set_of_nodes[set_of_nodes[nonleafnode].parent_node].child_node[group2_nodes.id] = group2_nodes.id;
        set_of_nodes[set_of_nodes[nonleafnode].parent_node].child_node.erase(nonleafnode);
        set_of_nodes.erase(nonleafnode);
        if (set_of_nodes[set_of_nodes[group1_nodes.id].parent_node].child_node.size() > 2)
        {
            divide_nonleafnode(set_of_nodes[group1_nodes.id].parent_node);
        }
    }
}

void divide_leafnode(int node, Data z)
{
    vector<int> group1_data_id, group2_data_id;
    int max_xdist_id1, max_xdist_id2, max_ydist_id1, max_ydist_id2;
    value_type max_xdist, max_ydist;
    value_type xdist, ydist;
    value_type x_axis_dist, y_axis_dist;
    value_type norm_xdist, norm_ydist;
    x_axis_dist = abs(set_of_nodes[node].MBR.x_max - set_of_nodes[node].MBR.x_min);
    y_axis_dist = abs(set_of_nodes[node].MBR.y_max - set_of_nodes[node].MBR.y_min);

    max_xdist = 0;
    max_ydist = 0;
    vector<int> init_data_id = set_of_nodes[node].data_id;

    value_type area1, area2, area3, all_area, calc_area, calc_area1, calc_area2;
    value_type max_area = 0;
    Range calc_MBR1, calc_MBR2;

    all_area = (abs(set_of_nodes[node].MBR.x_max - set_of_nodes[node].MBR.x_min)) * (abs(set_of_nodes[node].MBR.y_max - set_of_nodes[node].MBR.y_min));

    group1_data_id.clear();
    group2_data_id.clear();

    for (int i = 0; i < init_data_id.size(); i++)
    {
        for (int j = 0; j < init_data_id.size(); j++)
        {
            xdist = abs(objects[init_data_id[i] - 1].x - objects[init_data_id[j] - 1].x);
            ydist = abs(objects[init_data_id[i] - 1].y - objects[init_data_id[j] - 1].y);
            if (xdist > max_xdist)
            {
                max_xdist = xdist;
                max_xdist_id1 = init_data_id[i];
                max_xdist_id2 = init_data_id[j];
            }
            if (ydist > max_ydist)
            {
                max_ydist = ydist;
                max_ydist_id1 = init_data_id[i];
                max_ydist_id2 = init_data_id[j];
            }
        }
    }
    norm_xdist = max_xdist / x_axis_dist;
    norm_ydist = max_ydist / y_axis_dist;
    if (norm_xdist > norm_ydist)
    {
        group1_data_id.push_back(max_xdist_id1);
        group2_data_id.push_back(max_xdist_id2);
    }
    else
    {
        group1_data_id.push_back(max_ydist_id1);
        group2_data_id.push_back(max_ydist_id2);
    }
    auto itr1 = std::find(init_data_id.begin(), init_data_id.end(), group1_data_id[0]);
    if (itr1 != init_data_id.end())
    {
        init_data_id.erase(itr1);
    }
    auto itr2 = std::find(init_data_id.begin(), init_data_id.end(), group2_data_id[0]);
    if (itr2 != init_data_id.end())
    {
        init_data_id.erase(itr2);
    }
    int remaining_id;

    for (int i = 0; i < init_data_id.size(); i++)
    {
        calc_MBR1 = set_MBR(group1_data_id, init_data_id[i]);
        calc_MBR2 = set_MBR(group2_data_id, init_data_id[i]);
        area1 = (abs(calc_MBR1.x_max - calc_MBR1.x_min)) * (abs(calc_MBR1.y_max - calc_MBR1.y_min));
        area2 = (abs(calc_MBR2.x_max - calc_MBR2.x_min)) * (abs(calc_MBR2.y_max - calc_MBR2.y_min));
        calc_area1 = all_area - area1;
        calc_area2 = all_area - area2;
        if (calc_area1 <= calc_area2)
        {
            group1_data_id.push_back(init_data_id[i]);
        }
        if (calc_area1 > calc_area2)
        {
            group2_data_id.push_back(init_data_id[i]);
        }
        /*if (i % 2 == 0)
		{
			group1_data_id.push_back(init_data_id[i]);
		}
		else
		{
			group2_data_id.push_back(init_data_id[i]);
		}*/
    }
    Node group1_node, group2_node;

    group1_node.data_id = group1_data_id;
    group1_node.MBR = set_MBR(group1_node.data_id);
    group1_node.isleaf = true;
    group1_node.parent_node = set_of_nodes[node].parent_node;
    group1_node.depth = set_of_nodes[set_of_nodes[node].parent_node].depth + 1;
    group1_node.id = node_id;
    node_id++;
    group2_node.data_id = group2_data_id;
    group2_node.MBR = set_MBR(group2_node.data_id);
    group2_node.isleaf = true;
    group2_node.parent_node = set_of_nodes[node].parent_node;
    group2_node.depth = set_of_nodes[set_of_nodes[node].parent_node].depth + 1;
    group2_node.id = node_id;
    node_id++;
    if (group1_node.innerMBR(z) == true)
    {
        if (group2_node.innerMBR(z) == true)
        {
            if (group1_node.area_MBR() < group2_node.area_MBR())
            {
                group1_node.data_id.push_back(z.id);
            }
            else
            {
                group2_node.data_id.push_back(z.id);
            }
        }
        else
        {
            group1_node.data_id.push_back(z.id);
        }
    }
    else
    {
        if (group2_node.innerMBR(z) == true)
        {
            group2_node.data_id.push_back(z.id);
        }
        else
        {
            if (group1_node.increase_MBR(z) < group2_node.increase_MBR(z))
            {
                group1_node.data_id.push_back(z.id);
            }
            else
            {
                group2_node.data_id.push_back(z.id);
            }
        }
    }
    set_of_nodes[group1_node.id] = group1_node;
    set_of_nodes[group2_node.id] = group2_node;
    set_of_nodes[set_of_nodes[node].parent_node].child_node[group1_node.id] = group1_node.id;
    set_of_nodes[set_of_nodes[node].parent_node].child_node[group2_node.id] = group2_node.id;
    set_of_nodes[set_of_nodes[node].parent_node].child_node.erase(node);
    set_of_nodes.erase(node);
    if (set_of_nodes[set_of_nodes[group1_node.id].parent_node].child_node.size() > 2)
    {

        divide_nonleafnode(set_of_nodes[group1_node.id].parent_node);
    }
};

void insert(int node, Data z)
{
    if (set_of_nodes[node].isleaf == true) //葉ノードの時
    {
        bool flag = false;
        flag = set_of_nodes[node].ismax();
        if (flag == true)
        {
            divide_leafnode(node, z);
        }
        else
        {
            set_of_nodes[node].data_id.push_back(z.id);
            set_of_nodes[node].set_MBR(z);
        }
    }
    else //葉ノード以外の時
    {
        vector<int> node_id;
        for (auto itr = set_of_nodes[node].child_node.begin(); itr != set_of_nodes[node].child_node.end(); itr++)
        {
            Node temp = set_of_nodes[itr->second];
            if (temp.innerMBR(z) == true)
            {
                node_id.push_back(itr->first);
            }
        }
        if (node_id.empty() == true)
        {
            int insert_node_id;
            value_type min_area = 10000000;
            for (auto itr = set_of_nodes[node].child_node.begin(); itr != set_of_nodes[node].child_node.end(); itr++)
            {
                Node temp = set_of_nodes[itr->second];
                value_type area = temp.increase_MBR(z);
                if (area < min_area)
                {
                    min_area = area;
                    insert_node_id = itr->first;
                }
            }
            insert(insert_node_id, z);
        }
        else
        {
            value_type min_area = 100000000;
            value_type min_value;
            int insert_node_id;
            for (auto itr = node_id.begin(); itr != node_id.end(); itr++)
            {
                auto temp = set_of_nodes[node].child_node.begin();
                min_value = set_of_nodes[temp->second].area_MBR();
                if (min_value < min_area)
                {
                    min_area = min_value;
                    insert_node_id = temp->first;
                }
            }
            insert(insert_node_id, z);
        }
    }
};

int find_leafnode(int node, Data z)
{
    int x;
    int y;
    for (int i = 0; i < set_of_nodes[node].data_id.size(); i++)
    {
        x = set_of_nodes[node].data_id[i];
        if (set_of_nodes[node].data_id[i] == z.id)
        {
            y = node;
            break;
        }
        else
        {
            y = -1;
        }
    }
    return y;
}

int find_node(int node, Data z)
{
    int f = -1;
    for (auto itr = set_of_nodes[node].child_node.begin(); itr != set_of_nodes[node].child_node.end(); itr++)
    {
        if (set_of_nodes[itr->second].isleaf == true)
        {
            f = find_leafnode(itr->second, z);
            if (f != -1)
            {
                return f;
            }
        }
        else
        {
            if (set_of_nodes[itr->second].innerMBR(z) == true)
            {
                f = find_node(itr->second, z);
                if (f != -1)
                {
                    break;
                }
            }
        }
    }
    return f;
}
void get_all_node(int node)
{
    for (auto itr = set_of_nodes[node].child_node.begin(); itr != set_of_nodes[node].child_node.end(); itr++)
    {
        if (set_of_nodes[itr->second].isleaf == true)
        {
            child_node_id2.push_back(set_of_nodes[itr->second].id);
        }
        else
        {
            get_all_node(set_of_nodes[itr->second].id);
        }
    }
}

void get_leaf_node_id(Node node)
{
    if (node.isleaf == false)
    {
        for (auto itr = node.child_node.begin(); itr != node.child_node.end(); itr++)
        {
            get_leaf_node_id(set_of_nodes[itr->first]);
        }
    }
    if (node.isleaf == true)
    {
        child_node_id.push_back(node.id);
    }
};

value_type dist_nodes(Node node1, Node node2)
{
    value_type dist = 0;
    if (node1.MBR.x_max < node2.MBR.x_min && node1.MBR.y_min > node2.MBR.y_max)
    {
        dist = (sqrt(pow(abs(node1.MBR.x_max - node2.MBR.x_min), 2.0) + pow(abs(node1.MBR.y_min - node2.MBR.y_max), 2.0)));
    }
    if (node1.MBR.x_max < node2.MBR.x_min && node1.MBR.y_min < node2.MBR.y_max && node1.MBR.y_max > node2.MBR.y_min)
    {
        dist = (abs(node1.MBR.x_max - node2.MBR.x_min));
    }
    if (node1.MBR.x_max < node2.MBR.x_min && node1.MBR.y_max < node2.MBR.y_min)
    {
        dist = (sqrt(pow(abs(node1.MBR.x_max - node2.MBR.x_min), 2.0) + pow(abs(node1.MBR.y_max - node2.MBR.y_min), 2.0)));
    }
    if (node1.MBR.x_max > node2.MBR.x_min && node1.MBR.x_min < node2.MBR.x_max && node1.MBR.y_max < node2.MBR.y_min)
    {
        dist = (abs(node1.MBR.y_max - node2.MBR.y_min));
    }
    if (node1.MBR.x_min > node2.MBR.x_max && node1.MBR.y_max < node2.MBR.y_min)
    {
        dist = (sqrt(pow(abs(node1.MBR.x_max - node2.MBR.x_min), 2.0) + pow(abs(node1.MBR.y_max - node2.MBR.y_min), 2.0)));
    }
    if (node1.MBR.x_min > node2.MBR.x_max && node1.MBR.y_min < node2.MBR.y_max && node1.MBR.y_max > node2.MBR.y_min)
    {
        dist = (abs(node1.MBR.x_min - node2.MBR.x_max));
    }
    if (node1.MBR.x_min > node2.MBR.x_max && node1.MBR.y_min > node2.MBR.y_max)
    {
        dist = (sqrt(pow(abs(node1.MBR.x_min - node2.MBR.x_max), 2.0) + pow(abs(node1.MBR.y_min - node2.MBR.y_max), 2.0)));
    }
    if (node1.MBR.x_max > node2.MBR.x_min && node1.MBR.x_min < node2.MBR.x_max && node1.MBR.y_min > node2.MBR.y_max)
    {
        dist = (abs(node1.MBR.y_min - node2.MBR.y_max));
    }
    return dist;
}

value_type dist_data_to_node(Data data, Node node2)
{
    value_type dist = 0;
    if (data.x < node2.MBR.x_min && data.y > node2.MBR.y_max)
    {
        dist = (sqrt(pow(abs(data.x - node2.MBR.x_min), 2.0) + pow(abs(data.y - node2.MBR.y_max), 2.0)));
    }
    if (data.x < node2.MBR.x_min && data.y < node2.MBR.y_max && data.y > node2.MBR.y_min)
    {
        dist = (abs(data.x - node2.MBR.x_min));
    }
    if (data.x < node2.MBR.x_min && data.y < node2.MBR.y_min)
    {
        dist = (sqrt(pow(abs(data.x - node2.MBR.x_min), 2.0) + pow(abs(data.y - node2.MBR.y_min), 2.0)));
    }
    if (data.x > node2.MBR.x_min && data.x < node2.MBR.x_max && data.y < node2.MBR.y_min)
    {
        dist = (abs(data.y - node2.MBR.y_min));
    }
    if (data.x > node2.MBR.x_max && data.y < node2.MBR.y_min)
    {
        dist = (sqrt(pow(abs(data.x - node2.MBR.x_min), 2.0) + pow(abs(data.y - node2.MBR.y_min), 2.0)));
    }
    if (data.x > node2.MBR.x_max && data.y < node2.MBR.y_max && data.y > node2.MBR.y_min)
    {
        dist = (abs(data.x - node2.MBR.x_max));
    }
    if (data.x > node2.MBR.x_max && data.y > node2.MBR.y_max)
    {
        dist = (sqrt(pow(abs(data.x - node2.MBR.x_max), 2.0) + pow(abs(data.y - node2.MBR.y_max), 2.0)));
    }
    if (data.x > node2.MBR.x_min && data.x < node2.MBR.x_max && data.y > node2.MBR.y_max)
    {
        dist = (abs(data.y - node2.MBR.y_max));
    }
    if (data.x < node2.MBR.x_max && data.x > node2.MBR.x_min && data.y < node2.MBR.y_max && data.y > node2.MBR.y_min)
    {
        dist = 0;
    }
    return dist;
}

int share_parent_node(Node node1, Node node2)
{
    if (node1.parent_node == node2.parent_node)
    {
        return node1.parent_node;
    }
    else
    {

        share_parent_node(set_of_nodes[node1.parent_node], set_of_nodes[node2.parent_node]);
    }
}
