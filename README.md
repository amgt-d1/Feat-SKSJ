# Feat-SKSJ: an algorithm for top-k spatio-textual similarity join
## Data sources
- [Twitter](https://personal.ntu.edu.sg/gaocong/datacode.htm)
- [Places](https://archive.org/details/2011-08-SimpleGeo-CC0-Public-Spaces)

## Data File
- Twitter:  
    put "Tweet_sort" on /data.  
    file format (Delimiter : tab)  
        [Date, Object keywords(Delimiter : comma), Latitude, Longitude]  
    ```
    // exmple Tweet_sort
    2012-4-1 0:10:0	S,O,Bro,	38.4844	-75.8034
    2012-4-1 0:10:0	o,school,shit,Yay,	40.5207	-74.3624
    ```
    
- Places:  
    put "poi_place" on /data.  
    file format (Delimiter : tab)  
        [Object keywords(Delimiter : comma), Latitude, Longitude]  
    ```
    // exmple poi_place 
    Stoneham,theatrical,agency,	42.481499	-71.098735
    Hialeah,structural,engineer,	25.891515	-80.328017
    ```

## How to build
compiler: g++ 5.4.0  
command:
- proposed method: g++ -std=gnu++1y -o problem.out proposed.cpp set_data.cpp Rtree.cpp -O3
- baseline method: g++ -std=gnu++1y -o problem.out baseline.cpp set_data.cpp Rtree.cpp -O3
- naive method: g++ -std=gnu++1y -o problem.out lenear.cpp set_data.cpp Rtree.cpp -O3

## How to run
./problem.out

## How to tune parameters
### Set each parameter in class.h
```
// class.h  line 23
#define TOP_K 100
#define num 1000000
#define alpha 0.5
#define NODE_MAX 1200
#define L TOP_K
```
- TOP_K: The number of object pairs to output. The default is 100．
- num: The number of objects to retrive．The default is 1,000,000．
- alpha: The parameter to leverage the textual similarity and spatial similarity．The default is 0.5．
- NODE_MAX: The maximum number of objects contained in one leaf node. 
- L: The number of leaf nodes for similarity calclation in Init Threshold Calclation．

### Dataset setting
```
// proposed.cpp baseline.cpp lenear.cpp

set_data(objects2);	//Twitter
// set_data_place(objects2);	//Place
```
- Twitter: set_data(objects2);
- Places: set_data_place(objects2);

### Textual similarity setting
Change "score" in calc_textual_score (RTree.cpp line 23)
```
// RTree.cpp line 23
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
```
- Jaccard similarity: score = (double)num_key / (set_key.size());
- Cosine similarity: score = (double)num_key / sqrt(abs(z1.key.size() * z2.key.size()));
- Dice simialrity: score = (double)2 * num_key / (z1.key.size() + z2.key.size());

When executing the proposed method, change "jaccard" on proposed.cpp line 192 to 194
```
// proposed.cpp line 192
jaccard /= objects[itr->second.data_id[i] - 1].key.size() + objects[itr->second.data_id[j] - 1].key.size() - num_share_key; //Jaccard
// jaccard /= sqrt(objects[itr->second.data_id[i] - 1].key.size() * objects[itr->second.data_id[j] - 1].key.size());	//Cosine
// jaccard /= (objects[itr->second.data_id[i] - 1].key.size() + objects[itr->second.data_id[j] - 1].key.size());	//Dice
```
- Jaccard similarity: jaccard /= objects[itr->second.data_id[i] - 1].key.size() + objects[itr->second.data_id[j] - 1].key.size() - num_share_key;
- Cosine similarity: jaccard /= sqrt(objects[itr->second.data_id[i] - 1].key.size() * objects[itr->second.data_id[j] - 1].key.size());
- Dice similarity: jaccard /= (objects[itr->second.data_id[i] - 1].key.size() + objects[itr->second.data_id[j] - 1].key.size());

### Note
- This code is owned by one of my former master students Shohei Tsuruoka.
- If you have questions, contact `tsuruoka.shohei@ist.osaka-u.ac.jp`.
