#include "class.h"

//データ読み込み(Twitter)
void set_data(std::vector<Data> &z)
{
	char *ary[5];
	char *ary2[5];
	char *ctx;
	char key[MAX] = {0};
	char key2[40] = {0};
	char buf[MAX] = {0};
	unsigned int i = 0;
	unsigned int j = 0;
	int cons;
	string compare;

	cout << "Twitter" << endl;

	ifstream ifs("data/Tweet_sort");
	while (ifs.getline(buf, MAX))
	{
		if (i == num)
		{
			break;
		}
		z[i].id = i + 1;
		j = 0;
		ary2[0] = strtok(buf, "	");
		ary2[1] = strtok(NULL, "	");
		ary2[2] = strtok(NULL, "	");
		ary2[3] = strtok(NULL, "	");
		z[i].x = atof(ary2[2]);
		z[i].y = atof(ary2[3]);
		ary[0] = strtok(ary2[1], ",");
		memcpy(key2, ary[0], sizeof(key2));
		z[i].key.push_back(ary[0]);
		j++;
		while (true)
		{
			ary[1] = strtok(NULL, ",");
			if (ary[1] == NULL)
			{
				break;
			}
			cons = 0;
			compare = ary[1];
			for (int k = 0; k < z[i].key.size(); k++)
			{
				if (z[i].key[k] == compare)
				{
					cons = 1;
				}
			}
			if (cons == 0)
			{
				z[i].key.push_back(ary[1]);
				j++;
			}
		}
		i++;
	}
	ifs.close();
}

//データ読み込み(Places)
void set_data_place(std::vector<Data> &z)
{
	char *ary[5];
	char *ary2[5];
	char *ctx;
	char key[MAX] = {0};
	char key2[40] = {0};
	char buf[MAX] = {0};
	unsigned int i = 0;
	unsigned int j = 0;
	int cons;
	string compare;

	cout << "Place" << endl;

	ifstream ifs("data/poi_place");
	while (ifs.getline(buf, MAX))
	{
		if (i == num)
		{
			break;
		}
		z[i].id = i + 1;
		j = 0;
		ary2[0] = strtok(buf, "	");
		ary2[1] = strtok(NULL, "	");
		ary2[2] = strtok(NULL, "	");
		z[i].x = atof(ary2[1]);
		z[i].y = atof(ary2[2]);
		ary[0] = strtok(ary2[0], ",");
		memcpy(key2, ary[0], sizeof(key2));
		z[i].key.push_back(ary[0]);
		j++;
		while (true)
		{
			ary[1] = strtok(NULL, ",");
			if (ary[1] == NULL)
			{
				break;
			}
			cons = 0;
			compare = ary[1];
			for (int k = 0; k < z[i].key.size(); k++)
			{
				if (z[i].key[k] == compare)
				{
					cons = 1;
				}
			}
			if (cons == 0)
			{
				z[i].key.push_back(ary[1]);
				j++;
			}
		}
		i++;
	}
	ifs.close();
}


