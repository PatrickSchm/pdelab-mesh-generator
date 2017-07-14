std::map<int, int> readBulkPairs(std::string path) {

    std::string str;
    std::string vstr;
    if (!std::ifstream(path))
    {
        std::cout << "File: " << path << " not found!" << std::endl;
    }
    ifstream ifs (path);

    std::map<int, int> bulkPairs;
    while (std::getline(ifs, str)) {
        istringstream sv(str);
        int counter = 0;
        int mapIndx = -1;
        int pair = -1;
        while (std::getline(sv, vstr, ',')) {
            if (counter == 0) {
                mapIndx = stoi(vstr);
            } else {
                pair = stoi(vstr);
            }
            counter++;
        }
        bulkPairs[mapIndx] = pair;
    }
    return bulkPairs;
}

std::map<int, std::vector<int>> readBulkFracturePairs(std::string path) {

    std::string str;
    std::string vstr;
    if (!std::ifstream(path))
    {
        std::cout << "File: " << path << " not found!" << std::endl;
    }
    ifstream ifs (path);
    std::map<int, std::vector<int>> bulkPairs;
    while (std::getline(ifs, str)) {
        istringstream sv(str);
        int counter = 0;
        int mapIndx = -1;
        int pair = -1;
        std::vector<int> pair_dummy;
        while (std::getline(sv, vstr, ',')) {
            if (counter == 0) {
                mapIndx = stoi(vstr);
            } else {
                pair_dummy.push_back(stoi(vstr));
            }
            counter++;
        }
        bulkPairs[mapIndx] = pair_dummy;
    }
    return bulkPairs;
}