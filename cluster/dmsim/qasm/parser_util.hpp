#ifndef QASM_PARSER_UTIL
#define QASM_PARSER_UTIL

#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

/**************UTILITY FUNCTION DEFINIATIONS********************/
/**
 * @brief Split the given string based on the provided delim
 *
 * @param s
 * @param delim
 * @return vector<string>
 */
vector<string> split(const string &s, char delim);

/**
 * @brief Perform the division in case there is any in the expression
 *
 * @param str
 * @return double
 */
double get_value(std::string str);

/**
 * @brief Convert the expression to values
 *
 * @param str The string representation of the expression
 * @return double
 */
double get_param_value(string str);

/**
 * @brief Get the index of the target in the vector
 *
 * @param vec The vector that contains the target object
 * @param target Target to look for
 * @return int The found index of the target in the list
 */
int get_index(vector<string> vec, string target);

/**
 * @brief Print measurement outcomes
 *
 * @param counts Pointer to the counts dictionary.
 * @param repetition Number of shots performed.
 */
void print_counts(map<string, int> *counts, int repetition);

char *getCmdOption(char **begin, char **end, const std::string &option);
bool cmdOptionExists(char **begin, char **end, const std::string &option);

/************************** IMPLEMENTATION OF UTILITY FUNCTIONS **************************/

vector<string> split(const string &s, char delim)
{
    vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}

double get_value(std::string str)
{
    vector<string> numbers = split(str, '/');

    if (numbers.size() == 1)
    {
        return std::stod(numbers[0]);
    }
    else
    {
        return std::stod(numbers[0]) / std::stod(numbers[1]);
    }
}

double get_param_value(string str)
{
    // replace the pi with a decimal value
    string s = regex_replace(str, regex("PI"), "3.1415926");

    vector<string> vals = split(s, '*');

    if (vals.size() == 1)
    {
        return get_value(vals[0]);
    }
    else
    {
        return get_value(vals[0]) * get_value(vals[1]);
    }
}

int get_index(vector<string> vec, string target)
{
    for (int i = 0; i < vec.size(); i++)
    {
        if (vec[i] == target)
        {
            return i;
        }
    }
    return -1;
}

void print_counts(map<string, int> *counts, int repetition)
{
    assert(counts != NULL);
    printf("\n===============  Measurement (tests=%d) ================\n", repetition);

    for (auto const &[key, val] : *counts)
    {
        printf("\"%s\" : %d\n", key.c_str(), val);
        //printf("\"%s\" : %lf\n", key.c_str(), (double)val/(double)repetition);
    }
}

char *getCmdOption(char **begin, char **end, const string &option)
{
    char **itr = find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char **begin, char **end, const string &option)
{
    return find(begin, end, option) != end;
}

#endif
