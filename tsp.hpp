#ifndef KOMIWOJAZER_TSP_H
#define KOMIWOJAZER_TSP_H
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <functional>
#include <cfloat>
#include <memory>
#include <cstdlib>
#include <utility>
#include <cfloat>
#include <algorithm>
#define INF (NAN)

using VEC_int = std::vector<int>;
using VEC_double = std::vector<double>;
using MAT_double = std::vector<std::vector<double>>;
using u_int = long long unsigned int;
using VEC_doublePAIR = std::vector<std::pair<double,double>>;
using doublePAIR = std::pair<double,double>;
using PAIR_double_doublePAIR = std::pair<double,std::pair<double,double>>;
using VEC_PAIR_int_doublePAIR  = std::vector<std::pair<int,std::pair<double,double>>>;

class TSP_cost_matrix {
public:
    //TSP_cost_matrix() = default;
    TSP_cost_matrix(std::vector<std::vector<double>> cost_matrix_)
            : cost_matrix(cost_matrix_)
            , index_vector_row(create_row_index(cost_matrix_.size()))
            , index_vector_col(create_col_index(cost_matrix_[0].size()))
    {}
    //this->cost_matrix = cost_matrix;
    //this->index_vector_col = index_vector_col;
    //this->index_vector_row = index_vector_row;
    std::vector<int> create_row_index(long long unsigned int n);
    std::vector<int> create_col_index(long long unsigned int n);
    void reduce_all_rows();
    void reduce_all_cols();
    void verify_and_reduce_all_cols();
    void verify_and_reduce_all_rows();
    void if_2x2();
    void find_best_path();
    double lower_bound{0};
    std::vector<std::pair<double,double>> solution;
    std::vector<std::vector<double>> cost_matrix{{}};
    std::vector<int> index_vector_row;
    std::vector<int> index_vector_col;
    std::vector<double> final_solution{};
    ~TSP_cost_matrix() = default;
};

std::vector<int> tsp(const std::vector<std::vector<double>> &cost_matrix);
double get_forbidden_cost();

#endif //KOMIWOJAZER_TSP_H
