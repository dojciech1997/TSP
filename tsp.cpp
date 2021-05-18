#include "tsp.hpp"

VEC_int tsp(const MAT_double &my_cost_matrix) {
    TSP_cost_matrix MyMatrix(my_cost_matrix);
    u_int to_sort_uint = my_cost_matrix.size()-1;
    int to_sort = static_cast<int>(to_sort_uint);
    for (u_int i = 0; i < my_cost_matrix.size(); i++) {
        MyMatrix.verify_and_reduce_all_rows();
        MyMatrix.verify_and_reduce_all_cols();
        if (MyMatrix.cost_matrix.size() < 3) {
            MyMatrix.if_2x2();
            break;
        }
        MyMatrix.find_best_path();
    }
    MyMatrix.final_solution.push_back(MyMatrix.solution[0].first);
    int x = MyMatrix.solution[0].second;
    do {
        for (auto &e : MyMatrix.solution) {
            if (e != MyMatrix.solution[0]) {
                if (e.first == x) {
                    x = e.second;
                    MyMatrix.final_solution.push_back(e.first);
                    to_sort--;
                }
            }
        }
    }
    while (to_sort != 0);
    MyMatrix.final_solution.push_back(MyMatrix.solution[0].first);
    VEC_int final_solution_int;
    for (auto& e : MyMatrix.final_solution){
        int d = static_cast<int>(e);
        final_solution_int.push_back(d);
    }
    return final_solution_int;
}

VEC_int TSP_cost_matrix::create_row_index(u_int n){
    VEC_int v;
    for (u_int i = 1; i <= n; i++) {
        v.push_back(i);
    }
    return v;
}

VEC_int TSP_cost_matrix::create_col_index(u_int n){
    VEC_int v;
    for (u_int i = 1; i <= n; i++) {
        v.push_back(i);
    }
    return v;
}

double get_forbidden_cost() {
    return INF;
}

void TSP_cost_matrix::reduce_all_rows(){
    double minimum;
    for (u_int i=0; i<index_vector_row.size();i++){
        for (u_int e=0; e<index_vector_col.size();e++) {
            if (std::isnan(cost_matrix[i][e])) {
                continue;
            }
            minimum = cost_matrix[i][e];
            break;
        }
        for (u_int u=0; u<index_vector_col.size();u++) {
            if (std::isnan(cost_matrix[i][u])) {
                continue;
            }
            if (cost_matrix[i][u] < minimum) {
                minimum = cost_matrix[i][u];
            }
        }
        for (u_int z=0; z<index_vector_col.size();z++) {
            cost_matrix[i][z] -= minimum;
        }
        lower_bound += minimum;
    }
}
void TSP_cost_matrix::reduce_all_cols(){
    double minimum;
    for (u_int i=0; i<index_vector_col.size();i++){
        for (u_int e=0; e<index_vector_row.size();e++) {
            if (std::isnan(cost_matrix[e][i])) {
                continue;
            }
            minimum = cost_matrix[e][i];
            break;
        }
        for (u_int u=0; u<index_vector_row.size();u++) {
            if (std::isnan(cost_matrix[u][i])) {
                continue;
            }
            if (cost_matrix[u][i] < minimum) {
                minimum = cost_matrix[u][i];
            }
        }
        for (u_int z=0; z<index_vector_row.size();z++) {
            cost_matrix[z][i] -= minimum;
        }
        lower_bound += minimum;
    }
}
void TSP_cost_matrix::verify_and_reduce_all_cols(){
    int iszero = 0;
    for (u_int i = 0; i < index_vector_col.size(); i++){
        for (u_int e = 0; e < index_vector_row.size(); e++) {
            if (cost_matrix[e][i] == 0) {
                iszero++;
            }
        }
        if (iszero == 0) {
            reduce_all_cols();
            break;
        }
        iszero = 0;
    }
}

void TSP_cost_matrix::verify_and_reduce_all_rows(){
    int iszero = 0;
    for (u_int i = 0; i < index_vector_row.size(); i++){
        for (u_int e = 0; e < index_vector_col.size(); e++) {
            if (cost_matrix[i][e] == 0) {
                iszero++;
            }
        }
        if (iszero == 0) {
            reduce_all_rows();
            break;
        }
        iszero = 0;
    }
}

void TSP_cost_matrix::find_best_path() {
    VEC_doublePAIR opportunities;
    doublePAIR path;
    doublePAIR best_path;
    double minimum;
    double minimum1;
    double minimum2;
    double maximum;
    PAIR_double_doublePAIR minimum_pair;
    VEC_PAIR_int_doublePAIR vector_minimum;
    for (u_int i=0; i<index_vector_row.size();i++){
        for (u_int e=0; e<index_vector_col.size();e++) {
            if (cost_matrix[i][e] == 0) {
                path.first = index_vector_row[i];
                path.second = index_vector_col[e];
                opportunities.emplace_back(path);
            }
        }
    }
    for (auto& e : opportunities) {
        for (u_int a = 0; a < index_vector_row.size(); a++) {
            if (index_vector_row[a] == e.first) {
                for (u_int u = 0; u < index_vector_col.size(); u++) {
                    if (index_vector_col[u] == e.second) {
                        continue;
                    }
                    if (std::isnan(cost_matrix[a][u])) {
                        continue;
                    }
                    minimum1 = cost_matrix[a][u];
                    break;
                }
            }
        }
        for (u_int a = 0; a < index_vector_row.size(); a++) {
            if (index_vector_row[a] == e.first) {
                for (u_int i = 0; i < index_vector_col.size(); i++) {
                    if (index_vector_col[i] == e.second) {
                        continue;
                    }
                    if (std::isnan(cost_matrix[a][i])) {
                        continue;
                    }
                    if (cost_matrix[a][i] < minimum1) {
                        minimum1 = cost_matrix[a][i];
                    }
                }
            }
        }
        for (u_int a = 0; a < index_vector_col.size(); a++) {
            if (index_vector_col[a] == e.second) {
                for (u_int u = 0; u < index_vector_row.size(); u++) {
                    if (index_vector_col[u] == e.first) {
                        continue;
                    }
                    if (std::isnan(cost_matrix[u][a])) {
                        continue;
                    }
                    minimum2 = cost_matrix[u][a];
                    break;
                }
            }
        }
        for (u_int a = 0; a < index_vector_col.size(); a++) {
            if (index_vector_col[a] == e.second) {
                for (u_int i = 0; i < index_vector_row.size(); i++) {
                    if (index_vector_col[i] == e.first) {
                        continue;
                    }
                    if (std::isnan(cost_matrix[i][a])) {
                        continue;
                    }
                    if (cost_matrix[i][a] < minimum2) {
                        minimum2 = cost_matrix[i][a];
                    }
                }
            }
        }
        minimum = minimum1 + minimum2;
        minimum_pair.first = minimum;
        minimum_pair.second = e;
        vector_minimum.push_back(minimum_pair);
    }
    maximum = 0;
    for (auto& e : vector_minimum) {
        if (e.first > maximum) {
            maximum = e.first;
            best_path = e.second;
        }
    }
    for (u_int a = 0; a < index_vector_row.size(); a++) {
        for (u_int i = 0; i < index_vector_col.size(); i++) {
            if (index_vector_col[a] == best_path.first && index_vector_row[i] == best_path.second) {
                cost_matrix[i][a] = get_forbidden_cost();
            }
        }
    }
    for (u_int a = 0; a < index_vector_row.size(); a++) {
        for (u_int i = 0; i < index_vector_col.size(); i++) {
            if (index_vector_row[a] == best_path.first && index_vector_col[i] == best_path.second) {
                cost_matrix.erase(cost_matrix.begin() + a);
                for (u_int u = 0; u < index_vector_row.size()-1; u++) {
                    cost_matrix[u].erase(cost_matrix[u].begin() + i);
                }
                index_vector_row.erase(index_vector_row.begin() + a);
                index_vector_col.erase(index_vector_col.begin() + i);
                solution.push_back(best_path);
                break;
            }
        }
    }
}

void TSP_cost_matrix::if_2x2() {
    verify_and_reduce_all_rows();
    verify_and_reduce_all_cols();
    for (u_int i = 0; i < cost_matrix.size(); i++){
        for (u_int u = 0; u < cost_matrix[i].size(); u++){
            if (!std::isnan(cost_matrix[i][u])){
                solution.emplace_back(std::make_pair(index_vector_row[i], index_vector_col[u]));
                cost_matrix.erase(cost_matrix.begin()+i);
                cost_matrix[0].erase(cost_matrix[0].begin()+u);
                index_vector_col.erase(index_vector_col.begin()+u);
                index_vector_row.erase(index_vector_row.begin()+i);
                if (!std::isnan(cost_matrix[0][0])) {
                    solution.emplace_back(std::make_pair(index_vector_row[0], index_vector_col[0]));
                }
                break;
            }
        }
    }
}
